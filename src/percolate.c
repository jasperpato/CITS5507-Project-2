/*
 * CITS5507 HPC PROJECT 2
 * LATTICE PERCOLATION IN PARALLEL
 * 
 * Jasper Paterson 22736341
 * Allen Antony 22706998
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <util.h>

#include <omp.h>
#include <mpi.h>

#include "../include/constant.h"
#include "../include/stack.h"
#include "../include/site.h"
#include "../include/bond.h"
#include "../include/cluster.h"

/**
 * @return int start index of the region of the lattice allocated to this id.
 * Allocates n/n_threads+1 rows to the first n%n_threads ids and n/n_threads rows to the remaining ids.
 */
static int start_index(int n, int id, int num_ids)
{
  return id < n%num_ids ? n*id*(n/num_ids+1) : n*(n%num_ids)*(n/num_ids+1) + n*(id-n%num_ids)*(n/num_ids);
}

/**
 * @return int end index of the region of the lattice allocated to this id.
 */
static int end_index(int n, int id, int num_ids)
{
  return id < n%num_ids ? n*(id+1)*(n/num_ids+1) : n*(n%num_ids)*(n/num_ids+1) + n*(id+1-n%num_ids)*(n/num_ids);
}

/**
 * @return int number of nodes allocated to this id
 */
static int get_count(int n, int id, int num_ids)
{
  return id < n%num_ids ? n/num_ids+1 : n/num_ids;
}

/**
 * @brief Check if a site has at least one bond to a neighbour (and can therefore be the start of a cluster)
 *        Returns true even if neighbour is across the thread boundary (because the clusters need to be joined later).
 *        Only used in bond percolation.
 * @return short boolean true iff has at least one neighbour
 */
static short has_neighbours(Bond* b, int n, Site* s)
{
  // indices of bonds
  int ib[] = {
    s->r*n+s->c,          // left
    s->r*n+(s->c+n+1)%n,  // right
    s->r*n+s->c,          // top
    ((s->r+n+1)%n)*n+s->c // bottom
  };
  for(int i = 0; i < 4; ++i) {
    if((i<2 && b->h[ib[i]]) || (i>=2 && b->v[ib[i]])) return 1; 
  }
  return 0;
}

/** 
 * @return short boolean true iff site lies on a thread boundary (and therefore needs to maintain latest cluster info for possible future cluster joins)
 */
static short on_border(int n, int idx, int n_threads) {
  for(int i = 0; i < n_threads; ++i) {
    int start = start_index(n, i, n_threads); // region boundaries
    int end = end_index(n, i, n_threads);
    if((idx >= start && idx < start+n) || (idx >= end-n && idx < end)) return 1;
  }
  return 0;
}

/**
 * @brief get bottom neighbour to site s, intended for crossing the thread boundary in a cluster join.
 * @return site pointer to bottom neighbour if connected and separate cluster, else NULL
 */
static Site* bottom_neighbour(Site* a, Bond* b, int n, Site* s)
{
  int i = ((s->r+n+1)%n)*n+s->c; // index of bottom neighbour in a and bottom bond in b
  Site *nb = &a[i];
  if((!b && !nb->occupied) || (b && !b->v[i])) return NULL;
  if(s->cluster->id == nb->cluster->id) return NULL;
  return nb;
}

/**
 * @return Site** array of four neighbours' site pointers. Pointer is either a valid address
 *         if there is a connection to that neighbour, the neighbour is unseen and it is within thread bounds, otherwise NULL
 */
static void get_neighbours(Site* a, Bond* b, int n, int n_threads, Site* s, Site* nbs[], int start, int end)
{
  // indices of neighbours
  int is[] = {
    s->r*n+(s->c+n-1)%n,   // left
    s->r*n+(s->c+n+1)%n,   // right
    ((s->r+n-1)%n)*n+s->c, // top
    ((s->r+n+1)%n)*n+s->c  // bottom
  };
  // indices of bonds
  int ib[] = {
    s->r*n+s->c,          // left
    s->r*n+(s->c+n+1)%n,  // right
    s->r*n+s->c,          // top
    ((s->r+n+1)%n)*n+s->c // bottom
  };
  for(int i = 0; i < 4; ++i) {
    Site* nb = &a[is[i]];
    if(is[i] < start || is[i] >= end) nbs[i] = NULL; // out of thread bounds
    else if(!nb->seen && (
      (!b && nb->occupied) ||
      (b && ((i<2 && b->h[ib[i]]) || (i>=2 && b->v[ib[i]]))))
    ) {
      nb->seen = 1; // must mark site here for n=2 case when top neighbour == bottom neighbour
      nbs[i] = nb;
    } else nbs[i] = NULL;
  }
}

/**
 * @brief loop through neighbours and update cluster information
 */
static void DFS(Site* a, Bond* b, int n, int n_threads, Stack* st, int start, int end) {
  while(!is_empty(st)) {
    Site *s = pop(st);
    Site *nbs[4];
    get_neighbours(a, b, n, n_threads, s, nbs, start, end);
    Cluster* cl = s->cluster;
    for(int i = 0; i < 4; ++i) { // loop through connected, unseen neighbours
      if(!nbs[i]) continue; 
      Site* nb = nbs[i];
      nb->cluster = cl;
      ++cl->size;
      if(!cl->rows[nb->r]) cl->height++;
      if(!cl->cols[nb->c]) cl->width++;
      cl->rows[nb->r] = 1;
      cl->cols[nb->c] = 1;
      int idx = nb->r*n+nb->c;
      if(on_border(n, idx, n_threads)) cl->sites[cl->site_size++] = idx; // add index of neighbour to cluster's site index array
      add(st, nb);
    }
  }
}

/**
 * @brief Simulate percolation using DFS in a given region of the lattice.
 * @param a site array
 * @param b bond struct, either valid address if [-b] or NULL if [-s]
 */
static void percolate(Site* a, Bond* b, int n, int n_threads, CPArray* cpa, short tid)
{
  int start = start_index(n, tid, n_threads);
  int end = end_index(n, tid, n_threads);
  if(tid+1 == n_threads) end = n*n; // last region may be extended to reach end of lattice
  Stack* st = stack(n*n);
  for(int i = start; i < end; ++i) {
    Site *s = &a[i];
    if(!s->seen && ((!b && s->occupied) || (b && has_neighbours(b, n, s)))) { // if unseen and will form a cluster
      s->seen = 1;
      s->cluster = cluster(n, n_threads, i/n, i%n);
      Cluster *sc = s->cluster;
      int idx = s->r*n+s->c;
      if(on_border(n, idx, n_threads)) sc->sites[sc->site_size++] = idx;
      cpa->cls[cpa->size++] = sc; 
      add(st, s);
      DFS(a, b, n, n_threads, st, start, end);
    }   
  }
}

/**
 * @brief merges clusters along the bottom row of each thread region.
 *        Update cluster information and maintains a list of sites belonging to the cluster that are on a thread boundary.
 *        Only sites on a thread boundary need to have latest cluster information for future joins.
 */
static void join_clusters(Site* a, Bond* b, int n, int n_threads) {
  for(int i = 0; i < n_threads; ++i) {
    int row_end = end_index(n, i, n_threads);
    int row_start = row_end-n;
    for(int i = row_start; i < row_end; ++i) { // loop along bottom row of region
      Site *s = &a[i];
      Cluster *sc = s->cluster;
      if(!sc) continue;
      Site *nb = bottom_neighbour(a, b, n, s);
      if(!nb) continue;
      Cluster *nc = nb->cluster;
      // combine two clusters into sc
      sc->size += nc->size;
      for(int i = 0; i < n; ++i) {
        if(nc->rows[i]) {
          if(!sc->rows[i]) sc->height++;
          sc->rows[i] = 1;
        }
        if(nc->cols[i]) {
          if(!sc->cols[i]) sc->width++;
          sc->cols[i] = 1;
        }
      }
      for(int j = 0; j < nc->site_size; ++j) {
        int idx = nc->sites[j];
        sc->sites[sc->site_size++] = idx;
        if(idx != nb->r*n+nb->c) a[idx].cluster = sc; // don't overwrite neighbour until last
      }
      nc->id = -1; // mark as obsolete
      nb->cluster = sc; // now overwrite neighbour
    }
  }
}

/**
 * @brief scan each cluster and find num_clusters, max_cluster_size and row and col percolation.
 */
static void scan_clusters(CPArray* cpa, int n, int n_threads, int *num, int *max, short *rperc, short *cperc) {
  short rp = 0, cp = 0;
  int nm = 0, m = 0;
  for(int i = 0; i < n_threads; ++i) {
    for(int j = 0; j < cpa[i].size; ++j) {
      Cluster *cl = cpa[i].cls[j];
      if(cl->id == -1) continue; // this cluster was combined into another
      nm++;
      if(cl->size > m) m = cl->size;
      if(cl->height == n) rp = 1;
      if(cl->width == n) cp = 1;
    }
  }
  *rperc = rp;
  *cperc = cp;
  *num = nm;
  *max = m;
}

void print_params(short* a, Bond* b, int n, int n_threads, int num_workers, short site, char* fname, float p, int seed) {
  if(site) print_short_array(a, n);
  else print_bond(b, n);
  printf("\n%s\n%d CPU%s\n%d thread%s\n\nN: %d\n", site ? "Site" : "Bond", num_workers, num_workers > 1 ? "s" : "", n_threads, n_threads > 1 ? "s" : "", n);
  if(!fname) printf("P: %.2f\nS: %d\n", p, seed);
}

/**
 * USAGE: ./percolate [-s | -b] [-v] [-r SEED] [[-f LATTICE_FILENAME] | [N PROBABILITY]] [N_THREADS]
 * 
 * [-s | -b] site or bond percolation, default site
 * [-v] verbose
 * [-r SEED] number to seed srand, default time(NULL)
 * [-f LATTICE_FILENAME] file to scan lattice from
 * [N PROBABILITY] size of lattice and probability of site occupation or bond
 * [N_THREADS] number of threads to utilise
 */
int main(int argc, char *argv[])
{ 
  // all processes do command-line argument parsing
  double start = omp_get_wtime();

  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  short* a = NULL;
  Bond* b = NULL;

  // parse optional arguments
  short site = 1, verbose = 0;
  char* fname = NULL;
  unsigned int seed = time(NULL); // default unique seed

  int c;
  while ((c = getopt(argc, argv, "vbsf:r:")) != -1) {
    if(c == 'v') verbose = 1;              // silence printing
    else if(c == 'b') site = 0;            // bond
    else if(c == 'f') fname = optarg;      // scan lattice from file
    else if(c == 'r') seed = atoi(optarg); // seed rand with constant
  }
  // parse positional arguments
  int n;
  float p = -1.0;
  int n_threads = 1;

  if((fname && argc - optind < 1) || (!fname && argc - optind < 2)) {
    printf("Missing arguments.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); // exit(errno);
  }
  n = atoi(argv[optind++]);
  if(!fname) p = atof(argv[optind++]);
  if(argc - optind) n_threads = atoi(argv[optind]);

  if(n < 1 || n_threads < 1) {
    printf("Invalid arguments.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); // exit(errno);
  }
  // check n_threads
  int max_threads = omp_get_max_threads();
  if(n_threads > max_threads) n_threads = max_threads;
  if(n_threads > n) n_threads = n;
  omp_set_num_threads(n_threads);
  
  int max_clusters = n % 2 == 0 ? n*n/2 : (n-1)*(n-1)/2+1; // maximum number of size 1 clusters for a given n
  int num_workers = min(size, ceiling_divide(n, n_threads)); // utilise all threads first, then add nodes

  // master initialises lattice and sends to workers
  if(rank == MASTER) {
    srand(seed);
    if(site) {
      if(fname) a = file_short_array(fname, n);
      else a = short_array(n, p);
    }
    else {
      if(fname) b = file_bond(fname, n);
      else b = bond(n, p);
    }
    for(int r = 1; r < num_workers; ++r) {
      if(site) {
        MPI_Send(a, n*n, MPI_SHORT, r, TAG, MPI_COMM_WORLD);
        printf("Sent\n");
      }
      else {
        MPI_Send(b->v, n*n, MPI_SHORT, r, TAG, MPI_COMM_WORLD);
        MPI_Send(b->h, n*n, MPI_SHORT, r, TAG, MPI_COMM_WORLD);
      }
    }
    if(verbose) print_params(a, b, n, n_threads, num_workers, site, fname, p, seed);
  }
  // if(rank > MASTER && rank < num_workers) {
  //   printf("Rank %d\n", rank);
  //   if(site) {
  //     a = calloc(n*n, sizeof(short));
  //     MPI_Recv(a, n*n, MPI_SHORT, MASTER, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //   }
  //   else {
  //     b = calloc(1, sizeof(Bond));
  //     b->v = calloc(n*n, sizeof(short));
  //     b->h = calloc(n*n, sizeof(short));
  //     MPI_Recv(b->v, n*n, MPI_SHORT, MASTER, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //     MPI_Recv(b->h, n*n, MPI_SHORT, MASTER, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  //   }
  //   // if(verbose) print_params(a, b, n, n_threads, size, site, fname, p, seed);
  // }
  MPI_Finalize();
  exit(EXIT_SUCCESS);
}