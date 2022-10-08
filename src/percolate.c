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

static int get_start(int n, int n_rows, int id, int n_ids)
{
  return id < n_rows%n_ids ? n*id*(n_rows/n_ids+1) : n*(n_rows%n_ids)*(n_rows/n_ids+1) + n*(id-n_rows%n_ids)*(n_rows/n_ids);
}

static int get_n_rows(int n_rows, int id, int n_ids)
{
  return id < n_rows%n_ids ? n_rows/n_ids+1 : n_rows/n_ids;
}

static int max_clusters(int n, int n_rows)
{
  return n*n_rows; // overestimate for now
}

static short has_neighbours(Bond* b, int n, int idx)
{
  int r = idx/n, c = idx%n;
  int ib[] = {
    r*n+c,          // left
    r*n+(c+n+1)%n,  // right
    r*n+c,          // top
    ((r+n+1)%n)*n+c // bottom
  };
  for(int i = 0; i < 4; ++i) {
    if((i<2 && b->h[ib[i]]) || (i>=2 && b->v[ib[i]])) return 1; 
  }
  return 0;
}

static void get_neighbours(Site* sites, Bond* b, int nbs[], int idx, int n, int t_start, int t_end, int nt_rows)
{
  int r = idx/n, c = idx%n;
  // indices of neighbours
  int is[] = {
    r*n+(c+n-1)%n,   // left
    r*n+(c+n+1)%n,   // right
    ((r+n-1)%n)*n+c, // top
    ((r+n+1)%n)*n+c  // bottom
  };
  // indices of bonds
  int ib[] = {
    r*n+c,          // left
    r*n+(c+n+1)%n,  // right
    r*n+c,          // top
    ((r+n+1)%n)*n+c // bottom
  };
  for(int i = 0; i < 4; ++i) {
    int nbi = is[i], bi = ib[i];
    Site* nb = &sites[nbi];
    if(nbi < t_start || nbi >= t_end) nbs[i] = -1; // out of thread bounds
    else if(!nb->seen && (
      (!b && nb->occupied) ||
      (b && ((i<2 && b->h[bi]) || (i>=2 && b->v[bi]))))
    ) {
      nb->seen = 1; // must mark site here for n=2 case when top neighbour == bottom neighbour
      nbs[i] = nbi;
    } else nbs[i] = -1;
  }
}

static short on_border(int idx, int n, int nt_workers, int np_rows) {
  for(int tid = 0; tid < nt_workers; ++tid) {
    int start = get_start(n, np_rows, tid, nt_workers); // starting from zero
    int n_rows = get_n_rows(np_rows, tid, nt_workers);
    int end = start + n*n_rows;
    if((idx >= start && idx < start+n) || (idx >= end-n && idx < end)) return 1;
  }
  return 0;
}

static void DFS(Site* sites, Bond* b, Stack* st, int n, int t_start, int t_end, int nt_rows, int nt_workers, int np_rows) {
  while(!is_empty(st)) {
    int idx = pop(st);
    int nbs[4];
    get_neighbours(sites, b, nbs, idx, n, t_start, t_end, nt_rows);
    Cluster* cl = sites[idx].cluster;
    for(int i = 0; i < 4; ++i) { // loop through connected, unseen neighbours
      if(nbs[i] == -1) continue; 
      Site* nb = &sites[nbs[i]];
      nb->cluster = cl;
      ++cl->size;
      if(!cl->rows[nb->r]) cl->height++;
      if(!cl->cols[nb->c]) cl->width++;
      cl->rows[nb->r] = 1;
      cl->cols[nb->c] = 1;
      if(on_border(nbs[i], n, nt_workers, np_rows)) cl->sites[cl->site_size++] = nbs[i]; // add index of neighbour to cluster's site index array
      add(st, nbs[i]);
    }
  }
}

static void percolate(Site* sites, Bond* b, int n, int tid, int t_start, int nt_rows, int nt_workers, int np_rows, Cluster** clusters, int* n_clusters)
{
  int n_sites = n*nt_rows;
  int t_end = t_start + n_sites;
  Stack* st = stack(n_sites);
  for(int i = t_start; i < t_end; ++i) {
    Site *s = &sites[i];
    if(!s->seen && ((!b && s->occupied) || (b && has_neighbours(b, n, i)))) { // if unseen and will form a cluster
      s->seen = 1;
      s->cluster = cluster(n, nt_workers, s->r, s->c);
      Cluster *sc = s->cluster;
      if(on_border(i, n, nt_workers, np_rows)) sc->sites[sc->site_size++] = i;
      clusters[(*n_clusters)++] = sc; 
      add(st, i);
      DFS(sites, b, st, n, t_start, t_end, nt_rows, nt_workers, np_rows);
    }   
  }
}

void print_params(short* a, Bond* b, int n, int n_threads, int n_workers, short site, char* fname, float p, int seed) {
  if(site) print_short_array(a, n);
  else print_bond(b, n);
  char str[50];
  sprintf(str, "P: %.2f\nS: %d\n", p, seed);
  printf("\n%s\n%d CPU%s\n%d thread%s\n\nN: %d\n%s\n", site ? "Site" : "Bond", n_workers, n_workers > 1 ? "s" : "", n_threads, n_threads > 1 ? "s" : "", n, fname ? "" : str);
}

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
  unsigned int seed = time(NULL); // only used by master, could be inconsistent between processes

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
  int n_workers = min(size, ceiling_divide(n, n_threads)); // utilise all threads first, then add nodes

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
    for(int r = 1; r < n_workers; ++r) {
      if(site) {
        MPI_Send(a, n*n, MPI_SHORT, r, TAG, MPI_COMM_WORLD);
      }
      else {
        MPI_Send(b->v, n*n, MPI_SHORT, r, TAG, MPI_COMM_WORLD);
        MPI_Send(b->h, n*n, MPI_SHORT, r, TAG, MPI_COMM_WORLD);
      }
    }
    if(verbose) print_params(a, b, n, n_threads, n_workers, site, fname, p, seed);
  }
  if(rank > MASTER && rank < n_workers) {
    if(site) {
      a = calloc(n*n, sizeof(short));
      MPI_Recv(a, n*n, MPI_SHORT, MASTER, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else {
      b = calloc(1, sizeof(Bond));
      b->v = calloc(n*n, sizeof(short));
      b->h = calloc(n*n, sizeof(short));
      MPI_Recv(b->v, n*n, MPI_SHORT, MASTER, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(b->h, n*n, MPI_SHORT, MASTER, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }
  if(rank < n_workers) {
    int p_start = get_start(n, n, rank, n_workers);
    int np_rows = get_n_rows(n, rank, n_workers);
    int p_end = p_start + n*np_rows;
    int nt_workers = min(n_threads, np_rows);

    Site* sites = site_array(a, n);
    if(site) free(a); // occupation info now stored in sites

    Cluster*** t_clusters = calloc(nt_workers, sizeof(Cluster**));
    for(int i = 0; i < nt_workers; ++i) t_clusters[i] = calloc(max_clusters(n, np_rows), sizeof(Cluster*));
    int* nt_clusters = calloc(nt_workers, sizeof(int));

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      int t_start = p_start + get_start(n, np_rows, tid, nt_workers); // starting from zero
      int nt_rows = get_n_rows(np_rows, tid, nt_workers);
      if(tid < nt_workers) {
        percolate(sites, b, n, tid, t_start, nt_rows, nt_workers, np_rows, t_clusters[tid], &nt_clusters[tid]);
        // printf("Rank %d thread %d start %d end %d num rows %d num clusters %d\n", rank, tid, t_start, t_start+n*nt_rows, nt_rows, nt_clusters[tid]);
      }
      for(int i = 0; i < nt_clusters[tid]; ++i) printf("Rank %d thread %d cluster %d size %d\n", rank, tid, t_clusters[tid][i]->id, t_clusters[tid][i]->size);
    }
    // join clusters
    // send data to master
  }
  // master join clusters
  // print outputs
  MPI_Finalize();
  exit(EXIT_SUCCESS);
}