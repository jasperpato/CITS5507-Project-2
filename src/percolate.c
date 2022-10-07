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

static int get_start_index(int n, int n_rows, int id, int n_ids)
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

static void get_neighbours(Site* sites, Bond* b, int n, Site* s, Site* nbs[], int start_index, int end_index)
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
    Site* nb = &sites[is[i]];
    if(is[i] < start_index || is[i] >= end_index) nbs[i] = NULL; // out of thread bounds
    else if(!nb->seen && (
      (!b && nb->occupied) ||
      (b && ((i<2 && b->h[ib[i]]) || (i>=2 && b->v[ib[i]]))))
    ) {
      nb->seen = 1; // must mark site here for n=2 case when top neighbour == bottom neighbour
      nbs[i] = nb;
    } else nbs[i] = NULL;
  }
}

static short on_border(int idx, int n, int process_n_rows, int n_thread_workers) {
  for(int tid = 0; tid < n_thread_workers; ++tid) {
    int start_index = get_start_index(n, process_n_rows, tid, n_thread_workers); // starting from zero
    int n_rows = get_n_rows(process_n_rows, tid, n_thread_workers);
    int end_index = start_index + n*n_rows;
    if((idx >= start_index && idx < start_index+n) || (idx >= end_index-n && idx < end_index)) return 1;
  }
  return 0;
}

static void DFS(Site* sites, Bond* b, int n, int n_thread_workers, int process_n_rows, Stack* st, int start_index, int end_index) {
  while(!is_empty(st)) {
    Site *s = pop(st);
    Site *nbs[4];
    get_neighbours(sites, b, n, s, nbs, start_index, end_index);
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
      if(on_border(idx, n, process_n_rows, n_thread_workers)) cl->sites[cl->site_size++] = idx; // add index of neighbour to cluster's site index array
      add(st, nb);
    }
  }
}

static void percolate(Site* sites, Bond* b, int n, int start_index, int tid, int n_thread_workers, int thread_n_rows, int process_n_rows, Cluster** clusters, int* n_clusters)
{
  int n_sites = n*thread_n_rows;
  int end_index = start_index + n_sites;
  Stack* st = stack(n_sites);
  for(int i = start_index; i < end_index; ++i) {
    Site *s = &sites[i];
    if(!s->seen && ((!b && s->occupied) || (b && has_neighbours(b, n, s)))) { // if unseen and will form a cluster
      s->seen = 1;
      s->cluster = cluster(n, n_thread_workers, s->r, s->c);
      Cluster *sc = s->cluster;
      int idx = s->r*n+s->c;
      if(on_border(idx, n, process_n_rows, n_thread_workers)) sc->sites[sc->site_size++] = idx;
      clusters[(*n_clusters)++] = sc; 
      add(st, s);
      DFS(sites, b, n, n_thread_workers, process_n_rows, st, start_index, end_index);
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
    int process_start_index = get_start_index(n, n, rank, n_workers);
    int process_n_rows = get_n_rows(n, rank, n_workers);
    int process_end_index = process_start_index + n*process_n_rows;
    int n_thread_workers = min(n_threads, process_n_rows);

    Site* sites = site_array(a, n, process_start_index, process_end_index);
    if(site) free(a); // occupation info now stored in sites

    Cluster*** thread_clusters = calloc(n_thread_workers, sizeof(Cluster**));
    for(int i = 0; i < n_thread_workers; ++i) thread_clusters[i] = calloc(max_clusters(n, process_n_rows), sizeof(Cluster*));
    int* n_thread_clusters = calloc(n_thread_workers, sizeof(int));

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      int thread_start_index = get_start_index(n, process_n_rows, tid, n_thread_workers); // starting from zero
      int thread_n_rows = get_n_rows(process_n_rows, tid, n_thread_workers);
      if(tid < n_thread_workers) {
        percolate(sites, b, n, thread_start_index, tid, n_thread_workers, thread_n_rows, process_n_rows, thread_clusters[tid], &n_thread_clusters[tid]);
        printf("Rank %d thread %d num clusters %d\n", rank, tid, n_thread_clusters[tid]);
        for(int i = 0; i < n_thread_clusters[tid]; ++i) {
          printf("Rank %d thread %d cluster %d size %d\n", rank, tid, thread_clusters[tid][i]->id, thread_clusters[tid][i]->size);
        } 
      }
    }
    // join clusters
    // send data to master
  }
  // master join clusters
  // print outputs
  MPI_Finalize();
  exit(EXIT_SUCCESS);
}