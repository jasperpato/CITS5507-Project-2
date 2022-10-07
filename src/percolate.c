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

static int get_start_index(int n, int num_rows, int id, int num_ids)
{
  return id < num_rows%num_ids ? n*id*(num_rows/num_ids+1) : n*(num_rows%num_ids)*(num_rows/num_ids+1) + n*(id-num_rows%num_ids)*(num_rows/num_ids);
}

static int get_num_rows(int num_rows, int id, int num_ids)
{
  return id < num_rows%num_ids ? num_rows/num_ids+1 : num_rows/num_ids;
}

static int max_clusters(int n, int num_rows)
{
  return n*num_rows; // overestimate for now
}

// static void percolate(Site* a, Bond* b, int n, int count, int n_threads, Cluster* clusters, int* num_clusters, short tid)
// {
//   int start = start_index(count, tid, n_threads);
//   int end = end_index(count, tid, n_threads);
//   int sub_count = get_count(count, tid, n_threads);
//   Stack* st = stack(count);
//   for(int i = start; i < end; ++i) {
//     Site *s = &a[i];
//     if(!s->seen && ((!b && s->occupied) || (b && has_neighbours(b, n, s)))) { // if unseen and will form a cluster
//       s->seen = 1;
//       s->cluster = cluster(n, n_threads, s->r, s->c);
//       Cluster *sc = s->cluster;
//       int idx = s->r*n+s->c;
//       if(on_border(n, idx, n_threads)) sc->sites[sc->site_size++] = idx;
//       cpa->cls[cpa->size++] = sc; 
//       add(st, s);
//       DFS(a, b, n, n_threads, st, start, end);
//     }   
//   }
// }

void print_params(short* a, Bond* b, int n, int n_threads, int num_workers, short site, char* fname, float p, int seed) {
  if(site) print_short_array(a, n);
  else print_bond(b, n);
  char str[50];
  sprintf(str, "P: %.2f\nS: %d\n", p, seed);
  printf("\n%s\n%d CPU%s\n%d thread%s\n\nN: %d\n%s\n", site ? "Site" : "Bond", num_workers, num_workers > 1 ? "s" : "", n_threads, n_threads > 1 ? "s" : "", n, fname ? "" : str);
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
  
  // int max_clusters = n % 2 == 0 ? n*n/2 : (n-1)*(n-1)/2+1; // no longer relevant
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
      }
      else {
        MPI_Send(b->v, n*n, MPI_SHORT, r, TAG, MPI_COMM_WORLD);
        MPI_Send(b->h, n*n, MPI_SHORT, r, TAG, MPI_COMM_WORLD);
      }
    }
    if(verbose) print_params(a, b, n, n_threads, num_workers, site, fname, p, seed);
  }
  if(rank > MASTER && rank < num_workers) {
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
  if(rank < num_workers) {
    int process_start_index = get_start_index(n, n, rank, size);
    int process_num_rows = get_num_rows(n, rank, size);

    int num_thread_workers = min(n_threads, process_num_rows);

    Site* sites = site_array(a, n, process_start_index, process_num_rows);
    if(site) free(a); // occupation info now stored in sites

    Cluster*** thread_clusters = calloc(num_thread_workers, sizeof(Cluster**));
    for(int i = 0; i < num_thread_workers; ++i) thread_clusters[i] = calloc(max_clusters(n, process_num_rows), sizeof(Cluster*));
    int* num_thread_clusters = calloc(num_thread_workers, sizeof(int));

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      int thread_start_index = process_start_index + get_start_index(n, process_num_rows, tid, num_thread_workers);
      int thread_num_rows = get_num_rows(process_num_rows, tid, num_thread_workers);
      if(tid < num_thread_workers) {
        printf("Rank %d thread %d start index %d num rows %d\n", rank, tid, thread_start_index, thread_num_rows);
        // percolate(a, b, n, count, n_threads, &thread_clusters[tid], &num_thread_clusters[tid], tid);
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