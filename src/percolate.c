/*
 * CITS5507 HPC PROJECT 1
 * LATTICE PERCOLATION IN PARALLEL
 * 
 * Jasper Paterson 22736341
 * Allen Antony 22706998
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>
#include <errno.h>

#include "../include/stack.h"
#include "../include/site.h"
#include "../include/bond.h"
#include "../include/cluster.h"

/**
 * @return int start index of the region of the lattice allocated to this thread.
 * Allocates n/n_threads rows to the first n_threads-n%n_threads threads and n/n_threads+1 rows to the remaining threads.
 */
static int start_index(int n, int n_threads, int tid)
{
  return tid < n_threads - (n % n_threads) ? n*tid*(n/n_threads) : n*(n_threads-n%n_threads)*(n/n_threads) + n*(tid-(n_threads-n%n_threads))*(n/n_threads+1);
}

/**
 * @return int end index of the region of the lattice allocated to this thread.
 */
static int end_index(int n, int n_threads, int tid)
{
  return tid < n_threads - (n % n_threads) ? n*(tid+1)*(n/n_threads) : n*(n_threads-n%n_threads)*(n/n_threads) + n*((tid+1)-(n_threads-n%n_threads))*(n/n_threads+1);
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
    int start = start_index(n, n_threads, i); // region boundaries
    int end = end_index(n, n_threads, i);
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
  if(!s->cluster || !nb->cluster) {
    printf("Error if here.\n"); // both clusters should have been initialised
    return NULL;
  }
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
  int start = start_index(n, n_threads, tid);
  int end = end_index(n, n_threads, tid);
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
    int row_end = end_index(n, n_threads, i);
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

/**
 * USAGE: ./percolate [-s | -b] [-v] [-r SEED] [-p RESULTS_FILENAME] [[-f LATTICE_FILENAME] | [N PROBABILITY]] [N_THREADS]
 * 
 * [-s | -b] site or bond percolation, default site
 * [-v] silence printing
 * [-r SEED] number to seed srand, default time(NULL)
 * [-p RESULTS_FILENAME] file to append the results of the percolation
 * [-f LATTICE_FILENAME] file to scan lattice from
 * [N PROBABILITY] size of lattice and probability of site occupation or bond
 * [N_THREADS] number of threads to utilise
 * 
 */
int main(int argc, char *argv[])
{
  double start = omp_get_wtime();
  
  Site* a = NULL;
  Bond* b = NULL;
  
  // options
  short site = 1, verbose = 1;
  char* fname = NULL, *rname = NULL;
  unsigned int seed = time(NULL); // default unique seed

  // positional arguments
  int n;
  float p = -1.0;
  int n_threads = 1;
  
  int c;
  while ((c = getopt(argc, argv, "vbsf:p:r:")) != -1) {
    if(c == 'v') verbose = 0;   // silence printing
    else if(c == 'b') site = 0; // bond
    else if(c == 'f') {         // scan lattice from file
      if(!optarg) {
        printf("Error.\n");
        exit(errno);
      }
      fname = optarg;
    }
    else if(c == 'p') {         // write results to file
      if(!optarg) {
        printf("Error.\n");
        exit(errno);
      }
      rname = optarg;
    }
    else if(c == 'r') {         // seed rand with constant
      if(!optarg) {
        printf("Error.\n");
        exit(errno);
      }
      seed = atoi(optarg);
    }
  }
  if(fname) { // scan lattice from file
    if(argc - optind < 1) {
      printf("Invalid arguments.\n");
      exit(errno);
    }
    n = atoi(argv[optind++]);
    if(argc - optind > 0) n_threads = atoi(argv[optind]);
    if(site) {
      a = file_site_array(fname, n);
      if(!a) {
        printf("Error.\n");
        exit(errno);
      }
      print_site_array(a, n);
    } else {
      b = file_bond(fname, n);
      if(!b) {
        printf("Error.\n");
        exit(errno);
      }
      a = site_array(n, -1.0);
      if(!a) {
        printf("Memory error.\n");
        exit(errno);
      }
      print_bond(b, n);
    }
  } else { // initialise random lattice
    if(argc - optind < 2) {
      printf("Invalid arguments.\n");
      exit(errno);
    }
    n = atoi(argv[optind++]);
    p = atof(argv[optind++]);
    if(p > 1.0) p = 1.0;
    if(argc - optind > 0) n_threads = atoi(argv[optind]);
    
    srand(seed);
    if(site) {
      a = site_array(n, p);
      print_site_array(a, n);
    } else {
      b = bond(n, p);
      a = site_array(n, -1.0);
      print_bond(b, n);
    }
  }
  if(n < 1 || n_threads < 1) {
    printf("Invalid arguments.\n");
    exit(errno);
  }

  int max_clusters = n % 2 == 0 ? n*n/2 : (n-1)*(n-1)/2+1;

  int max_threads = omp_get_max_threads();
  if(n_threads > max_threads) n_threads = max_threads;
  if(n_threads > n) n_threads = n;
  omp_set_num_threads(n_threads);
  CPArray* cpa = cluster_array(n_threads, max_clusters); // each thread keeps an array of its cluster pointers 

  if(verbose) {
    printf("\n%s %d-Thread\n", site ? "Site" : "Bond", n_threads);
    printf("\nN: %d\n", n);
    if(p != -1.0) printf("P: %.2f\n", p);
    if(!fname) printf("S: %d\n", seed);
  }

  double init = omp_get_wtime();
  double init_time = init-start;
  if(verbose) printf("\n Init time: %9.6f\n", init_time);

  #pragma omp parallel
  {
    int num = omp_get_thread_num();
    percolate(a, b, n, n_threads, &cpa[num], num);
  }
  double pt = omp_get_wtime();
  double perc_time = pt-init;
  if(verbose) printf(" Perc time: %9.6f\n", perc_time);

  if(n_threads > 1) join_clusters(a, b, n, n_threads);
  double join = omp_get_wtime();
  double join_time = join-pt;
  if(verbose) printf(" Join time: %9.6f\n", join_time);

  int num = 0, max = 0;
  short rperc = 0, cperc = 0;
  scan_clusters(cpa, n, n_threads, &num, &max, &rperc, &cperc);
  double scan_time = omp_get_wtime()-join;
  if(verbose) printf(" Scan time: %9.6f\n", scan_time);

  double total_time = omp_get_wtime()-start;
  
  if(verbose) {
    printf("Total time: %9.6f\n", total_time);
    printf("\n   Num clusters: %d\n", num);
    printf("       Max size: %d\n", max);
    printf("Row percolation: %s\n", rperc ? "True" : "False");
    printf("Col percolation: %s\n\n", cperc ? "True" : "False");
  }

  if(rname && !fname) { // write results to results file
    FILE* f = fopen(rname, "a");
    if(!f) {
      printf("Error.\n");
      exit(errno);
    }
    fprintf(
      f,
      "%d,%f,%d,%d,%d,%d,%d,%d,%f,%f,%f,%f,%f\n",
      n, p, n_threads, seed, num, max, rperc, cperc, init_time, perc_time, join_time, scan_time, total_time
    );
    fclose(f);
  }
  exit(EXIT_SUCCESS);
}