/*
 * CITS5507 HPC PROJECT 2
 * LATTICE PERCOLATION USING MPI AND OPENMP
 * 
 * Jasper Paterson 22736341
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <util.h>
#include <string.h>

#include <omp.h>
#include <mpi.h>

#include "../include/constant.h"
#include "../include/stack.h"
#include "../include/site.h"
#include "../include/bond.h"
#include "../include/cluster.h"

/**
 * @return int start index of region for this id
 */
static int get_start(int n, int n_rows, int id, int n_ids)
{
  return id < n_rows%n_ids ? n*id*(n_rows/n_ids+1) : n*(n_rows%n_ids)*(n_rows/n_ids+1) + n*(id-n_rows%n_ids)*(n_rows/n_ids);
}

/**
 * @return int number of rows allocated to this id
 */
static int get_n_rows(int n_rows, int id, int n_ids)
{
  return id < n_rows%n_ids ? n_rows/n_ids+1 : n_rows/n_ids;
}

/**
 * @return int max number of clusters that can exist in region
 */
static int max_clusters(int n, int n_rows)
{
  return n_rows*(n/2);
}

/**
 * @return short 1 iff site has a bond to a neighbour
 */
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

/**
 * @brief fill nbs int array with valid neighbour index or -1
 */
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

/**
 * @brief traverse and update current cluster 
 */
static void DFS(Site* sites, Bond* b, Stack* st, int n, int t_start, int t_end, int nt_rows, int nt_workers, int p_start, int np_rows) {
  while(!is_empty(st)) {
    int idx = pop(st);
    int nbs[4];
    get_neighbours(sites, b, nbs, idx, n, t_start, t_end, nt_rows);
    int* cl = sites[idx].cluster;
    for(int i = 0; i < 4; ++i) { // loop through connected, unseen neighbours
      if(nbs[i] == -1) continue; 
      Site* nb = &sites[nbs[i]];
      nb->cluster = cl;
      cl[1]++;
      cl[2+nb->r] = 1;
      cl[2+n+nb->c] = 1;
      add(st, nbs[i]);
    }
  }
}

/**
 * @return short 1 iff cluster lies on thread border
*/
short on_border(Site *sites, int* sc, int n, int t_start, int t_end)
{
  for(int j = t_start; j < t_start+n; ++j) { // top row
    if(sites[j].cluster == sc) return 1;
  }
  for(int j = t_end-n; j < t_end; ++j) { // bottom row
    if(sites[j].cluster == sc) return 1;
  }
  return 0;
}

/**
 * @brief update num clusters, max size, rperc and cperc
*/
void update_stats(int* sc, int n, int* num, int* max, short* rperc, short* cperc)
{
  if(num) ++(*num);
  if(max && sc[1] > *max) *max = sc[1];
  int rs = 0, cs = 0;
  for(int j = 0; j < n; ++j) {
    if(rperc && sc[2+j]) ++rs;
    if(cperc && sc[2+n+j]) ++cs;
  }
  if(rperc && rs == n) *rperc = 1;
  if(cperc && cs == n) *cperc = 1;
}

/**
 * @brief loop through sites in region and search clusters
 */
static void percolate(Site* sites, Bond* b, int n, int tid, int t_start, int nt_rows, int nt_workers, int p_start, int np_rows, int** clusters, int* n_clusters, int* num, int* max, short* cperc)
{
  int n_sites = n*nt_rows;
  int t_end = t_start + n_sites;
  Stack* st = stack(n_sites);
  for(int i = t_start; i < t_end; ++i) {
    Site *s = &sites[i];
    if(!s->seen && ((!b && s->occupied) || (b && has_neighbours(b, n, i)))) { // if unseen and will form a cluster
      s->seen = 1;
      s->cluster = cluster(n, i);
      int *sc = s->cluster;
      add(st, i);
      DFS(sites, b, st, n, t_start, t_end, nt_rows, nt_workers, p_start, np_rows);
      update_stats(sc, n, num, max, NULL, cperc);
      if(on_border(sites, sc, n, t_start, t_end)) clusters[(*n_clusters)++] = sc; 
      else free(sc);
    }   
  }
  free_stack(st);
}

/**
 * @return int bottom neighbour index if connected else -1
 */
static int bottom_neighbour(Site* sites, Bond* b, int n, int i, int p_start, int np_rows)
{
  int r = i/n, c = i%n;
  int nbi = ((r+n+1)%n)*n+c; // index of bottom neighbour in a and bottom bond in b
  if(nbi < p_start || nbi >= p_start + n*np_rows) return -1; // in another process's region
  Site *s = &sites[i], *nb = &sites[nbi];
  if((!b && !nb->occupied) || (b && !b->v[nbi])) return -1;
  if(s->cluster[0] == nb->cluster[0]) return -1;
  return nbi;
}

/**
 * @brief if site j belonged to the old cluster, point it to the new cluster
 */
void update_cluster(Site *sites, int j, int nbi, int *nc, int *sc) {
  int *c = sites[j].cluster;
  if(j != nbi && c && c[0] == nc[0]) sites[j].cluster = sc;
}

/**
 * @brief loop over borders and connect clusters 
 */
static void join_clusters(Site* sites, Bond* b, int n, int n_workers, int start, int rows, int* num) {
  for(int w = 0; w < n_workers; ++w) {
    int s_start = start + get_start(n, rows, w, n_workers);
    int s_end = s_start + n*get_n_rows(rows, w, n_workers);
    for(int i = s_end-n; i < s_end; ++i) { // loop along bottom row of region
      Site *s = &sites[i];
      int *sc = s->cluster;
      if(!sc) continue;
      int nbi = bottom_neighbour(sites, b, n, i, start, rows);
      if(nbi == -1) continue;
      Site *nb = &sites[nbi];
      int *nc = nb->cluster;
      
      // combine two clusters into sc
      sc[1] += nc[1];
      for(int j = 0; j < n; ++j) {
        if(nc[2+j]) sc[2+j] = 1;
        if(nc[2+n+j]) sc[2+n+j] = 1;
      }
      // start and end of below region
      int s2_start = start + get_start(n, rows, (w+1)%n_workers, n_workers);
      int s2_end = s2_start + n*get_n_rows(rows, (w+1)%n_workers, n_workers);

      // loop along current boundary and joined cluster borders, update clusters
      for(int j = start; j < start+n; ++j) update_cluster(sites, j, nbi, nc, sc);
      for(int j = s_end-n; j < s_end; ++j) update_cluster(sites, j, nbi, nc, sc);
      for(int j = s2_start; j < s2_start+n; ++j) update_cluster(sites, j, nbi, nc, sc);
      for(int j = s2_end-n; j < s2_end; ++j) update_cluster(sites, j, nbi, nc, sc);

      nc[0] = -1; // mark as obsolete
      nb->cluster = sc; // now overwrite neighbour

      if(num) --(*num); // two clusters become one
    }
  }
}

/**
 * @brief pack data with the cluster id of each border site
 */
void copy_site_data(Site* sites, int i, int* data, int* di)
{
  int *c = sites[i].cluster;
  if(c) data[(*di)++] = c[0];
  else data[(*di)++] = -1;
}

/**
 * @brief pack unique cluster data into int array
 */
void copy_cluster_data(Site* sites, int n, int i, int* seen_cluster_ids, int* seen_index, int* data, int* di)
{
  int c_size = 2+2*n;
  int *c = sites[i].cluster;
  if(c && c[0] != -1 && !in_array(c[0], seen_cluster_ids, *seen_index)) {
    seen_cluster_ids[(*seen_index)++] = c[0];
    data[3]++;
    memcpy(&data[*di], c, c_size * sizeof(int));
    *di += c_size;
  }
}

/**
 * @brief Send all relevant information from worker to master
 */
void send_clusters(int rank, Site* sites, int n, int nt_workers, int** t_clusters, int* nt_clusters, int mc, int p_start, int p_end, int num, int max)
{
  int c_size = 2+2*n;
  int *data = calloc(4 + 2*n + c_size*(n+2), sizeof(int));
  int di = 4;
  
  data[0] = num;
  data[1] = max;
  for(int tid = 0; tid < nt_workers; ++tid) {
    for(int i = 0; i < nt_clusters[tid]; ++i) {
      int *c = (t_clusters + mc*tid)[i];
      update_stats(c, n, NULL, &data[1], (short*)&data[2], NULL);
    }
  }

  // copy site data along borders
  for(int i = p_start; i < p_start+n; ++i) copy_site_data(sites, i, data, &di); // top row
  for(int i = p_end-n; i < p_end; ++i) copy_site_data(sites, i, data, &di); // bottom row

  // seen clusters
  int seen_cluster_ids[n+2];
  int seen_index = 0;

  // copy cluster data along borders
  for(int i = p_start; i < p_start+n; ++i) copy_cluster_data(sites, n, i, seen_cluster_ids, &seen_index, data, &di); // top row
  for(int i = p_end-n; i < p_end; ++i) copy_cluster_data(sites, n, i, seen_cluster_ids, &seen_index, data, &di); // bottom row

  MPI_Send(data, 4+2*n, MPI_INT, MASTER, TAG, MPI_COMM_WORLD);
  MPI_Send(data + 4+2*n, c_size * data[3], MPI_INT, MASTER, TAG, MPI_COMM_WORLD);
  
  free(data);
}

/**
 * @brief receive cluster data from all workers and store in cluster array
 */
void recv_clusters(Site* sites, int n, int n_workers, int** p_clusters, int* np_clusters, int* num, int* max, short *cperc)
{
  int c_size = 2+2*n;
  int d_size = 0; // ints required to store all worker's clusters

  int* data = malloc((n_workers-1) * (4+2*n) * sizeof(int));
  
  for(int i = 0; i < n_workers-1; ++i) {
    int *d = data + i * (4+2*n);
    MPI_Recv(d, 4+2*n, MPI_INT, i+1, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    d_size += c_size * d[3];
  }

  *p_clusters = malloc(d_size * sizeof(int));

  for(int i = 0; i < n_workers-1; ++i) {
    int *d = data + i * (4+2*n);
    MPI_Recv(*p_clusters + c_size*(*np_clusters), c_size * d[3], MPI_INT, i+1, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    *np_clusters += d[3];

    // num clusters, max cluster size and cperc
    *num += d[0];
    if(d[1] > *max) *max = d[1];
    if(d[2]) *cperc = 1;

    // add site cluster pointers to sites
    int p_start = get_start(n, n, i+1, n_workers);
    int p_end = p_start + n*get_n_rows(n, i+1, n_workers);

    for(int j = p_start; j < p_start+n; ++j) { // top row
      for(int k = 0; k < *np_clusters; ++k) {
        if((*p_clusters + c_size*k)[0] == d[4+j-p_start]) {
          sites[j].cluster = *p_clusters + c_size*k;
          break;
        }
      }
    }

    for(int j = p_end-n; j < p_end; ++j) { // bottom row
      for(int k = 0; k < *np_clusters; ++k) {
        if((*p_clusters + c_size*k)[0] == d[4+j-p_end+2*n]) {
          sites[j].cluster = *p_clusters + c_size*k;
          break;
        }
      }
    }
  }
  free(data);
}

/**
 * @brief print the lattice and percolation parameters
 */
void print_params(short* a, Bond* b, int n, int n_threads, int n_workers, short site, char* fname, float p, int seed) {
  if(site) print_short_array(a, n);
  else print_bond(b, n);
  char str[50];
  sprintf(str, "P: %.2f\nS: %d\n", p, seed);
  printf("\n%s\n%d CPU%s\n%d thread%s\n\nN: %d\n%s\n", site ? "Site" : "Bond", n_workers, n_workers > 1 ? "s" : "", n_threads, n_threads > 1 ? "s" : "", n, fname ? "" : str);
  fflush(stdout);
}

int main(int argc, char *argv[])
{ 
  // all processes do command-line argument parsing
  double start_init = MPI_Wtime();

  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  short* a = NULL;
  Bond* b = NULL;

  // parse optional arguments
  short site = 1, verbose = 0;
  char *fname = NULL;
  unsigned int seed = time(NULL); // only used by master, could be inconsistent between processes

  int c;
  while ((c = getopt(argc, argv, "vbsf:r:o:")) != -1) {
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
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  n = atoi(argv[optind++]);
  if(!fname) p = atof(argv[optind++]);
  if(argc - optind) n_threads = atoi(argv[optind]);

  if(n < 1 || n_threads < 1) {
    printf("Invalid arguments.\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
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
      short* s = calloc(2*n*n, sizeof(short));
      b->v = s;
      b->h = s+n*n;
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

    // initialise cluster array for each thread
    int mc = max_clusters(n, get_n_rows(np_rows, 0, nt_workers));
    int** t_clusters = calloc(nt_workers * mc, sizeof(int*));
    int* nt_clusters = calloc(nt_workers, sizeof(int));

    // local stats
    int nums[nt_workers], maxs[nt_workers];
    for(int i = 0; i < nt_workers; ++i) { nums[i] = 0; maxs[i] = 0; }

    // process stats
    int num = 0, max = 0;
    short rperc = 0, cperc = 0;

    double start_perc = MPI_Wtime();

    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      int t_start = p_start + get_start(n, np_rows, tid, nt_workers);
      int nt_rows = get_n_rows(np_rows, tid, nt_workers);
      if(tid < nt_workers) percolate(sites, b, n, tid, t_start, nt_rows, nt_workers, p_start, np_rows, &t_clusters[mc*tid], &nt_clusters[tid], &nums[tid], &maxs[tid], &cperc);
    }

    // update process stats
    for(int i = 0; i < nt_workers; ++i) {
      num += nums[i];
      if(maxs[i] > max) max = maxs[i];
    }

    double start_tjoin = MPI_Wtime();
    if(nt_workers > 1) join_clusters(sites, b, n, nt_workers, p_start, np_rows, &num);

    if(rank > MASTER) send_clusters(rank, sites, n, nt_workers, t_clusters, nt_clusters, mc, p_start, p_end, num, max);

    else if(rank == MASTER) {
      int c_size = 2+2*n;

      double start_recv = MPI_Wtime();
      double start_pjoin = MPI_Wtime();
      double end_pjoin = MPI_Wtime();

      if(n_workers > 1) {

        // receive cluster data
        int* p_clusters = NULL; // malloc(c_size * (n_workers-1)*(n+2) * sizeof(int));
        int np_clusters = 0;
        recv_clusters(sites, n, n_workers, &p_clusters, &np_clusters, &num, &max, &cperc);

        // join process regions
        start_pjoin = MPI_Wtime();
        join_clusters(sites, b, n, n_workers, 0, n, &num);
        end_pjoin = MPI_Wtime();

        // scan workers' clusters
        for(int i = 0; i < np_clusters; ++i) {
          int *sc = &p_clusters[c_size*i];
          if(sc[0] != -1) update_stats(sc, n, NULL, &max, &rperc, &cperc);
        }
      }

      // scan master's clusters
      for(int tid = 0; tid < nt_workers; ++tid) {
        for(int i = 0; i < nt_clusters[tid]; ++i) {
          int *sc = (t_clusters + mc*tid)[i];
          if(sc[0] != -1) update_stats(sc, n, NULL, &max, &rperc, &cperc);
        }
      }

      double end = MPI_Wtime();

      if(verbose) {
        printf(" Init time %9.6f\n", start_perc-start_init);
        printf(" Perc time %9.6f\n", start_tjoin-start_perc);
        printf("Tjoin time %9.6f\n", start_recv-start_tjoin);
        printf(" Recv time %9.6f\n", start_pjoin-start_recv);
        printf("Pjoin time %9.6f\n", end_pjoin-start_pjoin);
        printf("Total time %9.6f\n", end-start_init);
        printf("\n");
        printf("   Num clusters: %d\n", num);
        printf("       Max size: %d\n", max);
        printf("Row percolation: %s\n", rperc ? "True" : "False");
        printf("Col percolation: %s\n", cperc ? "True" : "False");
      }
      else {
        printf(
          "%d,%f,%d,%d,%d,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f\n",
          n, p, n_workers, n_threads, seed, num, max, rperc, cperc, start_perc-start_init, start_tjoin-start_perc, start_recv-start_tjoin, start_pjoin-start_recv, end_pjoin-start_pjoin, end-start_init
        );
      }
    }
  }
  fflush(stdout);
  MPI_Finalize();
  exit(EXIT_SUCCESS);
}