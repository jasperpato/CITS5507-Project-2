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
  return n_rows*(n/2+1);
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

static void DFS(Site* sites, Bond* b, Stack* st, int n, int t_start, int t_end, int nt_rows, int nt_workers, int p_start, int np_rows) {
  while(!is_empty(st)) {
    int idx = pop(st);
    int nbs[4];
    get_neighbours(sites, b, nbs, idx, n, t_start, t_end, nt_rows);
    Cluster* cl = sites[idx].cluster;
    for(int i = 0; i < 4; ++i) { // loop through connected, unseen neighbours
      if(nbs[i] == -1) continue; 
      Site* nb = &sites[nbs[i]];
      nb->cluster = cl;
      cl->size++;
      if(!cl->rows[nb->r]) cl->height++;
      if(!cl->cols[nb->c]) cl->width++;
      cl->rows[nb->r] = 1;
      cl->cols[nb->c] = 1;
      add(st, nbs[i]);
    }
  }
}

static void percolate(Site* sites, Bond* b, int n, int tid, int t_start, int nt_rows, int nt_workers, int p_start, int np_rows, Cluster** clusters, int* n_clusters)
{
  int n_sites = n*nt_rows;
  int t_end = t_start + n_sites;
  Stack* st = stack(n_sites);
  for(int i = t_start; i < t_end; ++i) {
    Site *s = &sites[i];
    if(!s->seen && ((!b && s->occupied) || (b && has_neighbours(b, n, i)))) { // if unseen and will form a cluster
      s->seen = 1;
      s->cluster = cluster(n, nt_workers, i);
      Cluster *sc = s->cluster;
      clusters[(*n_clusters)++] = sc; 
      add(st, i);
      DFS(sites, b, st, n, t_start, t_end, nt_rows, nt_workers, p_start, np_rows);
    }   
  }
}

static int bottom_neighbour(Site* sites, Bond* b, int n, int i, int p_start, int np_rows)
{
  int r = i/n, c = i%n;
  int nbi = ((r+n+1)%n)*n+c; // index of bottom neighbour in a and bottom bond in b
  if(nbi < p_start || nbi >= p_start + n*np_rows) return -1; // in another process's region
  Site *s = &sites[i], *nb = &sites[nbi];
  if((!b && !nb->occupied) || (b && !b->v[nbi])) return -1;
  if(!s->cluster || !nb->cluster) {
    printf("Error if here s %p n %p.\n", s->cluster, nb->cluster); // both clusters should have been initialised
    return -1;
  }
  if(s->cluster->id == nb->cluster->id) return -1;
  return nbi;
}

static void join_clusters(Site* sites, Bond* b, int n, int n_workers, int start, int rows) {
  for(int w = 0; w < n_workers; ++w) {
    int s_start = start + get_start(n, rows, w, n_workers);
    int s_end = s_start + n*get_n_rows(rows, w, n_workers);
    for(int i = s_end-n; i < s_end; ++i) { // loop along bottom row of region
      Site *s = &sites[i];
      Cluster *sc = s->cluster;
      if(!sc) continue;
      int nbi = bottom_neighbour(sites, b, n, i, start, rows);
      if(nbi == -1) continue;
      Site *nb = &sites[nbi];
      Cluster *nc = nb->cluster;
      
      // combine two clusters into sc
      sc->size += nc->size;
      for(int j = 0; j < n; ++j) {
        if(nc->rows[j]) {
          if(!sc->rows[j]) sc->height++;
          sc->rows[j] = 1;
        }
        if(nc->cols[j]) {
          if(!sc->cols[j]) sc->width++;
          sc->cols[j] = 1;
        }
      }
      // loop along two borders of joined cluster and update cluster pointers if necessary
      for(int j = start; j < start+n; ++j) {
        Cluster *c = sites[j].cluster;
        if(j != nbi && c && c->id == nc->id) sites[j].cluster = sc;
      }
      int s2_end = start + get_start(n, rows, (w+1)%n_workers, n_workers) + n*get_n_rows(rows, (w+1)%n_workers, n_workers);
      for(int j = s2_end-n; j < s2_end; ++j) {
        Cluster *c = sites[j].cluster;
        if(j != nbi && c && c->id == nc->id) sites[j].cluster = sc;
      }
      nc->id = -1; // mark as obsolete
      nb->cluster = sc; // now overwrite neighbour
    }
  }
}

void send_clusters(int rank, Site* sites, int n, int nt_workers, Cluster*** t_clusters, int* nt_clusters, int p_start, int p_end)
{
  int p_stats[4] = {0,0,0,0}; // num clusters, max cluster size, col perc boolean, num border clusters
  for(int tid = 0; tid < nt_workers; ++tid) {
    for(int i = 0; i < nt_clusters[tid]; ++i) {
      Cluster *c = t_clusters[tid][i];
      if(c->id == -1) continue;
      p_stats[0]++;
      if(c->size > p_stats[1]) p_stats[1] = c->size;
      if(c->width == n) p_stats[2] = 1;
    }
  }
  int nc_attrs = 4 + 2*n; // number of ints that describes a cluster

  // group site data and then cluster data into one array
  int *data = calloc(2*n + nc_attrs*(n+2), sizeof(int));
  int di = 0;

  // seen clusters
  int *seen_cluster_ids = calloc(n+2, sizeof(int));
  int seen_index = 0;

  // gather site data
  for(int i = p_start; i < p_start+n; ++i) { // top row
    Cluster *c = sites[i].cluster;
    if(c) data[di++] = c->id;
    else data[di++] = -1;
  }
  for(int i = p_end-n; i < p_end; ++i) { // bottom row
    Cluster *c = sites[i].cluster;
    if(c) data[di++] = c->id;
    else data[di++] = -1;
  }

  // gather cluster data
  for(int i = p_start; i < p_start+n; ++i) { // top row
    Cluster *c = sites[i].cluster;
    if(c && c->id != -1 && in_array(c->id, seen_cluster_ids, seen_index)) {
      seen_cluster_ids[seen_index++] = c->id;
      p_stats[3]++;
      data[di++] = c->id;
      data[di++] = c->size;
      data[di++] = c->width;
      data[di++] = c->height;
      for(int k = 0; k < n; ++k) data[di++] = c->rows[k];
      for(int k = 0; k < n; ++k) data[di++] = c->cols[k];
    }
  }
  for(int i = p_end-n; i < p_end; ++i) { // bottom row
    Cluster *c = sites[i].cluster;
    if(c && c->id != -1 && in_array(c->id, seen_cluster_ids, seen_index)) {
      seen_cluster_ids[seen_index++] = c->id;
      p_stats[3]++;
      data[di++] = c->id;
      data[di++] = c->size;
      data[di++] = c->width;
      data[di++] = c->height;
      for(int k = 0; k < n; ++k) data[di++] = c->rows[k];
      for(int k = 0; k < n; ++k) data[di++] = c->cols[k];
    }
  }

  MPI_Send(p_stats, 4, MPI_INT, MASTER, TAG, MPI_COMM_WORLD);
  MPI_Send(data, 2*n + nc_attrs*p_stats[3], MPI_INT, MASTER, TAG, MPI_COMM_WORLD);
  free(data);
}

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
      int t_start = p_start + get_start(n, np_rows, tid, nt_workers);
      int nt_rows = get_n_rows(np_rows, tid, nt_workers);
      if(tid < nt_workers) percolate(sites, b, n, tid, t_start, nt_rows, nt_workers, p_start, np_rows, t_clusters[tid], &nt_clusters[tid]);
    }
    if(nt_workers > 1) join_clusters(sites, b, n, nt_workers, p_start, np_rows);
    if(rank > MASTER) send_clusters(rank, sites, n, nt_workers, t_clusters, nt_clusters, p_start, p_end);

    // receive cluster data
    else if(rank == MASTER && n_workers > 1) { 
      int nc_attrs = 4 + 2*n; // number of ints that describes a cluster
      int p_stats[n_workers-1][4]; // num clusters, max cluster size, col perc, num border clusters
      int **data = calloc(n_workers-1, sizeof(int*));
      for(int i = 0; i < n_workers-1; ++i) {
        MPI_Recv(p_stats[i], 4, MPI_INT, i+1, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int d_size = 2*n + nc_attrs*p_stats[i][3];
        data[i] = calloc(d_size, sizeof(int));
        MPI_Recv(data[i], d_size, MPI_INT, i+1, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      // process clusters
      Cluster*** p_clusters = calloc(n_workers, sizeof(Cluster**));
      int* np_clusters = calloc(n_workers, sizeof(int));
      for(int i = 0; i < n_workers; ++i) p_clusters[i] = calloc(max_clusters(n, get_n_rows(n, i, n_workers)), sizeof(Cluster*));
      printf("Hmm\n");
      // condense master clusters into array
      for(int tid = 0; tid < nt_workers; ++tid) {
        for(int i = 0; i < nt_clusters[tid]; ++i) {
          Cluster *c = t_clusters[tid][i];
          if(c->id != -1) p_clusters[0][np_clusters[0]++] = c;
        }
      }
      printf("Him\n");
      // convert worker data into clusters and add to array
      for(int i = 0; i < n_workers-1; ++i) {
        int d_size = 2*n + nc_attrs*p_stats[i][3];
        for(int j = 2*n; j < d_size; j+=nc_attrs) { // loop through process clusters
          Cluster *c = calloc(1, sizeof(Cluster));
          c->rows = calloc(n, sizeof(short)); c->cols = calloc(n, sizeof(short));
          c->id = data[i][j]; c->size = data[i][j+1]; c->width = data[i][j+2]; c->height = data[i][j+3];
          for(int k = 0; k < n; ++k) c->rows[k] = data[i][j+4+k];
          for(int k = 0; k < n; ++k) c->cols[k] = data[i][j+n+4+k];
          p_clusters[i+1][np_clusters[i+1]++] = c;
        }
      }
      printf("Hah\n");
      // add site cluster pointers to sites
      for(int i = 1; i < n_workers; ++i) {
        int p_start = get_start(n, n, i, n_workers);
        int p_end = p_start + n*get_n_rows(n, i, n_workers);
        for(int j = p_start; j < p_start+n; ++j) { // top row
          for(int k = 0; k < np_clusters[i]; ++k) {
            if(p_clusters[i][k]->id == data[i-1][j-p_start]) {
              sites[j].cluster = p_clusters[i][k];
              break;
            }
          }
        }
        for(int j = p_end-n; j < p_end; ++j) { // bottom row
          for(int k = 0; k < np_clusters[i]; ++k) {
            if(p_clusters[i][k]->id == data[i-1][j-p_end+2*n]) {
              sites[j].cluster = p_clusters[i][k];
              break;
            }
          }
        }
      }
      join_clusters(sites, b, n, n_workers, 0, n);
    }
    print_site_array(sites, n);
  }
  MPI_Finalize();
  exit(EXIT_SUCCESS);
}