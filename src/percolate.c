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

static short on_border(int idx, int n, int nt_workers, int p_start, int np_rows) {
  for(int tid = 0; tid < nt_workers; ++tid) {
    int t_start = p_start + get_start(n, np_rows, tid, nt_workers); // starting from zero
    int nt_rows = get_n_rows(np_rows, tid, nt_workers);
    int t_end = t_start + n*nt_rows;
    if((idx >= t_start && idx < t_start+n) || (idx >= t_end-n && idx < t_end)) return 1;
  }
  return 0;
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
      ++cl->size;
      if(!cl->rows[nb->r]) cl->height++;
      if(!cl->cols[nb->c]) cl->width++;
      cl->rows[nb->r] = 1;
      cl->cols[nb->c] = 1;
      if(on_border(nbs[i], n, nt_workers, p_start, np_rows)) cl->sites[cl->site_size++] = nbs[i]; // add index of neighbour to cluster's site index array
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
      if(on_border(i, n, nt_workers, p_start, np_rows)) sc->sites[sc->site_size++] = i;
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

static void join_clusters(Site* sites, Bond* b, int n, int nt_workers, int p_start, int np_rows) {
  for(int i = 0; i < nt_workers; ++i) {
    int r_end = p_start + get_start(n, np_rows, i, nt_workers) + n*get_n_rows(np_rows, i, nt_workers);
    int r_start = r_end-n;
    for(int i = r_start; i < r_end; ++i) { // loop along bottom row of region
      Site *s = &sites[i];
      Cluster *sc = s->cluster;
      if(!sc) continue;
      int nbi = bottom_neighbour(sites, b, n, i, p_start, np_rows);
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
      for(int j = 0; j < nc->site_size; ++j) {
        int idx = nc->sites[j];
        sc->sites[sc->site_size++] = idx;
        if(idx != nb->r*n+nb->c) sites[idx].cluster = sc; // don't overwrite neighbour until last
      }
      nc->id = -1; // mark as obsolete
      nb->cluster = sc; // now overwrite neighbour
    }
  }
}

void send_clusters(Site* sites, int n, int nt_workers, Cluster*** t_clusters, int* nt_clusters, int p_start, int p_end)
{
  int p_stats[3]; // num clusters, max cluster size, col perc boolean
  for(int tid = 0; tid < nt_workers; ++tid) {
    p_stats[0] += nt_clusters[tid];
    for(int i = 0; i < nt_clusters[tid]; ++i) {
      Cluster *c = t_clusters[tid][i];
      if(c->size > p_stats[1]) p_stats[1] = c->size;
      if(c->width == n) p_stats[2] = 1;
    }
  }
  // send current max cluster size and whether there is column percolation
  MPI_Send(p_stats, 2, MPI_INT, MASTER, TAG, MPI_COMM_WORLD);
  printf("Sent p_stats\n");

  int nborder_clusters = 0;
  int border_sites_size = 0; // keep track of total length of clusters' site arrays

  // cluster info
  int *cs = calloc(2+5*2*n, sizeof(int)); // first two ints are nborder_clusters and border_sites_size
  int j = 1;
  for(int i = p_start; i < p_end; ++i) { // top row
    Cluster *c = sites[i].cluster;
    if(c) {
      ++nborder_clusters;
      cs[j++] = c->id;
      cs[j++] = c->size;
      cs[j++] = c->width;
      cs[j++] = c->height;
      cs[j++] = c->site_size;
      border_sites_size += c->site_size;
    } else {
      for(int k = 0; k < 5; ++k) cs[j++] = -1; // no cluster at this site
    }
    if(i+1 == p_start+n) i = p_end-n; // jump to bottom row
  }
  cs[0] = nborder_clusters;
  cs[1] = border_sites_size;

  MPI_Send(cs, 2+5*2*n, MPI_INT, MASTER, TAG, MPI_COMM_WORLD);
  free(cs);
  printf("Sent cs\n");

  // row and col arrays
  int *rcs = calloc(2*n * nborder_clusters, sizeof(int));
  j = 0;
  for(int i = p_start; i < p_end; ++i) { // top row
    Cluster *c = sites[i].cluster;
    if(c) {
      for(int k = 0; k < n; ++k) rcs[j++] = c->rows[k];
      for(int k = 0; k < n; ++k) rcs[j++] = c->cols[k];
    }
    if(i+1 == p_start+n) i = p_end-n; // jump to bottom row
  }
  MPI_Send(rcs, 2*n * nborder_clusters, MPI_INT, MASTER, TAG, MPI_COMM_WORLD);
  free(rcs);
  printf("Sent rcs\n");

  // site arrays
  int *ss = calloc(border_sites_size, sizeof(int));
  j = 0;
  for(int i = p_start; i < p_end; ++i) { // top row
    Cluster *c = sites[i].cluster;
    if(c) {
      for(int k = 0; k < c->site_size; ++k) ss[j++] = c->sites[k];
    }
    if(i+1 == p_start+n) i = p_end-n; // jump to bottom row
  }
  MPI_Send(ss, border_sites_size, MPI_INT, MASTER, TAG, MPI_COMM_WORLD);
  free(ss);
  printf("Sent ss\n");
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
      int t_start = p_start + get_start(n, np_rows, tid, nt_workers);
      int nt_rows = get_n_rows(np_rows, tid, nt_workers);
      if(tid < nt_workers) percolate(sites, b, n, tid, t_start, nt_rows, nt_workers, p_start, np_rows, t_clusters[tid], &nt_clusters[tid]);
    }
    if(nt_workers > 1) join_clusters(sites, b, n, nt_workers, p_start, np_rows);
    if(rank > MASTER) send_clusters(sites, n, nt_workers, t_clusters, nt_clusters, p_start, p_end);
    if(rank == MASTER) { // receive cluster data
      if(n_workers > 1) {
        int p_stats[n_workers-1][3]; // num clusters, max cluster size, col perc bool
        int **cs = calloc(n_workers-1, sizeof(int*));
        int **rcs = calloc(n_workers-1, sizeof(int*));
        int **ss = calloc(n_workers-1, sizeof(int*));

        for(int i = 1; i < n_workers; ++i) {
          MPI_Recv(p_stats[i-1], 3, MPI_INT, i, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          cs[i-1] = calloc(2+5*2*n, sizeof(int)); // first two ints are nborder_clusters and border_sites_size
          MPI_Recv(cs[i-1], 2+5*2*n, MPI_INT, i, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          rcs[i-1] = calloc(2*n*cs[i][0], sizeof(int));
          MPI_Recv(rcs[i-1], 2*n*cs[i][0], MPI_INT, i, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          ss[i-1] = calloc(cs[i][1], sizeof(int));
          MPI_Recv(ss[i-1], cs[i][1], MPI_INT, i, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          printf("Rank %d nborder_clusters %d border_sites_size %d\n", i, cs[i][0], cs[i][1]);
        }

      }
    }
  }
  MPI_Finalize();
  exit(EXIT_SUCCESS);
}