/*
 * CITS5507 HPC PROJECT 1
 * LATTICE PERCOLATION IN PARALLEL
 * 
 * Jasper Paterson 22736341
 * Allen Antony 22706998
 */

#include "../include/cluster.h"

/** 
 * @return Cluster* pointer to a new cluster starting from site (r, c)  
 */
Cluster* cluster(int n, int n_threads, int r, int c) {
  Cluster *cl = (Cluster*)calloc(1, sizeof(Cluster));
  if(!cl) return NULL;
  cl->rows = (short*)calloc(n, sizeof(short));
  cl->cols = (short*)calloc(n, sizeof(short));
  if(!cl->rows || !cl->cols) return NULL;
  cl->id = r*n+c;
  cl->rows[r] = 1;
  cl->cols[c] = 1;
  cl->height = 1;
  cl->width = 1;
  cl->size = 1;
  cl->sites = calloc(2*n*n_threads, sizeof(int)); // only needs to store border sites
  if(!cl->sites) return NULL;
  return cl;
}

/** 
 * @return CPArray* an array of CPArrays, one for each thread
 */
CPArray* cluster_array(int n_threads, int num_clusters) {
  CPArray* cpa = calloc(n_threads, sizeof(CPArray));
  for(int i = 0; i < n_threads; ++i) cpa[i].cls = calloc(num_clusters, sizeof(Cluster*));
  return cpa;
}

void free_cluster(Cluster *cl) {
  free(cl->rows);
  free(cl->cols);
  free(cl->sites);
  free(cl);
}

void free_cparray(CPArray* cpa, int n_threads) {
  for(int i = 0; i < n_threads; ++i) free(cpa[i].cls);
  free(cpa);
}
