/*
 * CITS5507 HPC PROJECT 1
 * LATTICE PERCOLATION IN PARALLEL
 * 
 * Jasper Paterson 22736341
 * Allen Antony 22706998
 */

#include "../include/cluster.h"

/** 
 * @return Cluster* pointer to a new cluster starting from site i  
 */
Cluster* cluster(int n, int n_threads, int i) {
  Cluster *cl = (Cluster*)calloc(1, sizeof(Cluster));
  if(!cl) return NULL;
  cl->rows = (short*)calloc(n, sizeof(short));
  cl->cols = (short*)calloc(n, sizeof(short));
  if(!cl->rows || !cl->cols) return NULL;
  cl->id = i;
  cl->rows[i/n] = 1;
  cl->cols[i%n] = 1;
  cl->height = 1;
  cl->width = 1;
  cl->size = 1;
  return cl;
}

void free_cluster(Cluster *cl) {
  free(cl->rows);
  free(cl->cols);
  free(cl);
}
