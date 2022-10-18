/*
 * CITS5507 HPC PROJECT 2
 * LATTICE PERCOLATION USING MPI AND OPENMP
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
  short* s = calloc(2*n, sizeof(short));
  if(!cl || !s) return NULL;
  cl->rows = s; cl->cols = s+n;
  cl->id = i;
  cl->rows[i/n] = 1;
  cl->cols[i%n] = 1;
  cl->height = 1;
  cl->width = 1;
  cl->size = 1;
  return cl;
}
