/*
 * CITS5507 HPC PROJECT 2
 * LATTICE PERCOLATION USING MPI AND OPENMP
 * 
 * Jasper Paterson 22736341
 */

#include "../include/cluster.h"

/** 
 * @return int* pointer to a new cluster starting from site i
 */
int* cluster(int n, int i) {
  int* c = calloc(2+2*n, sizeof(int));
  if(!c) return NULL;
  c[0] = i; // id
  c[1] = 1; // size
  c[2+i/n] = 1; // rows
  c[2+n+i%n] = 1; // cols
  return c;
}
