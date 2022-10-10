/*
 * CITS5507 HPC PROJECT 1
 * LATTICE PERCOLATION IN PARALLEL
 * 
 * Jasper Paterson 22736341
 * Allen Antony 22706998
 */

#ifndef CLUSTER_H  
#define CLUSTER_H

#include <stdlib.h>
#include <mpi.h>

typedef struct Cluster {
  short *rows, *cols;
  int id, size, height, width, site_size; // id is initial site s index: s.r*n+s.c
  int* sites;
} Cluster;

/** 
 * @return Cluster* pointer to a new cluster starting from site (r, c)  
 */
Cluster* cluster(int, int, int);

void free_cluster(Cluster*);

#endif