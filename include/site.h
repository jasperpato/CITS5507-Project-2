/*
 * CITS5507 HPC PROJECT 1
 * LATTICE PERCOLATION IN PARALLEL
 * 
 * Jasper Paterson 22736341
 * Allen Antony 22706998
 */

#ifndef SITE_H
#define SITE_H

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#include "./constant.h"
#include "./util.h"
#include "./cluster.h"

typedef struct Site {
  int r, c;
  short occupied, seen;
  Cluster *cluster;
} Site;

/**
 * @param p probability of occupation. Negative p will skip the occupation step for performance (for bond percolation)
 * @return Site* pointer to site array
 */
Site* site_array(short*, int);
short* short_array(int, float);
short *file_short_array(char*, int);
void print_short_array(short*, int, int);
void print_site_array(Site*, int, int);
#endif