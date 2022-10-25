/*
 * CITS5507 HPC PROJECT 2
 * LATTICE PERCOLATION USING MPI AND OPENMP
 * 
 * Jasper Paterson 22736341
 */

#ifndef SITE_H
#define SITE_H

#include <stdlib.h>
#include <stdio.h>

#include "./constant.h"
#include "./util.h"
#include "./cluster.h"

typedef struct Site {
  int r, c;
  short occupied, seen;
  int* cluster;
} Site;

Site* site_array(short*, int);

short* short_array(int, float);

short *file_short_array(char*, int);

void print_short_array(short*, int);

void print_site_array(Site*, int);

#endif