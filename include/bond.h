/*
 * CITS5507 HPC PROJECT 2
 * LATTICE PERCOLATION USING MPI AND OPENMP
 * 
 * Jasper Paterson 22736341
 * Allen Antony 22706998
 */

#ifndef BOND_H
#define BOND_H

#include <stdlib.h>
#include <stdio.h>

#include "./constant.h"
#include "./util.h"
#include "./site.h"

typedef struct Bond {
  short *v, *h;
} Bond;

Bond* bond(int, float);

Bond* file_bond(char*, int);

void print_bond(Bond*, int);

#endif