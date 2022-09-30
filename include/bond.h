/*
 * CITS5507 HPC PROJECT 1
 * LATTICE PERCOLATION IN PARALLEL
 * 
 * Jasper Paterson 22736341
 * Allen Antony 22706998
 */

#ifndef BOND_H
#define BOND_H

#include <stdlib.h>
#include <stdio.h>

#include "./util.h"
#include "./site.h"

typedef struct Bond {
  short *v, *h;
} Bond;

/**
 * @param n size of lattice
 * @param p probability of a bond
 * @return Bond* pointer to a bond struct that contains vertical and horizontal bond information
 */
Bond* bond(int, float);

/**
 * @param n assumed lattice size
 * @return pointer to a bond struct scanned from a file
 */
Bond* file_bond(char*, int);

void print_bond(Bond*, int);

void free_bond(Bond*);

#endif