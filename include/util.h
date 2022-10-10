/*
 * CITS5507 HPC PROJECT 1
 * LATTICE PERCOLATION IN PARALLEL
 * 
 * Jasper Paterson 22736341
 * Allen Antony 22706998
 */

#ifndef UTIL_H
#define UTIL_H

/**
 * @return int number of digits in integer x, required for printing lattices
 */
int num_digits(int);
int min(int, int);
int ceiling_divide(int, int);
short in_array(int, int*, int);

#endif