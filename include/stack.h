/*
 * CITS5507 HPC PROJECT 2
 * LATTICE PERCOLATION USING MPI AND OPENMP
 * 
 * Jasper Paterson 22736341
 * Allen Antony 22706998
 */

#ifndef STACK_H
#define STACK_H 

#include <stdlib.h>

typedef struct Stack {
  int* stack;
  int first, last;
} Stack;

Stack* stack(int);

short is_empty(Stack*);

void add(Stack*, int);

int pop(Stack*);

void free_stack(Stack*);

#endif