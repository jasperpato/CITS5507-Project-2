/*
 * CITS5507 HPC PROJECT 1
 * LATTICE PERCOLATION IN PARALLEL
 * 
 * Jasper Paterson 22736341
 * Allen Antony 22706998
 */

#ifndef STACK_H
#define STACK_H 

#include "./site.h"

typedef struct Stack {
  Site** stack;
  int first, last;
} Stack;

/**
 * @return Stack* pointer to a stack data structure with capacity size
 */
Stack* stack(int);

short is_empty(Stack*);

void add(Stack*, Site*);

Site* pop(Stack*);

void free_stack(Stack*);

#endif