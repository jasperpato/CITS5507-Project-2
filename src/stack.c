/*
 * CITS5507 HPC PROJECT 2
 * LATTICE PERCOLATION USING MPI AND OPENMP
 * 
 * Jasper Paterson 22736341
 * Allen Antony 22706998
 */

#include "../include/stack.h"

/**
 * @return Stack* pointer to a stack data structure with capacity size
 */
Stack* stack(int size) {
  Stack* st = calloc(1, sizeof(Stack));
  if(!st) return NULL;
  st->stack = calloc(size, sizeof(int));
  if(!st->stack) return NULL;
  return st;
}

/**
 * @return short 1 iff stack is empty
 */
short is_empty(Stack* st) {
  return st->first == st->last;
}

/**
 * @brief add int s to stack
 */
void add(Stack* st, int s) {
  st->stack[st->last++] = s;
}

/**
 * @return int at top of stack
 */
int pop(Stack* st) {
  return st->stack[st->first++];
}

void free_stack(Stack* st)
{
  free(st->stack);
  free(st);
}