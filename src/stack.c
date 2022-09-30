/*
 * CITS5507 HPC PROJECT 1
 * LATTICE PERCOLATION IN PARALLEL
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
  st->stack = calloc(size, sizeof(Site*));
  if(!st->stack) return NULL;
  return st;
}

short is_empty(Stack* st) {
  return st->first == st->last;
}

void add(Stack* st, Site* s) {
  st->stack[st->last++] = s;
}

Site* pop(Stack* st) {
  return st->stack[st->first++];
}

void free_stack(Stack* st) {
  free(st->stack);
  free(st);
}