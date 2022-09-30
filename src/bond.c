/*
 * CITS5507 HPC PROJECT 1
 * LATTICE PERCOLATION IN PARALLEL
 * 
 * Jasper Paterson 22736341
 * Allen Antony 22706998
 */

#include "../include/bond.h"

/**
 * @param n size of lattice
 * @param p probability of a bond
 * @return Bond* pointer to a bond struct that contains vertical and horizontal bond information
 */
Bond* bond(int n, float p)
{
  Bond* b = calloc(1, sizeof(Bond));
  if(!b) return NULL;
  b->v = calloc(n*n, sizeof(short));
  b->h = calloc(n*n, sizeof(short));
  if(!b->v || !b->h) return NULL;
  for(int i = 0; i < n*n; ++i) {
    if((double)rand()/RAND_MAX < p) b->v[i] = 1;
    else b->v[i] = 0; 
    if((double)rand()/RAND_MAX < p) b->h[i] = 1;
    else b->h[i] = 0;
  }
  return b;
}

/**
 * @param n assumed lattice size
 * @return pointer to a bond struct scanned from a file
 */
Bond* file_bond(char* filename, int n)
{
  Bond* b = calloc(1, sizeof(Bond));
  if(!b) return NULL;
  b->v = calloc(n*n, sizeof(short));
  b->h = calloc(n*n, sizeof(short));
  if(!b->v || !b->h) return NULL;
  int ch, r = 0, c = 0, v_count = 0;
  short gap = 1, vert = 1;
  FILE* f = fopen(filename, "r");
  if(!f) return NULL;
  while((ch = getc(f)) != EOF) {
    if(v_count == n) break; // can have text underneath
    if(vert) {
      if(ch != '\n' && c == n) continue; // skip characters at end of line
      if(ch == ' ') {
        if(gap) gap = 0;
        else {
          gap = 1;
          ++c;
        }
        continue;
      } else if(ch == '\n') {
        vert = 0; c = 0; gap = 0;
      } else if(ch == '|') {
        b->v[r*n+c] = 1;
        gap = 1;
        ++c;
      }
    } else { // in a horizontal line
      if(ch != '\n' && c == n) continue; // skip characters at end of line
      if(ch == '\n') {
        vert = 1; gap = 1;
        c = 0;
        ++r; ++v_count;
      } else if(ch == 'O') {
        gap = 0;
      } else {
        b->h[r*n+c] = (ch == '-') ? 1 : 0;
        ++c;
      }
    }
  }
  return b;
}

void print_bond(Bond* b, int n)
{
  if(!b || n > 40 || n < 2) return;
  int s = num_digits(n-1);
  printf("\n");
  for(int i = 0; i < s; ++i) printf(" ");
  printf(" ");
  for(int c = 0; c < n; ++c) printf("\033[0;34m %*d\033[0;30m", s, c);
  printf("\n\n");
  for(int r = 0; r < n; ++r) {
    printf(" ");
    for(int i = 0; i < s; ++i) printf(" ");
    for(int c = 0; c < n; ++c) {
      for(int i = 0; i < s; ++i) printf(" ");
      printf(b->v[r*n+c] ? "\033[0;31m|\033[0;30m" : " ");
    }
    printf("\033[0;34m\n%*d \033[0;30m", s, r);
    for(int c = 0; c < n; ++c) {
      for(int i = 0; i < s; ++i) printf(b->h[r*n+c] ? "\033[0;31m-\033[0;30m" : " ");
      printf("O");
    }
    printf("\n");
  }
}

void free_bond(Bond* b) {
  free(b->v);
  free(b->h);
  free(b);
}