/*
 * CITS5507 HPC PROJECT 2
 * LATTICE PERCOLATION USING MPI AND OPENMP
 * 
 * Jasper Paterson 22736341
 */

#include "../include/site.h"

/**
 * @return short* square short array with length n and occupation probability p
 */
short* short_array(int n, float p)
{
  short *sites = calloc(n*n, sizeof(short));
  for(int i = 0; i < n*n; ++i) {
    if((double)rand()/RAND_MAX < p) sites[i] = 1;
  }
  return sites;
}

/** 
 * @return Site* array of site pointers based on short array
 */
Site* site_array(short* a, int n)
{
  Site* sites = calloc(n*n, sizeof(Site));
  if(!sites) return NULL;
  for(int i = 0; i < n*n; ++i) {
    sites[i].r = i/n;
    sites[i].c = i%n;
    if(a && a[i]) sites[i].occupied = 1;
  }
  return sites;
}

/** 
 * @return short* site lattice from file
 */
short* file_short_array(char* filename, int n) {
  short* s = calloc(n*n, sizeof(short));
  if(!s) return NULL;
  int ch, r = 0, c = 0;
  FILE* f = fopen(filename, "r");
  if(!f) return NULL;
  while((ch = getc(f)) != EOF) {
    if(ch == ' ' || (c == n && ch != '\n')) continue;
    else if(ch == '\n') {
      c = 0;
      ++r;
      if(r == n) break; // can have text underneath
      continue;
    }
    if(ch == 'X') s[r*n+c] = 1;
    ++c;
  }
  fclose(f);
  return s;
}

/**
 * @brief print site lattice 
 */
void print_short_array(short* a, int n)
{
  if(!a || n > PRINT_CUTOFF || n < 2) return;
  int s = num_digits(n-1);
  printf("\n ");
  for(int i = 0; i < s; ++i) printf(" ");
  for(int c = 0; c < n; ++c) printf("\033[0;34m %*d\033[0;30m", s, c);
  printf("\n\n");
  for(int r = 0; r < n; ++r) {
    printf("\033[0;34m%*d \033[0;30m", s, r);
    for(int c = 0; c < n; ++c) {
      for(int i = 0; i < s; ++i) printf(" ");
      if(a[r*n+c]) printf("\033[0;31mX\033[0;30m");
      else printf("O");
    }
    printf("\n");
  }
}

/**
 * @brief print site lattice 
 */
void print_site_array(Site* a, int n)
{
  if(!a || n > PRINT_CUTOFF || n < 2) return;
  int s = num_digits(n-1);
  printf("\n ");
  for(int i = 0; i < s; ++i) printf(" ");
  for(int c = 0; c < n; ++c) printf("\033[0;34m %*d\033[0;30m", s, c);
  printf("\n\n");
  for(int r = 0; r < n; ++r) {
    printf("\033[0;34m%*d \033[0;30m", s, r);
    for(int c = 0; c < n; ++c) {
      for(int i = 0; i < s; ++i) printf(" ");
      if(a[r*n+c].cluster) printf("\033[0;31m%3d\033[0;30m", a[r*n+c].cluster[0]);
      else printf(" O ");
    }
    printf("\n");
  }
}