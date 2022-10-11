/*
 * CITS5507 HPC PROJECT 1
 * LATTICE PERCOLATION IN PARALLEL
 * 
 * Jasper Paterson 22736341
 * Allen Antony 22706998
 */

#include "../include/site.h"

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

short* short_array(int n, float p)
{
  short *sites = calloc(n*n, sizeof(short));
  for(int i = 0; i < n*n; ++i) {
    if((double)rand()/RAND_MAX < p) sites[i] = 1;
  }
  return sites;
}

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
      printf("%s%3d%s", a[r*n+c].cluster ? "\033[0;31m" : "", a[r*n+c].cluster ? a[r*n+c].cluster->id : -1, a[r*n+c].cluster ? "\033[0;30m" : "");
    }
    printf("\n");
  }
}