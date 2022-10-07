/*
 * CITS5507 HPC PROJECT 1
 * LATTICE PERCOLATION IN PARALLEL
 * 
 * Jasper Paterson 22736341
 * Allen Antony 22706998
 */

#include "../include/site.h"

Site* site_array(short* a, int n, int start_row, int num_rows)
{
  Site* sites = calloc(n*num_rows, sizeof(Site));
  if(!sites) return NULL;
  for(int i = 0; i < num_rows*n; ++i) {
    sites[i].r = (start_row*n+i)/n;
    sites[i].c = (start_row*n+i)%n;
    if(a && a[start_row*n+i]) sites[i].occupied = 1;
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
  if(!a || n > 40 || n < 2) return;
  int s = num_digits(n-1);

  char str[n*n*s * 20 + 100]; // overestimate
  char *ptr = str;

  ptr += sprintf(ptr, "\n ");
  for(int i = 0; i < s; ++i) ptr += sprintf(ptr, " ");
  for(int c = 0; c < n; ++c) ptr += sprintf(ptr, "\033[0;34m %*d\033[0;30m", s, c);
  ptr += sprintf(ptr, "\n\n");
  for(int r = 0; r < n; ++r) {
    ptr += sprintf(ptr, "\033[0;34m%*d \033[0;30m", s, r);
    for(int c = 0; c < n; ++c) {
      for(int i = 0; i < s; ++i) ptr += sprintf(ptr, " ");
      if(a[r*n+c]) ptr += sprintf(ptr, "\033[0;31mX\033[0;30m");
      else ptr += sprintf(ptr, "O");
    }
    ptr += sprintf(ptr, "\n");
  }
  printf("%s", str);
}