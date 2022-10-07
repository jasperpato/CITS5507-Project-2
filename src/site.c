/*
 * CITS5507 HPC PROJECT 1
 * LATTICE PERCOLATION IN PARALLEL
 * 
 * Jasper Paterson 22736341
 * Allen Antony 22706998
 */

#include "../include/site.h"

Site* site_array(short* a, int n, int start_index, int end_index)
{
  Site* sites = calloc(end_index-start_index, sizeof(Site));
  if(!sites) return NULL;
  for(int i = start_index; i < end_index; ++i) {
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

void print_short_array(short* a, int n, int num_rows)
{
  if(!a || n > PRINT_CUTOFF || n < 2) return;
  int s = num_digits(n-1);

  char str[n*n*s * 20 + 100]; // overestimate
  char *ptr = str;

  ptr += sprintf(ptr, "\n ");
  for(int i = 0; i < s; ++i) ptr += sprintf(ptr, " ");
  for(int c = 0; c < n; ++c) ptr += sprintf(ptr, "\033[0;34m %*d\033[0;30m", s, c);
  ptr += sprintf(ptr, "\n\n");
  for(int r = 0; r < num_rows; ++r) {
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

void print_site_array(Site* a, int n, int num_rows)
{
  if(!a || n > PRINT_CUTOFF || n < 2) return;
  int s = num_digits(n-1);

  char str[n*n*s * 20 + 100]; // overestimate
  char *ptr = str;

  ptr += sprintf(ptr, "\n ");
  for(int i = 0; i < s; ++i) ptr += sprintf(ptr, " ");
  for(int c = 0; c < n; ++c) ptr += sprintf(ptr, " %*d", s, c);
  ptr += sprintf(ptr, "\n\n");
  for(int r = 0; r < num_rows; ++r) {
    ptr += sprintf(ptr, "%*d ", s, r);
    for(int c = 0; c < n; ++c) {
      for(int i = 0; i < s; ++i) ptr += sprintf(ptr, " ");
      if(a[r*n+c].occupied) ptr += sprintf(ptr, "X");
      else ptr += sprintf(ptr, "O");
    }
    ptr += sprintf(ptr, "\n");
  }
  printf("%s", str);
}