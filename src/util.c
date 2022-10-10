/*
 * CITS5507 HPC PROJECT 1
 * LATTICE PERCOLATION IN PARALLEL
 * 
 * Jasper Paterson 22736341
 * Allen Antony 22706998
 */

/**
 * @return int number of digits in integer x, required for printing lattices
 */
int num_digits(int x) {
  int c = 0;
  while(x != 0) {
    x /= 10;
    ++c;
  }
  return c;
}

int min(int a, int b) {
  return a <= b ? a : b;
}

int ceiling_divide(int a, int b) {
  return a / b + (a % b != 0);
}

short in_array(int el, int *a, int size)
{
  for(int i = 0; i < size; ++i) {
    if(el == a[i]) return 0;
  }
  return 1;
}