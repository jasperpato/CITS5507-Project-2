/*
 * CITS5507 HPC PROJECT 2
 * LATTICE PERCOLATION USING MPI AND OPENMP
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

/**
 * @return int min of a and b
 */
int min(int a, int b) {
  return a <= b ? a : b;
}
 /** 
  * @return int a divide b rounded up
  */
int ceiling_divide(int a, int b) {
  return a / b + (a % b != 0);
}

/**
 * @return short 1 iff el in array
 */
short in_array(int el, int *a, int size)
{
  for(int i = 0; i < size; ++i) {
    if(el == a[i]) return 1;
  }
  return 0;
}