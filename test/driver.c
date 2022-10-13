/*
 * CITS5507 HPC PROJECT 1
 * LATTICE PERCOLATION IN PARALLEL
 * 
 * Jasper Paterson 22736341
 * Allen Antony 22706998
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <errno.h>
#include <time.h>

#define ARG_LENGTH 40
#define ONAME "../test/results.csv"

/**
 * USAGE: ./driver [ >> RESULTS_FILE ] (pipe csv results from stdout to file)
 * 
 * Loops through n, p n_cpus and n_threads values and calls percolate program infinitely
 */
int main(int argc, char *argv[]) {

  int ns[] = {1000,2000,3000,4000,5000,6000}, n_size = 6;
  float ps[] = {0.2,0.4,0.8}; int p_size = 3;
  int nts[] = {1,4,8}, nt_size = 3;

  srand(time(NULL));
  while(1) {
    for(int n = 0; n < n_size; ++n) {
      for(int p = 0; p < p_size; ++p) {
        int seed = rand(); // same seed for each n_cpus and n_threads to make sure the results are identical
        for(int nt = 0; nt < nt_size; ++nt) {
          int pid = fork();
          if(pid == -1) exit(EXIT_FAILURE);
          else if(pid == 0) { // child
            // char *args[] = {"percolate", "-r", NULL, "-o", NULL, NULL, NULL, NULL, NULL};
            // args[2] = malloc(ARG_LENGTH*sizeof(char)); sprintf(args[2], "%d", seed); // same seed for all n_threads
            // args[4] = malloc(ARG_LENGTH*sizeof(char)); sprintf(args[4], "%s", ONAME);
            // args[5] = malloc(ARG_LENGTH*sizeof(char)); sprintf(args[5], "%d", ns[n]);
            // args[6] = malloc(ARG_LENGTH*sizeof(char)); sprintf(args[6], "%f", ps[p]);
            // args[7] = malloc(ARG_LENGTH*sizeof(char)); sprintf(args[7], "%d", nts[nt]);
            // execv("../src/percolate", args);

            char *args[] = {"--mpi=pmix", "../src/percolate", "-r", NULL, "-o", NULL, NULL, NULL, NULL, NULL};
            args[3] = malloc(ARG_LENGTH*sizeof(char)); sprintf(args[3], "%d", seed); // same seed for all n_threads
            args[5] = malloc(ARG_LENGTH*sizeof(char)); sprintf(args[5], "%s", ONAME);
            args[6] = malloc(ARG_LENGTH*sizeof(char)); sprintf(args[6], "%d", ns[n]);
            args[7] = malloc(ARG_LENGTH*sizeof(char)); sprintf(args[7], "%f", ps[p]);
            args[8] = malloc(ARG_LENGTH*sizeof(char)); sprintf(args[8], "%d", nts[nt]);
            execv("srun", args);
          }
          else { // parent
            int status;
            wait(&status); // wait for child
            if (!WIFEXITED(status)) printf("%s\n", strerror(status));
          }
        }
      }
    }
  }
  exit(EXIT_SUCCESS);
}