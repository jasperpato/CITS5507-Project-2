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

#define N_MIN 500
#define N_MAX 1000
#define N_STEP 500

#define P_MIN 0.2
#define P_MAX 1.0
#define P_STEP 0.2

#define P_RES (int)1e3

#define NT_MIN 2
#define NT_MAX 8
#define NT_STEP 2

#define ARG_LENGTH 40

#define ONAME "../test/results.csv"

/**
 * USAGE: ./driver [ >> RESULTS_FILE ] (pipe csv results from stdout to file)
 * 
 * Loops through n, p n_cpus and n_threads values and calls percolate program infinitely
 */
int main(int argc, char *argv[]) {
  srand(time(NULL));
  while(1) {
    for(int n = N_MIN; n <= N_MAX; n+=N_STEP) {
      for(int pi = (int)(P_MIN*P_RES); pi <= (int)(P_MAX*P_RES); pi+=(int)(P_STEP*P_RES)) {
        float p = (float)pi/P_RES;
        int seed = rand(); // same seed for each n_cpus and n_threads to make sure the results are identical
        for(int n_threads = NT_MIN; n_threads <= NT_MAX; n_threads+=NT_STEP) {
          int pid = fork();
          if(pid == -1) exit(EXIT_FAILURE);
          else if(pid == 0) { // child
            char *args[] = {"srun", "--mpi-pmix", "./percolate", "-r", NULL, "-o", NULL, NULL, NULL, NULL, NULL};
            args[4] = malloc(ARG_LENGTH*sizeof(char)); sprintf(args[4], "%d", seed); // same seed for all n_threads
            args[6] = malloc(ARG_LENGTH*sizeof(char)); sprintf(args[6], "%s", ONAME);
            args[7] = malloc(ARG_LENGTH*sizeof(char)); sprintf(args[6], "%d", n);
            args[8] = malloc(ARG_LENGTH*sizeof(char)); sprintf(args[7], "%f", p);
            args[9] = malloc(ARG_LENGTH*sizeof(char)); sprintf(args[9], "%d", n_threads);
            execv("../src/percolate", args);
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