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

#define N_MIN 100
#define N_MAX 3000
#define N_STEP 100

#define P_MIN 0.0
#define P_MAX 1.0
#define P_STEP 0.1

#define P_RES (int)1e3

#define NC_MIN 1
#define NC_MAX 4
#define NC_STEP 1

#define NT_MIN 1
#define NT_MAX 24
#define NT_STEP 1

#define ARG_LENGTH 40

#define RESULTS_FILE "../test/p2_results_kaya.csv"

/**
 * USAGE: ./driver [-k] [-f RESULTS_FILENAME] [-l LOG_FILENAME]
 * 
 * Loops through n, p and n_threads values and calls percolate program infinitely, which appends results to file.
 */
int main(int argc, char *argv[]) {
  char* resfile_name = RESULTS_FILE;
  char* logfile_name = NULL;
  int opt;
  while((opt = getopt(argc, argv, "f:l:")) != -1) {
    if(opt == 'f') resfile_name = optarg; 
    else if(opt == 'l') logfile_name = optarg;
    else {
      fprintf(stderr, "Usage: [-f RESULTS_FILENAME] [-l LOG_FILENAME]");
      exit(EXIT_FAILURE);
    }
  }
  srand(time(NULL));
  while(1) {
    for(int n = N_MIN; n <= N_MAX; n+=N_STEP) {
      for(int pi = (int)(P_MIN*P_RES); pi <= (int)(P_MAX*P_RES); pi+=(int)(P_STEP*P_RES)) {
        float p = (float)pi/P_RES;
        int seed = rand(); // same seed for each n_cpus and n_threads to make sure the results are identical
        for(int n_cpus = NC_MIN; n_cpus <= NC_MAX; n_cpus += NC_STEP) {
          for(int n_threads = NT_MIN; n_threads <= NT_MAX; n_threads+=NT_STEP) {
            int pid = fork();
            if(pid == -1) exit(EXIT_FAILURE);
            else if(pid == 0) { // child
              char *args[] = {"percolate", "-r", NULL, NULL, NULL, NULL, NULL, NULL};
              for(int i = 2; i < 7; ++i) args[i] = malloc(ARG_LENGTH*sizeof(char));
              sprintf(args[2], "%d", seed); // same seed for all n_threads
              sprintf(args[3], "%d", n);
              sprintf(args[4], "%f", p);
              sprintf(args[5], "%d", n_cpus);
              sprintf(args[6], "%d", n_threads);
              execv("../src/percolate", args);
            }
            else { // parent
              int status;
              wait(&status); // wait for child
              if (!WIFEXITED(status) && logfile_name != NULL) {
                FILE* logfile = fopen(logfile_name, "a");
                fprintf(logfile, "Error: %s with parameters N=%d, P=%f, NC=%d, NT=%d\n", strerror(status), n, p, n_cpus, n_threads);
                fclose(logfile);
              }
            }
          }
        }
      }
    }
  }
  exit(EXIT_SUCCESS);
}