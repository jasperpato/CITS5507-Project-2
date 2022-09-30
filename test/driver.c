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

#define NT_MIN 1
#define NT_MAX 4
#define NT_STEP 1

#define ARG_LENGTH 40

#define RESULTS_FILE "../test/results_4t.csv"


/**
 * USAGE: ./driver [-f RESULTS_FILENAME] [-l LOG_FILENAME]
 * 
 * Loops through n, p and n_threads values and calls percolate program infinitely, which appends results to file.
 */
int main(int argc, char *argv[]) {
  srand(time(NULL));
  int opt;
  char* resfile_name = NULL;
  char* logfile_name = NULL;

  while((opt = getopt(argc, argv, "f:l:")) != -1) {
    switch(opt) {
      case 'f':
        resfile_name = optarg;
        break;
      case 'l':
        logfile_name = optarg;
        break;
      default:
        fprintf(stderr, "Usage: [-f RESULTS_FILENAME] [-l LOG_FILENAME]");
        exit(EXIT_FAILURE);
    }
  }

  while(1) {
    for(int n = N_MIN; n <= N_MAX; n+=N_STEP) {
      for(int pi = (int)(P_MIN*P_RES); pi <= (int)(P_MAX*P_RES); pi+=(int)(P_STEP*P_RES)) {
        float p = (float)pi/P_RES;
        int seed = rand(); // same seed for each n_threads to make sure the results are identical
        for(int nt = NT_MIN; nt <= NT_MAX; nt+=NT_STEP) {
          int pid = fork();

          if(pid == -1) { 
            exit(EXIT_FAILURE);
          }

          else if(pid == 0) { // child
            char *args[] = {"percolate", "-v", "-p", NULL, "-r", NULL, NULL, NULL, NULL, NULL};
            args[3] = malloc(ARG_LENGTH*sizeof(char));
            for(int i = 5; i < 9; ++i) args[i] = malloc(ARG_LENGTH*sizeof(char));
            sprintf(args[3], "%s", resfile_name != NULL ? resfile_name : RESULTS_FILE);
            sprintf(args[5], "%d", seed); // same seed for all n_threads
            sprintf(args[6], "%d", n);
            sprintf(args[7], "%f", p);
            sprintf(args[8], "%d", nt);
            execv("../src/percolate", args);
            exit(EXIT_SUCCESS);
          }

          else {
            int status;
            wait(&status); // wait for child
            if (!WIFEXITED(status) && logfile_name != NULL) {
              FILE* logfile = fopen(logfile_name, "a");
              fprintf(logfile, "Error: %s with parameters N=%i, P=%f, NT=%i\n", strerror(status), n, p, nt);
              fclose(logfile);
            }
          }
        }
      }
    }
  }
  
  exit(EXIT_SUCCESS);
}