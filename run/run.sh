#!/bin/bash
#SBATCH --job-name=jp
#SBATCH --output=out.txt
#SBATCH --partition=cits5507
###SBATCH --nodes=4
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=59:59
#SBATCH --mem-per-cpu=4G

module load gcc/9.4.0 openmpi/4.0.5

#mpicc -fopenmp -o test_mpi test_mpi.c
#srun --mpi=pmix ./test_mpi

python3 run.py