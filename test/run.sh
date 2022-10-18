#!/bin/bash
#SBATCH --job-name=jp
#SBATCH --output=out.txt
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1

python3 run.py "$@"