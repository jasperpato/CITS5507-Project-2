#!/bin/bash
#SBATCH --job-name=jp
#SBATCH --output=out.txt

python3 run.py "$@"