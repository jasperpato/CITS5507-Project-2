import subprocess, random, sys # , time

s = '''#!/bin/bash
#SBATCH --job-name=jp
#SBATCH --output=/dev/null
#SBATCH --partition=cits5507
#SBATCH --nodes={}
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=59:59
#SBATCH --mem-per-cpu=4G

module load gcc/9.4.0 openmpi/4.0.5

python3 sub.py "$@" # >> results.csv

# cd ../src
# make
# cd ../run
# srun --mpi=pmix ../src/percolate -b -f ../lattice/bond20_2.txt 20 4
'''

ncs = [1,3,2,4]

random.seed()
r = random.randint(0, int(1e6))

for nc in ncs:
  with open('sub.sh', 'w') as f: f.write(s.format(nc))
  args = ['sbatch', 'sub.sh', str(r)]
  subprocess.run(args)
  # time.sleep(5)
