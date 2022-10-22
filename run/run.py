import subprocess, random, sys

s = '''#!/bin/bash
#SBATCH --job-name=jp
#SBATCH --output=/dev/null
#SBATCH --partition=cits5507
#SBATCH --nodes={}
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=59:59
#SBATCH --mem-per-cpu=4G
{}

module load gcc/9.4.0 openmpi/4.0.5

python3 sub.py "$@" >> results.csv

# cd ../src
# make
# cd ../run
# srun --mpi=pmix ../src/percolate -b -f ../lattice/bond20_2.txt 20 4
'''

ncs = [1,3,2,4]
nl = ['n010', 'n[011-013]', 'n[010-011]', 'n[010-013]']

random.seed()
r = random.randint(0, int(1e6))

for nc, n in zip(ncs, nl):
  with open('sub.sh', 'w') as f:
    f.write(s.format(nc, f'#SBATCH --nodelist={n}' if len(sys.argv)>1 else ''))
  args = ['sbatch', 'sub.sh', str(r)]
  subprocess.run(args)
