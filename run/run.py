import subprocess, random

s = '''#!/bin/bash
#SBATCH --job-name=jp
#SBATCH --partition=cits5507
#SBATCH --nodes={}
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=59:59
#SBATCH --mem-per-cpu=4G

module load gcc/9.4.0 openmpi/4.0.5

python3 sub.py "$@" >> results.csv
'''

ncs = [1,2,3,4]

loops = 1
random.seed()

l = 0
while not loops or l < loops:
  r = str(random.randint(0,10000))
  for nc in ncs:
    with open('sub.sh', 'w') as f: f.write(s.format(nc))
    args = ['sbatch', 'sub.sh', r]
    subprocess.run(args)
  l+=1