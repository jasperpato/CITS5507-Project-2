'''
CITS5507 HPC PROJECT 2
LATTICE PERCOLATION USING MPI AND OPENMP
 
Jasper Paterson 22736341

This script iterates through the lattice files and tests that all results are correct for a set of ncpus and nthreads.
The results files are tagged with their correct results.

Usage: python3 test.py scheduler
'''

import sys, subprocess, os

ncs = [1, 2, 3, 4]
nts = [1, 2, 4, 8]

s = f'''#!/bin/bash
#SBATCH --job-name=jp
#SBATCH --output=/dev/null
#SBATCH --partition=cits5507
#SBATCH --nodes={{}}
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task={min(nts)}
#SBATCH --mem-per-cpu=4G

module load gcc/9.4.0 openmpi/4.0.5

python3 test.py worker "$@"
'''

dir = '../lattice/'
file = '../test/test{}.txt'

def get_n(fname):
  n = ''
  for c in fname:
    if c in '1234567890': n+=c
    else:
      if len(n): break
  return n

if 'scheduler' in sys.argv:

  results = []
  if not os.path.exists('../test'): os.makedirs('../test')

  for nc in ncs:
    with open('sub.sh', 'w') as f: f.write(s.format(nc, nc))
    with open(file.format(nc), 'w') as f: pass
    args = ['sbatch', '--wait', 'sub.sh', str(nc)]
    out = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    with open(file.format(nc), 'r') as f: results.append(f.read())

  correct = ''
  for fname in sorted(os.listdir(dir)):
    path = dir + fname
    with open(path, 'r') as f:
      num, max = [int(x) for x in f.readlines()[-2:]]
      for nt in nts: correct += f'{num}, {max}\n'

  print('Correct' if len(set(results)) == 1 and results[0] == correct else 'Incorrect')

elif 'worker' in sys.argv:

   for fname in sorted(os.listdir(dir)):
    path = dir + fname
    for nt in nts:
      args = ['srun', '--mpi=pmix', '../src/percolate', '-s' if 'site' in path else '-b', '-f', path, get_n(fname), str(nt)]
      out = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      num, max = [int(x) for x in out.stdout.decode('utf-8').split(',')[5:7]]
      
      with open(file.format(sys.argv[2]), 'a') as f:
        f.write(f'{num}, {max}\n')