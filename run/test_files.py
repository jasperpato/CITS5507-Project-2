'''
CITS5507 HPC PROJECT 2
LATTICE PERCOLATION USING MPI AND OPENMP
 
Jasper Paterson 22736341

This script iterates through the lattice files and tests that all results are correct for a set of nnodes and nthreads.
The results files are tagged with their correct num_clusters and max_cluster.

Usage: python3 test_files.py
'''

import sys, subprocess, os

nns = [1, 2, 3, 4]
nts = [1, 2, 4, 8]

dir = '../lattice/'
file = '../test/test{}.txt'

s = f'''#!/bin/bash
#SBATCH --job-name=jp
#SBATCH --output={file} # /dev/null
#SBATCH --partition=cits5507
#SBATCH --nodes={{}}
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task={min(nts)}
#SBATCH --mem-per-cpu=4G

module load gcc/9.4.0 openmpi/4.0.5

python3 test_files.py worker "$@"
'''

def get_n(fname):
  n = ''
  for c in fname:
    if c in '1234567890': n+=c
    else:
      if len(n): break
  return n

if 'worker' not in sys.argv:

  results = []
  os.makedirs('../test', exist_ok=True)

  for nc in nns:
    with open('sub.sh', 'w') as f: f.write(s.format(nc, nc))
    args = ['sbatch', '--wait', 'sub.sh', str(nc)]
    out = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  
    with open(file.format(nc), 'r') as f: results.append(f.read())

  correct = ''
  for fname in sorted(os.listdir(dir)):
    path = dir + fname
    with open(path, 'r') as f:
      num, max = [int(x) for x in f.readlines()[-2:]]
      for nt in nts: correct += f'{num}, {max}\n'

  if len(results) == 0 or (len(set(results)) == 1 and results[0] == correct):
    print('Correct')
  else:
    print(f'Incorrect\n\nCorrect:\n{correct}\nGot:')
    for r in results: print(r)

elif 'worker' in sys.argv:

   for fname in sorted(os.listdir(dir)):
    path = dir + fname
    for nt in nts:
      args = ['srun', '--mpi=pmix', '../src/percolate', '-s' if 'site' in path else '-b', '-f', path, get_n(fname), str(nt)]
      out = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      num, max = [int(x) for x in out.stdout.decode('utf-8').split(',')[5:7]]
      print(f'{num}, {max}')