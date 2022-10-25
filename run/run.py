'''
CITS5507 HPC PROJECT 2
LATTICE PERCOLATION USING MPI AND OPENMP
 
Jasper Paterson 22736341

Runs the percolation program for each ncpus, nthreads, n, p, stores the results in a file for each ncpus.

Usage: python3 run.py scheduler
'''

import subprocess, random, sys, os

ns =  [29000] # range(500, 5000+1, 500)
ps =  [0.3] # [x/10 for x in range(5, 10+1)]

ncs = [4] # [1, 2, 3, 4]
nts = [28] # [1, 2, 4, 8]

loops = 1

s = f'''#!/bin/bash
#SBATCH --job-name=jp
#SBATCH --output=/dev/null
#SBATCH --partition=cits5507
#SBATCH --nodes={{}}
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task={min(nts)}
#SBATCH --mem-per-cpu=4G

module load gcc/9.4.0 openmpi/4.0.5

python3 run.py worker "$@"
'''

file = '../results/results{}.csv'

if not loops: exit()

if 'scheduler' in sys.argv:

  for nc in ncs:
    if '--overwrite' in sys.argv or not os.path.exists(file.format(nc)):
      with open(file.format(nc), 'w') as f: f.write('n,p,ncpus,nthreads,seed,num_clusters,max_cluster,rperc,cperc,init_time,perc_time,tjoin_time,recv_time,pjoin_time,total_time\n')

  random.seed()
  seeds = [str(random.randint(0, int(1e6))) for _ in range(loops)]

  for nc in ncs:
    with open('sub.sh', 'w') as f: f.write(s.format(nc))
    args = ['sbatch', 'sub.sh', str(nc), *seeds]
    subprocess.run(args)

elif 'worker' in sys.argv:

  seeds = sys.argv[3:3+loops]

  for i in range(loops):
    for n in ns:
      for p in ps:
        for nt in nts:
          args = ['srun', '--mpi=pmix', '../src/percolate', '-r', seeds[i], str(n), str(p), str(nt)]
          out = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
          with open(file.format(sys.argv[2]), 'a') as f: f.write(out.stdout.decode('utf-8'))