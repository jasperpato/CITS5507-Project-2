from subprocess import run

ns = [1000,2000]
ps = [0.3, 0.7]
ncs = [1, 2]
nts = [1]

if __name__ == '__main__':
  
  for n in ns:
    for p in ps:
      for nc in ncs:
        for nt in nts:
          args = ['srun', '--mpi=pmix', f'--nodes={nc}', f'--ntasks-per-node={nt}', '../src/percolate', str(n), str(p), str(nt)]
          print(run(args, capture_output=True).stdout.decode('utf-8'), end='')