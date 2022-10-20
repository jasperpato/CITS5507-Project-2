import subprocess, random

ns = range(500,4500,500)
ps = [0.25, 0.5, 0.75]
ncs = [1,2,3,4]
nts = [1,2,4]

LOOPS = 1

if __name__ == '__main__':
  l = 0
  random.seed()
  with open('results.csv', 'a') as f:
    while l < LOOPS:
      for n in ns:
        for p in ps:
          r = random.randint(0,10000)
          for nc in ncs:
            for nt in nts:
              args = ['srun', '--mpi=pmix', f'--nodes={nc}', f'--ntasks-per-node={nt}', '../src/percolate', '-r', str(r), str(n), str(p), str(nt)]
              # args = ['../src/percolate', '-r', str(r), str(n), str(p), str(nt)]
              out = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
              print(out.stdout.decode('utf-8'))
              # print(out.stderr.decode('utf-8'), end='')
              l += 1