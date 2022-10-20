import subprocess, random

ns = [500] # range(500,4500,500)
ps = [0.25] # , 0.5, 0.75]
ncs = [1,2,3,4]
nts = [2]

LOOPS = 1

if __name__ == '__main__':
  l = 0
  random.seed()
  # with open('results.csv', 'a') as f:
  while l < LOOPS:
    for n in ns:
      for p in ps:
        r = random.randint(0,10000)
        for nc in ncs:
          for nt in nts:
            args = ['srun', '--mpi=pmix', f'--nodes={nc}', '../src/percolate', '-r', str(r), str(n), str(p), str(nt)]
            # args = ['srun', '--mpi=pmix', '../src/percolate', '-r', str(r), str(n), str(p), str(nt)]
            # args = ['../src/percolate', '-r', str(r), str(n), str(p), str(nt)]
            out = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            print(out.stdout.decode('utf-8'), end='')
            # print(out.stderr.decode('utf-8'), end='')
            l += 1