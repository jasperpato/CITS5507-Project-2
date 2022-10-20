import subprocess, random, sys

ns = [500] # range(500,4500,500)
ps = [0.25] # , 0.5, 0.75]
nts = [2]

LOOPS = 1
NCPUS = sys.argv[1] if len(sys.argv) > 1 else 1

if __name__ == '__main__':
  l = 0
  random.seed()
  while not LOOPS or l < LOOPS:
    for n in ns:
      for p in ps:
        r = random.randint(0,10000)
        for nt in nts:
          args = ['srun', '--mpi=pmix', f'--nodes={NCPUS}', '../src/percolate', '-r', str(r), str(n), str(p), str(nt)]
          out = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
          print(out.stdout.decode('utf-8'), end='')
          l += 1