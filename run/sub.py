import subprocess, random, sys

ns = [500] # range(500,4500,500)
ps = [0.25, 0.5, 0.75]
nts = [1,2,4]

random.seed()
for n in ns:
  for p in ps:
    r = sys.argv[1] if len(sys.argv) > 1 else str(random.randint(0,10000))
    for nt in nts:
      args = ['srun', '--mpi=pmix', '../src/percolate', '-r', r, str(n), str(p), str(nt)]
      out = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      print(out.stdout.decode('utf-8'), end='')