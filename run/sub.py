import subprocess, random, sys

ns = [10000] # range(500,3500,500)
ps = [0.75] # , 0.5, 0.75]
nts = [1,2,4]

r = sys.argv[1]

for n in ns:
  for p in ps:
    for nt in nts:
      args = ['srun', '--mpi=pmix', '../src/percolate', '-r', r, str(n), str(p), str(nt)]
      out = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      print(out.stdout.decode('utf-8'), end='')