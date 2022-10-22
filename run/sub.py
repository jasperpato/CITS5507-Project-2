import subprocess, random, sys

file = 'results.csv'

ns = range(500,3500,500)
ps = [0.25, 0.5, 0.75]
nts = [1, 2, 4]

random.seed()


for n in ns:
  for p in ps:
    r = sys.argv[1] if len(sys.argv) > 1 else str(random.randint(0,int(1e6)))
    for nt in nts:
      args = ['srun', '--mpi=pmix', '../src/percolate', '-r', r, str(n), str(p), str(nt)]
      out = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      with open(file, 'a') as f: f.write(out.stdout.decode('utf-8'))
      # print(out.stderr.decode('utf-8'), end='')