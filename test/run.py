from subprocess import run
from sys import argv

if __name__ == '__main__':
  args = ['srun', '--mpi=pmix', '--nodes=2', '../src/percolate', *[a for a in argv[1:] if 'nodes' not in a]]
  print(args)
  print(run(args, capture_output=True).stdout.decode('utf-8'), end='')