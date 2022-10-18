from subprocess import run
from sys import argv

if __name__ == '__main__':
  print(run(['srun', '--mpi=pmix', *[a for a in argv[1:] if 'nodes' in a], '../src/percolate', *[a for a in argv[1:] if 'nodes' not in a]], capture_output=True).stdout.decode('utf-8'), end='')