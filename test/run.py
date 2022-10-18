from subprocess import run
from sys import argv

if __name__ == '__main__':
  print(run(['srun', '--mpi=pmix', '../src/percolate', *argv[1:]], capture_output=True).stdout.decode('utf-8'), end='')