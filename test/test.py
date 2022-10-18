from subprocess import run

if __name__ == '__main__':
  print(run(['srun', 'mpi-pmix', './percolate', '10', '0.5'], capture_output=True).stdout.decode('utf-8'))