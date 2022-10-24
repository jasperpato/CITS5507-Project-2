import os, subprocess

dir = '../lattice/'

ncs = [1]
nts= [1]

def get_n(fname):
  n = ''
  for c in fname[4:]:
    if c == '_': break
    n+=c
  return int(n)

def get_stats(fname):
  try:
    stats = [int(l.strip()) for l in open(fname, 'r').readlines()[-2:]]
    stats.append(('site' in fname))
    return tuple(stats)
  except: return 0, 0, False

if __name__ == '__main__':
  for fname in sorted(os.listdir(dir)):
    path = dir + fname
    num, max, site = get_stats(path)
    print(fname, num, max, site)
    for nc in ncs:
      for nt in nts:
        args = ['srun', '--partition=cits5507', '--tasks-per-node=1', '--cpus-per-task=8', '--mem-per-cpu=4G', f'--nodes={nc}', '--mpi=pmix', '../src/percolate', '-s' if site else '-b', '-f', path, str(get_n(fname)), str(nt)]
        out = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(out.stdout.decode('utf-8') + '\n')