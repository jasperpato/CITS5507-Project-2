import os, subprocess, csv

dir = '../lattice/'
file = '../results/results{}.csv'

ncs = [1, 2, 3, 4]

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

def test_files(ncs, nts):
  for fname in sorted(os.listdir(dir)):
    path = dir + fname
    num, max, site = get_stats(path)
    print(fname, num, max, site)
    for nc in ncs:
      for nt in nts:
        args = ['srun', '--partition=cits5507', '--tasks-per-node=1', '--cpus-per-task=8', '--mem-per-cpu=4G', f'--nodes={nc}', '--mpi=pmix', '../src/percolate', '-s' if site else '-b', '-f', path, str(get_n(fname)), str(nt)]
        out = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(out.stdout.decode('utf-8') + '\n')

def read(fname):
  rows = []
  with open(fname, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
      for k, v in row.items():
        if '.' in v: row[k] = float(v)
        else: row[k] = int(v)
      rows.append(row)
  return rows

def test_consistency(rows):
  d = {}
  for r in rows:
    t = tuple(r[p] for p in ('n', 'p', 'seed'))
    if t not in d: d[t] = []
    d[t].append(tuple(r[c] for c in ('num_clusters', 'max_cluster', 'rperc', 'cperc')))
  for t in d:
    s = set(d[t])
    if len(s) > 1: print(*s)

def min_data_points(rows):
  d = {}
  for r in rows:
    t = tuple(r[p] for p in ('n', 'p', 'ncpus', 'nthreads'))
    d[t] = d.get(t, 0) + 1
  m = min(d.values())
  for t, l in sorted(d.items()):
    if l == m: print(f'{t}: {l}')

if __name__ == '__main__':
  rows = []
  for nc in ncs: rows.extend(read(file.format(nc)))

  test_consistency(rows)
  min_data_points(rows)