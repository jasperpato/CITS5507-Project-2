'''
CITS5507 HPC PROJECT 2
LATTICE PERCOLATION USING MPI AND OPENMP
 
Jasper Paterson 22736341

Graph results.

Usage: python3 graph.py [-n N | -p P] [-t T_LIST] [-c C_LIST]
'''

import csv
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser

P_RES = 1e-3

ncs = [1, 2, 3, 4]
rfile = '../results/results{}.csv'
sfile = '../graph/graph_{}{}{}{}{}.png'

'''
Returns a list of dictionaries - one for each row in csvfile
'''
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

'''
Return data mapping x to mean total time
'''    
def average(data, verbose=False):
  lens = []
  for d in data.values():
    for x, times in d.items():
      lens.append(len(times))
      d[x] = np.mean(times) if times else 0
  if verbose: print(f'Min data points per parameter set: {min(lens) if lens else 0}')

'''
Gather data from rows that match the cval (cname could be n or p)
'''
def get_data(fname, cname, cval, group, ncpus, nthreads):
  results = []
  for nc in ncs: results.extend(read(fname.format(nc)))
  x = 'p' if cname == 'n' else 'n'
  data = {}
  for row in results:
    t = tuple(row[g] for g in group)
    if abs(cval-row[cname]) < P_RES and (not ncpus or (ncpus and row['ncpus'] in ncpus)) and (not nthreads or (nthreads and row['nthreads'] in nthreads)):
      if t not in data: data[t] = {} # map x to list of times
      if row[x] not in data[t]: data[t][row[x]] = []
      data[t][row[x]].append(row['total_time'])
  average(data)
  return data

'''
Graph all n_threads on one axis, keeping either n or p constant
'''
def graph(data, cname, cval, group, nsquared, save=False, file='', cs=None, ts=None):

  plt.figure(figsize=(8, 6))
  for t, d in data.items():
    zs = sorted(((k**2, v) for k, v in d.items()) if nsquared and cname == 'p' else d.items())
    xs, ys = [z[0] for z in zs], [z[1] for z in zs]
    plt.plot(xs, ys, label=str(t))

  plt.xlabel('Probability P' if cname == 'n' else ('Lattice size N*N' if nsquared else 'Lattice length N'))
  plt.ylabel('Mean total time (s)')
  plt.title(f"{'Probability P' if cname == 'p' else ('Lattice size N*N' if nsquared else 'Lattice length N')} = {cval}")
  
  # sort legend
  handles, labels = plt.gca().get_legend_handles_labels()
  labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
  plt.legend(handles, labels, title=str(group))
  
  if save:
    cval = f'0{int(cval*100)}' if cname == 'p' else cval
    c = '_c' + ''.join(map(lambda x: f'-{x}', cs)) if cs else ''
    t = '_t' + ''.join(map(lambda x: f'-{x}', ts)) if ts else ''
    plt.savefig(file.format(cname, cval, c, t, '_ns' if nsquared else ''))

  else: plt.show(block=True)

if __name__ == '__main__':
  a = ArgumentParser()
  a.add_argument('-n', type=int)
  a.add_argument('-p', type=float, default=0.4)
  a.add_argument('-t')
  a.add_argument('-c')
  a.add_argument('--all', action='store_true')
  a.add_argument('--n-squared', action='store_true')
  a.add_argument('--save', action='store_true')
  args = vars(a.parse_args())
  
  group = ('ncpus', 'nthreads')

  if args['all']:
    
    # all sets p nsquared
    for p in (0.2, 0.4, 0.6, 0.8):
      ncpus, nthreads = None, None
      data = get_data(rfile, 'p', p, group, ncpus, nthreads)
      graph(data, 'p', p, group, True, True, sfile, ncpus, nthreads)

    # all sets n
    for n in (1000, 2000, 4000, 5000):
      ncpus, nthreads = None, None
      data = get_data(rfile, 'n', n, group, ncpus, nthreads)
      graph(data, 'n', n, group, False, True, sfile, ncpus, nthreads)

  else:
    cname = 'n' if args['n'] else 'p'
    cval = args[cname]

    ncpus    = [int(c) for c in args['c'].split(',')] if args['c'] else None
    nthreads = [int(t) for t in args['t'].split(',')] if args['t'] else None

    data = get_data(rfile, cname, cval, group, ncpus, nthreads)
    graph(data, cname, cval, group, args['n_squared'], args['save'], sfile, ncpus, nthreads)

  
