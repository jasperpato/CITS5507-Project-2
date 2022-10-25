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
file = '../results/results{}.csv'

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
def average(data):
  lens = []
  for d in data.values():
    for x, times in d.items():
      lens.append(len(times))
      d[x] = np.mean(times) if times else 0
  print(f'Min data points per parameter set: {min(lens) if lens else 0}')

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
def graph(data, cname, cval, group, nsquared):

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
  
  plt.show(block=True)

if __name__ == '__main__':
  a = ArgumentParser()
  a.add_argument('-n', type=int)
  a.add_argument('-p', type=float, default=0.4)
  a.add_argument('-t')
  a.add_argument('-c')
  a.add_argument('--n-squared', action='store_true')
  args = vars(a.parse_args())
  
  cname = 'n' if args['n'] else 'p'
  cval = args[cname]
  group = ('ncpus', 'nthreads')

  ncpus    = [int(c) for c in args['c'].split(',')] if args['c'] else None
  nthreads = [int(t) for t in args['t'].split(',')] if args['t'] else None

  data = get_data(file, cname, cval, group, ncpus, nthreads)
  graph(data, cname, cval, group, args['n_squared'])

  
