'''
CITS5507 HPC PROJECT 1
LATTICE PERCOLATION IN PARALLEL

Jasper Paterson 22736341
Allen Antony 22706998

MUST BE UPDATED FOR N_CPUS

Reads results file and graphs all n_threads on a plot. Keeps either n or p constant as specified.

USAGE: python3 graph.py [--fname RESULTS_FILE] [-n N | -p P]

N keep lattice size constant at N
P keep probability constant at P
'''

import csv
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser

P_RES = 1e-3

RESULTS_FILE = 'results.csv'

'''
Returns a list of dictionaries - one for each row in csvfile
'''
def read_file(fname):
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
Remove total_times that are further than (stds * std) from mean
'''
def remove_outliers(data, stds):
  for d in data.values():
    for x, times in d.items():
      new = []
      mean, std = np.mean(times) if times else 0, np.std(times) if times else 0
      for time in times:
        if abs(time-mean) <= stds * std: new.append(time)
      d[x] = new

'''
Return data mapping x to mean total time
'''    
def average(data):
  for d in data.values():
    for x, times in d.items():
      d[x] = np.mean(times) if times else 0

'''
Gather data from rows that match the const (could be n or p)
'''
def get_data(fname, cname, cval, group, stds):
  results = read_file(fname)
  x = 'p' if cname == 'n' else 'n'
  data = {}
  for row in results:
    t = tuple(row[g] for g in group)
    if abs(cval-row[cname]) < P_RES:
      if t not in data: data[t] = {} # map x to list of times
      if row[x] not in data[t]: data[t][row[x]] = []
      data[t][row[x]].append(row['total_time'])
  remove_outliers(data, stds)
  average(data)
  return data

'''
Graph all n_threads on one axis, keeping either n or p constant
'''
def graph(data, cname, cval, group, nsquared):

  for t, d in data.items():
    zs = [(k**2, v) for k, v in d.items()] if nsquared and cname == 'p' else list(d.items())
    xs, ys = [z[0] for z in zs], [z[1] for z in zs]
    plt.plot(xs, ys, label=str(t))

  plt.xlabel('Probability P' if cname == 'n' else ('Lattice size N*N' if nsquared else 'Lattice length N'))
  plt.ylabel('Mean total time (s)')
  plt.title(f"{'Probability P' if cname == 'p' else ('Lattice size N*N' if nsquared else 'Lattice length N')} = {cval}")
  plt.legend(title=str(group))
  plt.show(block=True)

if __name__ == '__main__':
  a = ArgumentParser()
  a.add_argument('--fname', default=RESULTS_FILE)
  a.add_argument('-n', type=int)
  a.add_argument('-p', type=float, default=0.4)
  a.add_argument('-s', type=float, default=4)
  a.add_argument('--n-squared', action='store_true')
  args = vars(a.parse_args())
  
  cname = 'n' if args['n'] else 'p'
  cval = args[cname]
  group = ('ncpus', 'nthreads')

  data = get_data(args['fname'], cname, cval, group, args['s'])
  graph(data, cname, cval, group, args['n_squared'])

  
