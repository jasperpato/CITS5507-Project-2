'''
CITS5507 HPC PROJECT 1
LATTICE PERCOLATION IN PARALLEL

Jasper Paterson 22736341
Allen Antony 22706998
'''

'''
Reads results file and graphs all n_threads on a plot. Keeps either n or p constant as specified.

USAGE: python3 graph.py [-n N | -p P] [--fname RESULTS_FILE]

N keep lattice size constant at N
P keep probability constant at P
'''

import csv
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser

P_RES = 1e-3

'''
Returns a list of dictionaries - one for each row in csvfile
'''
def read_file(fname):
  rows = []
  with open(fname, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
      for i, (k, v) in enumerate(row.items()):
        if i == 1 or i > 7: row[k] = float(v)
        else: row[k] = int(v)
      rows.append(row)
  return rows

'''
Gather data from rows that match the const (could be n or p)
'''
def get_data(results, n_threads, const, c, time):
  x = 'p' if c == 'n' else 'n'
  data = tuple({} for _ in range(n_threads))
  for row in results:
    t = row['n_threads']-1
    if abs(const-row[c]) < P_RES:
      if row[x] not in data[t]: data[t][row[x]] = []
      data[t][row[x]].append(row[time])
  return data

'''
Remove total_times that are further than (stds * std) from mean
'''
def remove_outliers(data, n_threads, stds=4):
  count = 0
  for t in range(n_threads):
    for y, xs in data[t].items():
      new = []
      mean, std = np.mean(xs), np.std(xs)
      for i, x in enumerate(xs):
        if abs(x-mean) <= stds*std: new.append(x)
        else: count += 1
      data[t][y] = new
  print(f'Removed data points (stds={stds}): {count}')

'''
Replace lists of data with their mean to be graphed
'''    
def get_means(data, n_threads):
  lens, count, min_len, max_len = 0, 0, 0, 0
  
  m = tuple({} for _ in range(n_threads))
  for t in range(n_threads):
    for k, v in data[t].items():
      lens += len(v)
      if len(v) < min_len or min_len == 0: min_len = len(v)
      if len(v) > max_len: max_len = len(v)
      count += 1
      m[t][k] = np.mean(v)

  print(f'Minimum data points per parameter set: {min_len}')

  return m

'''
Graph all n_threads on one axis, keeping either n or p constant
'''
def graph(results, n_threads, const, c, time, n_squared=False):
  data = get_data(results, n_threads, const, c, time)
  remove_outliers(data, n_threads)
  m = get_means(data, n_threads)

  for t in range(n_threads):
    plt.plot([k**2 for k in m[t].keys()] if n_squared else m[t].keys(), m[t].values(), label=f'{t+1}')
  plt.xlabel('Probability P' if c == 'n' else ('Lattice size N*N' if n_squared else 'Lattice length N'))
  plt.ylabel('Mean total time (s)')
  plt.title(f"{'Probability P' if c == 'p' else ('Lattice size N*N' if n_squared else 'Lattice length N')} = {const}")
  plt.legend(title='Num threads')
  plt.show(block=True)

if __name__ == '__main__':
  a = ArgumentParser()
  a.add_argument('--fname', default='results_4t.csv')
  a.add_argument('--n-threads', default=4, type=int)
  a.add_argument('--time', default='total_time', help='Type of time to track. Options = {init_time, perc_time, join_time, scan_time, total_time}')
  a.add_argument('-n', type=int)
  a.add_argument('-p', type=float)
  a.add_argument('--n-squared', action='store_true')
  args = a.parse_args()

  results = read_file(args.fname)
  
  if(args.n): graph(results, args.n_threads, args.n, 'n', args.time)
  elif(args.p): graph(results, args.n_threads, args.p, 'p', args.time, args.n_squared)
