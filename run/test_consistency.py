'''
CITS5507 HPC PROJECT 2
LATTICE PERCOLATION USING MPI AND OPENMP
 
Jasper Paterson 22736341

Scan all results and make sure results are the same for all (n, p, seed) sets.

Usage: python3 test_consistency.py
'''

from graph import read

ncs = [1, 2, 3, 4]
file = '../results/results{}.csv'

if __name__ == '__main__':
  rows = []
  for nc in ncs: rows.extend(read(file.format(nc)))

  d = {}
  for r in rows:
    t = tuple(sorted(r[c] for c in ('n', 'p', 'seed')))
    if t not in d: d[t] = []
    d[t].append(tuple(sorted(r[c] for c in ('num_clusters', 'max_cluster', 'rperc', 'cperc'))))
  
  lens = tuple(len(r) for r in d.values())
  print('Avg:', round(sum(lens)/len(lens),1)) # average executions per (n, c, seed) set

  for t in d:
    if len(set(d[t])) > 1:
      print('Incorrect')
      for r in d[t]: print(r)
    # no output if all consistent
