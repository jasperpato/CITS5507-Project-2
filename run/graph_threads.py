'''
CITS5507 HPC PROJECT 2
LATTICE PERCOLATION USING MPI AND OPENMP
 
Jasper Paterson 22736341

Graph results, save as png.

Usage: python3 graph_threads.py
'''

from graph import read
import matplotlib.pyplot as plt
from pprint import pprint

rows = read('../results/threads4.csv')

d = {}
for r in rows:
  nt = r['nthreads']
  if nt not in d:
    d[nt] = []
  d[nt].append(r['total_time'])

pprint(d)

for nt in d:
  d[nt] = sum(d[nt]) / len(d[nt])

pprint(d)

plt.figure(figsize=(8, 6))
zs = sorted(d.items())
xs, ys = [z[0] for z in zs], [z[1] for z in zs]
plt.plot(xs, ys)
plt.xticks(range(8, 28+1, 2))

plt.xlabel('N_THREADS')
plt.ylabel('Mean total time (s)')
plt.title('N_NODES = 4, N = 5000, P = 0.3')

plt.savefig('../graph/threads.png')