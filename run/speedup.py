import sys
from graph import read
from pprint import pprint

ncs = [1, 2, 3, 4]
file = '../results/results{}.csv'

rows = []
for nc in ncs: rows.extend(read(file.format(nc)))

if len(sys.argv) < 3: exit()

n1 = tuple(int(x) for x in sys.argv[1].split(','))
n2 = tuple(int(x) for x in sys.argv[2].split(','))

d = {
  n1: {},
  n2: {}
}
for r in rows:
  t1 = tuple(r[c] for c in ('nnodes', 'nthreads'))
  if t1 not in d: continue
  t2 = tuple(r[c] for c in ('n', 'p'))
  if t2 not in d[t1]: d[t1][t2] = []
  d[t1][t2].append(r['total_time'])

for t1 in d:
  for t2 in d[t1]:
    d[t1][t2] = sum(d[t1][t2]) / len(d[t1][t2])

# pprint(d)

speedups = {}

for t2 in d[n1]:
  if t2[0] % 1000 == 0 and t2[1] == 0.3:
    speedups[t2] = round(d[n1][t2] / d[n2][t2], 3)

pprint(speedups)