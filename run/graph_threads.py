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

# plt.show(block=True)

plt.savefig('../graph/threads.png')