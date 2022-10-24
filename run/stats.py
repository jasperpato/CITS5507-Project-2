import csv

file = '../results/results{}.csv'

ncs = [1, 2, 3, 4]
params = ['n', 'p', 'ncpus', 'nthreads']

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

if __name__ == '__main__':
  d = {}
  rows = []
  for nc in ncs: rows.extend(read(file.format(nc)))
  for r in rows:
    t = tuple(r[p] for p in ('n', 'p', 'ncpus', 'nthreads'))
    if t not in d: d[t] = []
    d[t].append(r['total_time'])

  [print(x) for x in sorted(tuple((t, len(l)) for t, l in d.items()))]