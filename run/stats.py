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
    t = tuple(r[p] for p in params)
    d[t] = d.get(t, 0) + 1

  for t, c in sorted(d.items(), key=lambda x: (x[1], x[0])):
    print(f'{t}: {c}')