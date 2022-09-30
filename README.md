## CITS5507 HPC PROJECT 1

## LATTICE PERCOLATION IN PARALLEL

<br>
Jasper Paterson 22736341  <br>
Allen Antony 22706998  <br>
<br>
In this project we simulate percolation within a variable-length square lattice, and report the key details of the percolation: number of clusters, maximum cluster size, whether any cluster spans an entire row or column of the lattice, and time taken for the simulation.

### Usage

To run:

```console
cd src
make
./percolate [-s | -b] [-r SEED] [-p RESULTS_FILENAME] [[-f LATTICE_FILENAME] | [N PROBABILITY]] [N_THREADS]
```

- [-s | -b] site or bond percolation, default site
- [-v] silence printing
- [-r SEED] number to seed srand, default time(NULL)
- [-p RESULTS_FILENAME] file to append the results of the percolation
- [-f LATTICE_FILENAME] file to scan lattice from
- [N PROBABILITY] size of lattice and probability of site occupation or bond
- [N_THREADS] number of threads to utilise, default 1

Examples:

- Most common usage (N = 500, P = 0.4, N_THREADS= 1, then 2)

  ```console
  ./percolate -s 500 0.4

  ./percolate -b 500 0.4 2
  ```

- Load lattice from file

  ```console
  ./percolate -b -f ../lattice/bond20.txt 20 1

  ./percolate -b -f ../lattice/bond20.txt 20 4
  ```

- Same seed for identical results

  ```console
  ./percolate -r 999 10 0.4

  ./percolate -r 999 10 0.4 3
  ```

- Silence printing and append results to file
  ```console
  ./percolate -v -p ../test/results_4t.csv 100 0.1
  ```

To generate results:

```console
cd test
make
./driver
```

To graph results:

```console
cd test
python3 graph.py [-n N | -p P] [--fname RESULTS_FILE]
```

- [-n N] keep n constant at N and graph all N_THREADS
- [-p P] keep p constant at P and graph all N_THREADS
- [--fname RESULTS_FILE] results file name

Examples:

```console
python3 graph.py -n 2500

python3 graph.py -p 0.4
```
