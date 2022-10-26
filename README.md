## CITS5507 HPC PROJECT 2

## LATTICE PERCOLATION IN PARALLEL USING OPENMP AND MPI

<br>
Jasper Paterson 22736341
<br><br>
In this project I simulate percolation within a variable-length square lattice, and report the key details of the percolation: number of clusters, maximum cluster size, whether any cluster spans an entire row or column of the lattice, and time taken for the simulation.
<br><br>

### Usage

To execute a single percolation:

```console
cd src
bash make.sh
srun --mpi=pmix --partition=cits5507 --nodes=[N_NODES] --tasks-per-node=1 --cpus-per-task=[N_THREADS] ./percolate [-v] [-s | -b] [-r SEED] [[-f LATTICE_FILENAME N] | [N PROBABILITY]] [N_THREADS]
```

- [N_NODES] number of nodes to utilise in cluster
- [-v] verbose
- [-s | -b] site or bond percolation, default site
- [-r SEED] number to seed srand, default time(NULL)
- [-f LATTICE_FILENAME] file to scan lattice from
- [N PROBABILITY] size of lattice and probability of site occupation or bond
- [N_THREADS] number of threads to utilise, default 1

Examples:
```console
srun --mpi=pmix --partition=cits5507 --nodes=1 --tasks-per-node=1 --cpus-per-task=2 ./percolate 1000 0.5 2
```
```console
srun --mpi=pmix --partition=cits5507 --nodes=3 --tasks-per-node=1 --cpus-per-task=8 ./percolate -v -b -r 10 5000 0.5 8
```
```console
srun --mpi=pmix --partition=cits5507 --nodes=2 --tasks-per-node=1 --cpus-per-task=1 ./percolate -b -f ../lattice/bond20_2.txt 20
```
```console
srun --mpi=pmix --partition=cits5507 --nodes=2 --tasks-per-node=1 --cpus-per-task=1 ./percolate -v -f ../lattice/site12_2.txt 12
```

To generate results:

```console
cd run
python3 run.py
````

To test or graph results:
```console
cd run
python3 test_consistency.py
python3 test_files.py
python3 graph.py [-n N | -p P] [-t N_THREAD_LIST] [-c N_NODES_LIST] [--save]
```

