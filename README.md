# The openmp implementation of the GUPS
This repository holds the openmp implementation of the GUPS benchmarks.
The code in this repository is a copy version of the one in the 
[link](https://svn.mcs.anl.gov/repos/performance/benchmarks/randomaccess/hpcc/openmp).

# Configurations
## Number of threads
```bash
  $ export OMP_NUM_THREADS=16
```

## Working set size
 wall clock time format hh:mm:ss
 run on peroni with 16 threads
 |---------------+-----------------+----------------------------|
 | HPLMaxProcMem | wall clock time | max resident set size (kB) |
 |---------------+-----------------+----------------------------|
 |           2e7 |         3:55.86 | 264768                     |
 |           8e7 |        13:31.81 | 1,051,320                  |
 |           2e8 |        25:20.17 | 2,099,976                  |
 |---------------+-----------------+----------------------------|
