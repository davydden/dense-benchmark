#!/bin/bash
#   allocate 5 nodes with 40 CPU per node for 2 hour:
#PBS -l nodes=5:ppn=40,walltime=02:00:00
#   job name
#PBS -N DFT_dealii_cholesky
#   stdout and stderr files:
#PBS -o /home/woody/iwtm/iwtm108/deal.ii-mg-dft/_build_custom/benchmark/test.txt -e /home/woody/iwtm/iwtm108/deal.ii-mg-dft/_build_custom/benchmark/err.txt
#   first non-empty non-comment line ends PBS options

# submit with: qsub <name>.sh
cd /home/woody/iwtm/iwtm108/deal.ii-mg-dft/_build_custom/benchmark
module load openmpi/2.0.2-gcc

export OMP_NUM_THREADS=1

mpirun --npersocket 10 -np 1 cholesky 2>&1 | tee cholesky_1
mpirun --npersocket 10 -np 4 cholesky 2>&1 | tee cholesky_4
mpirun --npersocket 10 -np 9 cholesky 2>&1 | tee cholesky_9
mpirun --npersocket 10 -np 16 cholesky 2>&1 | tee cholesky_16
mpirun --npersocket 10 -np 25 cholesky 2>&1 | tee cholesky_25
mpirun --npersocket 10 -np 36 cholesky 2>&1 | tee cholesky_36
mpirun --npersocket 10 -np 49 cholesky 2>&1 | tee cholesky_49
mpirun --npersocket 10 -np 64 cholesky 2>&1 | tee cholesky_64
mpirun --npersocket 10 -np 81 cholesky 2>&1 | tee cholesky_81
mpirun --npersocket 10 -np 100 cholesky 2>&1 | tee cholesky_100
