#!/bin/bash -x
#PJM -g xg18i070 
#PJM --name "nbody.c"
#PJM -o "00stdout-%j.log"
#PJM -j
#PJM --norestart
#PJM -L rscgrp=regular-flat
#PJM -L node=1
#PJM -L elapse=02:00:00
#
module load intel
module load impi
module load mpi-fftw
ulimit -s unlimited
export OMP_NUM_THREADS=2
export OMP_STACKSIZE=512M
# Perform jobs
#time -p mpiexec.hydra -n 68 ./nbody.out > 00stdout.log 2>&1
time -p mpiexec.hydra -n 4 ./nbody.out > 00stdout.log 2>&1
