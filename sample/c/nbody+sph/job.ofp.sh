#!/bin/bash -x
#PJM -g xg18i070 
#PJM --name "nbodysph.c"
#PJM -o "00stdout-%j.log"
#PJM -j
#PJM --norestart
#PJM -L rscgrp=regular-flat
##PJM -L node=64
##PJM -L node=1
#PJM -L node=4
##PJM -L node=32
#PJM -L elapse=48:00:00
#
module load intel
module load impi
ulimit -s unlimited
export OMP_NUM_THREADS=2
export OMP_STACKSIZE=512M
mkdir -p result
# Perform jobs
#time -p mpiexec.hydra -n 2176 -prepend-rank ./nbodysph.out > 00stdout.log 2>&1
#time -p mpiexec.hydra -n 68 -prepend-rank ./nbodysph.out > 00stdout.log 2>&1
#time -p mpiexec.hydra -n 34 ./nbodysph.out > 00stdout.log 2>&1
time -p mpiexec.hydra -n 136 ./nbodysph.out > 00stdout.log 2>&1
#time -p mpiexec.hydra -n 1088 ./nbodysph.out > 00stdout.log 2>&1
#time -p ./ring_code.x > 00stdout.log 2>&1
