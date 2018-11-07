#!/bin/bash -x
#PJM -g xg18i070 
#PJM --name "nbody.cpp"
#PJM -o "00stdout-%j.log"
#PJM -j
#PJM --norestart
#PJM -L rscgrp=regular-flat
##PJM -L node=64
#PJM -L node=2
#PJM -L elapse=48:00:00
#PJM -m b,e
#PJM --mail-list daisuke.namekata@riken.jp
##PJM --stg-transfiles all
##PJM --stgin "./ring_code.x ./"
##PJM -s
##PJM --spath "stat-%j.log"
#
module load intel
module load impi
ulimit -s unlimited
export OMP_NUM_THREADS=4
#export OMP_NUM_THREADS=1
export OMP_STACKSIZE=512M
mkdir -p result
# Perform jobs
#time -p mpiexec.hydra -n 1088 -prepend-rank ./nbody.out > 00stdout.log 2>&1
time -p mpiexec.hydra -n 34 -prepend-rank ./nbody.out > 00stdout.log 2>&1
#time -p ./ring_code.x > 00stdout.log 2>&1
