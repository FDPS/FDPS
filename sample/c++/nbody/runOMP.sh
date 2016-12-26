#!/bin/bash

#SBATCH -p d024h
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH -J nbody
#SBATCH -o stdout.log
#SBATCH -e stderr.log

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
./nbody.out
RETCODE=$?
exit ${RETCODE}
