#!/bin/bash

#SBATCH -p d024h
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH -c 4
#SBATCH -J nbody
#SBATCH -o stdout.log
#SBATCH -e stderr.log

module load gnu/openmpi165
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
mpirun -mca btl ^openib -np `expr ${SLURM_NTASKS_PER_NODE} \* ${SLURM_JOB_NUM_NODES}` ./sph.out
RETCODE=$?
exit ${RETCODE}
