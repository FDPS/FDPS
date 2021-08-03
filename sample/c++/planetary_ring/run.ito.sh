#!/bin/bash

#PJM -L "rscunit=ito-a" 
#PJM -L "rscgrp=ito-m-dbg"
#PJM -L "vnode=16"
#PJM -L "vnode-core=36"
#PJM -L "elapse=1:00:00"
#PJM -j
#PJM -X

module load intel/2019.4

NUM_NODES=$PJM_VNODES
NUM_CORES=36
NUM_PROCS=64
NUM_THREADS=9

export I_MPI_PERHOST=`expr $NUM_CORES / $NUM_THREADS`
export I_MPI_FABRICS=shm:ofi
export I_MPI_PIN_DOMAIN=omp
export I_MPI_PIN_CELL=core

export OMP_NUM_THREADS=$NUM_THREADS
export KMP_STACKSIZE=8m
export KMP_AFFINITY=compact

export I_MPI_HYDRA_BOOTSTRAP=rsh
export I_MPI_HYDRA_BOOTSTRAP_EXEC=/bin/pjrsh
export I_MPI_HYDRA_HOST_FILE=${PJM_O_NODEINF}

mpiexec.hydra -n $NUM_PROCS ./ring_mono.out -N 5000000 -a 0.001 -e 2  -r -m -L 16 -t 0.3 -T 10.0 -n 128 -l 16  -o result_para/snp --flag_para_out --sat_mode 1 --n_smp 200 


# Intel Xeon Gold 6154 (Skylake-SP) (3.0 GHz (Turbo 3.7 GHz), 18 core）x 2 / node

# ito-q-dbg     1/4	 9	      42GB	1h	短時間ジョブ用、ノード内複数ジョブ共有可
# ito-ss-dbg    1	 36	      168GB	1h	短時間ジョブ用
# ito-s-dbg     4	 36x4	      168GBx4	1h	短時間ジョブ用
# ito-m-dbg     16	 36x16	      168Gx16	1h	短時間ジョブ用
# ito-l-dbg     64	 36x64	      168GBx64	1h	短時間ジョブ用（2018年9月より試験運用）
# ito-xl-dbg    128	 36x128       168GBx128	1h	短時間ジョブ用（2018年9月より試験運用）
# ito-xxl-dbg	256	 36x256       168GBx256 1h      短時間ジョブ用（2018年9月より試験運用）
# ito-single	1/36	 1	      4.6GB	168h	ノード内複数ジョブ共有可
# ito-qq	1/12	 3	      14GB	168h	ノード内複数ジョブ共有可
# ito-q		1/4	 9	      42GB	168h    ノード内複数ジョブ共有可
# ito-ss	1	 36	      168G	96h     シングルノードジョブ専用
# ito-s		4	 36x4	       168GB    48h    4ノードまで利用可能
# ito-m		16	 36x16	       168Gx16  24h    16ノードまで利用可能
# ito-l		64	 36x64	       168GBx4  12h    64ノードまで利用可能
# ito-xl	128	 36x128       168Gx128	6h     128ノードまで利用可能
# ito-xxl	256	 36x256       168Gx256	6h     256ノードまで利用可能

# 理論演算性能3,456 GFLOPS / node (倍精度)
# メモリDDR4 192 GB / node
# メモリ帯域幅255.9 GB/sec / node
