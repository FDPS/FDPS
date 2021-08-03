PIKG sample program of cold collapse N-body simulation

# overview
This sample codes generates interaction kernel between stars.
PIKG generates C++ header file "kernel.hpp" from kernel.pikg.
See kernel_orig.hpp which is the hand-written source code of interaction kernel which pikg is expected to generate.

# how to compile and run
make
./nbody.out

# for AVX2 mode
use_avx2=yes make
./nbody.out
# for AVX-512 mode
use_avx512=yes make
./nbody.out
# for ARM SVE mode
use_arm_sve=yes make
./nbody.out
# for CUDA mode
use_cuda=yes make
./nbody.out
