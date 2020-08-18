PIKG sample program of cold collapse N-body simulation

# overview
This sample codes generates interaction kernel between stars.
PIKG generates C++ file "kernel.cpp" and C header file "kernel.h" from kernel.pikg.
User code "main.c" includes "kernel.h" and call kernel function "calc_gravity" which is defined inside "kernel.cpp".

# how to compile and run
make
./nbody.out
