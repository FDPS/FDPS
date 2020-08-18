PIKG sample program of molecular dynamics simulation of Lennard-Jones fluid

# overview
This sample codes generates interaction kernel between Lennard-Jones fluids under periodic boundary condition.
PIKG generates C++ header file "kernel.hpp" from kernel.pikg.
See kernel_orig.hpp which is the hand-written source code of interaction kernel which pikg is expected to generate.

# how to compile and run
make
./main.out

# expected result
If the total energy is conserved, the program shows "test passed" at the end of run.
See ./exact which is the detailed log of succcessful run
