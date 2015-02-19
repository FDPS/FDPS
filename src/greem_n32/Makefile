PS_PATH  = -I../src_parallel
PS_PATH += -I../src_parallel/particle_mesh

#CCC = time mpiFCCpx
#CCC = time FCCpx
#CFLAGS = -Kfast,openmp -NRtrap
#CFLAGS = -Kfast,openmp
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CFLAGS += -Xg
#CFLAGS += -Nfjcex
#CFLAGS += -Krestp=all
#CFLAGS = -Kfast,openmp,restp=all

CCC = mpicxx-openmpi-gcc49
#CCC = mpicxx-mpich-gcc49
#CCC = g++

CFLAGS = -O3 -ffast-math -funroll-loops -DMPICH_IGNORE_CXX_SEEK #-Wall
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#TOWDIM = PARTICLE_SIMULATOR_TOW_DIMENSION

INCLUDE_FFTW = -I/Users/tanikawa/mylib/fftw-3.3.4/include
LIB_FFTW = -L/Users/tanikawa/mylib/fftw-3.3.4/lib -lfftw3f -lfftw3f_mpi
LIB_PM   = -L../src_parallel/particle_mesh -lpm

CFLAGS += $(INCLUDE_FFTW)
CLIBS += $(LIB_FFTW) $(LIB_PM)

OBJ = main.o
HDR = myheader.hpp
ETC = Makefile

run: $(OBJ) $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -o run $(OBJ) $(CLIBS)

main.o: main.cc $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -c main.cc

clean:
	rm -f *~
distclean:
	make clean
	rm -f *.o
allclean:
	make distclean
	rm -f run
