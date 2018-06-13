#======================================================================
#   Numerical Libraries and Compilers
#======================================================================

FDPS_LOC = ../../../../../
FDPS_INC = -I$(FDPS_LOC)/src -I../../inc

CXXFLAGS_COMMON = -std=c++11 -O3 -ffast-math -funroll-loops $(FDPS_INC)

CUDAPATH?=/usr/local/cuda
NVCC=$(CUDAPATH)/bin/nvcc
NVCCFLAGS= -I$(CUDAPATH)/include -I$(CUDAPATH)/samples/common/inc -L$(CUDAPATH)/lib64 $(FDPS_INC) -Xcompiler "-O3"

LDFLAGS+=-L$(CUDAPATH)/lib64  -lcudart

# [1] Serial
CXX = g++
CXXFLAGS = $(CXXFLAGS_COMMON)
#CXXFLAGS+= -DPARTICLE_SIMULATOR_DEBUG_PRINT
# [2] OpenMP
#CXX = g++
#CXXFLAGS = $(CXXFLAGS_COMMON) -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
# [3] MPI
#CXX = mpicxx
#CXXFLAGS = $(CXXFLAGS_COMMON) -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -L$(MPI_PATH)/lib
# [4] MPI + OpenMP
#CXX = mpicxx
#CXXFLAGS = $(CXXFLAGS_COMMON) -DPARTICLE_SIMULATOR_MPI_PARALLEL -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

#----------------------------------------------------------------------
#   Source and Object files
#----------------------------------------------------------------------
SRC   = main.cpp
HED   = user_defined_class.h kernel.h

all:	gpu

gpu:	Makefile $(SRC) $(HED) kernel.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SRC) kernel.o -o gpu.out

reuse:	Makefile $(SRC) $(HED) kernel.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(SRC) kernel.o  -DENABLE_REUSE -o reuse.out

clean:
	rm -f *.out *.o

distclean: clean
	rm -rf result

kernel.o:	kernel.h kernel.cu
	$(NVCC) $(NVCCFLAGS) kernel.cu -c -o kernel.o

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
