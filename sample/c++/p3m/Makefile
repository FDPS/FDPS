#======================================================================
#   Numerical Libraries and Compilers
#======================================================================
FFTW_LOC = /nfshome/namekata/gnu
FFTW_INC = -I$(FFTW_LOC)/include
FFTW_LIB = -L$(FFTW_LOC)/lib -lfftw3_mpi -lfftw3_omp -lfftw3 -lm
FDPS_LOC = ../../../
FDPS_INC = -I$(FDPS_LOC)/src -I$(FDPS_LOC)/src/particle_mesh
FDPS_LIB = -L$(FDPS_LOC)/src/particle_mesh -lpm
LDFLAGS = $(FFTW_LIB) $(FDPS_LIB) 

CXXFLAGS_COMMON = -std=c++11 -O3 -ffast-math -funroll-loops $(FFTW_INC) $(FDPS_INC)
#CXXFLAGS_COMMON = -std=c++11 -O0 -Wall -Wextra -ftrapv -fexceptions -g3 $(FFTW_INC) $(FDPS_INC)
# [1] Serial
#CXX = g++
#CXXFLAGS = $(CXXFLAGS_COMMON) -DPARTICLE_SIMULATOR_DEBUG_PRINT
# [2] OpenMP
#CXX = g++
#CXXFLAGS = $(CXXFLAGS_COMMON) -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
# [3] MPI
#CXX = mpicxx
#CXXFLAGS = $(CXXFLAGS_COMMON) -DPARTICLE_SIMULATOR_MPI_PARALLEL 
# [4] MPI + OpenMP
CXX = mpicxx
CXXFLAGS = $(CXXFLAGS_COMMON) -DPARTICLE_SIMULATOR_MPI_PARALLEL -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp

#----------------------------------------------------------------------
#   WORKDIR
#----------------------------------------------------------------------
WORKDIR = ./work

#----------------------------------------------------------------------
#   Source and Object files
#----------------------------------------------------------------------
SRCS	= main.cpp

TARGET = p3m.x

$(TARGET): Makefile $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)
	mv $(TARGET) $(WORKDIR)

clean:
	rm -f $(TARGET) *.o

