#======================================================================
#   Numerical Libraries and Compilers
#======================================================================
FFTW_LOC = /work4a/local/fftw
FFTW_INC = -I$(FFTW_LOC)/include
FFTW_LIB = -L$(FFTW_LOC)/lib -lfftw3f_mpi -lfftw3f_omp -lfftw3f -lm
#FFTW_LIB = -L$(FFTW_LOC)/lib -lfftw3_mpi -lfftw3_omp -lfftw3 -lm
FDPS_LOC = ../../../
PM_ROOT  = $(FDPS_LOC)/src/particle_mesh
FDPS_INC = -I$(FDPS_LOC)/src -I$(PM_ROOT)
FDPS_LIB = -L$(PM_ROOT) -lpm
LDFLAGS = $(FDPS_LIB) $(FFTW_LIB) 

CXXFLAGS_COMMON = -std=c++11 -O3 -ffast-math -funroll-loops $(FFTW_INC) $(FDPS_INC)
# [1] MPI only
CXX = mpicxx
CXXFLAGS = -DPARTICLE_SIMULATOR_MPI_PARALLEL 
# [2] MPI + OpenMP
#CXX = mpicxx
#CXXFLAGS = -DPARTICLE_SIMULATOR_MPI_PARALLEL -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

CXXFLAGS += $(CXXFLAGS_COMMON)

#----------------------------------------------------------------------
#   WORKDIR
#----------------------------------------------------------------------
WORKDIR = ./work

#----------------------------------------------------------------------
#   Source and Object files
#----------------------------------------------------------------------
SRCS	= main.cpp

TARGET = p3m.x

$(TARGET): Makefile $(PM_ROOT)/libpm.a $(WORKDIR) $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)
	mv $(TARGET) $(WORKDIR)

$(PM_ROOT)/libpm.a: $(PM_ROOT)/param_fdps.h
	cd $(PM_ROOT) && $(MAKE) allclean libpm.a

$(WORKDIR):
	mkdir -p $(WORKDIR)

clean:
	rm -f $(TARGET) *.o

distclean: clean
	rm -rf $(WORKDIR)

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
