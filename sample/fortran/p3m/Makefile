#-----------------------------------------------------------------------
#   Configuration of compile
#-----------------------------------------------------------------------
# Variables that must be specified by users
# (i) Variables related to FDPS
FDPS_LOC = ../../../
PM_ROOT  = $(FDPS_LOC)/src/particle_mesh
FDPS_INC = -I$(FDPS_LOC)/src -I$(PM_ROOT)
FDPS_LIB = -L$(PM_ROOT) -lpm
FDPS_FTN_MOD_DIR = $(FDPS_LOC)/src/fortran_interface/modules
FDPS_FTN_IF_GENERATOR = $(FDPS_LOC)/scripts/gen_ftn_if.py

# (ii) Variables related to FFTW
FFTW_LOC = /work4a/local/fftw
FFTW_INC = -I$(FFTW_LOC)/include
FFTW_LIB = -L$(FFTW_LOC)/lib -lfftw3f_mpi -lfftw3f_omp -lfftw3f -lm

# (iii) Variables to specify compilers and compile options
FC=mpif90
CXX=mpic++
MPI_LOC  = /usr/lib64/openmpi
MPI_LIB  = -L$(MPI_LOC)/lib -lmpi -lmpi_f77 -lmpi_f90
# [Option 1] w/o optimization
#FCFLAGS = -std=f2003 -O0 -Wall
#CXXFLAGS = -Wall -Wextra -ftrapv -fexceptions -g3 $(FDPS_INC)
# [Option 2] w/ optimization 
FCFLAGS = -std=f2003 -O3 -ffast-math -funroll-loops -finline-functions
CXXFLAGS = -O3 -ffast-math -funroll-loops $(FFTW_INC) $(FDPS_INC)
LDFLAGS = -lgfortran $(FDPS_LIB) $(FFTW_LIB) $(MPI_LIB) 
# Activate PM module
FCFLAGS  += -DPARTICLE_SIMULATOR_USE_PM_MODULE
CXXFLAGS += -DPARTICLE_SIMULATOR_USE_PM_MODULE
# OpenMP options
#FCFLAGS  += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
# MPI options
FCFLAGS  += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

#-----------------------------------------------------------------------
#   Source files
#-----------------------------------------------------------------------
%.o : %.F90
	$(FC) $(FCFLAGS) -c $<
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<

SRC_USER_DEFINED_TYPE = user_defined.F90
SRC_USER = f_main.F90
SRC_FDPS_MOD = $(wildcard $(FDPS_FTN_MOD_DIR)/*.F90)
SRC_FTN = $(SRC_FDPS_MOD) \
	  $(SRC_USER_DEFINED_TYPE) \
	  FDPS_module.F90 \
	  $(SRC_USER)
SRC_CXX = FDPS_ftn_if.cpp \
	  FDPS_Manipulators.cpp \
	  main.cpp

OBJ_USER_DEFINED_TYPE	= $(SRC_USER_DEFINED_TYPE:F90=o)
OBJ_USER		= $(SRC_USER:F90=o)
OBJ_FDPS_MOD		= $(notdir $(SRC_FDPS_MOD:F90=o))
OBJ_FTN			= $(notdir $(SRC_FTN:F90=o))
OBJ_CXX			= $(SRC_CXX:cpp=o)
OBJ			= $(OBJ_FTN) $(OBJ_CXX)

VPATH = $(FDPS_FTN_MOD_DIR)
TARGET = p3m.out

$(TARGET): $(OBJ) $(PM_ROOT)/libpm.a
	$(CXX) $(OBJ) $(CXXFLAGS) -o $(TARGET) $(LDFLAGS)

$(PM_ROOT)/libpm.a: $(PM_ROOT)/param_fdps.h
	cd $(PM_ROOT) && $(MAKE) allclean libpm.a

$(SRC_CXX) FDPS_module.F90: $(SRC_USER_DEFINED_TYPE)
	$(FDPS_FTN_IF_GENERATOR) $(SRC_USER_DEFINED_TYPE) --output ./ 

FDPS_super_particle.o: FDPS_vector.o FDPS_matrix.o

$(OBJ_USER_DEFINED_TYPE): $(OBJ_FDPS_MOD)

FDPS_module.o: $(OBJ_USER_DEFINED_TYPE)

$(OBJ_USER): $(OBJ_USER_DEFINED_TYPE) FDPS_module.o

clean:
	rm -f *.o *.s *.mod $(TARGET) *.dat

distclean: clean
	rm -f $(SRC_CXX) FDPS_Manipulators.h  FDPS_module.F90 user_defined.hpp
	cd $(PM_ROOT) && $(MAKE) allclean

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
