#-----------------------------------------------------------------------
#   Configuration of compile
#-----------------------------------------------------------------------
# Variables that must be specified by users
# (i) Variables related to FDPS
FDPS_LOC = ../../../
PM_ROOT  = $(FDPS_LOC)/src/particle_mesh
FDPS_INC = -I$(FDPS_LOC)/src -I$(PM_ROOT)
FDPS_LIB = -L$(PM_ROOT) -lpm
FDPS_INC += -I$(FDPS_LOC)/src/c_interface/headers
FDPS_C_IF_GENERATOR = $(FDPS_LOC)/scripts/gen_c_if.py

# (ii) Variables related to FFTW
FFTW_LOC = /work4/local/fftw
FFTW_INC = -I$(FFTW_LOC)/include
FFTW_LIB = -L$(FFTW_LOC)/lib -lfftw3f_mpi -lfftw3f_omp -lfftw3f -lm

# (iii) Variables to specify compilers and compile options
CC=mpicc
CXX=mpic++
# [Option 1] w/o optimization
#CFLAGS   = -Wall -Wextra -ftrapv -fexceptions -g3 $(FDPS_INC)
#CXXFLAGS = -Wall -Wextra -ftrapv -fexceptions -g3 $(FDPS_INC)
# [Option 2] w/ optimization 
CFLAGS   = -O3 -ffast-math -funroll-loops -finline-functions $(FDPS_INC)
CXXFLAGS = -O3 -ffast-math -funroll-loops $(FFTW_INC) $(FDPS_INC)
LDFLAGS  = $(FFTW_LIB) $(FDPS_LIB)
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
%.o : %.c
	$(CC) $(CFLAGS) -c $<
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<

HDR_USER_DEFINED_TYPE = user_defined.h
SRC_USER = user_defined.c \
	   c_main.c
SRC_CXX = FDPS_ftn_if.cpp \
	  FDPS_Manipulators.cpp \
	  main.cpp

OBJ_USER = $(SRC_USER:c=o)
OBJ_CXX	 = $(SRC_CXX:cpp=o)
OBJ	 = $(OBJ_USER) $(OBJ_CXX)

TARGET = p3m.out

$(TARGET): $(OBJ) $(PM_ROOT)/libpm.a
	$(CXX) $(OBJ) $(CXXFLAGS) -o $(TARGET) $(LDFLAGS)

$(PM_ROOT)/libpm.a: $(PM_ROOT)/param_fdps.h
	cd $(PM_ROOT) && $(MAKE) allclean libpm.a

FDPS_c_if.h $(SRC_CXX): $(HDR_USER_DEFINED_TYPE) Makefile
	$(FDPS_C_IF_GENERATOR) user_defined.h --output ./

$(OBJ_USER): FDPS_c_if.h 

clean:
	rm -f *.o *.s $(TARGET) *.dat

distclean: clean
	rm -f $(SRC_CXX) FDPS_c_if.h FDPS_Manipulators.h user_defined.hpp
	cd $(PM_ROOT) && $(MAKE) allclean

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
