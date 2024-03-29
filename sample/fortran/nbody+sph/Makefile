#-----------------------------------------------------------------------
#   Configuration of compile
#-----------------------------------------------------------------------
# Variables that must be specified by users
# (i) Variables related to FDPS
FDPS_LOC = ../../../
FDPS_INC = -I$(FDPS_LOC)/src 
FDPS_FTN_MOD_DIR = $(FDPS_LOC)/src/fortran_interface/modules
FDPS_FTN_IF_GENERATOR = $(FDPS_LOC)/scripts/gen_ftn_if.py

# (iii) Variables to specify compilers and compile options
# Serial or OpenMP cases
FC=gfortran
CXX=g++
# MPI case
#FC=mpif90
#CXX=mpicxx
#FC=mpif90-openmpi-gcc7
#CXX=mpicxx-openmpi-gcc7
# [Option 1] w/o optimization
#FCFLAGS = -std=f2003 -O0 -Wall
#CXXFLAGS = -std=c++11 -Wall -Wextra -ftrapv -fexceptions -g3 $(FDPS_INC)
# [Option 2] w/ optimization 
FCFLAGS = -std=f2003 -O3 -ffast-math -funroll-loops -finline-functions
CXXFLAGS = -std=c++17 -O3 -ffast-math -funroll-loops $(FDPS_INC)
LDFLAGS = -lgfortran 
# OpenMP options
#FCFLAGS  += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
# MPI options
#FCFLAGS  += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

use_phantom_grape_x86 = yes

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

# Simulation control macros 
FCFLAGS += -DENABLE_VARIABLE_SMOOTHING_LENGTH
FCFLAGS += -DUSE_ENTROPY
FCFLAGS += -DUSE_BALSARA_SWITCH
#FCFLAGS += -DUSE_PRESCR_OF_THOMAS_COUCHMAN_1992
FCFLAGS += -DISOTHERMAL_EOS

FCFLAGS += -DINITIAL_CONDITION=0 
#FCFLAGS += -DINITIAL_CONDITION=1 
#FCFLAGS += -DINITIAL_CONDITION=2 
#FCFLAGS += -DINITIAL_CONDITION=3 

#CXXFLAGS += -DREAD_DATA_WITH_BYTESWAP

ifeq ($(use_phantom_grape_x86),yes)
PG_ROOT = $(FDPS_LOC)/src//phantom_grape_x86/G5/newton/libpg5
SRC_PG  = $(PG_ROOT)/phantom_grape_g5_x86.F90
OBJ_PG  = $(notdir $(SRC_PG:F90=o))
FCFLAGS  += -DENABLE_PHANTOM_GRAPE_X86
LDFLAGS  += -L$(PG_ROOT) -lpg5
PG_BUILD = cd $(PG_ROOT) && $(MAKE) distclean libpg5.a
PG_CLEAN = cd $(PG_ROOT) && $(MAKE) distclean
else
PG_ROOT  = 
SRC_PG   = 
OBJ_PG   = 
PG_BUILD =
PG_CLEAN = 
endif

#-----------------------------------------------------------------------
#   Source files
#-----------------------------------------------------------------------
%.o : %.F90
	$(FC) $(FCFLAGS) -c $<
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<

SRC_UDT_DEP = mathematical_constants.F90
SRC_USER_DEFINED_TYPE = user_defined.F90
SRC_USER_FTN = f_main.F90 \
	       physical_constants.F90 \
	       ic.F90 \
	       leapfrog.F90
SRC_USER_CXX = tipsy_file_reader.cpp
SRC_FDPS_MOD = $(wildcard $(FDPS_FTN_MOD_DIR)/*.F90)
SRC_FTN = $(SRC_FDPS_MOD) \
	  $(SRC_UDT_DEP) \
	  $(SRC_USER_DEFINED_TYPE) \
	  FDPS_module.F90 \
	  $(SRC_USER_FTN)
SRC_CXX = FDPS_ftn_if.cpp \
	  FDPS_Manipulators.cpp \
	  main.cpp \
	  $(SRC_USER_CXX)

OBJ_UDT_DEP		= $(SRC_UDT_DEP:F90=o)
OBJ_USER_DEFINED_TYPE	= $(SRC_USER_DEFINED_TYPE:F90=o)
OBJ_USER_FTN		= $(SRC_USER_FTN:F90=o)
OBJ_USER_CXX		= $(SRC_USER_CXX:cpp=o)
OBJ_FDPS_MOD		= $(notdir $(SRC_FDPS_MOD:F90=o))
OBJ_FTN			= $(notdir $(SRC_FTN:F90=o))
OBJ_CXX			= $(SRC_CXX:cpp=o) 
OBJ			= $(OBJ_FTN) $(OBJ_CXX) $(OBJ_PG)

VPATH = $(FDPS_FTN_MOD_DIR) $(PG_ROOT)
TARGET = nbodysph.out

$(TARGET): $(OBJ) result
	$(PG_BUILD)
	$(CXX) $(OBJ) $(CXXFLAGS) -o $(TARGET) $(LDFLAGS)

result:
	mkdir -p result

$(SRC_USER_DEFINED_TYPE) $(SRC_USER_FTN): macro_defs.h

$(SRC_CXX) FDPS_module.F90: $(SRC_USER_DEFINED_TYPE)
	$(FDPS_FTN_IF_GENERATOR) $(SRC_USER_DEFINED_TYPE) --output ./ 

FDPS_super_particle.o: FDPS_vector.o FDPS_matrix.o

$(OBJ_USER_DEFINED_TYPE): $(OBJ_FDPS_MOD) $(OBJ_UDT_DEP) $(OBJ_PG)

FDPS_module.o: $(OBJ_USER_DEFINED_TYPE)

$(OBJ_USER_FTN): $(OBJ_USER_DEFINED_TYPE) FDPS_module.o

ic.o: physical_constants.o tipsy_file_reader.o

clean:
	rm -f *.o *.s *.mod $(TARGET) *.dat

distclean: clean
	$(PG_CLEAN)
	rm -f FDPS_ftn_if.cpp FDPS_Manipulators.cpp main.cpp FDPS_Manipulators.h FDPS_module.F90 user_defined.hpp 
	rm -rf result

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
