#-----------------------------------------------------------------------
#   Configuration of compile
#-----------------------------------------------------------------------
# Variables that must be specified by users
# (i) Variables related to FDPS
FDPS_LOC = ../../../
FDPS_INC = -I$(FDPS_LOC)/src 
FDPS_INC += -I$(FDPS_LOC)/src/c_interface/headers
FDPS_C_IF_GENERATOR = $(FDPS_LOC)/scripts/gen_c_if.py

# (ii) Variables to specify compilers and compile options
# Serial or OpenMP cases
CC=gcc
CXX=g++
# MPI case
#CC=mpicc
#CXX=mpic++
# [Option 1] w/o optimization
#CFLAGS   = -Wall -Wextra -ftrapv -fexceptions -g3 $(FDPS_INC)
#CXXFLAGS = -Wall -Wextra -ftrapv -fexceptions -g3 $(FDPS_INC)
# [Option 2] w/ optimization 
CFLAGS = -O3 -ffast-math -funroll-loops $(FDPS_INC)
CXXFLAGS = -O3 -ffast-math -funroll-loops $(FDPS_INC)
LDFLAGS = 
# OpenMP options
#CFLAGS   += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
# MPI options
#CFLAGS   += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

#use_phantom_grape_x86 = yes

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

# Simulation control macros 
CFLAGS += -DENABLE_VARIABLE_SMOOTHING_LENGTH
CFLAGS += -DUSE_ENTROPY
CFLAGS += -DUSE_BALSARA_SWITCH
#CFLAGS += -DUSE_PRESCR_OF_THOMAS_COUCHMAN_1992
CFLAGS += -DISOTHERMAL_EOS

CFLAGS += -DINITIAL_CONDITION=0 
#CFLAGS += -DINITIAL_CONDITION=1 
#CFLAGS += -DINITIAL_CONDITION=2 
#CFLAGS += -DINITIAL_CONDITION=3 

ifeq ($(use_phantom_grape_x86),yes)
PG_ROOT = $(FDPS_LOC)/src//phantom_grape_x86/G5/newton/libpg5
CFLAGS  += -DENABLE_PHANTOM_GRAPE_X86
LDFLAGS  += -L$(PG_ROOT) -lpg5
PG_BUILD = cd $(PG_ROOT) && $(MAKE) distclean libpg5.a
PG_CLEAN = cd $(PG_ROOT) && $(MAKE) distclean
else
PG_ROOT  = 
PG_BUILD =
PG_CLEAN = 
endif

#-----------------------------------------------------------------------
#   Source files
#-----------------------------------------------------------------------
%.o : %.c
	$(CC) $(CFLAGS) -c $<
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<

HDR_USER_DEFINED_TYPE = user_defined.h
SRC_USER = user_defined.c \
	   mathematical_constants.c \
	   physical_constants.c \
	   leapfrog.c \
	   ic.c \
	   c_main.c 
SRC_CXX = FDPS_ftn_if.cpp \
	  FDPS_Manipulators.cpp \
	  main.cpp

OBJ_USER = $(SRC_USER:c=o)
OBJ_CXX	 = $(SRC_CXX:cpp=o)
OBJ	 = $(OBJ_USER) $(OBJ_CXX)

TARGET = nbodysph.out

$(TARGET): $(OBJ) result
	$(PG_BUILD)
	$(CXX) $(OBJ) $(CXXFLAGS) -o $(TARGET) $(LDFLAGS)

result:
	mkdir -p result

FDPS_c_if.h $(SRC_CXX): $(HDR_USER_DEFINED_TYPE) Makefile
	$(FDPS_C_IF_GENERATOR) user_defined.h --output ./

$(OBJ_USER): FDPS_c_if.h 

clean:
	rm -f *.o *.s $(TARGET) *.dat

distclean: clean
	$(PG_CLEAN)
	rm -f $(SRC_CXX) FDPS_c_if.h FDPS_Manipulators.h user_defined.hpp 
	rm -rf result

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
