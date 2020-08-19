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
#CXXFLAGS = -std=c++11 -Wall -Wextra -ftrapv -fexceptions -g3 $(FDPS_INC)
# [Option 2] w/ optimization 
CFLAGS = -O3 -ffast-math -funroll-loops $(FDPS_INC)
CXXFLAGS = -std=c++11 -O3 -ffast-math -funroll-loops $(FDPS_INC)
LDFLAGS = 
# OpenMP options
#CFLAGS   += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
# MPI options
#CFLAGS   += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL


#use_pikg_x86=yes
PIKG_ROOT = $(FDPS_LOC)/pikg/
PIKG = $(PIKG_ROOT)/bin/pikg
ifeq ($(use_pikg_x86),yes)
INC += -I$(PIKG_ROOT)/inc
CFLAGS+= -DUSE_PIKG_KERNEL -I$(PIKG_ROOT)/inc
CXXFLAGS+= -DUSE_PIKG_KERNEL -DPIKG_USE_FDPS_VECTOR -I$(PIKG_ROOT)/inc

# reference option
CONVERSION_TYPE=reference

# AVX2 options
#CONVERSION_TYPE = AVX2
#CXXFLAGS+= -mavx2 -mfma -ffast-math

# AVX-512 options
#CONVERSION_TYPE = AVX-512
#CXXFLAGS+= -mavx512f -mavx512dq -ffast-math

PIKG_FLAGS= --conversion-type $(CONVERSION_TYPE)
OBJ_PIKG = kernel_epep.o 
PIKG_BUILD = make $(OBJ_PIKG)
else
PIKG_BUILD =
endif

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
TARGET = nbody.out

$(TARGET): $(OBJ_PIKG) $(OBJ) result
	$(PIKG_BUILD)
	$(CXX) $(OBJ) $(OBJ_PIKG) $(CXXFLAGS) -o $(TARGET) $(LDFLAGS)

result:
	mkdir -p result

$(OBJ_PIKG): $(SRC_CXX)

FDPS_c_if.h $(SRC_CXX): $(HDR_USER_DEFINED_TYPE) Makefile
	$(FDPS_C_IF_GENERATOR) user_defined.h --output ./

kernel_epep.cpp:	kernel_epep.pikg
	$(PIKG) $(PIKG_FLAGS) -i kernel_epep.pikg --c-interface kernel_epep.h --epi-name epi_grav --epj-name epj_grav --force-name force_grav --kernel-name pikg_calc_grav_ep_ep -o $@

$(OBJ_USER): FDPS_c_if.h 

clean:
	rm -f *.o *.s $(TARGET) *.dat

distclean: clean
	rm -f $(SRC_CXX) FDPS_c_if.h FDPS_Manipulators.h user_defined.hpp kernel_pikg.cpp kernel_epep.h kernel_epsp.h kernel_epep.cpp kernel_epsp.cpp
	rm -rf result

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
