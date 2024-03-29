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
CXXFLAGS = -std=c++17 -O3 -ffast-math -funroll-loops $(FDPS_INC)
LDFLAGS = 
# OpenMP options
#CFLAGS   += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
# MPI options
#CFLAGS   += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

#use_pikg_x86=yes

#use_pikg_x86=yes
PIKG_ROOT = $(FDPS_LOC)/pikg/
PIKG = $(PIKG_ROOT)/bin/pikg

ifeq ($(use_pikg_x86),yes)
PIKG_ROOT = $(FDPS_LOC)/pikg
PIKG_HDR_DIR = $(PIKG_ROOT)/inc
PIKG = $(PIKG_ROOT)/bin/pikg
CFLAGS += -DUSE_PIKG_KERNEL -I$(PIKG_HDR_DIR)
CXXFLAGS += -DUSE_PIKG_KERNEL -DPIKG_USE_FDPS_VECTOR -I$(PIKG_HDR_DIR)

# Set convertion type
# (i) Reference mode
CONVERSION_TYPE=reference

# (ii) AVX2 mode
#CONVERSION_TYPE = AVX2
#CXXFLAGS += -mavx2 -mfma -ffast-math
# (iii) AVX-512 mode
#CONVERSION_TYPE = AVX-512
#CXXFLAGS += -mavx512f -mavx512dq -ffast-math
# (iv) ARM SVE mode
#CONVERSION_TYPE = A64FX
#CXXFLAGS += -Kfast

EPI = epi_grav
EPJ = epj_grav
FORCE = force_grav
SRC_CXX_PIKG = kernel_epep.cpp
HDR_C_PIKG = kernel_epep.h 

else
PIKG_HDR_DIR =
PIKG =
CONVERSION_TYPE =
EPI =
EPJ =
FORCE =
SRC_CXX_PIKG =
HDR_C_PIKG =
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
SRC_CXX_FDPS = FDPS_ftn_if.cpp \
	       FDPS_Manipulators.cpp \
	       main.cpp
SRC_CXX_ALL = $(SRC_CXX_FDPS) $(SRC_CXX_PIKG)

OBJ_USER	= $(SRC_USER:c=o)
OBJ_CXX_FDPS	= $(SRC_CXX_FDPS:cpp=o)
OBJ_CXX_PIKG	= $(SRC_CXX_PIKG:cpp=o)
OBJ_CXX_ALL	= $(OBJ_CXX_FDPS) $(OBJ_CXX_PIKG)
OBJ		= $(OBJ_USER) $(OBJ_CXX_ALL)


TARGET = nbody.out

$(TARGET): $(OBJ) result
	$(CXX) $(OBJ) $(CXXFLAGS) -o $(TARGET) $(LDFLAGS)

result:
	mkdir -p result


$(OBJ_USER): FDPS_c_if.h  $(HDR_C_PIKG)

FDPS_c_if.h $(SRC_CXX_FDPS): $(HDR_USER_DEFINED_TYPE) Makefile
	$(FDPS_C_IF_GENERATOR) $(HDR_USER_DEFINED_TYPE) --output ./

$(SRC_CXX_PIKG) $(HDR_C_PIKG): kernel_epep.pikg
	$(PIKG) --conversion-type $(CONVERSION_TYPE) --epi-name $(EPI) --epj-name $(EPJ) --force-name $(FORCE) --kernel-name "pikg_calc_grav_ep_ep" --c-interface $(HDR_C_PIKG) -o $(SRC_CXX_PIKG) -i $<

clean:
	rm -f *.o *.s $(TARGET) *.dat

distclean: clean
	rm -f $(SRC_CXX_FDPS) FDPS_c_if.h FDPS_Manipulators.h user_defined.hpp
	rm -f $(SRC_CXX_PIKG) $(HDR_C_PIKG)
	rm -rf result

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
