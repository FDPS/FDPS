#-----------------------------------------------------------------------
#   Configuration of compile
#-----------------------------------------------------------------------
# Variables that must be specified by users
# (i) Variables related to FDPS
FDPS_LOC = ../../../
FDPS_INC = -I$(FDPS_LOC)/src 
FDPS_FTN_MOD_DIR = $(FDPS_LOC)/src/fortran_interface/modules
FDPS_FTN_IF_GENERATOR = $(FDPS_LOC)/scripts/gen_ftn_if.py
#FDPS_FTN_MOD_DIR = $(FDPS_LOC)/src/fortran_interface-v1/modules
#FDPS_FTN_IF_GENERATOR = $(FDPS_LOC)/scripts/gen_ftn_if-v1.py

# (ii) Variables to specify compilers and compile options
# Serial or OpenMP cases
FC=gfortran
CXX=g++
# MPI case
#FC=mpif90
#CXX=mpicxx
# [Option 1] w/o optimization
#FCFLAGS = -std=f2003 -O0 -Wall
#CXXFLAGS = -std=c++17 -Wall -Wextra -ftrapv -fexceptions -g3 $(FDPS_INC)
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

#use_phantom_grape_x86 = yes
#use_pikg_x86 = yes

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

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

ifeq ($(use_pikg_x86),yes)
PIKG_ROOT = $(FDPS_LOC)/pikg
PIKG_FTN_MOD_DIR = $(PIKG_ROOT)/src/fortran_interface/modules
PIKG = $(PIKG_ROOT)/bin/pikg
FCFLAGS += -DENABLE_PIKG_KERNEL_X86
CXXFLAGS += -DPIKG_USE_FDPS_VECTOR -I$(PIKG_ROOT)/inc

# Set conversion type
# (i) Reference mode
CONVERSION_TYPE = reference
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
SRC_CXX_PIKG = pikg_kernel_ep_ep.cpp
SRC_FTN_PIKG = pikg_module_ep_ep.F90
else
PIKG_FTN_MOD_DIR =
PIKG =
CONVERSION_TYPE =
EPI =
EPJ =
FORCE =
SRC_CXX_PIKG = 
SRC_FTN_PIKG = 
endif

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
SRC_PIKG_MOD = $(wildcard $(PIKG_FTN_MOD_DIR)/*.F90)
SRC_FTN_ALL = $(SRC_FDPS_MOD) \
	      $(SRC_PIKG_MOD) \
	      $(SRC_USER_DEFINED_TYPE) \
	      FDPS_module.F90 \
	      $(SRC_FTN_PIKG) \
	      $(SRC_USER)
SRC_CXX_FDPS = FDPS_ftn_if.cpp \
	       FDPS_Manipulators.cpp \
	       main.cpp
SRC_CXX_ALL = $(SRC_CXX_FDPS) $(SRC_CXX_PIKG)

OBJ_USER_DEFINED_TYPE	= $(SRC_USER_DEFINED_TYPE:F90=o)
OBJ_USER		= $(SRC_USER:F90=o)
OBJ_FDPS_MOD		= $(notdir $(SRC_FDPS_MOD:F90=o))
OBJ_PIKG_MOD		= $(notdir $(SRC_PIKG_MOD:F90=o))
OBJ_FTN_PIKG		= $(SRC_FTN_PIKG:F90=o)
OBJ_FTN_ALL		= $(notdir $(SRC_FTN_ALL:F90=o))
OBJ_CXX_FDPS            = $(SRC_CXX_FDPS:cpp=o)
OBJ_CXX_PIKG		= $(SRC_CXX_PIKG:cpp=o)
OBJ_CXX_ALL		= $(SRC_CXX_ALL:cpp=o)
OBJ_PIKG		= $(OBJ_CXX_PIKG) $(OBJ_FTN_PIKG)
OBJ			= $(OBJ_FTN_ALL) $(OBJ_CXX_ALL) $(OBJ_PG)

VPATH = $(FDPS_FTN_MOD_DIR) $(PIKG_FTN_MOD_DIR) $(PG_ROOT)
TARGET = nbody.out

$(TARGET): $(OBJ) result
	$(PG_BUILD)
	$(CXX) $(OBJ) $(CXXFLAGS) -o $(TARGET) $(LDFLAGS)

result:
	mkdir -p result

$(SRC_CXX_FDPS) FDPS_module.F90: $(SRC_USER_DEFINED_TYPE)
	$(FDPS_FTN_IF_GENERATOR) $(SRC_USER_DEFINED_TYPE) --output ./ 

FDPS_super_particle.o: FDPS_vector.o FDPS_matrix.o

$(OBJ_USER_DEFINED_TYPE): $(OBJ_FDPS_MOD) $(OBJ_PIKG_MOD) $(OBJ_PIKG) $(OBJ_PG)

FDPS_module.o: $(OBJ_USER_DEFINED_TYPE)

$(OBJ_USER): $(OBJ_USER_DEFINED_TYPE) FDPS_module.o $(OBJ_PIKG)

$(OBJ_PIKG): $(SRC_CXX_PIKG) $(SRC_CXX_FDPS)

$(SRC_CXX_PIKG) $(SRC_FTN_PIKG): kernel_ep_ep.pikg
	$(PIKG) --conversion-type $(CONVERSION_TYPE) --epi-name $(EPI) --epj-name $(EPJ) --force-name $(FORCE) --kernel-name "pikg_calc_grav_ep_ep" --fortran-interface $(basename $(SRC_FTN_PIKG)) --output $(SRC_CXX_PIKG) -i $<

clean:
	rm -f *.o *.s *.mod $(TARGET) *.dat

distclean: clean
	$(PG_CLEAN)
	rm -f $(SRC_CXX_FDPS) FDPS_Manipulators.h  FDPS_module.F90 user_defined.hpp 
	rm -f $(SRC_CXX_PIKG) $(SRC_FTN_PIKG)
	rm -rf result

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
