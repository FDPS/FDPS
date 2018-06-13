CC := mpicc
CXX := mpicxx
CFLAGS = -O3 
# Note that in some environment, the options -ffast-math funroll-loops
# results in an irreproducible error.
CFLAGS += -Wall
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
CFLAGS += -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
# fdps-autotest-set-vars-1st
# (DO NOT CHANGE THE LINE ABOVE AND ITS POSITION MUST BE
#  JUST BEFORE THE LINE BELOW)
CXXFLAGS := $(CFLAGS)
CXXFLAGS += -std=c++11

all: treepm

TREEPM_OBJ = treepm.o make_directory.o timing.o
TREEPM_DEP := $(TREEPM_OBJ) treepm.hpp cosmology.hpp run_param.hpp

FDPS_ROOT=../../../
PM_ROOT=$(FDPS_ROOT)/src/particle_mesh

CXXFLAGS += -I$(FDPS_ROOT)/src
CXXFLAGS += -I$(PM_ROOT)

FFTW_ROOT = /work4a/local/fftw
FFTW_LIBS = -L$(FFTW_ROOT)/lib -lfftw3f_mpi -lfftw3f

CXXFLAGS += -I$(FFTW_ROOT)/include
CFLAGS += -I$(FFTW_ROOT)/include
CFLAGS += -I$(FDPS_ROOT)/src
CLIBS +=  -L$(PM_ROOT) -lpm $(FFTW_LIBS)


$(PM_ROOT)/libpm.a: $(PM_ROOT)/param_fdps.h
	cd $(PM_ROOT) && $(MAKE) allclean libpm.a


EXTERNAL_LIB = $(PM_ROOT)/param_fdps.h $(PM_ROOT)/libpm.a

use_phantom_grape_x86 = no
# fdps-autotest-set-vars-2nd
# (DO NOT CHANGE THE LINE ABOVE AND ITS POSITION MUST BE
#  AFTER THE DEFINITION OF use_phantom_grape_x86)

ifeq ($(use_phantom_grape_x86),yes)

PHANTOM_ROOT = $(FDPS_ROOT)/src/phantom_grape_x86/G5/table/
CXXFLAGS  += -DENABLE_PHANTOM_GRAPE_X86
CXXFLAGS += -I$(PHANTOM_ROOT)

PG_CC = $(CXX)
PG_CFLAGS = -DSFT_FOR_PP -DSFT_FOR_PM
PG_CFLAGS += -O3
#PG_CFLAGS += -ffast-math -funroll-loops -Wall

gravity_kernel_avx2.o:$(PHANTOM_ROOT)/gravity_kernel_avx2.c
	$(PG_CC) $(PG_CFLAGS) -c -o $@ $^
phantom_g5.o:$(PHANTOM_ROOT)/phantom_g5.c
	$(PG_CC) $(PG_CFLAGS) -c -o $@ $^
pg5_table.o:$(PHANTOM_ROOT)/pg5_table.c
	$(PG_CC) $(PG_CFLAGS) -c -o $@ $^
libpg5.a:gravity_kernel_avx2.o phantom_g5.o pg5_table.o
	ar -r $@ $^
	ranlib $@
CLIBS   += -L./ -lpg5
EXTERNAL_LIB += libpg5.a

endif

treepm.o: treepm.hpp cosmology.hpp run_param.hpp

treepm:  $(TREEPM_DEP) $(EXTERNAL_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $(TREEPM_OBJ) $(CLIBS) -lm 

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $*.cpp -o $*.o

clean: 
	-rm -rf treepm *.o *.a
	cd $(PM_ROOT) && $(MAKE) allclean

distclean: clean
	-rm -rf *~

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
