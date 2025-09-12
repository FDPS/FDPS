CXX := mpicxx
CXXFLAGS = -std=c++17 -O3 
# Note that in some environment, the options -ffast-math funroll-loops
# results in an irreproducible error.
CXXFLAGS += -Wall
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
CXXFLAGS += -fopenmp
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
# fdps-autotest-set-vars-1st
# (DO NOT CHANGE THE LINE ABOVE AND ITS POSITION MUST BE
#  JUST BEFORE THE LINE BELOW)
#  CXXFLAGS += -DPARTICLE_SIMULATOR_USE_64BIT_KEY


all: pmmm.out

FDPS_ROOT=../../../../
PMM_ROOT=$(FDPS_ROOT)/src/particle_mesh_multipole

CXXFLAGS += -I$(FDPS_ROOT)/src
CXXFLAGS += -I$(PMM_ROOT)

FFTW_ROOT = /usr/lib/x86_64-linux-gnu
FFTW_LIBS = -L$(FFTW_ROOT)/lib -lfftw3_mpi -lfftw3_omp -lfftw3
CXXFLAGS += -I$(FFTW_ROOT)/include
CLIBS +=  $(FFTW_LIBS)

DEP_FILES = *.cpp
DEP_FILES += $(filter-out kernel_pikg.hpp, $(wildcard *.hpp))

#use_pikg = yes
ifeq ($(use_pikg),yes)
PIKG_ROOT = $(FDPS_ROOT)/pikg
PIKG = $(PIKG_ROOT)/bin/pikg
CXXFLAGS += -DENABLE_PIKG
CXXFLAGS += -DPIKG_USE_FDPS_VECTOR
CXXFLAGS += -I$(PIKG_ROOT)/inc
DEP_FILES += kernel_pikg.hpp
CONVERSION_TYPE = AVX2
CXXFLAGS += -mavx2 -mfma -ffast-math
PIKG_FLAGS = --conversion-type $(CONVERSION_TYPE)
endif

pmmm.out: pmmm.cpp $(DEP_FILES)
	$(CXX) $(CXXFLAGS) -o $@ pmmm.cpp $(CLIBS) -lm

kernel_pikg.hpp: kernel_epep.pikg
	$(PIKG) $(PIKG_FLAGS) --epi-name Epi --epj-name Epj --force-name Force --kernel-name CalcGravityEpEp -i kernel_epep.pikg -o $@

clean: 
	-rm -rf pmmm.out kernel_pikg.hpp *.o *.a

distclean: clean
	-rm -rf *~ 

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
