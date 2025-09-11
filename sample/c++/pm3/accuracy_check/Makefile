CC := mpicc
CXX := mpicxx
MAKE := make
FDPS_ROOT=../../../../

#PM3_ROOT=$(FDPS_ROOT)/src/particle_mesh_multipole
FFTW_ROOT = /usr/lib/x86_64-linux-gnu
CFLAGS = -std=c++17 -O3 
CFLAGS += -I$(FFTW_ROOT)/include
CFLAGS += -I$(FDPS_ROOT)/src
CFLAGS += -I$(PM3_ROOT)/
CFLAGS += -Wall
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
CFLAGS += -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

#CFLAGS += -DPS_DEBUG_PRINT_setRootCell
#CFLAGS += -DPS_DEBUG_PRINT_SetLocalEssentialTreeToGlobalTreeLong
##CFLAGS += -DPARTICLE_SIMULATOR_DEBUG_PRINT
#CFLAGS += -DPS_DEBUG_PRINT_makeGlobalTreeOnly
##CFLAGS += -DSANITY_CHECK_REALLOCATABLE_ARRAY
#CFLAGS += -DPS_DEBUG_PRINT_linkCellGlobalTreeOnly
#CFLAGS += -DPS_DEBUG_PRINT_LinkCell

FFTW_LIBS = -L$(FFTW_ROOT) -lfftw3_mpi -lfftw3_omp -lfftw3
CLIBS +=  $(FFTW_LIBS)

PROGRAM = main.out
all: ${PROGRAM}

#use_pm = yes
ifeq ($(use_pm),yes)
PM_ROOT=$(FDPS_ROOT)/src/particle_mesh
$(PM_ROOT)/libpm.a: $(PM_ROOT)/param_fdps.h
	cd $(PM_ROOT) && $(MAKE) allclean libpm.a
CFLAGS += -I$(PM_ROOT)/
CFLAGS += -DUSE_PM
CLIBS += -L$(PM_ROOT) -lpm
DEP_FILES += $(PM_ROOT)/libpm.a
endif

#use_pikg = yes
ifeq ($(use_pikg),yes)
PIKG_ROOT = $(FDPS_ROOT)/pikg
PIKG = $(PIKG_ROOT)/bin/pikg
CFLAGS += -DENABLE_PIKG
CFLAGS += -DPIKG_USE_FDPS_VECTOR
CFLAGS += -I$(PIKG_ROOT)/inc
DEP_FILES += kernel_pikg.hpp
CONVERSION_TYPE = AVX2
CFLAGS += -mavx2 -mfma -ffast-math
PIKG_FLAGS = --conversion-type $(CONVERSION_TYPE)
endif

#CFLAGS += -DPARTICLE_SIMULATOR_DEBUG_PRINT_PMMM_CALC_FORCE_ALL
#MEMCHK += -fsanitize=address -g

main.out: main.cpp  $(DEP_FILES)
	$(CXX) $(MEMCHK) $(CFLAGS) -o $@ $^ $(CLIBS)

test_multipolemoment.out: test_multipolemoment.cpp  $(DEP_FILES)
	$(CXX) -DDEV_CODE $(MEMCHK) $(CFLAGS) -o $@ $^ $(CLIBS)

kernel_pikg.hpp: kernel_epep.pikg
	$(PIKG) $(PIKG_FLAGS) --epi-name Epi --epj-name Epj --force-name Force --kernel-name CalcGravityEpEp -i kernel_epep.pikg -o $@

clean:
	rm -f *.out

