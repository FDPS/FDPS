SRC_FDPS_DIR = ../../../src/
INC = -I$(SRC_FDPS_DIR)
CXX = mpicxx
#CXX = g++
CXXFLAGS = -O3 -std=c++17
CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CXXFLAGS += -Wall

SRC = main.cpp

use_pikg_x86 = yes

ifeq ($(use_pikg_x86),yes)
PIKG_ROOT = ../../../pikg
PIKG = $(PIKG_ROOT)/bin/pikg
INC += -I$(PIKG_ROOT)/inc
CXXFLAGS += -DUSE_PIKG_KERNEL -DPIKG_USE_FDPS_VECTOR
#CXXFLAGS += -DKERNEL_64BIT

# reference option
#CONVERSION_TYPE=reference

# AVX2 options
#CONVERSION_TYPE = AVX2
#CXXFLAGS+= -mavx2 -mfma -ffast-math

# AVX-512 options
CONVERSION_TYPE = AVX-512
CXXFLAGS+= -mavx512f -mavx512dq -ffast-math

PIKG_FLAGS= --conversion-type $(CONVERSION_TYPE)
SRC += user_defined_kernel.hpp kernel_ep.hpp kernel_sp.hpp kernel_sp_quad.hpp
endif

ring_mono.out: $(SRC)
	$(CXX) $(INC) $(CXXFLAGS)  $< -o $@

ring_quad.out: $(SRC)
	$(CXX) $(INC) $(CXXFLAGS)  -DQUAD $< -o $@

user_defined_kernel.hpp: kernel_ep.hpp kernel_sp.hpp kernel_sp_quad.hpp

kernel_ep.hpp: kernel_ep.pikg
	$(PIKG) $(PIKG_FLAGS) --epi-name Epi0 --epj-name Epj0 --force-name Force0 --class-file user_defined_kernel.hpp --kernel-name CalcForceEpEpImpl -i $< -o $@

kernel_sp.hpp: kernel_sp.pikg
	$(PIKG) $(PIKG_FLAGS) --epi-name Epi1 --epj-name Epj1 --force-name Force1 --class-file user_defined_kernel.hpp --kernel-name CalcForceEpSpImpl -i $< -o $@

kernel_sp_quad.hpp: kernel_sp_quad.pikg
	$(PIKG) $(PIKG_FLAGS) --epi-name Epi2 --epj-name Epj2 --force-name Force2 --class-file user_defined_kernel.hpp --kernel-name CalcForceEpSpQuadImpl -i $< -o $@

clean:
	rm *.out kernel_ep.hpp kernel_sp.hpp
