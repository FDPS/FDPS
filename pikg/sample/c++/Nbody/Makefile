CC=g++
CCFLAGS= -O2 

PIKG_ROOT=../../../
PIKG=$(PIKG_ROOT)/bin/pikg
CONVERSION_TYPE=reference
CCFLAGS += -I$(PIKG_ROOT)/inc

PIKG_OPTION= --conversion-type $(CONVERSION_TYPE) --epi-name Particle --epj-name Particle --force-name Particle --class-file particle.hpp

SRC=nbody.cpp

ifeq ($(use_avx2),yes)
CONVERSION_TYPE=AVX2
CCFLAGS += -mavx2 -mfma -ffast-math
HED = kernel.hpp
endif

ifeq ($(use_avx512),yes)
CONVERSION_TYPE=AVX-512
CCFLAGS += -mavx512f -mavx512dq -ffast-math
HED = kernel.hpp
endif

ifeq ($(use_arm_sve),yes)
CONVERSION_TYPE=A64FX
CCFLAGS += -Kfast
HED = kernel.hpp
endif

ifeq ($(use_cuda),yes)
CONVERSION_TYPE=CUDA
SRC+= kernel.o
CC=nvcc
NVCC=nvcc
CCFLAGS+= -DUSE_CUDA_KERNEL
NVCCFLAGS= $(CCFLAGS)
endif

ifeq ($(CONVERSION_TYPE),reference)
HED = kernel.hpp
endif

all:	nbody.out
kernel.cu:	kernel.pikg
	$(PIKG) $(PIKG_OPTION) --output $@ -i $<
kernel.o:	kernel.cu
	$(NVCC) $(NVCCFLAGS) $< -c -o $@
kernel.hpp:	kernel.pikg
	$(PIKG) $(PIKG_OPTION) --output $@ -i $<
nbody.out:	$(SRC) $(HED)
	$(CC) $(CCFLAGS) $(SRC) -o $@
test:	nbody.out
	./nbody.out
clean:
	rm -r result nbody.out kernel.hpp kernel.cu kernel.o
