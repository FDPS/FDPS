CC=g++
CCFLAGS= -O2 

PIKG_ROOT=../../../
PIKG=$(PIKG_ROOT)/bin/pikg
CONVERSION_TYPE=reference
CCFLAGS += -I$(PIKG_ROOT)/inc

PIKG_OPTION= --conversion-type $(CONVERSION_TYPE) --epi-name Particle --epj-name Particle --force-name Particle

# AVX2 mode
#CONVERSION_TYPE=AVX2 
#CCFLAGS += -mavx2 -mfma -ffast-math

# AVX-512 mode
#CONVERSION_TYPE=AVX-512
#CCFLAGS += -mavx512f -mavx512dq -ffast-math

# ARM SVE mode
#CONVERSION_TYPE=A64FX
#CCFLAGS += -Kfast

all:	nbody.out
kernel.hpp:	kernel.pikg
	$(PIKG) $(PIKG_OPTION) --output $@ -i $<
nbody.out:	nbody.cpp kernel.hpp Makefile
	$(CC) $(CCFLAGS) nbody.cpp -o $@
test:	nbody.out
	./nbody.out
clean:
	rm -r result nbody.out kernel.hpp
