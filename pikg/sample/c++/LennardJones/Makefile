CC=g++
CCFLAGS= -O2 -std=c++11

PIKG_ROOT=../../../
PIKG=$(PIKG_ROOT)/bin/pikg
CONVERSION_TYPE=reference

CCFLAGS += -I$(PIKG_ROOT)/inc

PIKG_OPTION= --conversion-type $(CONVERSION_TYPE) --epi-name EP --epj-name EP --force-name Force --class-file main.cpp
# AVX2 mode
ifeq ($(use_avx2),yes)
CONVERSION_TYPE=AVX2 
CCFLAGS += -mavx2 -mfma -ffast-math
endif
# AVX-512 mode
ifeq ($(use_avx512),yes)
CONVERSION_TYPE=AVX-512
CCFLAGS += -mavx512f -mavx512dq -ffast-math
endif
# ARM SVE mode
ifeq ($(use_arm_sve),yes)
CONVERSION_TYPE=A64FX
CCFLAGS += -Kfast
endif

all:	main.out
kernel.hpp:	kernel.pikg
	$(PIKG) $(PIKG_OPTION) --output $@ -i $<
main.out:	main.cpp kernel.hpp Makefile
	$(CC) $(CCFLAGS) main.cpp -o $@
test:	main.out
	./main.out
clean:
		rm kernel.hpp main.out
