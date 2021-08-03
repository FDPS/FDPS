CC=g++
CCFLAGS= -O2 -std=c++11

PIKG_ROOT=../../
PIKG=$(PIKG_ROOT)/bin/pikg
CONVERSION_TYPE=reference

CCFLAGS += -I$(PIKG_ROOT)/inc

PIKG_OPTION= --epi-name Ep --epj-name Ep --force-name Force --class-file main.cpp

all:	main.out
kernel.hpp:	kernel.pikg
	$(PIKG) $(PIKG_OPTION) --output $@ -i $<
main.out:	kernel.hpp
	$(CC) $(CCFLAGS) main.cpp -o $@
ref:
	$(PIKG) $(PIKG_OPTION) --conversion-type reference --output kernel.hpp -i kernel.pikg
	$(CC) $(CCFLAGS) main.cpp -o main.out
	./main.out

avx2:
	$(PIKG) $(PIKG_OPTION) --conversion-type AVX2 --output kernel.hpp -i kernel.pikg
	$(CC) $(CCFLAGS) -mavx2 -mfma -DSIMD_KERNEL main.cpp -o main.out
	./main.out
avx512:
	$(PIKG) $(PIKG_OPTION) --conversion-type AVX-512 --output kernel.hpp -i kernel.pikg
	$(CC) $(CCFLAGS) -mavx512f -mavx512dq -DSIMD_KERNEL main.cpp -o main.out
	./main.out
clean:
		rm kernel.hpp main.out
