CC=gcc
CFLAGS= -O2 -I$(PIKG_ROOT)/inc

CXX=g++
CXXFLAGS= -O2 -I$(PIKG_ROOT)/inc

PIKG_ROOT=../../../
PIKG=$(PIKG_ROOT)/bin/pikg
CONVERSION_TYPE=reference
CFLAGS += 

PIKG_OPTION= --conversion-type $(CONVERSION_TYPE) --epi-name particle --epj-name particle --force-name particle --c-interface --kernel-name calc_gravity --class-file user_defined.h

# AVX2 mode
#CONVERSION_TYPE=AVX2
#CXXFLAGS += -mavx2 -mfma -ffast-math

# AVX-512 mode
#CONVERSION_TYPE=AVX-512
#CXXFLAGS += -mavx512f -mavx512dq -ffast-math

# ARM SVE mode
#CONVERSION_TYPE=A64FX
#CXXFLAGS += -Kfast

SRC_C=nbody.c
SRC_CXX=kernel.cpp
OBJ=$(SRC_CXX:.cpp=.o) $(SRC_C:.c=.o)

.SECONDARY:	$(SRC_CXX)
%.o : %.c
	$(CC) $(CFLAGS) -c $<
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<
%.cpp : %.pikg
	$(PIKG) $(PIKG_OPTION) --output $@ -i $<

all:	nbody.out

nbody.out:	$(OBJ) Makefile
	$(CXX) $(CCFLAGS) $(OBJ) -o $@
clean:
	rm -r nbody.out kernel.cpp kernel.h *.o
