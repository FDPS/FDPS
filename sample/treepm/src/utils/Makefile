CC=gcc
CFLAGS=-std=c99 -O3 -D__FFTW2__ 

FFTW_ROOT=/usr/local/fftw-2
CFLAGS += -I$(FFTW_ROOT)/include
CLIBS = -L$(FFTW_ROOT)/lib -lsrfftw -lsfftw

POWER_SPEC_OBJ = power_spec.o input_data.o cosmology.o calc_power_spec.o
POWER_SPEC_DEP := $(POWER_SPEC_OBJ)

CONCAT_DATA_OBJ = concat_data.o input_data.o output_data.o
CONCAT_DATA_DEP = $(CONCAT_DATA_OBJ)

all: power_spec concat_data

power_spec: $(POWER_SPEC_DEP)
	$(CC) -o $@ $(POWER_SPEC_OBJ) $(CLIBS) -lm

concat_data: $(CONCAT_DATA_DEP)
	$(CC) -o $@ $(CONCAT_DATA_OBJ) $(CLIBS) -lm

clean:
	\rm -rf *.o

distclean: clean
	\rm -rf *~ power_spec

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $*.cpp -o $*.o
