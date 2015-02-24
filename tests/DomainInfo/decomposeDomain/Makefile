#PS_PATH  = -I../../../src_parallel

#CCC = mpicxx-openmpi-gcc49

#CFLAGS = -O3 -ffast-math -funroll-loops -DMPICH_IGNORE_CXX_SEEK #-Wall
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

ALL = testf64 testf32
EXE = mainf64 mainf32
HDR = ../../basic_particle.hpp ../check_domain_info.hpp
ETC = Makefile

all: $(ALL)

testf64: mainf64
	@error=`mpirun-openmpi-gcc49 -np 8 ./mainf64 >& /dev/null; echo $$?`; \
	if test $$error -eq 0 ; then \
			echo "SUCCESS: DomainInfo::decomposeDomain: F64"; exit $$error; \
	elif test $$error -eq 1 ; then \
			echo "ERROR[1]: DomainInfo::decomposeDomain: F64"; exit $$error; \
	fi \

testf32: mainf32
	@error=`mpirun-openmpi-gcc49 -np 8 ./mainf32 >& /dev/null; echo $$?`; \
	if test $$error -eq 0 ; then \
			echo "SUCCESS: DomainInfo::decomposeDomain: F32"; exit $$error; \
	elif test $$error -eq 1 ; then \
			echo "ERROR[1]: DomainInfo::decomposeDomain: F32"; exit $$error; \
	fi \

mainf64: mainf64.cpp $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -o $@ mainf64.cpp $(CLIBS)

mainf32: mainf32.cpp $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -o $@ mainf32.cpp $(CLIBS)

clean:
	rm -f *~

distclean:
	make clean
	rm -f *.o

allclean:
	make distclean
	rm -f $(EXE)
