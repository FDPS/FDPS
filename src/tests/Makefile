PS_PATH  = -I../../../src_parallel

CCC = mpicxx-openmpi-gcc49
CFLAGS = -O3 -ffast-math -funroll-loops -DMPICH_IGNORE_CXX_SEEK #-Wall
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

check:
	make -C DomainInfo CCC="$(CCC)" CFLAGS="$(CFLAGS)" PS_PATH="$(PS_PATH)"
	make -C ParticleSystem CCC="$(CCC)" CFLAGS="$(CFLAGS)" PS_PATH="$(PS_PATH)"

clean:
	rm -f *~
	make -C DomainInfo clean
	make -C ParticleSystem clean

distclean:
	rm -f *~
	make -C DomainInfo distclean
	make -C ParticleSystem distclean

allclean:
	rm -f *~
	make -C DomainInfo allclean
	make -C ParticleSystem allclean

