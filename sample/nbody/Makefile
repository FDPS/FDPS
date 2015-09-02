PS_PATH = -I../../src/

CC = g++
#CC = time CC
#CC = mpicxx
#CFLAGS = -O3 -Wall
CFLAGS = -O3 -ffast-math -funroll-loops -Wall
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

#use_phantom_grape_x86 = yes

ifeq ($(use_phantom_grape_x86),yes)
PS_PATH += -I../../src/phantom_grape_x86/G5/newton/libpg5/
CFLAGS  += -DENABLE_PHANTOM_GRAPE_X86
CLIBS   += -L../../src/phantom_grape_x86/G5/newton/libpg5/ -lpg5
endif

PROGRAM = nbody.out

nbody.out:nbody.cpp 
	$(CC) $(MULEXP) $(PS_PATH) $(CFLAGS) -o $@ $< $(CLIBS)

clean:	
	rm -f $(PROGRAM)

distclean:	
	rm -f $(PROGRAM) *~
