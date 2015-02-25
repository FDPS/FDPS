PS_PATH = -I../../src/

CC = g++
#CC = mpicxx
CFLAGS = -O3 -ffast-math -funroll-loops
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

MULEXP = -DMONOPOLE

SRC = nbody.cpp
PROGRAM = $(SRC:%.cpp=%.out)

$(PROGRAM):$(SRC)
	$(CC) $(MULEXP) $(PS_PATH) $(CFLAGS) -o $@ $<

clean:
	rm $(PROGRAM)

