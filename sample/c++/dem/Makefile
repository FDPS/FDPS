INC = -I ../../../src
CC = g++
#CC = mpic++
CFLAGS = -O3 -ffast-math -funroll-loops -std=c++17

#FDPS FLAGS
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
CFLAGS += -DPARTICLE_SIMULATOR_ALL_64BIT_PRECISION
CFLAGS += -DFAST_COMM_FOR_2ND_EXCHANGE
## FDPS 2D MODE
#CFLAGS += -DPARTICLE_SIMULATOR_TWO_DIMENSION

CPPOBJS = main.cpp
PROGRAM = a.out

.PHONY:	clean all

all:	$(CPPOBJS)
	@echo "Linking object files..."
	@$(CC) $(CFLAGS) $(CPPOBJS) -o $(PROGRAM) $(INC)
	@echo "Link Success! [$(PROGRAM)]"

%.o:	%.cpp
	@echo "Building $< ..."
	@$(CC) -c $< $(CFLAGS) $(INC)
	@echo "[$< OK]"

clean:
	-rm *.out *.o

