PS_PATH = -I../../src/
#PS_PATH = -I../../fdps/fdps/src/
#PS_PATH = -I../../sandbox/iwasawa_box/FDPS-master/src/

CC = g++
#CC = mpicxx
CFLAGS = -O3 -ffast-math -funroll-loops
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

#MULEXP = -DMONOPOLE

SRC = vdwtest.cpp
PROGRAM = $(SRC:%.cpp=%.out)

$(PROGRAM):$(SRC)
	$(CC) $(MULEXP) $(PS_PATH) $(CFLAGS) -o $@ $<

a.out: test.cpp
	$(CC) $(MULEXP) $(PS_PATH) $(CFLAGS) test.cpp

clean:
	rm $(PROGRAM)

