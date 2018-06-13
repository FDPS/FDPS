PS_PATH = -I ../../../src/

CC = time g++
#CC = mpic++-mpich-gcc49
CFLAGS = -O3
#CFLAGS += -ffast-math
#CFLAGS += -funroll-loops 
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
#CFLAGS += -fopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CFLAGS += #-Wall -Wformat=2 -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wfloat-equal -Wpointer-arith

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
CPPHDRS = $(wildcard *.h)
PROGRAM = sph.out

.PHONY:	clean all

all:	$(CPPOBJS) $(CPPHDRS)
	@echo "Linking object files..."
	@$(CC) $(CFLAGS) $(WARNINGS) $(CPPOBJS) -o $(PROGRAM) $(LIBS) $(PS_PATH)
	@echo "Link Success! [$(PROGRAM)]"

%.o:	%.cpp $(CPPHDRS)
	@echo "Bulding $< ..."
	@$(CC) -c $< $(CFLAGS) $(WARNINGS) $(PS_PATH)
	@echo "[$< OK]"

clean:
	-rm *.out *.o

distclean: clean
	rm -rf result

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
