PS_PATH = -I ../fdps/src/src_parallel/
#PS_PATH = -I ./header/

#CPPC = mpiFCCpx
#FLAGS = -Kfast,openmp
#FLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL
#FLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#FLAGS += -Xg
#FLAGS += -Nfjcex
#FLAGS += -Krestp=all
#FLAGS = -Kfast,openmp,restp=all

CPPC = mpicxx-openmpi-gcc49
#CPPC = g++
FLAGS = -O3 -ffast-math -funroll-loops #-Wall -Wformat=2 -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wfloat-equal -Wpointer-arith #-Weffc++
FLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
FLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

CPPOBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
CPPHDRS = $(wildcard *.h)
PSHDRS  = $(wildcard ./header/*.hpp)
PROGRAM = a.out

.PHONY:	clean

$(PROGRAM):	$(CPPOBJS) $(CPPHDRS) $(PSHDRS)
	$(CPPC) $(FLAGS) $(WARNINGS) $(CPPOBJS) -o $(PROGRAM) $(LIBS) $(PS_PATH)

%.o:	%.cpp $(CPPHDRS)
	$(CPPC) -c $< $(FLAGS) $(WARNINGS) $(PS_PATH)


clean:
	-rm *.out *.o

