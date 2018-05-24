CC = mpicxx
#CC := time mpiFCCpx

CFLAGS = -O3 -ffast-math -funroll-loops -DMPICH_IGNORE_CXX_SEEK #-Wall
#CFLAGS = -O3 -Kfast

INCLUDE_FFTW = -I/work4a/local/fftw/include

FFTW3_IDENTIFIER = FFTW3_PARALLEL

CFLAGS += $(INCLUDE_FFTW)
CFLAGS += -D$(FFTW3_IDENTIFIER)
CFLAGS += -DVERBOSE_MODE
CFLAGS += -DParallel_IO
CFLAGS += -DRESTART2
CFLAGS += -DEXCHANGE_COMM_NONBLOCKING
CFLAGS += -DUSING_MPI_PARTICLE
CFLAGS += -DBOUNDARY_COMM_NONBLOCKING
CFLAGS += -DBOUNDARY_COMM_NONBLOCKING_OVERLAP
CFLAGS += -DGRAPE_OFF
#CFLAGS += -DFFTW_DOUBLE

OBJ = pm_parallel.o tools.o io.o evolve.o decomposition.o \
	ptobmp.o pp2.o
HDR = gadget_param.h param.h qsort_omp.h treepm_header.h \
	openmp_param.h pm_parallel.h tools.h pp.h particle_mesh_class.hpp \
	param_fdps.h

libpm.a : $(OBJ) $(HDR)
	$(AR) rv libpm.a $(OBJ); \
	ranlib libpm.a; \
	sh make_header.sh $(FFTW3_IDENTIFIER) > particle_mesh.hpp

pm_parallel.o: pm_parallel.cpp $(HDR)
	$(CC) $(CFLAGS) -c pm_parallel.cpp

tools.o: tools.cpp $(HDR)
	$(CC) $(CFLAGS) -c tools.cpp

io.o: io.cpp $(HDR)
	$(CC) $(CFLAGS) -c io.cpp

evolve.o: evolve.cpp $(HDR)
	$(CC) $(CFLAGS) -c evolve.cpp

decomposition.o: decomposition.cpp $(HDR)
	$(CC) $(CFLAGS) -c decomposition.cpp

ptobmp.o: ptobmp.cpp $(HDR)
	$(CC) $(CFLAGS) -c ptobmp.cpp

pp2.o: pp2.cpp $(HDR)
	$(CC) $(CFLAGS) -c pp2.cpp

clean:
	rm -f *~
distclean:
	make clean
	rm -f *.o
allclean:
	make distclean
	rm -f *.a
