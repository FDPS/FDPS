CC=gcc
CFLAGS=-mavx -O3 -Wall -ffast-math -funroll-loops

#enable_openmp = yes
#enable_avx2 = yes

ifeq ($(enable_openmp), yes)
CFLAGS += -fopenmp -DENABLE_OPENMP
endif

LIBPG5OBJ = gravity.o phantom_g5.o rsqrt.o pg5_fortran.o
ifeq ($(enable_avx2), yes)
LIBPG5OBJ = gravity_avx2.o phantom_g5.o rsqrt.o pg5_fortran.o
endif

libpg5.a: $(LIBPG5OBJ)
	ar -r $@ $^
	ranlib $@

clean:
	rm -rf *.o
distclean: clean
	rm -rf *~
	rm -rf libpg5.a
