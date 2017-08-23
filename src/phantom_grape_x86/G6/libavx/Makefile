CC = gcc -march=core-avx2 -O3 -ffast-math -funroll-loops -fopenmp -Wall
AR = ar

HDR = gp6util.h gravity.h avx_type.h timeprof.h avx.h
OBJ = phantom_g6.o gravity.o timeprof.o
ETC = Makefile

#enable_avx2 = yes
ifeq ($(enable_avx2),yes)
OBJ = phantom_g6.o gravity_avx2.o timeprof.o
endif

libg6avx.a : $(OBJ)
	$(AR) rv $@ $^
	ranlib $@

clean:
	rm -f *.o *.s *~ libg6avx.a
