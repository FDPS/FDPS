CC = gcc -march=core-avx2 -O3 -ffast-math -funroll-loops -fopenmp
CFLAGS += -DNORMAL
#CFLAGS += -DONLYNNB
#CFLAGS += -DONLYNEIGHBOUR

enable_avx2 = yes

TARGET = banana01

OBJ    = $(TARGET).o io.o calc_energy.o hermite.o \
	timestep.o timeprof.o
HEADER = $(TARGET).h particle.h io.h calc_energy.h hermite.h \
	globaltime.h timestep.h step.h timeprof.h
ETC = Makefile

LIBDIR = ./libavx/
LDFLGS = -L$(LIBDIR) -lg6avx
LIBFLE = $(LIBDIR)/libg6avx.a
INCDIR = $(LIBDIR)
ICFLGS = -I$(INCDIR)

$(TARGET) : $(OBJ) $(HEADER) $(ETC) $(LIBFLE)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJ) -lm $(LDFLGS)
$(TARGET).o : $(TARGET).c $(HEADER) $(ETC)
	$(CC) $(CFLAGS) -c $(TARGET).c $(ICFLGS)
io.o : io.c $(HEADER) $(ETC)
	$(CC) $(CFLAGS) -c io.c
calc_energy.o : calc_energy.c $(HEADER) $(ETC)
	$(CC) $(CFLAGS) -c calc_energy.c
hermite.o : hermite.c $(HEADER) $(ETC)
	$(CC) $(CFLAGS) -c hermite.c $(ICFLGS)
timestep.o : timestep.c $(HEADER) $(ETC)
	$(CC) $(CFLAGS) -c timestep.c
timeprof.o : timeprof.c $(HEADER) $(ETC)
	$(CC) $(CFLAGS) -c timeprof.c

$(LIBFLE):
	make -C libavx CC="$(CC)" enable_avx2=$(enable_avx2)

$(PHANTOM).s : $(PHANTOM).c $(HEADER) $(ETC)
	$(CC) -S $(CFLAGS) -c $(PHANTOM).c

clean :
	rm -f *.o *.s *~ $(TARGET)
	make -C libavx clean
