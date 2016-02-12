CC = cc -march=core-avx2
CFLAGS = -O3 -ffast-math -funroll-loops

#enable_openmp=yes
enable_avx2=yes
#enable_symmetric=yes
#enable_symmetric=no

ifeq ($(enable_openmp), yes)
LDFLAGS += -fopenmp
CFLAGS += -fopenmp -DENABLE_OPENMP
endif

ifeq ($(enable_symmetric), yes)
CFLAGS += -DSYMMETRIC
endif

OBJS = direct.o cpu.o libpg5/libpg5.a

direct: $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS) -lm
direct.o: direct.c
	$(CC) $(CFLAGS) -c direct.c
cpu.o: cpu.c
	$(CC) $(CFLAGS) -c cpu.c
libpg5/libpg5.a: #libpg5/libpg5.a
	make -C libpg5 CC="$(CC)" CFLAGS="$(CFLAGS)" enable_openmp=$(enable_openmp) enable_avx2=$(enable_avx2) libpg5.a

clean:
	rm -rf *.o
	make -C libpg5 clean
distclean: clean
	rm -rf direct *~
	make -C libpg5 distclean
