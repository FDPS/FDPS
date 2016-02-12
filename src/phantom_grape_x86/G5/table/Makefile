#CC = gcc
#CC=gcc-4.8.1 -march=core-avx2
#CFLAGS = -O3 -ffast-math -funroll-loops -Wall
OBJS = gravity_kernel.o phantom_g5.o pg5_table.o

ifeq ($(enable_avx2), yes)
OBJS = gravity_kernel_avx2.o phantom_g5.o pg5_table.o
endif

ifeq ($(enable_ntdr), yes)
OBJS = gravity_kernel_avx2_ntdr.o phantom_g5.o pg5_table.o
endif

#all: $(OBJS)
#
#asm: gravity_kernel.s
#
#clean:
#	rm *.o *.s *.a *.swp a.out -f
#
#distclean: clean
#	rm -f *~
#
#.c.o:
#	$(CC) $(CFLAGS) -c $^
#
#.c.s:
#	$(CC) $(CFLAGS) -S $^
#
#libpg5.a: $(OBJS)
#	ar -r $@ $^
#	ranlib $@

libpg5.a: $(OBJS)
	ar -r $@ $^
	ranlib $@

clean:
	rm -rf *.o

distclean: clean
	rm -rf *~
	rm -rf libpg5.a
