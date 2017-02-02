#CC = icc
#CFLAGS = -std=c99 -O3 -openmp -Wall -Wno-missing-declarations -fast -xavx -vec-report1

#CC = gcc
ifeq ($(COMPILER), intel)
	CC = icc 
else
	CC = gcc
endif

#CC = /opt/intel/bin/icc

CFLAGS = -std=c99 -O3 -fopenmp -g

INCLUDES = -I .
LIBS = -lxml2 -lm
TEST = test

all: hpg-sw hpg-sw-bench

hpg-sw: smith_waterman.o sse.o macros.o $(TEST)/hpg_sw.h
	$(CC) $(CFLAGS) smith_waterman.o macros.o sse.o $(TEST)/hpg_sw.c -o bin/hpg-sw $(INCLUDES) $(LIBS)

hpg-sw-bench: smith_waterman.o sse.o emboss.o macros.o $(TEST)/hpg_sw_bench.h
	$(CC) $(CFLAGS) smith_waterman.o macros.o emboss.o sse.o $(TEST)/hpg_sw_bench.c -o bin/hpg-sw-bench $(INCLUDES) $(LIBS)

macros.o: macros.c macros.h
	$(CC) $(CFLAGS) -c macros.c $(INCLUDES)

emboss.o: emboss.c emboss.h
	$(CC) $(CFLAGS) -c emboss.c $(INCLUDES)

sse.o: sse.c sse.h
	$(CC) $(CFLAGS) -c sse.c $(INCLUDES)

smith_waterman.o: smith_waterman.c smith_waterman.h
	$(CC) $(CFLAGS) -c smith_waterman.c $(INCLUDES)

clean:
	-rm -rf main *.o bin/hpg-sw bin/hpg-sw-bench




