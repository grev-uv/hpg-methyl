all: csa_run

#CC = gcc43
CC = gcc

#CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -msse4.2 -O9 -fPIC
CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -O9 -fPIC -std=c99

csa_run: run_queries.c csa.a interface.c interface.h
	$(CC) $(CFLAGS) -o csa_run run_queries.c interface.c ../csa.a

