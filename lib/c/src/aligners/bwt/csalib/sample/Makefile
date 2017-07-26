all: leftfrequent rightfrequent

#CC = gcc43
CC = gcc

#CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -msse4.2 -O9 -fPIC -I..
CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -O9 -fPIC -std=c99 -I..

leftfrequent: leftfrequent.c ../csa.h ../csa.a
	$(CC)  $(CFLAGS) -o leftfrequent leftfrequent.c ../csa.a
rightfrequent: rightfrequent.c ../csa.h ../csa.a
	$(CC)  $(CFLAGS) -o rightfrequent rightfrequent.c ../csa.a
