COMMONS_DIR = ../commons/

CC = gcc
CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -O3 -std=c99 -fPIC
DEPEN = queue.o utils.o sais.o dbwt.o
OBJECTS = $(DEPEN) main.o

.PHONY: clean all

all: ../bin/dbwt

dbwt.a: $(DEPEN)
	ar rc dbwt.a $(DEPEN)

../bin/dbwt: main.c dbwt.a commons
	$(CC) $(CFLAGS) -o $@ main.c $(COMMONS_DIR)/commons.a dbwt.a

commons:
	make -C $(COMMONS_DIR) commons.a

clean:
	rm -f *~ \#*\# .*.swp *.o *.a
