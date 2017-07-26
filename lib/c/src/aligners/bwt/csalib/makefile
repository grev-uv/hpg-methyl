COMMONS_DIR = ../commons/

#CC = gcc43
CC = gcc

#CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -msse4.2 -O9 -fPIC
CFLAGS = -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -O9 -fPIC -std=c99

CFILES = csa.c mmap.c psi1.c diskbuf.c lf_dna.c lf_wt.c densearray.c huffman.c comparray.c psi2.c densearray2.c sparsearray2.c sparsearray.c cst.c lf_bit.c
OFILES = csa.o mmap.o psi1.o diskbuf.o lf_dna.o lf_wt.o densearray.o huffman.o comparray.o psi2.o densearray2.o sparsearray2.o sparsearray.o cst.o lf_bit.o
HFILES = typedef.h csa6.h mman.h psi1.h diskbuf.h lf_dna.h lf_wt.h densearray.h huffman.h comparray.h psi2.h densearray2.h sparsearray2.h sparsearray.h lf_bit.h cst.h

all: ../bin/mkcsa ../bin/csa ../bin/unbwt ../bin/mkcst ../bin/approx

csalib.a: $(OFILES)
	ar rc csalib.a $(OFILES)

../bin/mkcsa: mkarray2.c csalib.a commons
	$(CC) $(CFLAGS) -o $@ mkarray2.c csalib.a $(COMMONS_DIR)/commons.a
../bin/csa: test_csa.c csalib.a commons
	$(CC) $(CFLAGS) -o $@ test_csa.c csalib.a $(COMMONS_DIR)/commons.a
../bin/unbwt: unbwt.c csalib.a commons
	$(CC) $(CFLAGS) -o $@ unbwt.c csalib.a $(COMMONS_DIR)/commons.a
../bin/mkcst: mkcst.c
	$(CC) $(CFLAGS) -o $@ mkcst.c
../bin/approx: test_approx.c approx.c csalib.a commons
	$(CC) $(CFLAGS) -o $@ test_approx.c approx.c csalib.a $(COMMONS_DIR)/commons.a

commons:
	make -C $(COMMONS_DIR) commons.a

clean:
	rm -f *~ \#*\# .*.swp *.o *.a
