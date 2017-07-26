CC = nvcc
CFLAGS = --compiler-options -Wall,-fPIC,-m64 -O9 -Xptxas -v -arch=sm_13 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -DFM_COMP_64 #-DCSALIB_SEARCH #-DVERBOSE_DBG

DEPEN = gpu.o

.PHONY: clean

gpu.a: $(DEPEN)
	ar rc gpu.a $(DEPEN)

gpu.o: gpu.cu gpu.cuh
	$(CC) $(CFLAGS) -c gpu.cu
	
clean:
	rm -f *~ \#*\# .*.swp *.o *.a
