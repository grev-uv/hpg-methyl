#Debian distro location
#CUDA_LOCATION=/usr/lib/nvidia-cuda-toolkit

CC = gcc
#NVCC = nvcc
CFLAGS = -Wall -m64 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 -O9 -std=c99 -fPIC -fopenmp -DFM_COMP_64 #-DCSALIB_SEARCH #-I $(CUDA_LOCATION)/include/ -L /usr/local/cuda/lib64/ -lcudart -msse4.2 #-DVERBOSE_DBG
#NVCCFLAGS = --compiler-options -Wall,-fopenmp,-m64 -O3 -Xptxas -v -arch=sm_13 -DFM_COMP_32 -D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64 #-DVERBOSE_DBG

DEPEN = io.o csafm.o search.o preprocess.o results.o

.PHONY: clean

search.a: $(DEPEN)
	ar rc search.a $(DEPEN)

clean:
	rm -f *~ \#*\# .*.swp *.o *.a
