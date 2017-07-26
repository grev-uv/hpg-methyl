#CC = gcc
#CC = /opt/intel/bin/icc

ifeq ($(COMPILER), intel)
	CC = icc 
else
	CC = gcc
endif


CFLAGS =  -O3 -fopenmp -std=c99 -D_GNU_SOURCE -DVECTOR_O_64BIT_COMPRESSION -DVERBOSE -g
#CFLAGS = -DVERBOSE_DBG -Wall
#CFLAGS = -Wall -pg

NVCC = nvcc
NVCCFLAGS = -g -G -Xptxas -v -arch=sm_20 -DVECTOR_O_64BIT_COMPRESSION -DVERBOSE
#NVCCFLAGS = -g -G -Xptxas -v  -arch=sm_12 -DVECTOR_O_COMPRESSION
#NVCCFLAGS = -g -G -Xptxas -v

# Directories' path
ifdef LIBS_HOME
	LIBS_ROOT = $(LIBS_HOME)
else
	LIBS_ROOT = $(PWD)/../../..
endif

## This param must point to samtools project, 'make' is assumed
SAMTOOLS_DIR = ~/appl/bioinfo-c/libs/ext/samtools-0.1.18

CONTAINERS_DIR = $(LIBS_ROOT)/common-libs/containers
COMMONS_DIR = $(LIBS_ROOT)/common-libs/commons
FASTQ_DIR = $(LIBS_ROOT)/bioinfo-libs/bioformats/fastq
BAM_DIR = $(LIBS_ROOT)/bioinfo-libs/bioformats/bam-sam

CONTAINERS_DIR = $(LIBS_ROOT)/common-libs/containers
COMMONS_DIR = $(LIBS_ROOT)/common-libs/commons
BIOFORMATS_DIR = $(LIBS_ROOT)/bioinfo-libs/bioformats

# Include path and libraries
INCLUDES = -I . -I $(LIBS_ROOT) -I $(LIBS_ROOT)/common-libs -I $(LIBS_ROOT)/bioinfo-libs -I $(SAMTOOLS_DIR)
LIBS = -L $(SAMTOOLS_DIR) -lcprops -fopenmp -lcurl -Wl,-Bsymbolic-functions -lm -lbam -lz

INCLUDES_STATIC = -I . -I $(LIBS_ROOT) -I $(LIBS_ROOT)/../include -I $(LIBS_ROOT)/common-libs -I $(LIBS_ROOT)/bioinfo-libs -I $(SAMTOOLS_DIR)
LIBS_STATIC = -L $(SAMTOOLS_DIR)  -L $(LIBS_ROOT) -L $(LIBS_ROOT)/cprops -lcprops -fopenmp -lm -lbam -lz

#CUDA_HOME = /usr/local/cuda
CUDA_HOME = /opt/cuda/4.1/cuda
CUDA_LIBS = -L $(CUDA_HOME)/lib64 -lcudart

# Source and object files
MAIN_SOURCES = bwt.c BW_io.c BW_search.c BW_preprocess.c
MAIN_OBJECTS = bwt.o BW_io.o BW_search.o BW_preprocess.o

GPU_SOURCES = bwt_gpu.c
GPU_OBJECTS = bwt_gpu.o

CUDA_SOURCES = gpu.cu BW_gpu.cu
CUDA_OBJECTS = gpu.o BW_gpu.o

MISC_OBJECTS = $(FASTQ_DIR)/fastq_file.o $(FASTQ_DIR)/fastq_read.o $(FASTQ_DIR)/fastq_batch.o $(FASTQ_DIR)/fastq_batch_reader.o $(BAM_DIR)/alignment.o $(CONTAINERS_DIR)/linked_list.o $(CONTAINERS_DIR)/array_list.o $(CONTAINERS_DIR)/list.o $(COMMONS_DIR)/log.o $(COMMONS_DIR)/system_utils.o $(COMMONS_DIR)/string_utils.o

TEST = test


#LIBS_DIR = ../../../libs
#COMMONS_DIR = $(LIBS_DIR)/commons
#COMMONS_CUDA_DIR = $(LIBS_DIR)/commons-cuda
#CONTAINERS_DIR = $(LIBS_DIR)/containers
#FASTQ_DIR = $(LIBS_DIR)/bioformats/fastq
#SAM_DIR = $(LIBS_DIR)/bioformats/bam-sam
#SAMTOOLS_DIR = ../../../hpg-aligner/samtools-0.1.18/
#SAMTOOLS_DIR = /home/jtarraga/appl/bioinfo-c-old/ngs-hpc-pipeline/bam-hpc-tools/samtools-0.1.18
#SAMTOOLS_DIR = ../../../samtools-0.1.18/
#CPROPS_DIR = ../../../libcprops-0.1.12/

#LIB = -L$(SAM_DIR) -L$(SAMTOOLS_DIR) -L$(CPROPS_DIR) -lpthread -lcprops -lbam -lm -lz -fopenmp 


#CINCLUDES = -I. -I$(COMMONS_DIR) -I$(COMMONS_CUDA_DIR) -I$(CONTAINERS_DIR) -I$(FASTQ_DIR) -I$(CPROPS_DIR) -I$(SAMTOOLS_DIR) -I$(SAM_DIR) 
#CUINCLUDES = -I. -I$(COMMONS_DIR) -I$(COMMONS_CUDA_DIR) -I$(CONTAINERS_DIR) -I$(FASTQ_DIR) -I$(CPROPS_DIR) -I$(SAMTOOLS_DIR) -I$(SAM_DIR) -I/opt/cuda/include

# Targets
all: test-preprocess test-search-cpu

test-preprocess: compile-static $(TEST)/test-preprocess.c
	$(CC) $(CFLAGS) $(MAIN_OBJECTS) $(MISC_OBJECTS) $(TEST)/test-preprocess.c -o bin/test-preprocess $(INCLUDES) $(LIBS_STATIC)

test-search-cpu: compile-static $(TEST)/test-search-cpu.c
	$(CC) $(CFLAGS) $(MAIN_OBJECTS) $(MISC_OBJECTS) $(TEST)/test-search-cpu.c -o bin/test-search-cpu $(INCLUDES) $(LIBS_STATIC)


test-search-batch-cpu: compile-static $(TEST)/test-search-batch-cpu.c
	$(CC) $(CFLAGS) $(MAIN_OBJECTS) $(MISC_OBJECTS) $(TEST)/test-search-batch-cpu.c -o bin/test-search-batch-cpu $(INCLUDES) $(LIBS_STATIC)

test-search-gpu: compile-static compile-gpu $(TEST)/test-search-gpu.c
	$(CC) $(CFLAGS) $(CUDA_OBJECTS) $(GPU_OBJECTS) $(MAIN_OBJECTS) $(MISC_OBJECTS) $(TEST)/test-search-gpu.c -o bin/test-search-gpu $(INCLUDES) $(LIBS_STATIC) $(CUDA_LIBS)

test-search-seeds-gpu: compile-static compile-gpu $(TEST)/test-search-gpu.c
	$(CC) $(CFLAGS) $(CUDA_OBJECTS) $(GPU_OBJECTS) $(MAIN_OBJECTS) $(MISC_OBJECTS) $(TEST)/test-search-seeds-gpu.c -o bin/test-search-seeds-gpu $(INCLUDES) $(LIBS_STATIC) $(CUDA_LIBS)

compile-gpu: compile-cuda
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c $(GPU_SOURCES) $(INCLUDES_STATIC)

compile-cuda: 
	$(NVCC) $(NVCCFLAGS) -c $(CUDA_SOURCES) $(INCLUDES_STATIC)


compile: compile-dependencies
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c $(MAIN_SOURCES) $(INCLUDES)

compile-static: compile-dependencies-static
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c $(MAIN_SOURCES) $(INCLUDES_STATIC)

compile-dependencies:
	cd $(CONTAINERS_DIR) && make COMPILER=$(COMPILER) compile && \
	cd $(COMMONS_DIR) && make COMPILIER=$(COMPILER) compile && \
	cd $(BIOFORMATS_DIR)/fastq && make COMPILER=$(COMPILER) compile && \
	cd $(BIOFORMATS_DIR)/bam-sam && make COMPILER=$(COMPILER) compile

compile-dependencies-static:
	cd $(CONTAINERS_DIR) && make COMPILER=$(compiler) compile && \
	cd $(COMMONS_DIR) && make COMPILER=$(COMPILER) compile && \
	cd $(BIOFORMATS_DIR)/fastq && make COMPILER=$(COMPILER) compile && \
	cd $(BIOFORMATS_DIR)/bam-sam && make COMPILER=$(COMPILER) compile
clean:
	-rm -f *~ \#*\# *.o $(ALL) 

