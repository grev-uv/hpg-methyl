CC = gcc
CFLAGS = -std=c99 -O3 -D_GNU_SOURCE
CFLAGS_DEBUG = -std=c99 -g

# Directories' path
ifdef LIBS_HOME
	LIBS_ROOT = $(LIBS_HOME)
else
	LIBS_ROOT = $(PWD)/../../../..
endif

SRC_DIR = ..
COMMONS_DIR = $(LIBS_ROOT)/common-libs/commons
CONTAINERS_DIR = $(LIBS_ROOT)/common-libs/containers

INCLUDES = -I $(SRC_DIR) -I $(LIBS_HOME) -I $(LIBS_ROOT) -I $(LIBS_ROOT)/common-libs -I $(LIBS_ROOT)/bioinfo-libs
LIBS = -lpthread -lcheck -lcprops -fopenmp 


#SRC_OBJECTS = $(SRC_DIR)/vcf_filters.c $(SRC_DIR)/vcf_file.c $(SRC_DIR)/vcf_read.c $(SRC_DIR)/vcf_reader.c $(SRC_DIR)/vcf_write.c $(SRC_DIR)/vcf_batch.c $(SRC_DIR)/vcf_util.c
SRC_OBJECTS = $(SRC_DIR)/fastq_*.c

MISC_OBJECTS = $(CONTAINERS_DIR)/array_list.o $(CONTAINERS_DIR)/list.o $(COMMONS_DIR)/string_utils.o $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/log.o

# Test executables
FILTER_OBJECTS = test_fastq_filter.c $(SRC_OBJECTS) $(MISC_OBJECTS)
#STATS_OBJECTS = test_stats.c $(SRC_OBJECTS) $(REGION_TABLE_OBJECTS) $(GFF_OBJECTS) $(MISC_OBJECTS)


all: array_list.o list.o string_utils.o file_utils.o log.o

filter.test: $(FILTER_OBJECTS)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o $@ $(FILTER_OBJECTS) $(INCLUDES) $(LIBS)
#	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -o $@ $(FILTER_OBJECTS) $(INCLUDES) $(LIBS)

stats.test: $(STATS_OBJECTS)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o $@ $(STATS_OBJECTS) $(INCLUDES) $(LIBS)

array_list.o:
	cd $(CONTAINERS_DIR) && make array_list.o

list.o:
	cd $(CONTAINERS_DIR) && make list.o
	
string_utils.o:
	cd $(COMMONS_DIR) && make string_utils.o

file_utils.o:
	cd $(COMMONS_DIR) && make file_utils.o

log.o:
	cd $(COMMONS_DIR) && make log.o
		
clean:
	rm *.o
	rm *.test

