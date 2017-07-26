CC = gcc
CFLAGS = -std=c99 -O3
CFLAGS_DEBUG = -std=c99 -g

# Directories' path
ifdef LIBS_HOME
	LIBS_ROOT = $(LIBS_HOME)
else
	LIBS_ROOT = $(PWD)/../..
endif

SRC_DIR = ..
CONTAINERS_DIR = $(LIBS_ROOT)/common-libs/containers
COMMONS_DIR = $(LIBS_ROOT)/common-libs/commons
REGION_DIR = $(LIBS_ROOT)/bioinfo-libs/bioformats/features/region
GFF_DIR = $(LIBS_ROOT)/bioinfo-libs/bioformats/gff

INCLUDES = -I $(SRC_DIR) -I $(LIBS_HOME) -I $(LIBS_HOME)/common-libs -I $(LIBS_HOME)/bioinfo-libs
LIBS = -lpthread -lcheck -lcprops -fopenmp


#SRC_OBJECTS = $(SRC_DIR)/vcf_filters.c $(SRC_DIR)/vcf_file.c $(SRC_DIR)/vcf_read.c $(SRC_DIR)/vcf_reader.c $(SRC_DIR)/vcf_write.c $(SRC_DIR)/vcf_batch.c $(SRC_DIR)/vcf_util.c
SRC_OBJECTS = $(SRC_DIR)/vcf_*.c
GFF_OBJECTS = $(GFF_DIR)/gff_file.c $(GFF_DIR)/gff_read.c $(GFF_DIR)/gff_reader.c $(GFF_DIR)/gff_write.c $(GFF_DIR)/gff_batch.c
REGION_TABLE_OBJECTS = $(REGION_DIR)/region_table.o $(REGION_DIR)/region_table_utils.o $(REGION_DIR)/region.o

MISC_OBJECTS = $(CONTAINERS_DIR)/list.o $(COMMONS_DIR)/string_utils.o $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/log.o

# Test executables
FILTER_OBJECTS = test_vcf_filters.c $(SRC_OBJECTS) $(REGION_TABLE_OBJECTS) $(GFF_OBJECTS) $(MISC_OBJECTS)
REGION_PARSE_OBJECTS = test_parse_regions.c $(GFF_OBJECTS) $(REGION_TABLE_OBJECTS) $(MISC_OBJECTS) $(SRC_DIR)/vcf_util.c
STATS_OBJECTS = test_stats.c $(SRC_OBJECTS) $(REGION_TABLE_OBJECTS) $(GFF_OBJECTS) $(MISC_OBJECTS)


all: list.o string_utils.o region.o region_table.o

filter.test: $(FILTER_OBJECTS)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o $@ $(FILTER_OBJECTS) $(INCLUDES) $(LIBS)
#	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -o $@ $(FILTER_OBJECTS) $(INCLUDES) $(LIBS)

parse_region.test: $(REGION_PARSE_OBJECTS)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o $@ $(REGION_PARSE_OBJECTS) $(INCLUDES) $(LIBS)

stats.test: $(STATS_OBJECTS)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o $@ $(STATS_OBJECTS) $(INCLUDES) $(LIBS)

list.o:
	cd $(CONTAINERS_DIR) && make list.o

string_utils.o:
	cd $(COMMONS_DIR) && make string_utils.o

region.o:
	cd $(REGION_DIR) && make LIBS_HOME=$(LIBS_HOME) region.o

region_table.o:
	cd $(REGION_DIR) && make LIBS_HOME=$(LIBS_HOME) region_table.o

region_table_utils.o:
	cd $(REGION_DIR) && make LIBS_HOME=$(LIBS_HOME) region_table_utils.o

clean:
	rm list.o
	rm string_utils.o
	rm region.o
	rm region_table.o
	rm *.test

