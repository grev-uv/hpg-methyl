CC = gcc
CFLAGS = -std=c99 -O3 -D_GNU_SOURCE -D_XOPEN_SOURCE=600
CFLAGS_DEBUG = -std=c99 -g -D_GNU_SOURCE

# Directories' paths
ifdef LIBS_HOME
	LIBS_ROOT = $(LIBS_HOME)
else
	LIBS_ROOT = $(PWD)/../../..
endif

CONTAINERS_DIR = $(LIBS_ROOT)/common-libs/containers
COMMONS_DIR = $(LIBS_ROOT)/common-libs/commons
BIOFORMATS_DIR = $(LIBS_ROOT)/bioinfo-libs/bioformats
REGION_DIR = $(BIOFORMATS_DIR)/features/region

# Include path and libraries
INCLUDES = -I . -I $(LIBS_ROOT) -I $(LIBS_ROOT)/common-libs -I $(LIBS_ROOT)/bioinfo-libs
LIBS = -lcprops -fopenmp -lcurl -Wl,-Bsymbolic-functions -lm

INCLUDES_STATIC = -I . -I $(LIBS_ROOT) -I $(LIBS_ROOT)/../include -I $(LIBS_ROOT)/common-libs -I $(LIBS_ROOT)/bioinfo-libs
LIBS_STATIC = -L $(LIBS_ROOT) -lcprops -fopenmp -lm

# Source and object files
# MAIN_OBJECTS = vcf_file.c vcf_file_structure.c vcf_reader.c vcf_ragel.c vcf_write.c vcf_batch.c vcf_filters.c vcf_stats.c vcf_util.c
MAIN_OBJECTS = vcf_file.c vcf_file_structure.c vcf_reader.c vcf_ragel.c vcf_write.c vcf_filters.c vcf_stats.c vcf_util.c
MISC_OBJECTS = $(CONTAINERS_DIR)/array_list.o $(CONTAINERS_DIR)/list.o $(COMMONS_DIR)/http_utils.o $(COMMONS_DIR)/string_utils.o $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/log.o
REGION_TABLE_OBJECTS = $(REGION_DIR)/region_table.o $(REGION_DIR)/region_table_utils.o $(REGION_DIR)/region.o
GFF_OBJECTS = $(BIOFORMATS_DIR)/gff/*.o
VCF_OBJECTS = ./vcf_*.o

VCF_BASIC_OBJECTS = main_basic.c $(MAIN_OBJECTS) $(MISC_OBJECTS)
VCF_BATCH_OBJECTS = main_batch.c $(VCF_OBJECTS) $(MISC_OBJECTS) $(REGION_TABLE_OBJECTS) $(GFF_OBJECTS)
VCF_FILTER_OBJECTS = main_filter.c $(MAIN_OBJECTS) $(MISC_OBJECTS) $(REGION_TABLE_OBJECTS) $(GFF_OBJECTS)


# Targets
all: compile-dependencies compile

compile: compile-dependencies
	$(CC) $(CFLAGS) -c $(MAIN_OBJECTS) $(INCLUDES) $(LIBS)

compile-static: compile-dependencies-static
	$(CC) $(CFLAGS) -c $(MAIN_OBJECTS) $(INCLUDES_STATIC) $(LIBS_STATIC)

ragel:
	ragel -o vcf_ragel.c -C -G2 vcf.ragel

ragel_graphviz:
	ragel -o vcf_ragel.dot -V -C -G2 vcf.ragel

main_batch: compile-dependencies compile
	$(CC) $(CFLAGS) -o $@ $(VCF_BATCH_OBJECTS) $(INCLUDES) $(LIBS)

main_filter: compile-dependencies ragel compile $(VCF_FILTER_OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(VCF_FILTER_OBJECTS) $(INCLUDES) $(LIBS)

compile-dependencies:
	cd $(CONTAINERS_DIR) && make compile && \
	cd $(COMMONS_DIR) && make compile && \
	cd $(BIOFORMATS_DIR) && make region gff
#	cd $(REGION_DIR) && make region.o && \
	cd $(BIOFORMATS_DIR)/gff && make compile

compile-dependencies-static:
	cd $(CONTAINERS_DIR) && make compile && \
	cd $(COMMONS_DIR) && make compile && \
	cd $(BIOFORMATS_DIR) && make region gff
#	cd $(REGION_DIR) && make region.o && \
	cd $(BIOFORMATS_DIR)/gff && make compile-static

clean:
	rm -f *.o
	rm main_batch
	rm main_filter

