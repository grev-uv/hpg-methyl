CC = gcc
CFLAGS = -std=c99 -O3
CFLAGS_DEBUG = -std=c99 -g

# Directories' path
ifdef LIBS_HOME
	LIBS_ROOT = $(LIBS_HOME)
else
	LIBS_ROOT = $(PWD)/../../..
endif

CONTAINERS_DIR = $(LIBS_ROOT)/common-libs/containers
COMMONS_DIR = $(LIBS_ROOT)/common-libs/commons
REGION_DIR = $(LIBS_ROOT)/bioinfo-libs/bioformats/features/region

# Inclusion paths and libraries
#INCLUDES = -I $(CONTAINERS_DIR) -I $(COMMONS_DIR) -I $(REGION_DIR) -I $(LIBS_ROOT)/bioformats/vcf/
INCLUDES = -I . -I $(LIBS_ROOT) -I $(LIBS_ROOT)/common-libs -I $(LIBS_ROOT)/bioinfo-libs
LIBS = -lcprops -fopenmp -lcurl -Wl,-Bsymbolic-functions

# INCLUDES_STATIC = -I $(CONTAINERS_DIR) -I $(COMMONS_DIR) -I $(REGION_DIR) -I $(LIBS_ROOT)/bioformats/vcf/ -I $(LIBS_ROOT)/../include
INCLUDES_STATIC = -I . -I $(LIBS_ROOT) -I $(LIBS_ROOT)/../include -I $(LIBS_ROOT)/common-libs -I $(LIBS_ROOT)/bioinfo-libs
LIBS_STATIC = -L $(LIBS_ROOT) -lcprops -fopenmp -lcurl -Wl,-Bsymbolic-functions

# Source and object files
MAIN_OBJECTS = bed_file.c bed_read.c bed_reader.c bed_write.c bed_batch.c
MISC_OBJECTS = $(CONTAINERS_DIR)/list.o $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/http_utils.o $(COMMONS_DIR)/log.o $(COMMONS_DIR)/string_utils.o
REGION_TABLE_OBJECTS = $(CONTAINERS_DIR)/region_table.o $(REGION_DIR)/region.o
GFF_BATCH_OBJECTS = main_batch.c $(MAIN_OBJECTS) $(MISC_OBJECTS) $(REGION_TABLE_OBJECTS)


# Targets
all: compile-dependencies compile

compile: compile-dependencies
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -c $(MAIN_OBJECTS) $(INCLUDES) $(LIBS)

compile-static: compile-dependencies-static
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -c $(MAIN_OBJECTS) $(INCLUDES_STATIC) $(LIBS_STATIC)

ragel:
	ragel -o bed_reader.c -C -G2 bed.ragel

ragel_graphviz:
	ragel -o bed_ragel.dot -V -C -G2 bed.ragel

main_batch: compile-dependencies compile ragel $(GFF_BATCH_OBJECTS)
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o $@ $(GFF_BATCH_OBJECTS) $(INCLUDES) $(LIBS)
#	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -o $@ $(GFF_BATCH_OBJECTS) $(INCLUDES) $(LIBS)

compile-dependencies:
	cd $(COMMONS_DIR) && make compile && \
	cd $(CONTAINERS_DIR) && make compile  && \
	cd $(REGION_DIR) && make region.o

compile-dependencies-static:
	cd $(COMMONS_DIR) && make compile && \
	cd $(CONTAINERS_DIR) && make compile && \
	cd $(REGION_DIR) && make region.o

clean:
	rm -f *.o
	rm -f main_batch

