CC = gcc
CFLAGS = -std=c99 -O3 -D_XOPEN_SOURCE=600
CFLAGS_DEBUG = -std=c99 -g -D_XOPEN_SOURCE=600

# Directories' path
ifdef LIBS_HOME
	LIBS_ROOT = $(LIBS_HOME)
else
	LIBS_ROOT = $(PWD)/../../..
endif

CONTAINERS_DIR = $(LIBS_ROOT)/common-libs/containers
COMMONS_DIR = $(LIBS_ROOT)/common-libs/commons
BIOFORMATS_DIR = $(LIBS_ROOT)/bioinfo-libs/bioformats

# Inclusion paths and libraries
# INCLUDES = -I $(CONTAINERS_DIR) -I $(COMMONS_DIR) -I $(BIOFORMATS_DIR) -I $(BIOFORMATS_DIR)/vcf/
INCLUDES = -I . -I $(LIBS_ROOT) -I $(LIBS_ROOT)/common-libs -I $(LIBS_ROOT)/bioinfo-libs
LIBS = -lcprops -fopenmp

# INCLUDES_STATIC = -I $(CONTAINERS_DIR) -I $(COMMONS_DIR) -I $(BIOFORMATS_DIR) -I $(LIBS_ROOT)/bioformats/vcf/ -I $(LIBS_ROOT)/../include
INCLUDES_STATIC = -I . -I $(LIBS_ROOT) -I $(LIBS_ROOT)/../include -I $(LIBS_ROOT)/common-libs -I $(LIBS_ROOT)/bioinfo-libs
LIBS_STATIC = -L $(LIBS_ROOT) -lcprops -fopenmp

# Source and object files
MAIN_OBJECTS = ped_file.c ped_read.c ped_reader.c ped_write.c ped_batch.c ped_error.c ped_util.c
MISC_OBJECTS = $(CONTAINERS_DIR)/list.o $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/log.o $(COMMONS_DIR)/string_utils.o \
$(BIOFORMATS_DIR)/family/family.o $(BIOFORMATS_DIR)/family/checks_family.o
PED_BATCH_OBJECTS = main_batch.c $(MAIN_OBJECTS) $(MISC_OBJECTS)


# Targets
all: compile-dependencies compile

compile: compile-dependencies
	$(CC) $(CFLAGS_DEBUG) -c $(MAIN_OBJECTS) $(INCLUDES) $(LIBS)

compile-static: compile-dependencies-static
	$(CC) $(CFLAGS) -c $(MAIN_OBJECTS) $(INCLUDES_STATIC) $(LIBS_STATIC)

ragel:
	ragel -o ped_reader.c -C -G2 ped.ragel

ragel_graphviz:
	ragel -o ped_ragel.dot -V -C -G2 ped.ragel

main_batch: compile-dependencies compile ragel $(PED_BATCH_OBJECTS)
#	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -o $@ $(PED_BATCH_OBJECTS) $(INCLUDES) $(LIBS)
	$(CC) $(CFLAGS_DEBUG) -D_XOPEN_SOURCE=600 -o $@ $(PED_BATCH_OBJECTS) $(INCLUDES) $(LIBS)

compile-dependencies:
	cd $(COMMONS_DIR) && make compile && \
	cd $(CONTAINERS_DIR) && make compile && \
	cd $(BIOFORMATS_DIR) && make family

compile-dependencies-static:
	cd $(COMMONS_DIR) && make compile && \
	cd $(CONTAINERS_DIR) && make compile && \
	cd $(BIOFORMATS_DIR) && make family

clean:
	rm -f *.o
	rm main_batch

