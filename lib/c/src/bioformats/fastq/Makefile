#CC = /opt/intel/bin/icc
#CC = gcc
ifeq ($(compiler), intel)
        CC = /opt/intel/bin/icc
else
        CC = gcc
        compiler=gcc
endif

CFLAGS = -std=c99 -O3 -D_GNU_SOURCE -g
CFLAGS_DEBUG = -std=c99 -g -D_GNU_SOURCE

# Directories' paths
ifdef LIBS_HOME
	LIBS_ROOT = $(LIBS_HOME)
else
	LIBS_ROOT = $(PWD)/../../..
endif

COMMONS_DIR = $(LIBS_ROOT)/common-libs/commons
CONTAINERS_DIR = $(LIBS_ROOT)/common-libs/containers

# Include path and libraries
INCLUDES = -I . -I $(LIBS_ROOT) -I $(LIBS_ROOT)/common-libs -I $(LIBS_ROOT)/bioinfo-libs
LIBS = -L/usr/lib/x86_64-linux-gnu -lconfig -fopenmp -Wl,-Bsymbolic-functions -lxml2 -lcurl -lcprops

# Source and object files
MAIN_OBJECTS = fastq_batch.c fastq_batch_reader.c fastq_file.c fastq_gzfile.c fastq_filter.c fastq_read.c fastq_stats.c
#MISC_OBJECTS = $(CONTAINERS_DIR)/array_list.o $(CONTAINERS_DIR)/list.o $(COMMONS_DIR)/http_utils.o $(COMMONS_DIR)/string_utils.o $(COMMONS_DIR)/file_utils.o $(COMMONS_DIR)/log.o
MAIN_TEST = main.c $(MAIN_OBJECTS)
# Targets
all: compile-dependencies compile

compile: compile-dependencies
	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600 -c $(MAIN_OBJECTS) $(INCLUDES) $(LIBS)

compile-dependencies:
	cd $(CONTAINERS_DIR) && make compiler=$(compiler) && \
	cd $(COMMONS_DIR) && make compiler=$(compiler)

main: compile-dependencies $(MAIN_TEST)
	$(CC) $(CFLAGS_DEBUG)  -c main.c $(MAIN_OBJECTS) $(INCLUDES) $(LIBS)
	$(CC) $(CFLAGS_DEBUG) -o $@ *.o $(INCLUDES) $(COMMONS_DIR)/*.o $(CONTAINERS_DIR)/array_list.o $(CONTAINERS_DIR)/list.o $(LIBS)
#	$(CC) $(CFLAGS) -D_XOPEN_SOURCE=600  $(MAIN_TEST) $(INCLUDES) $(LIBS) -o main-test
	
clean:
	rm -f *.o
