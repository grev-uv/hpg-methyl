CC = gcc
CFLAGS = -std=c99 -O3 -D_BSD_SOURCE -D_XOPEN_SOURCE=600

ifdef LIBS_HOME
	LIBS_ROOT = $(LIBS_HOME)
else
	LIBS_ROOT = $(PWD)/../../../..
endif

# SRC_DIR = $(PWD)
# COMMONS_DIR = $(LIBS_ROOT)/commons

INCLUDES = -I $(PWD) -I $(LIBS_ROOT)/common-libs -I $(LIBS_ROOT)/bioinfo-libs -I $(LIBS_ROOT)
LIBS = -lcprops -lcurl -Wl,-Bsymbolic-functions
LIBS_TEST = -lcheck

#MISC_FILES = $(LIBS_ROOT)/common-libs/commons/string_utils.o

all: dna_map_region.o region.o region_table.o region_table_utils.o rna_map_region.o

dna_map_region.o: compile-dependencies dna_map_region.c
	$(CC) $(CFLAGS) -c dna_map_region.c $(MISC_FILES) $(INCLUDES) $(LIBS)

region.o: compile-dependencies region.c
	$(CC) $(CFLAGS) -c region.c $(MISC_FILES) $(INCLUDES) $(LIBS)

region_table.o: compile-dependencies region_table.c
	$(CC) $(CFLAGS) -c region_table.c $(MISC_FILES) $(INCLUDES) $(LIBS)

region_table_utils.o: compile-dependencies region_table_utils.c
	$(CC) $(CFLAGS) -c region_table_utils.c $(MISC_FILES) $(INCLUDES) $(LIBS)

rna_map_region.o: compile-dependencies rna_map_region.c
	$(CC) $(CFLAGS) -c rna_map_region.c $(MISC_FILES) $(INCLUDES) $(LIBS)

compile-dependencies:
	cd $(LIBS_ROOT)/common-libs/commons && make compile

clean:
	rm -f dna_map_region.o
	rm -f region.o
	rm -f region_table.o
	rm -f region_table_utils.o
	rm -f rna_map_region.o
