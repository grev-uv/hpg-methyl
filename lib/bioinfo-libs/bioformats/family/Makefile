CC = gcc
CFLAGS = -std=c99 -O3 -D_GNU_SOURCE -D_XOPEN_SOURCE=600
CFLAGS_DEBUG = -std=c99 -g -D_GNU_SOURCE -D_XOPEN_SOURCE=600

# Directories' paths
ifdef LIBS_HOME
	LIBS_ROOT = $(LIBS_HOME)
else
	LIBS_ROOT = $(PWD)/../../..
endif

COMMONS_DIR = $(LIBS_ROOT)/common-libs
BIOINFO_DIR = $(LIBS_ROOT)/bioinfo-libs

# Include path and libraries
INCLUDES = -I . -I $(COMMONS_DIR) -I $(BIOINFO_DIR)
LIBS = -fopenmp -Wl,-Bsymbolic-functions

# Source and object files
MAIN_OBJECTS = family.c checks_family.c 
#MISC_OBJECTS = $(COMMONS_DIR)/commons/string_utils.o $(COMMONS_DIR)/commons/log.o

# Targets
#all: compile-dependencies compile
all: compile

#compile: compile-dependencies
compile:
	$(CC) $(CFLAGS) -c $(MAIN_OBJECTS) $(MISC_OBJECTS) $(INCLUDES) $(LIBS)

compile-dependencies:
	cd $(COMMONS_DIR) && make

clean:
	rm -f *.o
