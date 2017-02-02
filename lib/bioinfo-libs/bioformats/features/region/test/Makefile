
CC = gcc
CFLAGS = -std=c99 -O3

# Directories' path
ifdef LIBS_HOME
	LIBS_ROOT = $(LIBS_HOME)
else
	LIBS_ROOT = $(PWD)/../..
endif

SRC_DIR = ..
COMMONS_DIR = $(LIBS_ROOT)/common-libs/commons

INCLUDES = -I $(SRC_DIR) -I $(LIBS_ROOT)/common-libs -I $(COMMONS_DIR)
LIBS = -lcheck -lcprops 

TEST_REGION_OBJS = test_region.c $(SRC_DIR)/region.o $(COMMONS_DIR)/string_utils.o

all: string_utils.o region.o test_region

test_region: $(TEST_REGION_OBJS)
	$(CC) $(CFLAGS) -g -o $@ $(TEST_REGION_OBJS) $(INCLUDES) $(LIBS)

region.o:
	cd .. && make region.o

string_utils.o:
	cd $(COMMONS_DIR) && make string_utils.o

clean:
	rm test_region

