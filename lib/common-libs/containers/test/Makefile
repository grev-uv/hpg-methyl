CC = gcc
CFLAGS = -std=c99 -O3
CFLAGS_DEBUG = -std=c99 -g

ROOT = ../..
SRC_DIR = ..
COMMONS_DIR = $(ROOT)/commons

INCLUDES = -I $(SRC_DIR) -I $(ROOT)
LIBS = -lm -lcheck -lcprops

TEST_ARRAY_LIST_OBJS = test_array_list.c $(SRC_DIR)/array_list.o $(COMMONS_DIR)/string_utils.o $(COMMONS_DIR)/log.o

all: test_array_list

test_array_list: $(TEST_ARRAY_LIST_OBJS)
	$(CC) $(CFLAGS_DEBUG) -g -D_BSD_SOURCE -o $@ $(TEST_ARRAY_LIST_OBJS) $(INCLUDES) $(LIBS)

array_list.o:
	cd .. && make array_list.o

string_utils.o:
	cd .. && make string_utils.o

log.o:
	cd .. && make log.o

clean:
	test ! -f test_array_list || rm test_array_list

