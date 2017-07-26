CC = gcc
CFLAGS = -static -std=c99 -O3 -D_XOPEN_SOURCE=600
CFLAGS_DEBUG = -static -std=c99 -g -D_XOPEN_SOURCE=600

SRC_DIR = ..
ROOT = ~/appl/bioinfo-c

INCLUDES = -I $(SRC_DIR)
LIBS = -lgsl -lcheck -lpthread -lm

OBJECTS = $(SRC_DIR)/array_utils.o $(SRC_DIR)/p_adjust.o $(SRC_DIR)/R_utils.o
# $(SRC_DIR)/math_utils.


all: test-array-utils

test-array-utils: array_utils.o
	$(CC) $(CFLAGS) -o $@ test_array_utils.c $(OBJECTS) $(INCLUDES) $(LIBS)

array_utils.o:
	cd $(SRC_DIR) && make array_utils.o

p_adjust.o: p_adjust.h p_adjust.c
	cd $(SRC_DIR) && make p_adjust.o

R_utils.o: R_utils.h R_utils.c
	cd $(SRC_DIR) && make R_utils.o

#ls_fit.o:
#	$(CC) $(CFLAGS) -c ls_fit.c $(INCLUDES) $(LIBS)

#math_utils.o:
#	cd $(SRC_DIR) && make math_utils.o

#matrix.o:
#	$(CC) $(CFLAGS) -c matrix.c $(INCLUDES) $(LIBS)

#stats.o:
#	$(CC) $(CFLAGS) -c stats.c$(INCLUDES) $(LIBS)

clean:
	rm -f *.o
	rm -f $(SRC_DIR)/*.o
