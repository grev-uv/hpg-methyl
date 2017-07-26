BIN = ../bin
LIB = .
SRC = ..

CC = gcc
CFLAGS = -Wall -O3 -std=c99
#CFLAGS_DEBUG = -Wall -pg

INCLUDES = -I$(SRC) -I/usr/include/libxml2
LIBS = -lcurl -Wl,-Bsymbolic-functions -lcheck -lxml2

all: test-file-utils test-http-utils test-log test-string-utils test-result


test-file-utils: 
	$(CC) $(CFLAGS) -o $@ test_file_utils.c $(SRC)/string_utils.o $(SRC)/file_utils.o $(INCLUDES) 

test-http-utils:
	$(CC) $(CFLAGS) -o $@ test_http_utils.c $(SRC)/http_utils.o $(INCLUDES) $(LIBS)

test-log:
	$(CC) $(CFLAGS) -o $@ test_log.c $(SRC)/string_utils.o $(SRC)/log.o $(INCLUDES) 

test-string-utils:
	$(CC) $(CFLAGS) -o $@ test_string_utils.c $(SRC)/string_utils.o $(INCLUDES) 

test-result: result_file_test.c ../result.h ../result.c
	$(CC) $(CFLAGS) -o $@ result_file_test.c $(SRC)/result.o $(INCLUDES) $(LIBS)

clean:
	rm -f test-file-utils
	rm -f test-http-utils
	rm -f test-log
	rm -f test-string-utils
	rm -f test-result