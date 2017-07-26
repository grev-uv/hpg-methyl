CC = gcc
CFLAGS = -static -std=c99 -O3 -Wall
CFLAGS_DEBUG = -static -std=c99 -g

LIBS_DIR = $(PWD)/..
COMMONS_LIB = $(LIBS_DIR)/common-libs/commons

INCLUDES = -I$(PWD) -I$(PWD)/algebra -I$(PWD)/data -I$(PWD)/stats -I$(COMMONS_LIB)
LIBS = -lgsl

OBJECTS = anova.o anova_result.o array_utils.o correlation.o correlation_matrix.o cox.o fisher.o fitting.o math_utils.o matrix.o p_adjust.o statistics.o stats_inference.o R_utils.o ttest.o ttest_result.o


all: $(OBJECTS)

compile:
	$(CC) $(CFLAGS) -c *.c algebra/*.c data/*.c stats/*.c $(INCLUDES) $(LIBS)

algebra:
	$(CC) $(CFLAGS) -c algebra/*.c $(INCLUDES) $(LIBS)

data:
	$(CC) $(CFLAGS) -c data/*.c $(INCLUDES) $(LIBS)

stats:
	$(CC) $(CFLAGS) -c stats/*.c $(INCLUDES) $(LIBS)

anova.o:
	$(CC) $(CFLAGS) -c stats/anova.c $(INCLUDES) $(LIBS)

anova_result.o:
	$(CC) $(CFLAGS) -c stats/anova_result.c $(INCLUDES) $(LIBS)

array_utils.o:
	$(CC) $(CFLAGS) -c data/array_utils.c $(INCLUDES) $(LIBS)

correlation.o:
	$(CC) $(CFLAGS) -c stats/correlation.c $(INCLUDES) $(LIBS)

correlation_matrix.o:
	$(CC) $(CFLAGS) -c stats/correlation_matrix.c $(INCLUDES) $(LIBS)

cox.o:
	$(CC) $(CFLAGS) -c stats/cox.c $(INCLUDES) $(LIBS)

fitting.o:
	$(CC) $(CFLAGS) -c algebra/fitting.c $(INCLUDES) $(LIBS)

fisher.o:
	$(CC) $(CFLAGS) -c stats/fisher.c $(INCLUDES) $(LIBS)

math_utils.o:
	$(CC) $(CFLAGS) -c math_utils.h $(INCLUDES) $(LIBS)

matrix.o:
	$(CC) $(CFLAGS) -c data/matrix.c $(INCLUDES) $(LIBS)

p_adjust.o:
	$(CC) $(CFLAGS) -c stats/p_adjust.c $(INCLUDES) $(LIBS)

R_utils.o:
	$(CC) $(CFLAGS) -c data/R_utils.c $(INCLUDES) $(LIBS)

statistics.o:
	$(CC) $(CFLAGS) -c stats/statistics.c $(INCLUDES) $(LIBS)

stats_inference.o:
	$(CC) $(CFLAGS) -c stats/stats_inference.c $(INCLUDES) $(LIBS)

ttest.o:
	$(CC) $(CFLAGS) -c stats/ttest.c $(INCLUDES) $(LIBS)

ttest_result.o:
	$(CC) $(CFLAGS) -c stats/ttest_result.c $(INCLUDES) $(LIBS)


clean:
	rm -f *.o
