
#ifndef COMMONS_H
#define COMMONS_H

#include <stdio.h>
#include <sys/time.h>

#include "log.h"

#define PHRED33 	33
#define PHRED64		64


//====================================================================================
//  commons.h
//
//  commons structures and prototypes
//====================================================================================

extern int number_of_batchs;

extern int time_flag;

extern double read_time;
extern struct timeval t1_read, t2_read;

extern double gpu_time;
extern struct timeval t1_gpu, t2_gpu;

extern double cpu_time;
extern struct timeval t1_cpu, t2_cpu;

extern double kmers_time;
extern struct timeval t1_kmers, t2_kmers;

extern double result_time;
extern struct timeval t1_result, t2_result;

extern double reporting_time;
extern struct timeval t1_reporting, t2_reporting;

extern double write_time;
extern struct timeval t1_write, t2_write;

extern double mean_reads_per_batch;
extern double mean_batch_size;

//------------------------------------------------------------------------------------

int count_file_lines(char* filename, int max_line_length);
char* trim(char* input);
int is_a_number(char* string);

//------------------------------------------------------------------------------------

#define timevars() struct timeval t1, t2;

#define min(x,y)	(((x) < (y)) ? (x) : (y))

#define max(x,y)	(((x) > (y)) ? (x) : (y))

#define tic(msg)				\
  printf("--->> " msg " \n");			\
  fflush(stdout);				\
  gettimeofday(&t1, NULL);

#define toc()								                       \
  gettimeofday(&t2, NULL);						                       \
  printf("<<--- completo en %.0f usecs\n", (t2.tv_sec-t1.tv_sec)*1e6+(t2.tv_usec-t1.tv_usec)); \
  fflush(stdout);

#define start_timer(start_time)				\
  gettimeofday(&start_time, NULL);

#define stop_timer(start_time, end_time, timer)					                       \
  gettimeofday(&end_time, NULL);						                       \
  timer = timer + (end_time.tv_sec-start_time.tv_sec)*1e6 + (end_time.tv_usec-start_time.tv_usec);

#endif /* COMMONS_H */
