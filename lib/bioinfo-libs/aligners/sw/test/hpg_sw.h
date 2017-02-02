#ifndef HPG_SW_H
#define HPG_SW_H

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "macros.h"
#include "sse.h"
#include "smith_waterman.h"
//#include "avx.h"

void run_sse(char *q_filename, char *r_filename, 
	     float gap_open, float gap_extend, char *matrix_filename,
	     int batch_size, int num_threads, char *out_filename);

void display_usage(char *msg);
int file_exists(const char *filename);

#endif // HPG_SW_H
