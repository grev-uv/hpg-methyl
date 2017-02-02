#ifndef FASTQ_STATS_H
#define FASTQ_STATS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fastq_read.h"
#include "containers/array_list.h"


/**
* @brief Fastq read statistics
*
* Structure for storing statistics from a fastq read
*/
typedef struct fastq_read_stats {
    int length;
    int Ns;
    int N_out_quality;

    int num_A;
    int num_T;
    int num_C;
    int num_G;

    float quality_average;	/**< Mean quality average. */
} fastq_read_stats_t;

typedef struct fastq_read_stats_options {
	int min_qual;
	int max_qual;

	int num_threads;

} fastq_read_stats_options_t;


fastq_read_stats_t *fastq_read_stats_new();

void fastq_read_stats_init(fastq_read_stats_t *fq_read_stats);

void fastq_read_stats_free(fastq_read_stats_t *read_stats);

void fastq_reads_stats_free(array_list_t *read_stats);


fastq_read_stats_options_t *fastq_read_stats_options_new(int min_qual, int max_qual, int num_threads);

void fastq_read_stats_options_free(fastq_read_stats_options_t *read_stats);


int fastq_read_stats(fastq_read_t *fq_read, fastq_read_stats_options_t *fq_read_stats_options, fastq_read_stats_t *read_stats);

int fastq_reads_stats(array_list_t *fq_reads, fastq_read_stats_options_t *fq_read_stats_options, array_list_t *reads_stats);

void fastq_read_stats_print(fastq_read_stats_t *fq_read_stats);



//size_t fastq_stats_se(array_list_t *fq_reads, array_list_t *reads_stats);

//size_t fastq_stats_pe(array_list_t *fq_reads, array_list_t *reads_stats);

//int fastq_stats(fastq_read_t *fq_read, fastq_read_stats_t *read_stats);


#endif	/*  FASTQ_STATS_H  */
