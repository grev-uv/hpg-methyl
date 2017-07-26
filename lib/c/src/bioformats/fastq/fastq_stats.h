#ifndef FASTQ_STATS_H
#define FASTQ_STATS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fastq_read.h"
#include "containers/array_list.h"


//------------------------------------------------------------------------

#define QUALITY_PHRED33_VALUE  33
#define QUALITY_PHRED33_NAME  "phred33"

#define QUALITY_PHRED64_VALUE  64
#define QUALITY_PHRED64_NAME  "phred64"

//------------------------------------------------------------------------



#define NUM_KMERS    1024

//------------------------------------------------------------------------

typedef struct kmer {
  int id;		      /**< Kmer id */
  char string[6];	      /**< String representation of the kmer (e.g.: AATCG) */
  size_t counter;	      /**< kmer total counter */

  size_t counter_by_pos_size;
  size_t *counter_by_pos;     /**< kmer counter by nucleotide position */
} kmer_t;

int kmers_sort(const void* k1, const void* k2);
char* kmers_string(int index, char* kmer);

//------------------------------------------------------------------------

typedef struct fastq_read_stats_options {
  //  int min_qual;
  //  int max_qual;
  int kmers_on;
  int num_threads;
} fastq_read_stats_options_t;

//------------------------------------------------------------------------

/**
* @brief Fastq read statistics
*
* Structure for storing statistics from a fastq read
*/
typedef struct fastq_read_stats {
  int length;	         /**< Read length. */

  int num_A;
  int num_T;
  int num_C;
  int num_G;
  int num_N;
  //int N_out_quality;

  float quality_average; /**< Mean quality average. */

  int kmers_on;
  kmer_t kmers[NUM_KMERS]; /**< kmers information. */
} fastq_read_stats_t;

fastq_read_stats_t *fastq_read_stats_new();

void fastq_read_stats_init(fastq_read_stats_t *fq_read_stats);

void fastq_read_stats_free(fastq_read_stats_t *read_stats);

void fastq_reads_stats_free(array_list_t *read_stats);


fastq_read_stats_options_t *fastq_read_stats_options_new(int kmers, int num_threads);

void fastq_read_stats_options_free(fastq_read_stats_options_t *read_stats);


int fastq_read_stats(fastq_read_t *fq_read, fastq_read_stats_options_t *fq_read_stats_options, fastq_read_stats_t *read_stats);

int fastq_reads_stats(array_list_t *fq_reads, fastq_read_stats_options_t *fq_read_stats_options, array_list_t *reads_stats);

void fastq_read_stats_print(fastq_read_stats_t *fq_read_stats);

#ifdef __cplusplus
}
#endif

#endif	/*  FASTQ_STATS_H  */
