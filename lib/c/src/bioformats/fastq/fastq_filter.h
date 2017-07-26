/*
 * fastq_filter.h
 *
 *  Created on: Sep 27, 2012
 *      Author: imedina
 */

#ifndef FASTQ_FILTER_H
#define FASTQ_FILTER_H

#ifdef __cplusplus
extern "C" {
#endif

#include "containers/array_list.h"

#include "fastq_filter.h"

/**

 * @brief A filter selects a subcollection of reads which fulfill some condition.
 * @details A filter selects a subcollection of reads which fulfill some condition.
 * It is mandatory to provide the array list of reads to filter and a array list to store in
 * the reads that failed the filter's test.
 *
 */
typedef struct fastq_filter_options {
  int min_read_length;
  int max_read_length;

  int min_read_quality;
  int max_read_quality;
  int max_out_of_quality;

  int left_length;
  int min_left_quality;
  int max_left_quality;

  int right_length;
  int min_right_quality;
  int max_right_quality;

  int max_N;
} fastq_filter_options_t;


fastq_filter_options_t *fastq_filter_options_new(int min_read_length, int max_read_length,
						 int min_read_quality, int max_read_quality,
						 int max_out_of_quality, int left_length,
						 int min_left_quality, int max_left_quality,
						 int right_length, int min_right_quality,
						 int max_right_quality, int max_N);


void fastq_filter_options_free(fastq_filter_options_t *options);

array_list_t *fastq_filter(array_list_t *reads, array_list_t *passed, 
			   array_list_t *failed, fastq_filter_options_t *filter_options);

#ifdef __cplusplus
}
#endif

#endif /* FASTQ_FILTER_H_ */
