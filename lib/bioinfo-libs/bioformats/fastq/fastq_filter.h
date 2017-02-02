/*
 * fastq_filter.h
 *
 *  Created on: Sep 27, 2012
 *      Author: imedina
 */

#ifndef FASTQ_FILTER_H
#define FASTQ_FILTER_H

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
	int min_length;
	int max_length;

	float min_quality;
	float max_quality;

	int max_Ns;
	int max_N_out_quality;

} fastq_filter_options_t;


fastq_filter_options_t *fastq_filter_options_new(int min_length, int max_lentgh, float min_quality, float max_quality, int max_Ns, int max_N_out_quality);

void fastq_filter_options_free(fastq_filter_options_t *options);

array_list_t *fastq_filter(array_list_t *reads, array_list_t *passed, array_list_t *failed, fastq_filter_options_t *filter_options);


#endif /* FASTQ_FILTER_H_ */
