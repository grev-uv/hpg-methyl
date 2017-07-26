/*
 * fastq_edit.h
 *
 *  Created on: May 23, 2013
 *      Author: jtarraga
 *              
 */

#ifndef FASTQ_EDIT_H
#define FASTQ_EDIT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "containers/array_list.h"

#include "fastq_filter.h"

//------------------------------------------------------------------------

#define PHRED_33_TO_64    1
#define PHRED_64_TO_33    2

//------------------------------------------------------------------------

/**
 *
 */
typedef struct fastq_edit_options {
  // to left trim
  int left_length;
  int min_left_quality;
  int max_left_quality;

  // to right trim
  int right_length;
  int min_right_quality;
  int max_right_quality;

  // convert to N those nucleotides with quality
  int min_N_quality;
  int max_N_quality;

  // convert quality
  int convert_quality;
} fastq_edit_options_t;


fastq_edit_options_t *fastq_edit_options_new(int left_length, int min_left_quality, 
					     int max_left_quality, int right_length, 
					     int min_right_quality, int max_right_quality, 
					     int min_N_quality, int max_N_quality,
					     int convert_quality);


void fastq_edit_options_free(fastq_edit_options_t *options);

int fastq_edit(array_list_t *reads, fastq_edit_options_t *options);

#ifdef __cplusplus
}
#endif

#endif /* FASTQ_EDIT_H_ */
