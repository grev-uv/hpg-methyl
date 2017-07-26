#include <stdlib.h>

#include "fastq_read.h"
#include "fastq_filter.h"
#include "fastq_stats.h"

//------------------------------------------------------------------------

fastq_filter_options_t *fastq_filter_options_new(int min_read_length, int max_read_length,
						 int min_read_quality, int max_read_quality,
						 int max_out_of_quality, int left_length,
						 int min_left_quality, int max_left_quality,
						 int right_length, int min_right_quality,
						 int max_right_quality, int max_N) {


  fastq_filter_options_t *b = (fastq_filter_options_t *) malloc(sizeof(fastq_filter_options_t));

  b->min_read_length = min_read_length;
  b->max_read_length = max_read_length;
  
  b->min_read_quality = min_read_quality;
  b->max_read_quality = max_read_quality;
  b->max_out_of_quality = max_out_of_quality;

  b->left_length = left_length;
  b->min_left_quality = min_left_quality;
  b->max_left_quality = max_left_quality;

  b->right_length = right_length;
  b->min_right_quality = min_right_quality;
  b->max_right_quality = max_right_quality;

  b->max_N = max_N;

  return b;
}

//------------------------------------------------------------------------

void fastq_filter_options_free(fastq_filter_options_t *b) {
  if (b) {
    free(b);
  }
}

//------------------------------------------------------------------------

array_list_t *fastq_filter(array_list_t *reads, array_list_t *passed, 
			   array_list_t *failed, fastq_filter_options_t *options) {
  char *sequence, *quality;
  fastq_read_t *read;

  
  int left_qual, acc_left_qual;
  int right_qual, acc_right_qual, right_start;
  int q, qual_on, qual, acc_qual;
  int out_of, N_on, num_N;

  int read_length, num_items = array_list_size(reads);

  //  #pragma omp parallel for schedule(dynamic, 100000)
  for (size_t i = 0; i < num_items; i++) {
    read = array_list_get(i, reads);
    sequence = read->sequence;
    quality = read->quality;
    read_length = read->length;

    // min. and max. read length
    //    printf("read_length = %i, (min, max) = (%i, %i)\n", read_length, options->min_read_length, options->max_read_length);
    if (read_length < options->min_read_length || read_length > options->max_read_length) {
      //      printf("\t...failed\n");
      array_list_insert(read, failed);
      continue;
    }

    right_start = read_length - options->right_length - 1;

    // compute basic statistics
    out_of = 0;
    acc_qual = 0; qual = 0;
    acc_left_qual = 0; left_qual = 0;
    acc_right_qual = 0; right_qual = 0;
    num_N = 0;
    
    for (size_t j = 0; j < read_length; j++) {
      switch (sequence[j]) {
      case 'N':	num_N++;
	break;
      default:	break;
      }
      q = quality[j];
      acc_qual += q;
      if (j < options->left_length) acc_left_qual += q;
      if (j > right_start) acc_right_qual += q;
      if (q < options->min_read_quality || q > options->max_read_quality) {
	out_of++;
      }
    }

    // min. and max. read quality
    qual = round(1.0f * acc_qual / read_length);
    //    printf("qual = %i, (min, max) = (%i, %i)\n", qual, options->min_read_quality, options->max_read_quality);
    if (qual < options->min_read_quality || qual > options->max_read_quality) {
      //      printf("\t...failed\n");
      array_list_insert(read, failed);
      continue;
    }

    // max. N
    //    printf("num_N = %i, (max) = (%i)\n", num_N, options->max_N);
    if (num_N > options->max_N) {
      //      printf("\t...failed\n");
      array_list_insert(read, failed);
      continue;
    }

    // max. nucleotides out of quality
    //    printf("out_of = %i, (max) = (%i)\n", out_of, options->max_out_of_quality);
    if (out_of > options->max_out_of_quality) {
      //      printf("\t...failed\n");
      array_list_insert(read, failed);
      continue;
    }

    // min. and max. left quality
    if (options->left_length) {
      left_qual = round(1.0f * acc_left_qual / options->left_length);
      //      printf("left_qual = %i, (min, max) = (%i, %i)\n", left_qual, options->min_left_quality, options->max_left_quality);
      if (left_qual < options->min_left_quality || left_qual > options->max_left_quality) {
	//	printf("\t...failed\n");
	array_list_insert(read, failed);
	continue;
      }
    }

    // min. and max. right quality
    if (options->right_length) {
      right_qual = round(1.0f * acc_right_qual / options->right_length);
      //      printf("right_qual = %i, (min, max) = (%i, %i)\n", right_qual, options->min_right_quality, options->max_right_quality);
      if (right_qual < options->min_right_quality || right_qual > options->max_right_quality) {
	//	printf("\t...failed\n");
	array_list_insert(read, failed);
	continue;
      }
    }

    // passed all filters
    //    printf("----> PASSED !!!!\n");
    array_list_insert(read, passed);
  }

  return passed;
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
