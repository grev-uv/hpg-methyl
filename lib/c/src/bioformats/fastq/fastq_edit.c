#include <stdlib.h>

#include "fastq_read.h"
#include "fastq_edit.h"

//------------------------------------------------------------------------

fastq_edit_options_t *fastq_edit_options_new(int left_length, int min_left_quality, 
					     int max_left_quality, int right_length, 
					     int min_right_quality, int max_right_quality, 
					     int min_N_quality, int max_N_quality,
					     int convert_quality) {

  fastq_edit_options_t *b = (fastq_edit_options_t *) malloc(sizeof(fastq_edit_options_t));


  b->left_length = left_length;
  b->min_left_quality = min_left_quality;
  b->max_left_quality = max_left_quality;

  b->right_length = right_length;
  b->min_right_quality = min_right_quality;
  b->max_right_quality = max_right_quality;

  b->min_N_quality = min_N_quality;
  b->max_N_quality = max_N_quality;

  b->convert_quality = convert_quality;

  return b;
}

//------------------------------------------------------------------------

void fastq_edit_options_free(fastq_edit_options_t *b) {
  if (b) {
    free(b);
  }
}

//------------------------------------------------------------------------

int fastq_edit(array_list_t *reads, fastq_edit_options_t *options) {
  char *sequence, *quality;
  fastq_read_t *read;

  size_t num_edited = 0;

  int left_qual, acc_left_qual;
  int right_qual, acc_right_qual, right_start;
  int len, start_trim, end_trim, k;
    /*qual_on, qual, acc_qual;

  int out_of, N_on, num_N;
  */

  int read_length, num_items = array_list_size(reads);

  //  #pragma omp parallel for schedule(dynamic, 100000)
  for (size_t i = 0; i < num_items; i++) {
    read = array_list_get(i, reads);
    sequence = read->sequence;
    quality = read->quality;
    read_length = read->length;

    // init
    start_trim = 0;
    end_trim = read_length - 1;

    // trimming left ?
    len = options->left_length;
    if (len > 0) {
      acc_left_qual = 0;
      for (size_t j = 0; j < len; j++) {
	acc_left_qual += quality[j];
      }
      left_qual = round(1.0f * acc_left_qual / len);
      if (left_qual < options->min_left_quality || left_qual > options->max_left_quality) {
	start_trim = len;
      }
    }

    // trimming right ?
    len = options->right_length;
    if (len > 0) {
      acc_right_qual = 0;
      end_trim = read_length - len - 1;
      for (size_t j = end_trim; j < read_length; j++) {
	acc_right_qual += quality[j];
      }
      right_qual = round(1.0f * acc_right_qual / len);
      if (right_qual >= options->min_right_quality && right_qual <= options->max_right_quality) {
	end_trim = read_length - 1;
      }
    }

    if (start_trim != 0 || end_trim != read_length - 1) {

      // we must trim
      num_edited++;
      if (start_trim != 0) {
	k = 0;
	for (size_t j = start_trim; j <= end_trim; j++) {
	  sequence[k] = sequence[j];
	  quality[k] = quality[j];
	  k++;
	}
	read->length = k;
      } else {
	read->length = end_trim;
      }
      sequence[read->length] = 0;
      quality[read->length] = 0;
    }
  }

  return num_edited;
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
