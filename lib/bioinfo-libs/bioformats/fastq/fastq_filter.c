#include <stdlib.h>

#include "fastq_read.h"
#include "fastq_filter.h"
#include "fastq_stats.h"

fastq_filter_options_t *fastq_filter_options_new(int min_length, int max_lentgh, float min_quality, 
						 float max_quality, int max_Ns, int max_N_out_quality) {
  fastq_filter_options_t *options = (fastq_filter_options_t *) malloc(sizeof(fastq_filter_options_t));
  
  options->min_length = min_length;
  options->max_length = max_lentgh;
  options->min_quality = min_quality;
  options->max_quality = max_quality;
  options->max_Ns = max_Ns;
  options->max_N_out_quality = max_N_out_quality;
  
  return options;
}

void fastq_filter_options_free(fastq_filter_options_t *options) {
  if (options != NULL) {
    free(options);
  }
}

array_list_t *fastq_filter(array_list_t *reads, array_list_t *passed, array_list_t *failed, 
			   fastq_filter_options_t *options) {
  fastq_read_t *read;
  
  fastq_read_stats_t *fq_read_stats = fastq_read_stats_new();
  fastq_read_stats_options_t *stats_options = fastq_read_stats_options_new(options->min_length, 
									   options->max_length, 4);
  
  #pragma omp parallel for schedule(dynamic, 100000)
  for (size_t i=0; i<reads->size; i++) {
    //		fastq_read_stats_init(fq_read_stats);
    read = array_list_get(i, reads);
    if (read->length >= options->min_length && read->length <= options->max_length) {
      fastq_read_stats(read, stats_options, fq_read_stats);
      //			fastq_read_stats_print(fq_read_stats);
      if (fq_read_stats->quality_average >= options->min_quality && 
	  fq_read_stats->quality_average <= options->max_quality && 
	  fq_read_stats->Ns < options->max_Ns) {
	array_list_insert(read, passed);
      } else {
	array_list_insert(read, failed);
      }
    } else {
      // get read stats
      array_list_insert(read, failed);
    }
  }
  fastq_read_stats_free(fq_read_stats);
  fastq_read_stats_options_free(stats_options);
  
  return passed;
}
