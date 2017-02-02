#include "fastq_stats.h"

fastq_read_stats_t *fastq_read_stats_new() {
  fastq_read_stats_t *fq_read_stats = (fastq_read_stats_t *)malloc(sizeof(fastq_read_stats_t));

  fastq_read_stats_init(fq_read_stats);
  
  return fq_read_stats;
}

void fastq_read_stats_init(fastq_read_stats_t *fq_read_stats) {
  if(fq_read_stats != NULL) {
    fq_read_stats->length = 0;
    fq_read_stats->quality_average = 0.0;
    fq_read_stats->Ns = 0;
    fq_read_stats->N_out_quality = 0;
    
    fq_read_stats->num_A = 0;
    fq_read_stats->num_T = 0;
    fq_read_stats->num_C = 0;
    fq_read_stats->num_G = 0;
	}
}

void fastq_read_stats_free(fastq_read_stats_t *fq_read_stats) {
  if(fq_read_stats != NULL) {
    free(fq_read_stats);
  }
}

void fastq_reads_stats_free(array_list_t *fq_reads_stats) {
  if(fq_reads_stats != NULL) {
    array_list_free(fq_reads_stats, free);
  }
}


fastq_read_stats_options_t *fastq_read_stats_options_new(int min_qual, int max_qual, int num_threads) {
  fastq_read_stats_options_t *read_stats_options = (fastq_read_stats_options_t *)malloc(sizeof(fastq_read_stats_options_t));
  
  read_stats_options->min_qual = min_qual;
  read_stats_options->max_qual = max_qual;
  read_stats_options->num_threads = num_threads;
  
  return read_stats_options;
}

void fastq_read_stats_options_free(fastq_read_stats_options_t *read_stats_options) {
	if(read_stats_options != NULL) {
		free(read_stats_options);
	}
}


int fastq_read_stats(fastq_read_t *fq_read, 
		     fastq_read_stats_options_t *fq_read_stats_options, 
		     fastq_read_stats_t *fq_read_stats) {
  if (fq_read != NULL) {
    fastq_read_stats_init(fq_read_stats);
    fq_read_stats->length = fq_read->length;
    for (size_t i=0; i<fq_read->length; i++) {
      if (fq_read->sequence[i] == 'N' || fq_read->sequence[i] == 'n') {
	fq_read_stats->Ns++;
      }
      if (fq_read->quality[i] < fq_read_stats_options->min_qual || 
	  fq_read->quality[i] > fq_read_stats_options->max_qual) {
	fq_read_stats->N_out_quality++;
      }
      switch (fq_read->sequence[i]) {
      case 'A':	fq_read_stats->num_A++;
	break;
      case 'T':	fq_read_stats->num_T++;
	break;
      case 'C':	fq_read_stats->num_C++;
	break;
      case 'G':	fq_read_stats->num_G++;
	break;
      default:	break;
      }
      fq_read_stats->quality_average += fq_read->quality[i];
    }
    fq_read_stats->quality_average /= (float) fq_read->length;
    return 0;
  }
  return 1;
}

int fastq_reads_stats(array_list_t *fq_reads, 
		      fastq_read_stats_options_t *fq_read_stats_options, 
		      array_list_t *reads_stats) {

  fastq_read_t *fq_read;
  fastq_read_stats_t *fq_read_stats;

  if (fq_reads != NULL) {

    size_t num_reads = array_list_size(fq_reads);
    for (size_t r = 0; r < num_reads; r++) {
      fq_read = array_list_get(r, fq_reads);
      if (fq_read != NULL) {
	fq_read_stats = fastq_read_stats_new();
	fq_read_stats->length = fq_read->length;
	for (size_t i = 0; i < fq_read->length; i++) {
	  if (fq_read->sequence[i] == 'N' || fq_read->sequence[i] == 'n') {
	    fq_read_stats->Ns++;
	  }
	  if (fq_read->quality[i] < fq_read_stats_options->min_qual || 
	      fq_read->quality[i] > fq_read_stats_options->max_qual) {
	    fq_read_stats->N_out_quality++;
	  }
	  switch (fq_read->sequence[i]) {
	  case 'A':	fq_read_stats->num_A++;
	    break;
	  case 'T':	fq_read_stats->num_T++;
	    break;
	  case 'C':	fq_read_stats->num_C++;
	    break;
	  case 'G':	fq_read_stats->num_G++;
	    break;
	  default:	break;
	  }
	  fq_read_stats->quality_average += fq_read->quality[i];
	}
	fq_read_stats->quality_average /= (float) fq_read->length;
	array_list_insert((void *) fq_read_stats, reads_stats);
      }
    }
  }
  return reads_stats->size;
}

void fastq_read_stats_print(fastq_read_stats_t *fq_read_stats) {
  printf("Quality avg.: %f, A,T,C,G: %i,%i,%i,%i, Ns: %i\n", 
	 fq_read_stats->quality_average, fq_read_stats->num_A, fq_read_stats->num_T, 
	 fq_read_stats->num_C, fq_read_stats->num_G, fq_read_stats->Ns);
}

