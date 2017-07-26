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
    //    fq_read_stats->N_out_quality = 0;
    
    fq_read_stats->num_A = 0;
    fq_read_stats->num_T = 0;
    fq_read_stats->num_C = 0;
    fq_read_stats->num_G = 0;
    fq_read_stats->num_N = 0;

    if (fq_read_stats->kmers_on) {
      fq_read_stats->kmers_on = 0;
      for(int i = 0; i < NUM_KMERS; i++) {
	fq_read_stats->kmers[i].counter = 0;
	fq_read_stats->kmers[i].counter_by_pos_size = 0;
      }
    }
  }
}

void fastq_read_stats_free(fastq_read_stats_t *fq_read_stats) {
  if(fq_read_stats != NULL) {
    if (fq_read_stats->kmers_on) {
      for(int i = 0; i < NUM_KMERS; i++) {
	if (fq_read_stats->kmers[i].counter_by_pos_size > 0) {
	  free(fq_read_stats->kmers[i].counter_by_pos);
	}
      }
    }
    free(fq_read_stats);
  }
}

void fastq_reads_stats_free(array_list_t *fq_reads_stats) {
  if(fq_reads_stats != NULL) {
    array_list_free(fq_reads_stats, (void *) fastq_read_stats_free);
  }
}


fastq_read_stats_options_t *fastq_read_stats_options_new(int kmers_on, int num_threads) {
  fastq_read_stats_options_t *read_stats_options = (fastq_read_stats_options_t *)malloc(sizeof(fastq_read_stats_options_t));
  
  read_stats_options->kmers_on = kmers_on;
  read_stats_options->num_threads = num_threads;
  
  return read_stats_options;
}

void fastq_read_stats_options_free(fastq_read_stats_options_t *read_stats_options) {
  if (read_stats_options != NULL) {
    free(read_stats_options);
  }
}


int fastq_read_stats(fastq_read_t *fq_read, 
		     fastq_read_stats_options_t *fq_read_stats_options, 
		     fastq_read_stats_t *fq_read_stats) {
  if (fq_read != NULL) {
    kmer_t *kmer;
    char *kmer_str;
    int kmer_index, kmers_on = fq_read_stats_options->kmers_on;
    
    int read_length = fq_read->length;

    fq_read_stats->length = read_length;
    fq_read_stats->kmers_on = kmers_on;

    for (size_t i = 0; i < read_length; i++) {
      if (kmers_on && (i <= read_length - 5)) {
	// kmers handling
	kmer_str = &fq_read->sequence[i];
	if ((kmer_str[0] != 78) && (kmer_str[1] != 78) && (kmer_str[2] != 78) && 
	    (kmer_str[3] != 78) && (kmer_str[4] != 78)) {
	  kmer_index = 0;
	  kmer_index += ((kmer_str[4]  >> 1) & 3);
	  kmer_index += (((kmer_str[3] >> 1) & 3) << 2);
	  kmer_index += (((kmer_str[2] >> 1) & 3) << 4);
	  kmer_index += (((kmer_str[1] >> 1) & 3) << 6);
	  kmer_index += (((kmer_str[0] >> 1) & 3) << 8);
	  
	  kmer = &fq_read_stats->kmers[kmer_index];
	  kmer->counter++;
	  kmer->id = kmer_index;
	  kmer->string[0] = kmer_str[0];
	  kmer->string[1] = kmer_str[1];
	  kmer->string[2] = kmer_str[2];
	  kmer->string[3] = kmer_str[3];
	  kmer->string[4] = kmer_str[4];
	  kmer->string[5] = '\0';
	  if (kmer->counter_by_pos_size == 0) {
	    kmer->counter_by_pos_size = read_length;
	    kmer->counter_by_pos = (size_t *) calloc(read_length, sizeof(size_t));
	  }
	  kmer->counter_by_pos[i]++;
	}
      }
      
      // nucleotide counters
      switch (fq_read->sequence[i]) {
      case 'A':	fq_read_stats->num_A++;
	break;
      case 'T':	fq_read_stats->num_T++;
	break;
      case 'C':	fq_read_stats->num_C++;
	break;
      case 'G':	fq_read_stats->num_G++;
	break;
      case 'N':	fq_read_stats->num_N++;
	break;
      default:	break;
      }
      
      // quality accumulator
      fq_read_stats->quality_average += fq_read->quality[i];
    }
    // mean quality
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
    kmer_t *kmer;
    char *kmer_str;
    int kmer_index, kmers_on = fq_read_stats_options->kmers_on;
    
    int read_length;

    size_t num_reads = array_list_size(fq_reads);
    for (size_t r = 0; r < num_reads; r++) {
      fq_read = array_list_get(r, fq_reads);
      if (fq_read != NULL) {
	read_length = fq_read->length;

	fq_read_stats = fastq_read_stats_new();
	fq_read_stats->length = read_length;
	fq_read_stats->kmers_on = kmers_on;

	for (size_t i = 0; i < read_length; i++) {
	  if (kmers_on && (i <= read_length - 5)) {
	    // kmers handling
	    kmer_str = &fq_read->sequence[i];
	    if ((kmer_str[0] != 78) && (kmer_str[1] != 78) && (kmer_str[2] != 78) && 
		(kmer_str[3] != 78) && (kmer_str[4] != 78)) {
	      kmer_index = 0;
	      kmer_index += ((kmer_str[4]  >> 1) & 3);
	      kmer_index += (((kmer_str[3] >> 1) & 3) << 2);
	      kmer_index += (((kmer_str[2] >> 1) & 3) << 4);
	      kmer_index += (((kmer_str[1] >> 1) & 3) << 6);
	      kmer_index += (((kmer_str[0] >> 1) & 3) << 8);

	      kmer = &fq_read_stats->kmers[kmer_index];
	      kmer->counter++;
	      kmer->id = kmer_index;
	      kmer->string[0] = kmer_str[0];
	      kmer->string[1] = kmer_str[1];
	      kmer->string[2] = kmer_str[2];
	      kmer->string[3] = kmer_str[3];
	      kmer->string[4] = kmer_str[4];
	      kmer->string[5] = '\0';
	      if (kmer->counter_by_pos_size == 0) {
		kmer->counter_by_pos_size = read_length;
		kmer->counter_by_pos = (size_t *) calloc(read_length, sizeof(size_t));
	      }
	      kmer->counter_by_pos[i]++;
	    }
	  }
	
	  // nucleotide counters
	  switch (fq_read->sequence[i]) {
	  case 'A':	fq_read_stats->num_A++;
	    break;
	  case 'T':	fq_read_stats->num_T++;
	    break;
	  case 'C':	fq_read_stats->num_C++;
	    break;
	  case 'G':	fq_read_stats->num_G++;
	    break;
	  case 'N':	fq_read_stats->num_N++;
	    break;
	  default:	break;
	  }

	  // quality accumulator
	  fq_read_stats->quality_average += fq_read->quality[i];
	}
	// mean quality
	fq_read_stats->quality_average /= (float) fq_read->length;

	// finally insert stats in the list
	array_list_insert((void *) fq_read_stats, reads_stats);
      }
    }
  }
  return reads_stats->size;
}

void fastq_read_stats_print(fastq_read_stats_t *fq_read_stats) {
  printf("Quality avg.: %f, A,T,C,G: %i,%i,%i,%i, N: %i\n", 
	 fq_read_stats->quality_average, fq_read_stats->num_A, fq_read_stats->num_T, 
	 fq_read_stats->num_C, fq_read_stats->num_G, fq_read_stats->num_N);
}

char* kmers_string(int index, char* kmer) {
    int rest, quotient;

    quotient = index;
    kmer[5] = '\0';

    for (int i = 0; i < 5; i++) {
        rest = quotient % 4;
        quotient /= 4;

        switch (rest) {
            case 0:
                kmer[4-i] = 'A';
                break;
            case 1:
                kmer[4-i] = 'C';
                break;
            case 2:
                kmer[4-i] = 'T';
                break;
            case 3:
                kmer[4-i] = 'G';
                break;
        }
    }

    return kmer;
}

int kmers_sort(const void* k1, const void* k2) {
    kmer_t* kmer1 = (kmer_t*) k1;
    kmer_t* kmer2 = (kmer_t*) k2;

    return (kmer2->counter - kmer1->counter);
}

