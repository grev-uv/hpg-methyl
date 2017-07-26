#ifndef BAM_STATS_H
#define BAM_STATS_H

/*
 * bam_stats.h
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <string.h>

#include "containers/array_list.h"
#include "bioformats/features/region/region_table.h"
#include "samtools/bam.h"

//------------------------------------------------------------------------

typedef struct bam_stats {
  // mapped
  int mapped;

  // strand
  int strand;
  
  // number of errors
  int num_errors;
  
  // cigar handling: number of indels and length
  int num_indels;
  int indels_length;

  // quality
  int quality;

  // unique alignment
  int unique_alignment;
  
  // handling pairs
  int single_end;
  int unmapped_pair_1;
  int unmapped_pair_2;
  int mapped_pair_1;
  int mapped_pair_2;
  int isize;

  // mapping length
  int seq_length;

  // nucleotide content
  int num_As;
  int num_Cs;
  int num_Gs;
  int num_Ts;
  int num_Ns;
  int num_GCs;
} bam_stats_t;

bam_stats_t *bam_stats_new();
void bam_stats_free(bam_stats_t *p);

//------------------------------------------------------------------------

typedef struct bam_stats_options {
  region_table_t *region_table;
  char **sequence_labels;
} bam_stats_options_t;

bam_stats_options_t *bam_stats_options_new(region_table_t *region_table,  
					   char **sequence_labels);
void bam_stats_options_free(bam_stats_options_t *p);


//------------------------------------------------------------------------

bam_stats_t *bam1_stats(bam1_t *bam1, bam_stats_options_t *opts);
int bam1s_stats(array_list_t *bam1s, bam_stats_options_t *opts,
		array_list_t *bam1s_stats);

//------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

#endif // end of BAM_STATS_H

//------------------------------------------------------------------------
//------------------------------------------------------------------------
