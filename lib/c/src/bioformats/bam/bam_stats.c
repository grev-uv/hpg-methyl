/*
 * bam_stats.c
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#include "bam_stats.h"

//------------------------------------------------------------------------
// bam_stats_t
//------------------------------------------------------------------------

bam_stats_t *bam_stats_new() {
  bam_stats_t *p = (bam_stats_t *) calloc(1, sizeof(bam_stats_t));

  // mapped
  p->mapped = 0;

  // strand
  p->strand = 0;
  
  // number of errors
  p->num_errors = 0;
  
  // cigar handling: number of indels and length
  p->num_indels = 0;
  p->indels_length = 0;

  // quality
  p->quality = 0;

  // unique alignment
  p->unique_alignment = 0;

  // handling pairs
  p->single_end = 0;
  p->unmapped_pair_1 = 0;
  p->unmapped_pair_2 = 0;
  p->mapped_pair_1 = 0;
  p->mapped_pair_2 = 0;
  p->isize = 0;

  // mapping length
  p->seq_length = 0;

  // nucleotide content
  p->num_As = 0;
  p->num_Cs = 0;
  p->num_Gs = 0;
  p->num_Ts = 0;
  p->num_Ns = 0;
  p->num_GCs = 0;

  return p;
}

void bam_stats_free(bam_stats_t *p) {
  if (p) {
    free(p);
  }
}

//------------------------------------------------------------------------
// bam_stats_options_t
//------------------------------------------------------------------------

bam_stats_options_t *bam_stats_options_new(region_table_t *region_table,  
					   char **sequence_labels) {
  bam_stats_options_t *p = (bam_stats_options_t *) calloc(1, sizeof(bam_stats_options_t));

  p->region_table = region_table;
  p->sequence_labels = sequence_labels;

  return p;
}

void bam_stats_options_free(bam_stats_options_t *p) {
  if (p) {
    free(p);
  }
}


//------------------------------------------------------------------------
// bam1_stats
//------------------------------------------------------------------------

bam_stats_t *bam1_stats(bam1_t *bam1, bam_stats_options_t *opts) {
  
  bam_stats_t *bam_stats = NULL;
  uint32_t bam_flag = (uint32_t) bam1->core.flag;

  if (bam_flag & BAM_FUNMAP) {
    // not mapped, then return
    bam_stats = bam_stats_new();
    bam_stats->mapped = 0;
    return bam_stats;
  }

  if (opts->region_table) {
    region_t region;
    region.chromosome = opts->sequence_labels[bam1->core.tid];
    region.start_position = bam1->core.pos;
    region.end_position = region.start_position + bam1->core.l_qseq;
    region.strand = NULL;
    region.type = NULL;
    
    if (find_region(&region, opts->region_table)) {
      bam_stats = bam_stats_new();
    } else {
      return NULL;
    }
  } else {
    bam_stats = bam_stats_new();
  }

  // mapped !!
  bam_stats->mapped = 1;
  
  bam_stats->strand = (int) ((bam_flag & BAM_FREVERSE) > 0);
  
  // number of errors
  bam_stats->num_errors = bam_aux2i(bam_aux_get(bam1, "NM"));
  
  // cigar handling: number of indels and length
  uint32_t cigar_int, *cigar = bam1_cigar(bam1);
  int num_cigar_ops = (int) bam1->core.n_cigar; 
  for (int j = 0; j < num_cigar_ops; j++) {
    cigar_int = cigar[j];
    switch (cigar_int & BAM_CIGAR_MASK) {
    case BAM_CINS:  //I: insertion to the reference
    case BAM_CDEL:  //D: deletion from the reference
      bam_stats->num_indels++;
      bam_stats->indels_length += (cigar_int >> BAM_CIGAR_SHIFT);
      break;
    }
  }

  // quality
  bam_stats->quality = bam1->core.qual;

  // unique alignment
  if (!(bam_flag & BAM_FSECONDARY)) {
    bam_stats->unique_alignment = 1;
  }

  // handling pairs
  bam_stats->single_end = 1;
  if (bam_flag & BAM_FPAIRED) {
    bam_stats->single_end = 0;    
    if (bam_flag & BAM_FUNMAP) {
      if (bam_flag & BAM_FREAD1) {
	bam_stats->unmapped_pair_1 = 1;
      } else {
	bam_stats->unmapped_pair_2 = 1;
      }
    } else {
      if (bam_flag & BAM_FREAD1) {
	bam_stats->mapped_pair_1 = 1;
      } else {
	bam_stats->mapped_pair_2 = 1;
      }
    }
    
    if (!(bam_flag & BAM_FUNMAP) && !(bam_flag & BAM_FMUNMAP) && (bam_flag & BAM_FPROPER_PAIR)) { 
      bam_stats->isize = abs(bam1->core.isize);
    }
  }

  // mapping length
  char *bam_seq = (char *)bam1_seq(bam1);
  int seq_len = bam1->core.l_qseq;
  bam_stats->seq_length = seq_len;

  // nucleotide content
  for (int i = 0; i < seq_len; i++) {
    switch (bam1_seqi(bam_seq, i)) {
    case 1:
      bam_stats->num_As++;
      break;
    case 2:
      bam_stats->num_Cs++;
      break;
    case 4:
      bam_stats->num_Gs++;
      break;
    case 8:
      bam_stats->num_Ts++;
      break;
    case 15:
      bam_stats->num_Ns++;
      break;
    }
  }
  bam_stats->num_GCs = bam_stats->num_Gs + bam_stats->num_Cs;

  return bam_stats;
}

//------------------------------------------------------------------------
// bam1s_stats
//------------------------------------------------------------------------

int bam1s_stats(array_list_t *bam1s,
	       bam_stats_options_t *opts,
	       array_list_t *bam1s_stats) {

  bam1_t *bam1;
  bam_stats_t *stats;

  size_t num_items = array_list_size(bam1s);
  for (int i = 0; i < num_items; i++) {
    bam1 = array_list_get(i, bam1s);
    stats = bam1_stats(bam1, opts);
    array_list_insert(stats, bam1s_stats);
  }
  return array_list_size(bam1s_stats);
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
