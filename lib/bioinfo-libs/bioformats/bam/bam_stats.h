#ifndef BAM_STATS_H
#define BAM_STATS_H

/*
 * bam_stats.h
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#include <stdlib.h>
#include <string.h>

//#include "argtable2.h"
//#include "libconfig.h"
//#include "commons/log.h"
//#include "commons/system_utils.h"
//#include "commons/file_utils.h"

#include "commons/workflow_scheduler.h"
#include "commons/sqlite/sqlite3.h"
#include "containers/khash.h"
#include "containers/array_list.h"
#include "bioformats/features/region/region_table.h"
#include "bioformats/db/db_utils.h"
#include "bioformats/bam/bam_file.h"
#include "bioformats/bam/bam_db.h"


//------------------------------------------------------------------------

#define NUM_ERRORS_STATS 20
#define QUALITY_STATS   256

//------------------------------------------------------------------------

typedef struct bam_stats_input {
  int num_threads;
  int batch_size;
  region_table_t *region_table;
  char *in_filename;
  void *db;
  void *hash;
} bam_stats_input_t;

bam_stats_input_t *bam_stats_input_new(char *in_filename, region_table_t *region_table,
				       int num_threads, int batch_size, void *db, void *hash);
void bam_stats_input_free(bam_stats_input_t *input);

//------------------------------------------------------------------------

typedef struct bam_stats_output {
  // global statistics
  size_t ref_length;
  int num_sequences;
  int single_end;

  size_t num_reads;
  size_t num_unique_alignments;
  size_t num_mapped_reads;
  size_t num_unmapped_reads;
  size_t num_mapped_reads_1;
  size_t num_unmapped_reads_1;
  size_t num_mapped_reads_2;
  size_t num_unmapped_reads_2;

  size_t min_alignment_length;
  size_t max_alignment_length;

  // stats per strand
  size_t num_unique_alignments_strand[2];
  size_t num_mapped_reads_strand[2];

  // errors stats
  size_t num_indels;
  size_t indels_acc;
  size_t num_errors[NUM_ERRORS_STATS + 1];

  // nucleotide content
  size_t num_nucleotides;
  size_t num_As;
  size_t num_Cs;
  size_t num_Ts;
  size_t num_Gs;  
  size_t num_Ns;  
  size_t GC_content[100];

  // insert
  size_t min_insert_size;
  size_t max_insert_size;
  size_t insert_size_acc;

  // quality
  size_t min_quality;
  size_t max_quality;
  size_t quality_acc;
  size_t quality[QUALITY_STATS];

  // coverage (depth): global and per chromosome
  size_t unmapped_nts;
  double depth;
  char** sequence_labels;
  int **sequence_depths_per_nt;
  size_t* sequence_lengths;
  double* depth_per_sequence;

  //  khash_t(32) *gc_hash;
} bam_stats_output_t;

bam_stats_output_t *bam_stats_output_new();
void bam_stats_output_free(bam_stats_output_t *output);

//------------------------------------------------------------------------

void bam_stats(bam_stats_input_t *input, bam_stats_output_t *output);

//------------------------------------------------------------------------

#endif // end of BAM_STATS_H

//------------------------------------------------------------------------
//------------------------------------------------------------------------
