#ifndef BAM_FILTER_H
#define BAM_FILTER_H

/*
 * bam_filter.h
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

#include "commons/sqlite/sqlite3.h"
#include "commons/workflow_scheduler.h"
#include "containers/array_list.h"
#include "bioformats/bam/bam_file.h"
#include "bioformats/features/region/region_table.h"


//------------------------------------------------------------------------

typedef struct bam_filter_input {
  int num_threads;
  int batch_size;

  int by_mapped;
  int by_unmapped;
  int by_proper_pairs;
  int by_unique;

  int by_num_errors;
  int min_num_errors;
  int max_num_errors;

  int by_quality;
  int min_quality;
  int max_quality;

  int by_length;
  int min_length;
  int max_length;

  region_table_t *region_table;
  char *in_filename;
  char *out_dirname;


} bam_filter_input_t;

bam_filter_input_t *bam_filter_input_new(char *in_filename, char *out_dirname,
					 int by_mapped, int by_unmapped, 
					 int by_proper_pairs, int by_unique,
					 int by_num_errors, int min_num_errors, int max_num_errors,
					 int by_quality, int min_quality, int max_quality,
					 int by_length, int min_length, int max_length,
					 region_table_t *region_table,
					 int num_threads,int batch_size);

void bam_filter_input_free(bam_filter_input_t *input);

//------------------------------------------------------------------------

void bam_filter(bam_filter_input_t *input);

//------------------------------------------------------------------------

#endif // end of BAM_FILTER_H

//------------------------------------------------------------------------
//------------------------------------------------------------------------
