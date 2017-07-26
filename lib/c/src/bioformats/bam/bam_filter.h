#ifndef BAM_FILTER_H
#define BAM_FILTER_H

/*
 * bam_filter.h
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <string.h>

//#include "argtable2.h"
//#include "libconfig.h"
//#include "commons/log.h"
//#include "commons/system_utils.h"
//#include "commons/file_utils.h"

#include "sqlite/sqlite3.h"
#include "commons/workflow_scheduler.h"
#include "containers/array_list.h"
#include "bioformats/bam/bam_file.h"
#include "bioformats/features/region/region_table.h"


//------------------------------------------------------------------------

typedef struct bam_filter_options {
  int proper_pairs;
  int unique;

  int min_num_errors;
  int max_num_errors;

  int min_quality;
  int max_quality;

  int min_length;
  int max_length;

  region_table_t *region_table;
} bam_filter_options_t;

bam_filter_options_t *bam_filter_options_new(int unique, int proper_pairs,
					     int min_length, int max_length,
					     int min_quality, int max_quality,
					     int min_num_errors, int max_num_errors,
					     region_table_t *region_table);

void bam_filter_options_free(bam_filter_options_t *p);

//------------------------------------------------------------------------

void bam_filter(array_list_t *bam1s, array_list_t *passed_bam1s,
		array_list_t *failed_bam1s, bam_filter_options_t *opts);

//------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

#endif // end of BAM_FILTER_H

//------------------------------------------------------------------------
//------------------------------------------------------------------------
