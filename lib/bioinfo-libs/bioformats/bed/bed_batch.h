#ifndef BED_BATCH_H
#define BED_BATCH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <assert.h>

#include <containers/array_list.h>

#include "bed_file_structure.h"

//====================================================================================
//  bed_batch.h
//
//  bed_batch_t structures and prototypes
//====================================================================================

/**
 * Struct which represents a batch of BED records whose fields have been 
 * already loaded.
 */
typedef struct bed_batch {
    array_list_t *records;      /**< Records in the block */
    char *text;                 /**< Input buffer with the data for the records */
} bed_batch_t;

bed_batch_t* bed_batch_new(size_t size);

void bed_batch_free(bed_batch_t *bed_batch);

void add_record_to_bed_batch(bed_record_t *record, bed_batch_t *bed_batch);

int bed_batch_is_empty(bed_batch_t *bed_batch);

int bed_batch_is_full(bed_batch_t *bed_batch);

void bed_batch_print(FILE *fd, bed_batch_t *bed_batch);

#endif
