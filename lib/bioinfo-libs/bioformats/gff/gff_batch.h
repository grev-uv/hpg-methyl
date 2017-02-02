#ifndef GFF_BATCH_H
#define GFF_BATCH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <assert.h>

#include <containers/array_list.h>

#include "gff_file_structure.h"

//====================================================================================
//  gff_batch.h
//
//  gff_batch_t structures and prototypes
//====================================================================================

/**
 * Struct which represents a batch of GFF records whose fields have been 
 * already loaded.
 */
typedef struct gff_batch {
    array_list_t *records;      /**< Records in the block */
    char *text;                 /**< Input buffer with the data for the records */
} gff_batch_t;

gff_batch_t* gff_batch_new(size_t size);

void gff_batch_free(gff_batch_t *gff_batch);

void add_record_to_gff_batch(gff_record_t *record, gff_batch_t *gff_batch);

int gff_batch_is_empty(gff_batch_t *gff_batch);

int gff_batch_is_full(gff_batch_t *gff_batch);

void gff_batch_print(FILE *fd, gff_batch_t *gff_batch);

#endif
