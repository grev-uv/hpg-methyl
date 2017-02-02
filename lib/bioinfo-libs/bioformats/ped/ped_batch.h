#ifndef PED_BATCH_H
#define PED_BATCH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include <containers/list.h>

#include "ped_file_structure.h"

//====================================================================================
//  ped_batch.h
//
//  ped_batch_t structures and prototypes
//====================================================================================

/**
 * Struct which represents a batch of PED records whose fields have been 
 * already loaded.
 */
typedef list_t ped_batch_t;
// typedef struct ped_batch {
//        size_t num_records;	// Number of records read
//        size_t size;             // Max buffer size (can be greater than num_records)
//        ped_record_t **records;	// Records read
// } ped_batch_t;

ped_batch_t* ped_batch_new(size_t size);

void ped_batch_free(ped_batch_t *ped_batch);

void add_record_to_ped_batch(ped_record_t *record, ped_batch_t *ped_batch);

int ped_batch_is_empty(ped_batch_t *ped_batch);

int ped_batch_is_full(ped_batch_t *ped_batch);

void ped_batch_print(FILE *fd, ped_batch_t *ped_batch);

#endif
