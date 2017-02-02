#ifndef BED_WRITE_H
#define BED_WRITE_H

#include <stdio.h>
#include <string.h>

#include <containers/linked_list.h>
#include <containers/cprops/linked_list.h>

#include "bed_file_structure.h"
#include "bed_batch.h"

//====================================================================================
//  bed_write.h
//
//  bed writing functions prototypes
//====================================================================================

int bed_write_to_file(bed_file_t *bed_file, FILE *fd);

/* ************ Header management functions **********************/

void write_bed_header_entry(bed_header_entry_t *entry, FILE *fd);

/* ************ Record management functions **********************/

void write_bed_batch(bed_batch_t *bed_batch, FILE *fd);

void write_bed_record(bed_record_t* bed_record, FILE *fd);


#endif 
