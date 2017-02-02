#ifndef GFF_WRITE_H
#define GFF_WRITE_H

#include <stdio.h>
#include <string.h>

#include <containers/linked_list.h>
#include <containers/cprops/linked_list.h>

#include "gff_file_structure.h"
#include "gff_batch.h"

//====================================================================================
//  gff_write.h
//
//  gff writing functions prototypes
//====================================================================================

int gff_write_to_file(gff_file_t *gff_file, FILE *fd);

/* ************ Header management functions **********************/

void write_gff_header_entry(gff_header_entry_t *entry, FILE *fd);

/* ************ Record management functions **********************/

void write_gff_batch(gff_batch_t *gff_batch, FILE *fd);

void write_gff_record(gff_record_t* gff_record, FILE *fd);


#endif 
