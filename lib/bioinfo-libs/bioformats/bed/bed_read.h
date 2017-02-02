#ifndef BED_READ_H
#define BED_READ_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>

#include <commons/log.h>

#include "bed_file_structure.h"

//====================================================================================
//  bed_read.h
//
//  bed reading functions prototypes
//====================================================================================



/* ************ Header management functions **********************/

bed_header_entry_t* bed_header_entry_new();

void set_bed_header_entry_text(char *text, size_t length, bed_header_entry_t *entry);

/* ************ Record management functions **********************/

bed_record_t* bed_record_new();

void set_bed_record_chromosome(char* chromosome, size_t length, bed_record_t* bed_record);

void set_bed_record_start(unsigned long start, bed_record_t* bed_record);

void set_bed_record_end(unsigned long end, bed_record_t* bed_record);

void set_bed_record_name(char* name, size_t length, bed_record_t* bed_record);

void set_bed_record_score(float score, bed_record_t* bed_record);

void set_bed_record_strand(char strand, bed_record_t* bed_record);

void set_bed_record_thickstart(unsigned long thickstart, bed_record_t* bed_record);

void set_bed_record_thickend(unsigned long thickend, bed_record_t* bed_record);

void set_bed_record_itemrgb(char* item_rgb, size_t length, bed_record_t* bed_record);

void set_bed_record_blockcount(int block_count, bed_record_t* bed_record);

void set_bed_record_blocksizes(char* block_sizes, size_t length, bed_record_t* bed_record);

void set_bed_record_blockstarts(char* block_starts, size_t length, bed_record_t* bed_record);

#endif
