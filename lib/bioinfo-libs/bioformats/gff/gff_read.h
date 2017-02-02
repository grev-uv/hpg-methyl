#ifndef GFF_READ_H
#define GFF_READ_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>

#include <commons/log.h>

#include "gff_file_structure.h"

//====================================================================================
//  gff_read.h
//
//  gff reading functions prototypes
//====================================================================================



/* ************ Header management functions **********************/

gff_header_entry_t* gff_header_entry_new();

void set_gff_header_entry_text(char *text, size_t length, gff_header_entry_t *entry);

/* ************ Record management functions **********************/

gff_record_t* gff_record_new();

void set_gff_record_sequence(char* sequence, size_t length, gff_record_t* gff_record);

void set_gff_record_source(char* source, size_t length, gff_record_t* gff_record);

void set_gff_record_feature(char* feature, size_t length, gff_record_t* gff_record);

void set_gff_record_start(unsigned long start, gff_record_t* gff_record);

void set_gff_record_end(unsigned long end, gff_record_t* gff_record);

void set_gff_record_score(float score, gff_record_t* gff_record);

void set_gff_record_strand(char strand, gff_record_t* gff_record);

void set_gff_record_frame(char frame, gff_record_t* gff_record);

void set_gff_record_attribute(char* attribute, size_t length, gff_record_t* gff_record);

#endif
