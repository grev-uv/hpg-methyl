#ifndef RNA_MAP_REGION_H
#define RNA_MAP_REGION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>

#include <commons/string_utils.h>

typedef struct rna_map_region {
	uint32_t start_position;
	uint32_t end_position;
	char *chromosome;
	unsigned char strand:1;
	unsigned char mapped:1;
	unsigned char rwmapped:1;
	unsigned char spliced:1;
	unsigned char rwspliced:1;
	unsigned char rep;
	// char* cigar;
} rna_map_region_t;


void rna_map_region_free(void * region);

void rna_print_region(rna_map_region_t* reg);

#endif	/*  RNA_MAP_REGION_H	*/
