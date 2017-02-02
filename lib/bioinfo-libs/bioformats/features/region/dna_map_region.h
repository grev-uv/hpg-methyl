#ifndef DNA_MAP_REGION_H
#define DNA_MAP_REGION_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>

#include <commons/string_utils.h>

typedef struct dna_map_region {
	uint32_t start_position;
	uint32_t end_position;
	char *chromosome;
	unsigned char strand:1;
	unsigned char mapped:1;
	unsigned char rwmapped:1;
	// unsigned char rep:6;
	// char* cigar;
} dna_map_region_t;


void dna_map_region_free(void * region);

int dna_map_region_equal_soft(dna_map_region_t *region1,dna_map_region_t* region2);

int dna_map_region_equal_hard(dna_map_region_t *region1,dna_map_region_t* region2);

void dna_print_region(dna_map_region_t* region);

void dna_fprint_region(FILE *f,dna_map_region_t* region);

#endif	/*  DNA_MAP_REGION_H	*/
