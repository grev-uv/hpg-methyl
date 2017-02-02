#ifndef GENOME_H
#define GENOME_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>

#include "commons/log.h"
#include "containers/cprops/hashtable.h"

//====================================================================================
//  structures and prototypes
//====================================================================================

typedef struct genome {
  size_t genome_length;
  unsigned int num_chromosomes;
  size_t *chr_name_length;
  size_t *chr_size;
  size_t *chr_offset;
  unsigned char* X;
  char **chr_name;
  char **code_table;
} genome_t;

genome_t* genome_new(char* sequence_filename, char* directory);

void genome_free(genome_t* genome_p);

void genome_read_sequence(char* sequence, unsigned int strand, char* chromosome, unsigned long int* start_p, unsigned long int* end_p, genome_t* genome_p);

void genome_read_sequence_by_chr_index(char* sequence, unsigned int strand, unsigned int chr, size_t* start_p, size_t* end_p, genome_t* genome_p);

void generate_codes();

unsigned char *load_binary_dna(char *dna_binary_filename, size_t *size);

char **load_array_codes();

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

#endif  // GENOME_H
