#ifndef _SEARCH_IO_
#define _SEARCH_IO_

#include <string.h>
#include <stdbool.h>

#include "../bwt_commons.h"

typedef struct {
  uint8_t *vector;
  uint64_t n;
  uint64_t dollar; //Position ending with the $ symbol (the first in the reference)
} ref_vector;

void read_ref_vector(ref_vector *vector, const char *directory, const char *name);
void save_ref_vector(ref_vector *vector, const char *directory, const char *name);

//Data structure for chromosome or exome separation positions in the reference
#define INDEX_EXOME 240000
#define IDMAX 500

typedef struct {
  int max_chromosomes;
  char *chromosome;
  uintmax_t *start;
  uintmax_t *end;
  uintmax_t *offset;
  uintmax_t size;
} exome;

exome *exome_new();
void exome_free(exome *ex);

void load_exome_file(exome *ex, const char *directory);
void save_exome_file(exome *ex, bool reverse, const char *directory);

void save_config(char *nucleotide, bool duplicate_strand, const char *directory);
void read_config(char *nucleotide, bool *duplicate_strand, const char *directory);

void encode_reference(ref_vector *X, exome *ex, const char *ref_path, bwt_config_t *bwt_config);
bool nextFASTAToken(FILE *queries_file, char *uncoded, uint8_t *coded, uintmax_t *nquery, bwt_config_t *bwt_config);

static inline uintmax_t binsearch(uintmax_t *array, uintmax_t size, uintmax_t key) {

  if( !array ) return 0;

  uintmax_t *p = array;
  uintmax_t w;

  while( size > 0 ) {

    w=size/2;

    if ( p[w] <= key ) {
      p+=w+1;
      size-=w+1;
    } else {
      size=w;
    }

  }

  return p - array;
}

#endif

/*
  This file is part of GNU BWT Aligner.

  Copyright José Salavert Torres, Ignacio Blanquer Espert
  Universitat Politècnica of València
  - Institut de Instrumentació per a la Imatge Molecular (GRyCAP)
  with partial support of Bull (S.A) Informatique
  2010 - 2013

  GNU BWT Aligner is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  GNU BWT Aligner is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with GNU BWT Aligner.  If not, see <http://www.gnu.org/licenses/>.
*/
