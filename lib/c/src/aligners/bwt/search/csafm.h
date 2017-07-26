#ifndef _SEARCH_CSAFM_
#define _SEARCH_CSAFM_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <stdbool.h>

#include "../bwt_commons.h"

#if   defined SA_64

typedef uint64_t SA_TYPE;
typedef int64_t S_SA_TYPE;

#elif defined SA_32

typedef uint32_t SA_TYPE;
typedef int32_t S_SA_TYPE;

#elif defined SA_16

typedef uint16_t SA_TYPE;
typedef int16_t S_SA_TYPE;

#elif defined SA_8

typedef uint8_t SA_TYPE;
typedef int8_t S_SA_TYPE;

#else

typedef uint32_t SA_TYPE;
typedef int32_t S_SA_TYPE;

#endif

#if   defined FM_COMP_64

typedef uint64_t FM_COMP_TYPE;
#define FM_COMP_VALUE 64

#elif defined FM_COMP_32

typedef uint32_t FM_COMP_TYPE;
#define FM_COMP_VALUE 32

#endif

typedef struct {

  SA_TYPE **desp; // nA
  SA_TYPE siz; //Real number of columns of the uncompressed matrix
  SA_TYPE n_desp;
  SA_TYPE m_desp;

#if defined FM_COMP_32 || FM_COMP_64
  FM_COMP_TYPE **count; // nA
  SA_TYPE n_count;
  SA_TYPE m_count;
#endif

} comp_matrix;

typedef struct {

  SA_TYPE *vector;
  SA_TYPE n;

} vector;

typedef struct {

  SA_TYPE *vector;

  SA_TYPE siz; //Real size of the uncompressed vector
  SA_TYPE n;
  SA_TYPE ratio;

} comp_vector;

void reverse_strand_C(vector *r_C, vector *s_C, vector *r_C1, vector *s_C1, bwt_config_t config);
void reverse_strand_O(comp_matrix *r_O, comp_matrix *s_O, bwt_config_t config);

void read_vector(vector *vector, const char *directory, const char *name);
void read_comp_vector(comp_vector *vector, const char *directory, const char *name);
void read_comp_matrix(comp_matrix *matrix, const char *directory, const char *name);

void save_vector(vector *vector, const char *directory, const char *name);
void save_comp_vector(comp_vector *vector, const char *directory, const char *name);
void save_comp_matrix(comp_matrix *matrix, const char *directory, const char *name);

void free_comp_matrix(comp_matrix *reverse, comp_matrix *strand);

#ifdef VERBOSE_DBG

#define print_comp_matrix_count(M,n,m,siz)\
{\
	SA_TYPE bit, bbb;\
	printf("Matrix " #M ":\n");\
	for (SA_TYPE i_=0; i_<((SA_TYPE) (n)); i_++) {\
		printf("%ju: ", (uintmax_t) i_);\
		for (SA_TYPE j_=0; j_<((SA_TYPE) (siz)); j_++) {\
			bbb = j_ / FM_COMP_VALUE;\
			bit  = j_ % FM_COMP_VALUE;\
			printf("%ju ", (uintmax_t) (((M)[i_][bbb] >> bit) & 1));\
			if (bit+1 == FM_COMP_VALUE) printf("\t");\
		}\
		printf("\n");\
	}\
}

#if defined FM_COMP_32 || FM_COMP_64

#define print_comp_matrix(O)\
{\
  print_comp_matrix_count((O).count, (O).n_count, (O).m_count, (O).siz);\
  print_matrix((O).desp, (O).n_desp, (O).m_desp);\
}

#else

#define print_comp_matrix(O)\
  print_matrix((O).desp, (O).n_desp, (O).m_desp);
 
#endif

#else

#define print_comp_matrix_count(M,n,m,siz);
#define print_comp_matrix(O);

#endif

#ifdef __SSE4_2__
#include <smmintrin.h>
#endif

#if   defined FM_COMP_32

#ifdef __SSE4_2__
#define popcount(x) _mm_popcnt_u32(x)
#else
static inline unsigned int popcount(uint32_t i) {
  i = i - ((i >> 1) & 0x55555555);
  i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
  return (((i + (i >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24;
}
#endif

#elif defined FM_COMP_64

#ifdef __SSE4_2__
#define popcount(x) _mm_popcnt_u64(x)
#else
static inline unsigned int popcount(uint64_t i) {
  i = i - ((i >> 1) & 0x5555555555555555);
  i = (i & 0x3333333333333333) + ((i >> 2) & 0x3333333333333333);
  return (((i + (i >> 4)) & 0xF0F0F0F0F0F0F0F) * 0x101010101010101) >> 56;
}
#endif

#endif

static inline SA_TYPE get_O(SA_TYPE n, SA_TYPE m, comp_matrix *O) {

#if defined FM_COMP_32 || FM_COMP_64

  SA_TYPE pos, desp;
  pos  = m / FM_COMP_VALUE;
  desp = m % FM_COMP_VALUE;
  return O->desp[n][pos] + popcount( O->count[n][pos] << (FM_COMP_VALUE - (desp + 1)) );

#else

  return O->desp[n][m];

#endif

}

static inline uint8_t get_B_from_O(SA_TYPE m, comp_matrix *O) {

#if defined FM_COMP_32 || FM_COMP_64

  m++;

  SA_TYPE pos, desp;
  pos  = m / FM_COMP_VALUE;
  desp = m % FM_COMP_VALUE;

  for(uint8_t i=0; i<O->n_count; i++) {
    if ( ( (O->count[i][pos] >> desp) & ((FM_COMP_TYPE) 1)) != 0) return i;
  }

#else

  for(uint8_t i=0; i<O->n_desp; i++)
    if ( (O->desp[i][m] < O->desp[i][m+1]) ) return i;

#endif

  return (uint8_t) -1;

}

static inline SA_TYPE getScompValue(SA_TYPE m, comp_vector *Scomp, vector *C, comp_matrix *O) {

  SA_TYPE i,j;
  uint8_t b_aux;
  
  i=m; j=0;
  
  while (i % Scomp->ratio) {
    
    b_aux = get_B_from_O(i, O);
    
    if (b_aux == (uint8_t) -1) {
      
      i=0;

    } else {

      i = C->vector[b_aux] + get_O(b_aux, i+1, O);

    }
    
    j++;

  }

  return (Scomp->vector[i / Scomp->ratio] + j) % (O->siz-1);

}

static inline SA_TYPE getRcompValue(SA_TYPE m, comp_vector *Rcomp, vector *C, comp_matrix *O) {

	SA_TYPE i, j, k;
	uint8_t b_aux;

	i = (Rcomp->ratio - (m % Rcomp->ratio)) % Rcomp->ratio;
	k = m + i;

	if (k < Rcomp->siz) {
		j = Rcomp->vector[k / Rcomp->ratio];
	} else {
		j = Rcomp->vector[0];
		i = Rcomp->siz - m;
	}

	while (i) {

		b_aux = get_B_from_O(j, O);

		if (b_aux == (uint8_t) -1) {

			j=0;

		} else {

			j = C->vector[b_aux] + get_O(b_aux, j+1, O);

		}

		i--;

	}

	return j;

}

#endif
