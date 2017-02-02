/*
    bwa_gpu a set of tools which allow short sequence alignment using the Burrows-Wheeler transform    usign both CPU and GPU approaches.
    Copyright (C) 2011  Jose Salavert Torres, Ignacio Blanquer Espert,
                        Andres Tomas Dominguez, Vicente Hernandez Garcia,
	 		Ignacio Medina Castello, Joaquin Tarraga Gimenez,
			Joaquin Dopazo Blazquez

    Contact e-mail: josator@fiv.upv.es, iblanque@dsic.upv.es, atomas@dsic.upv.es,
                    vhernand@dsic.upv.es, imedina@cipf.es, jtarraga@cipf.es,
                    jdopazo@cipf.es

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _BW_IO_
#define _BW_IO_

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<limits.h>
#include<stddef.h>
#include<malloc.h>
#include<ctype.h>
#include<errno.h>


#include "commons/string_utils.h"

#define MAX_MISMATCHES  20

#define INDEX_EXOME 1024
#define IDMAX 200

#define MAXLINE     200
#define MAXLINECOMP  50
#define MAX_READ_GPU 256000

#define TRUE            1
#define FALSE           0

//TODO: Usar estas macros para clarificar el codigo.
#define MATCH  0
#define DELETION  1
#define MISMATCH  2
#define INSERTION 3

//const char alph_rep[] ={'A','C','G','T'};

#ifdef VERBOSE_DBG
#define printUIntMatrix(M,n,m)			          \
  printf("Matrix " #M ":\n");			          \
  for (int i_=0; i_<((int) (n)); i_++) {		  \
    printf("%d: ", i_);					  \
    for (int j_=0; j_<((int) (m)); j_++) {		  \
      printf("%u ", (M)[i_][j_]);     		          \
    }                                                     \
    printf("\n");				          \
  }
#else
#define printUIntMatrix(M,n,m);
#endif

#ifdef VERBOSE_DBG
#define print32bitUIntMatrix(M,n,m,siz)					\
  {									\
    unsigned int bit, b32;						\
    printf("Matrix " #M ":\n");						\
    for (int i_=0; i_<((int) (n)); i_++) {				\
      printf("%d: ", i_);						\
      for (int j_=0; j_<((int) (siz)); j_++) {				\
	b32 = j_ / 32;							\
	bit  = j_ % 32;							\
	printf("%u ", ((M)[i_][b32] >> bit) & 1);			\
	if (bit+1 == 32) printf("\t");					\
      }									\
      printf("\n");							\
    }									\
  }
#else
#define print32bitUIntMatrix(M,n,m,siz);
#endif

#ifdef VERBOSE_DBG
#define print64bitUIntMatrix(M,n,m,siz)					\
  {									\
    unsigned int bit, b64;						\
    printf("Matrix " #M ":\n");						\
    for (int i_=0; i_<((int) (n)); i_++) {				\
      printf("%d: ", i_);						\
      for (int j_=0; j_<((int) (siz)); j_++) {				\
	b64 = j_ / 64;							\
	bit  = j_ % 64;							\
	printf("%llu ", ((M)[i_][b64] >> bit) & 1);			\
	if (bit+1 == 64) printf("\t");					\
      }									\
      printf("\n");							\
    }									\
  }
#else
#define print64bitUIntMatrix(M,n,m,siz);
#endif

#ifdef VERBOSE_DBG

#if   defined VECTOR_O_32BIT_COMPRESSION

#define printCompMatrix(O)					      \
  print32bitUIntMatrix((O).count, (O).n_count, (O).m_count, (O).siz); \
  printUIntMatrix((O).desp, (O).n_desp, (O).m_desp);

#elif defined VECTOR_O_64BIT_COMPRESSION

#define printCompMatrix(O)					      \
  print64bitUIntMatrix((O).count, (O).n_count, (O).m_count, (O).siz); \
  printUIntMatrix((O).desp, (O).n_desp, (O).m_desp);

#else

#define printCompMatrix(O)				\
  printUIntMatrix((O).desp, (O).n_desp, (O).m_desp);
  
#endif

#else
#define printCompMatrix(O);
#endif

#ifdef VERBOSE_DBG
#define printUIntVector(V,n)				\
  printf("Vector " #V ":\n");				\
  for (size_t i_=0; i_<((size_t) (n)); i_++) {		\
    if(((int)(V)[i_])==-1)				\
      printf("- ");				        \
    else                                                \
      printf("%u ", (V)[i_]);				\
  }                                                     \
  printf("\n");
#else
#define printUIntVector(V,n);
#endif

#ifdef VERBOSE_DBG
#define printIntVector(V,n)				\
  printf("Vector " #V ":\n");				\
  for (size_t i_=0; i_<((size_t) (n)); i_++) {		\
      printf("%d ", (V)[i_]);				\
  }                                                     \
  printf("\n");
#else
#define printIntVector(V,n);
#endif

#ifdef VERBOSE_DBG
#define printString(S) printf("String " #S ":\n%s\n", S);
#else
#define printString(S);
#endif

#define checkMalloc(D, path)						 \
  if ((D)==NULL) {							 \
    fprintf(stderr, "Data structure " #D " in %s is too large", (path)); \
    exit(1);								 \
  }

//--------------------------------------------------------------------------------------
// In case the has to be opened in write mode, the user launching the aligner application
// must have enough permissions to access the destination directory with R+W permissions.
// Now the error output logs the corresponding file access error, which in that case
// would be "Access denied"
// - Date: 11 / 11 / 2016
// - Who: Cesar
#define checkFileOpen(fp, path)						\
  if (!(fp)) {								\
    fprintf(stderr, "Error opening file %s: %s\n", (path), strerror(errno));		\
    exit(1);								\
  }

#define checkFileRead(err, nmemb, path)					\
  if ((err) != (nmemb)) {						\
    fprintf(stderr, "Error reading file '%s'\n", (path));		\
    exit(1);								\
  }

#define checkFileWrite(err, nmemb, path)			\
  if ((err) != (nmemb)) {					\
    fprintf(stderr, "Error writing file '%s'\n", (path));	\
    exit(1);							\
  }

#define timevars() struct timeval t1, t2;

#define tic_bwt(msg)				   \
  printf(">> " msg "\n");			   \
  fflush(stdout);				   \
  gettimeofday(&t1, NULL);

#define toc_bwt()					\
  gettimeofday(&t2, NULL);						                 \
  printf("<< Finished %.0f usecs\n", (t2.tv_sec-t1.tv_sec)*1e6+(t2.tv_usec-t1.tv_usec)); \
  fflush(stdout);

typedef struct {

  size_t siz;

  unsigned int *desp[4]; // nA
  size_t n_desp;
  size_t m_desp;

#if   defined VECTOR_O_32BIT_COMPRESSION

  unsigned int *count[4]; // nA
  size_t n_count;
  size_t m_count;

#elif defined VECTOR_O_64BIT_COMPRESSION

  unsigned long long *count[4];  // nA
  size_t n_count;
  size_t m_count;

#endif

} comp_matrix;

typedef struct {

  unsigned int *vector;
  size_t n;

} vector;

typedef struct {

  unsigned int *vector;
  size_t n;
  size_t siz;
  size_t ratio;

} comp_vector;

typedef struct {

  char *vector;
  size_t n;

} byte_vector;

typedef struct {
  
  char chromosome[INDEX_EXOME*IDMAX];
  unsigned int start[INDEX_EXOME];
  unsigned int end[INDEX_EXOME];
  unsigned int offset[INDEX_EXOME];
  /*
  char *chromosome;
  unsigned int *start;
  unsigned int *end;
  unsigned int *offset;*/
  unsigned int size;
} exome;

typedef struct {
  size_t k, l;
  int start, pos, end;
  char dir;
  char num_mismatches;
  char err_kind[MAX_MISMATCHES];
  char err_base[MAX_MISMATCHES];
  int err_pos[MAX_MISMATCHES];
  unsigned int read_index;
} result;

typedef struct {
  result *list;
  unsigned int num_results;
  unsigned int max_results;
  unsigned int read_index;
} results_list;


//--------------------------------------------------------------------------------------
// In case the build is being compiled with all optimizations disabled (ie: debug build)
// inline functions must be forward declared to enable the compiler to find the appropriate
// simbols to link into them on all the translation units.
// - Date: 14 / 11 / 2016
// - Who: Cesar
#ifdef __GNUC__
#ifdef __NO_INLINE__

static void new_results_list(results_list*, unsigned int);
static void init_result(result*, char);
static void bound_result(result*, int, int);
static void change_result(result*, size_t, size_t, int);
static void add_mismatch(result*, char, char, int);
static void modify_last_mismatch3(result*, char, char, int);
static void modify_last_mismatch2(result*, char, char);
static void modify_last_mismatch1(result*, char);
static void copy_result(result*, result*);
static void add_result(result*, results_list*);
static unsigned int getScompValue(size_t, comp_vector*, vector*, comp_matrix*);
static unsigned int getRcompValue(size_t, comp_vector*, vector*, comp_matrix*);
static char getBfromO(size_t, comp_matrix*);
static unsigned int getOcompValue(size_t, size_t, comp_matrix*);

#if defined VECTOR_O_32BIT_COMPRESSION
static unsigned int bitcount(unsigned int);
#elif defined VECTOR_O_64BIT_COMPRESSION
static unsigned int bitcount(unsigned long long);
#endif /* VECTOR_O_XBIT_COMPRESSION */

#endif /* __NO_INLINE__ */
#endif /* __GNUC__ */
//--------------------------------------------------------------------------------------


inline void new_results_list(results_list *r_list, unsigned int max_results) {

  //printf("Size of result %lu\n", sizeof(result)*CHAR_BIT);

  r_list->list = (result *) malloc(max_results * sizeof(result));
  checkMalloc(r_list->list, "new_result_list");

  r_list->num_results = 0;
  r_list->max_results = max_results;

}

inline void init_result(result *r, char dir) {
  r->dir = dir;
  r->num_mismatches = 0;
}

inline void bound_result(result *r, int start, int end) {
  r->start = start;
  r->end = end;
}

inline void change_result(result *r, size_t k, size_t l, int pos) {
  r->k = k;
  r->l = l;
  r->pos = pos;
}

inline void add_mismatch(result *r, char err_kind, char base, int position) {

  if (r->num_mismatches < MAX_MISMATCHES) {

    size_t pos_mismatch = r->num_mismatches;

    r->err_kind[pos_mismatch] = err_kind;
    r->err_base[pos_mismatch] = base;
    r->err_pos[pos_mismatch] = position;

    r->num_mismatches++;

    //  } else {

    //    fprintf(stderr, "ERROR: Number of allowed mismatches exceeded: %d\n", r->num_mismatches);
    //exit(1);

  }

}

inline void modify_last_mismatch3(result *r, char err_kind, char base, int position) {

  size_t pos_mismatch = r->num_mismatches-1; // Last mismatch

  r->err_kind[pos_mismatch] = err_kind;
  r->err_base[pos_mismatch] = base;
  r->err_pos[pos_mismatch] = position;

}

inline void modify_last_mismatch2(result *r, char err_kind, char base) {

  size_t pos_mismatch = r->num_mismatches-1; // Last mismatch

  r->err_kind[pos_mismatch] = err_kind;
  r->err_base[pos_mismatch] = base;

}

inline void modify_last_mismatch1(result *r, char err_kind) {

  r->err_kind[r->num_mismatches-1] = err_kind;

}

inline void copy_result(result *dest, result *orig) {

  init_result(dest, orig->dir);
  bound_result(dest, orig->start, orig->end);
  change_result(dest, orig->k, orig->l, orig->pos);

  int i;
  for(i=0; i<orig->num_mismatches; i++) {
    add_mismatch(dest, orig->err_kind[i], orig->err_base[i], orig->err_pos[i]);
  }

}

inline void add_result(result *orig, results_list *r_list) {
  //printf("\tCALL to add result, %i-%i\n", 
  //	 r_list->num_results, r_list->max_results);
  if (r_list->num_results < r_list->max_results) {
    //printf("\t------------>ADD RESULT k=%lu, l=%lu, start=%i, pos=%i, end=%i\n",
    //   orig->k, orig->l, orig->start, orig->pos,orig->end);
    result *dest;
    dest = &r_list->list[r_list->num_results];
    r_list->num_results++;

    copy_result(dest, orig);

    //  } else {

    //    fprintf(stderr, "ERROR: Number of allowed results exceeded: %u\n", r_list->max_results);
    //exit(1);

  }

}

#if defined VECTOR_O_32BIT_COMPRESSION

inline unsigned int bitcount(unsigned int i) {

  //Parallel binary bit add
  i = i - ((i >> 1) & 0x55555555);
  i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
  return (((i + (i >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24;

}

inline unsigned int getOcompValue(size_t n, size_t m, comp_matrix *O) {

  size_t pos, desp;

  pos  = m / 32;
  desp = m % 32;
  
  return O->desp[n][pos] + bitcount( O->count[n][pos] << (32 - (desp + 1)) );

}

#elif defined VECTOR_O_64BIT_COMPRESSION

inline unsigned int bitcount(unsigned long long i) { //TODO: Probar la llamada en C y la de GCC para ver si usan la instrucción máquina.

  //Parallel binary bit add
  i = i - ((i >> 1) & 0x5555555555555555);
  i = (i & 0x3333333333333333) + ((i >> 2) & 0x3333333333333333);
  return (((i + (i >> 4)) & 0xF0F0F0F0F0F0F0F) * 0x101010101010101) >> 56;
}

inline unsigned int getOcompValue(size_t n, size_t m, comp_matrix *O) {

  size_t pos, desp;

  pos  = m / 64;
  desp = m % 64;

  //printf("****************n = %lu\tm = %lu \n", n, m);

  return O->desp[n][pos] + bitcount( O->count[n][pos] << (64 - (desp + 1)) );

}

#endif

inline char getBfromO(size_t m, comp_matrix *O) {

#if   defined VECTOR_O_32BIT_COMPRESSION

  m++;

  size_t pos, desp;
  pos  = m / 32;
  desp = m % 32;

  size_t i;
  for(i=0; i<O->n_count; i++) {
    if ( ( (O->count[i][pos] >> desp) & 1) != 0) return i;
  }

#elif defined VECTOR_O_64BIT_COMPRESSION

  m++;

  size_t pos, desp;
  pos  = m / 64;
  desp = m % 64;

  size_t i;
  for(i=0; i<O->n_count; i++) {
    if ( ( (O->count[i][pos] >> desp) & 1) != 0) return i;
  }

#else

  for(size_t i=0; i<O->n_desp; i++)
    if ( (O->desp[i][m] < O->desp[i][m+1]) ) return i;

#endif

  return -1;

}

inline unsigned int getScompValue(size_t m, comp_vector *Scomp, vector *C, comp_matrix *O) {

  size_t i,j;
  
  i=m; j=0;
  
  char b_aux;
  
  while (i % Scomp->ratio) {
    
    b_aux = getBfromO(i, O);
    
    if ((int) b_aux == -1) {
      
      i=0;

    } else {

#if defined VECTOR_O_32BIT_COMPRESSION || VECTOR_O_64BIT_COMPRESSION
      i = C->vector[(int)b_aux] +  getOcompValue((int)b_aux, i+1, O);
#else
      i = C->vector[(int)b_aux] + O->desp[(int)b_aux][i+1/*0 is -1*/];
#endif
      
    }
    
    j++;

  }

  return (Scomp->vector[i / Scomp->ratio] + j) % (O->siz-1);

}

inline unsigned int getScompValueB(size_t m, comp_vector *Scomp, vector *C, comp_matrix *O, byte_vector *B) {
  
  size_t i, j;
  
  i=m; j=0;
  
  char b_aux;
  
  while (i % Scomp->ratio) {
    
    b_aux = B->vector[i];
    
    if ((int) b_aux == -1) {

      i=0;

    } else {

#if defined VECTOR_O_32BIT_COMPRESSION || VECTOR_O_64BIT_COMPRESSION
      i = C->vector[(int)b_aux] + getOcompValue((int)b_aux, i+1, O);
#else
      i = C->vector[(int)b_aux] + O->desp[(int)b_aux][i+1/*0 is -1*/];
#endif

    }

    j++;

  }

  return (Scomp->vector[i / Scomp->ratio] + j) % (O->siz-1);

}

inline unsigned int getRcompValue(size_t m, comp_vector *Rcomp, vector *C, comp_matrix *O) {

  size_t i, j, k;

  i = (Rcomp->ratio - (m % Rcomp->ratio)) % Rcomp->ratio;
  k = m + i;

  if (k < Rcomp->siz) {
    j = Rcomp->vector[k / Rcomp->ratio];
  } else {
    j = Rcomp->vector[0];
    i = Rcomp->siz - m;
  }

  char b_aux;

  while (i) {

    b_aux = getBfromO(j, O);

    if ((int) b_aux == -1) {

      j=0;

    } else {

#if defined VECTOR_O_32BIT_COMPRESSION || VECTOR_O_64BIT_COMPRESSION
      j = C->vector[(int)b_aux] + getOcompValue((int)b_aux, j+1, O);
#else
      j = C->vector[(int)b_aux] + O->desp[(int)b_aux][j+1/*0 is -1*/];
#endif

    }

    i--;

  }

  return j;

}

inline unsigned int getRcompValueB(size_t m, comp_vector *Rcomp, vector *C, comp_matrix *O, byte_vector *B) {

  size_t i, j, k;

  i = (Rcomp->ratio - (m % Rcomp->ratio)) % Rcomp->ratio;
  k = m + i;

  if(k < Rcomp->siz) {
    j = Rcomp->vector[k / Rcomp->ratio];
  } else {
    j = Rcomp->vector[0];
    i = Rcomp->siz - m;
  }

  char b_aux;

  while (i) {

    b_aux = B->vector[j];

    if ((int) b_aux == -1) {

      j=0;

    } else {

#if defined VECTOR_O_32BIT_COMPRESSION || VECTOR_O_64BIT_COMPRESSION
      j = C->vector[(int)b_aux] + getOcompValue((int)b_aux, j+1, O);
#else
      j = C->vector[(int)b_aux] + O->desp[(int)b_aux][j+1/*0 is -1*/];
#endif

    }

    i--;

  }

  return j;

}

void reverseStrandC(vector *r_C, vector *s_C, vector *r_C1, vector *s_C1);
void reverseStrandO(comp_matrix *r_O, comp_matrix *s_O);
void freeCompMatrix(comp_matrix *matrix);

void readUIntVector(vector *vector, const char *directory, const char *name);
void readUIntCompVector(comp_vector *vector, const char *directory, const char *name);
void readCharVector(byte_vector *vector, const char *directory, const char *name);
void readCompMatrix(comp_matrix *matrix, const char *directory, const char *name);

void saveUIntVector(vector *vector, const char *directory, const char *name);
void saveUIntCompVector(comp_vector *vector, const char *directory, const char *name);
void saveCharVector(byte_vector *vector, const char *directory, const char *name);
void saveCompMatrix(comp_matrix *matrix, const char *directory, const char *name);

void initReplaceTable(const char *str);

char *replaceBases(char *uncoded, char *coded, size_t length);

int nextFASTAToken(FILE *queries_file, char *uncoded, char *coded, unsigned int *nquery, char *compressed, unsigned int*ncompress);

void load_reference(byte_vector *X, int duplicate, exome *ex, const char *path);

void load_exome_file(exome *ex, const char *directory);
void save_exome_file(exome *ex, const char *directory);

void revstring(char *X, size_t nX);

size_t comp4basesInChar(char *X, size_t nX, char *Y);
unsigned int binsearch(unsigned int *array, unsigned int size, size_t key);

void initialize_init_mask();
int write_results(results_list *r_list, exome* ex, comp_vector *S, comp_vector *Si, vector *C, comp_matrix *O, comp_matrix *Oi, char *mapping, int nW, int type, FILE *fp);

#endif
