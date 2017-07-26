#ifndef BWT_COMMONS_H
#define BWT_COMMONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include <stdbool.h>

#define MAXLINE      1000
#define MAXLINECOMP  250

//The meaning of Insertion and Deletion may be swapped to the meaning of the SAM file format
#define MATCH	  0
#define DELETION  1
#define MISMATCH  2
#define INSERTION 3

typedef struct {

	uint8_t *table;
  char *nucleotides;
  char *reverse;

	uint8_t nA;
  uint8_t AA, CC, GG, TT;

	bool inverse_sa;
  bool duplicate_strand;

} bwt_config_t;

#ifndef min
#define min(x,y) (((x)<(y))?(x):(y))
#endif

#ifndef max
#define max(x,y) (((x)>(y))?(x):(y))
#endif

#define check_syntax(argc, n, params)\
  if ((argc)!=(n)) {\
    fprintf(stderr, "Syntax:\n\t%s\n", (params));\
    exit(EXIT_FAILURE);\
  }

#define check_malloc(D, path)\
  if ((D)==NULL) {\
    fprintf(stderr, "Data structure " #D " in %s is too large\n", (path));\
    exit(EXIT_FAILURE);\
  }

#define check_file_open(fp, path)\
  if (!(fp)) {\
    fprintf(stderr, "Error opening file: %s\n", (path));\
    exit(EXIT_FAILURE);\
  }

#define check_file_read(err, nmemb, path)\
  if ((err) != (nmemb)) {\
    fprintf(stderr, "Error reading file '%s'\n", (path));\
    exit(EXIT_FAILURE);\
  }

#define check_file_write(err, nmemb, path)\
  if ((err) != (nmemb)) {\
    fprintf(stderr, "Error writing file '%s'\n", (path));\
    exit(EXIT_FAILURE);\
  }

#ifdef VERBOSE_DBG

#define print_matrix(M,n,m)\
{\
  printf("Matrix " #M ":\n");\
  for (uintmax_t i_=0; i_<((uintmax_t) (n)); i_++) {\
    printf("%ju: ", (uintmax_t) i_);\
    for (uintmax_t j_=0; j_<((uintmax_t) (m)); j_++) {\
      printf("%ju ", (uintmax_t) (M)[i_][j_]);\
    }\
    printf("\n");\
  }\
}

#define print_vector(V,n)\
{\
  printf("Vector " #V ":\n");\
  for (uintmax_t i_=0; i_<((uintmax_t) (n)); i_++) {\
		printf("%ju ", (uintmax_t) (V)[i_]);\
  }\
  printf("\n");\
}

#define print_string(S)\
	printf("String " #S ":\n%s\n", S);

#else

#define print_matrix(M,n,m);
#define print_vector(V,n);
#define print_string(S);

#endif

void *mymalloc(size_t n);
void *myrealloc(void *ptr, size_t next, size_t last);
void myfree(void *p, size_t s);
void report_mem(const char *s);

extern size_t cur_alloc, max_alloc;


//-----------------------------------------------------------------------------

void bwt_init_replace_table(bwt_config_t *config, char *nucleotides);
void bwt_free_replace_table(bwt_config_t *config);
/**
 *  @brief Encodes a sequence of plain nucleotides
 *  @param dest pointer to destination with encoded nucleotides
 *  @param src pointer to char with plain nucleotides
 *  @param length length of the nucleotide sequence
 * 
 *  Encodes a sequence of plain nucleotides
 */
void encode_bases(uint8_t* dest, char* src, uintmax_t length, uint8_t *table);

/**
 *  @brief Decodes a sequence of plain nucleotides
 *  @param dest pointer to destination with plain nucleotides
 *  @param src pointer to char with encoded nucleotides
 *  @param length length of the nucleotide sequence
 * 
 *  Encodes a sequence of plain nucleotides
 */
void decode_bases(char* dest, uint8_t* src, uintmax_t length, char *rev_table);

void revstring(uint8_t *X, uintmax_t nX);
void revstrand(uint8_t *X, uintmax_t nX, char *reverse);

void duplicate_reverse(uint8_t *X, uintmax_t nX, char *reverse);

#endif
