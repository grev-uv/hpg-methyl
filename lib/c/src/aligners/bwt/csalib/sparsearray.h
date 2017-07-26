/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _SPARSEARRAY_H_
#define _SPARSEARRAY_H_

#include "densearray.h"

//#define LOWCHAR 0

//#define SPARSEARRAY_RANK1 (1<<0)
//#define SPARSEARRAY_RANK0 (1<<1)
//#define SPARSEARRAY_SELECT1 (1<<2)
//#define SPARSEARRAY_SELECT0 (1<<3)

#if 0
#define logD 6

#define PBS (sizeof(bitvec_t)*8)
#define D (1<<logD)
#define logM 5
#define M (1<<logM)
#define logP 8
#define P (1<<logP)
#define logLL 16    // size of word
#define LL (1<<logLL)
//#define logLLL 7
#define logLLL 5
//#define LLL 128
//#define LLL 32
#define LLL (1<<logLLL)
//#define logL 10
//#define logL (logLL-3)
#define logL (logLL-1-5)
#define L (1<<logL)
#endif

#ifndef min
 #define min(x,y) ((x)<(y)?(x):(y))
#endif

typedef struct {
  i64 n,m,d;
  i64 size;
  bitvec_t *hi;
#if LOWCHAR
  byte *low;
#else
  bitvec_t *low;
#endif
  densearray *sd0,*sd1;
  int opt;

} sparsearray;

int sparsearray_construct(sparsearray *sa, i64 n, bitvec_t *buf, int opt);
int sparsearray_construct_init(sparsearray *sa, i64 n, i64 m);
int sparsearray_construct_set(sparsearray *sa, i64 i, i64 x);
int sparsearray_construct_end(sparsearray *sa, int opt);
i64 sparsearray_write(sparsearray *sa, FILE *f);
void sparsearray_read(sparsearray *sa, uchar **map);

i64 sparsearray_select(sparsearray *sa, i64 i);
i64 sparsearray_rank(sparsearray *sa, i64 i);



#endif // _SPARSEARRAY_H_
