/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _DENSEARRAY2_H_
#define _DENSEARRAY2_H_

#include "typedef.h"
#include "densearray.h"

#ifdef __SSE4_2__
#include <smmintrin.h>
#define POPCOUNT(x) _mm_popcnt_u64(x)
#else
#define POPCOUNT(x) popcount(x)
#endif

//#define SDARRAY_RANK1 (1<<0)
//#define SDARRAY_RANK0 (1<<1)
//#define SDARRAY_SELECT1 (1<<2)
//#define SDARRAY_SELECT0 (1<<3)
//#define SDARRAY_NOBUF (1<<4)


#ifndef _BITVEC_T_
#define _BITVEC_T_
typedef u64 bitvec_t;
#endif

typedef struct {
  word *r; // SB でのrankの値
  word *s; // selectの値
  bitvec_t *b; // ベクトル
} densearray_sb;

i64 densearray_sb_construct(densearray_sb *da, i64 n, bitvec_t *buf, int opt);
i64 densearray_sb_write(densearray_sb *da, i64 n, int opt, FILE *f);
i64 densearray_sb_read(densearray_sb *da, i64 n, int opt, uchar **map);
i64 densearray_sb_select(densearray_sb *da, i64 i,int f);
i64 densearray_sb_rank(densearray_sb *da, i64 i);
i64 densearray_sb_rank0(densearray_sb *da, i64 i);
int densearray_sb_getbit(densearray_sb *da, i64 i);
i64 densearray_sb_rank_bit(densearray_sb *da, i64 i, int *c);

#if 0
#ifndef _MYTIMESTRUCT_
#define _MYTIMESTRUCT_
typedef struct timeb mytimestruct;
void mygettime(mytimestruct *t);
double mylaptime(mytimestruct *before,mytimestruct *after);
#endif
#endif

#endif // _DENSEARRAY_H_
