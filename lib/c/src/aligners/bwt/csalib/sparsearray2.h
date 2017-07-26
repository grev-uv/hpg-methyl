/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _SPARSEARRAY2_H_
#define _SPARSEARRAY2_H_

#include "densearray2.h"
#include "mman.h"

#define logD 6

#define PBS (sizeof(bitvec_t)*8)
//#define D (1<<logD)

#ifndef min
 #define min(x,y) ((x)<(y)?(x):(y))
#endif

typedef struct {
  i64 n;
  i64 m;
  i64 d;
  bitvec_t *hi;
  bitvec_t *low;
  densearray_sb *sd0,*sd1;
} sparsearray_sb;

i64 sparsearray_sb_construct(sparsearray_sb *sa, i64 n, bitvec_t *buf, int opt);
int sparsearray_sb_construct_init(sparsearray_sb *sa, i64 n, i64 m);
int sparsearray_sb_construct_set(sparsearray_sb *sa, i64 i, i64 x);
i64 sparsearray_sb_construct_end(sparsearray_sb *sa, int opt);
i64 sparsearray_sb_write(sparsearray_sb *sa, int opt, FILE *f);
i64 sparsearray_sb_read(sparsearray_sb *sa, int opt, uchar **map);

i64 sparsearray_sb_select(sparsearray_sb *sa, i64 i);
i64 sparsearray_sb_rank(sparsearray_sb *sa, i64 i);
i64 sparsearray_sb_rank0(sparsearray_sb *sa, i64 i);
int sparsearray_sb_getbit(sparsearray_sb *sa, i64 i);
i64 sparsearray_sb_rank_bit(sparsearray_sb *sa, i64 i, int *c);

typedef struct {
  i64 n,m,d,k;
  i64 size;
//  i64 *s;
  uchar *s;
  sparsearray_sb *sa;
  bitvec_t *hi, *low;
  int opt;

} sparsearray4;

i64 sparsearray4_construct(sparsearray4 *sa, i64 n, bitvec_t *buf, int opt);
int sparsearray4_construct_init(sparsearray4 *sa, i64 n, i64 m);
int sparsearray4_construct_set(sparsearray4 *sa, i64 i, i64 x);
i64 sparsearray4_construct_end(sparsearray4 *sa, ushort L, int opt);
i64 sparsearray4_write(sparsearray4 *sa, FILE *f);
i64 sparsearray4_read(sparsearray4 *sa, uchar **map);

i64 sparsearray4_select(sparsearray4 *sa, i64 i);
i64 sparsearray4_rank(sparsearray4 *sa, i64 i);
i64 sparsearray4_rank0(sparsearray4 *sa, i64 i);
int sparsearray4_getbit(sparsearray4 *sa, i64 i);
i64 sparsearray4_rank_and_bit(sparsearray4 *sa, i64 i, int *c);



#endif // _SPARSEARRAY2_H_
