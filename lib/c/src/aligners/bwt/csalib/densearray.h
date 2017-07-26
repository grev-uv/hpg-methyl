/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _DENSEARRAY_H_
#define _DENSEARRAY_H_

#include "csa.h"
#include "typedef.h"
#include "mman.h"

#ifdef __SSE4_2__
#include <smmintrin.h>
#define POPCOUNT(x) _mm_popcnt_u64(x)
#else
#define POPCOUNT(x) popcount(x)
#endif

#define SDARRAY_RANK1 (1<<0)
//#define SDARRAY_RANK0 (1<<1)
#define SDARRAY_RR (1<<1)
#define SDARRAY_SELECT1 (1<<2)
#define SDARRAY_SELECT0 (1<<3)
#define SDARRAY_LENHUF (1<<4)
#define SDARRAY_NOBUF (1<<5)
#define SDARRAY_COMPRANK (1<<6)
#define SDARRAY_COMPPTR (1<<7)
#define SDARRAY_SUC (1<<5)


#ifndef _BITVEC_T_
#define _BITVEC_T_
typedef u64 bitvec_t;
#endif

typedef struct {
  i64 n; // ベクトルの長さ
  i64 m; // 1の数
  i64 k; // 整数のバイト数
  i64 size; // 索引サイズ (ベクトルは含まない)
  bitvec_t *buf; // ベクトル
  int opt; // サポートする操作
  i64 ml,ms; // selectの表の要素数
  i64 rrr;

// for select
  dword *lp; // 32bit なので注意
  dword *sl;
  word *ss;
  i64 *p;

//for rank
//  dword *rl; // 32bit なので注意
  uchar *rl;
  word *rs;

} densearray;


typedef struct {
  i64 n; // ベクトルの長さ
  i64 m; // 1の数
  i64 size; // 索引サイズ (ベクトルは含まない)
  bitvec_t *buf; // ベクトル
  word *rs; // small blockごとの 1 の数
  dword *rl; // large blockごとの 1 の数
  i64 leaf_ofs;

  int opt; // サポートする操作

} densearray2;


void densearray_construct(densearray *da, i64 n, bitvec_t *buf, int opt);
int densearray_construct_init(densearray *da, i64 n);
int densearray_construct_set(densearray *da, i64 i, int x);
i64 densearray_construct_end(densearray *da, ushort L, int opt);
void densearray_make_selecttbl(void);

int densearray_getbit(densearray *da, i64 i);
i64 densearray_rank(densearray *da, i64 i);
i64 densearray_rank0(densearray *da, i64 i);
i64 densearray_rank_and_bit(densearray *da, i64 i, int *c);
i64 densearray_select(densearray *da, i64 i,int f);
i64 densearray_select2(densearray *da, i64 i,int f);
i64 densearray_write(densearray *da, FILE *f);
void densearray_read(densearray *da, uchar **map);

void densearray2_construct(densearray2 *da, i64 n, bitvec_t *buf, int opt);
i64 densearray2_rank(densearray2 *da, i64 i);
i64 densearray2_select(densearray2 *da, i64 i,int f);



#if 0
#ifndef _MYTIMESTRUCT_
#define _MYTIMESTRUCT_
typedef struct timeb mytimestruct;
void mygettime(mytimestruct *t);
double mylaptime(mytimestruct *before,mytimestruct *after);
#endif
#endif
#endif // _DENSEARRAY_H_
