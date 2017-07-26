/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _COMPARRAY_H_
#define _COMPARRAY_H_

#include "typedef.h"
#include "sparsearray.h"
#include "huffman.h"

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
//#define SDARRAY_LENHUF (1<<4)

#ifndef _BITVEC_T_
#define _BITVEC_T_
typedef u64 bitvec_t;
#endif


typedef struct {
  i64 n; // ベクトルの長さ
  i64 m; // 1の数
  i64 size; // 索引サイズ (ベクトル含む)
  bitvec_t *buf; // ベクトル
  i64 bufsize; // ベクトルのサイズ (bitvec_tの数)
  uchar opt; // オプション (1の数をハフマン符号で符号化するか)
  ushort L;

  uchar k1; // pl のサイズ
  uchar k2; // rl のサイズ

  Huffman *h;

// pointer to buf
  sparsearray *psa;
  uchar *pl; // bufのオフセット．ビット単位．
  word *ps; // large block 内での buf のオフセット．ビット単位

//for rank
  sparsearray *rsa;
  uchar *rl;
  word *rs;

} comparray;


void comparray_construct(comparray *da, i64 n, bitvec_t *buf, ushort L, int opt);
int comparray_construct_init(comparray *da, i64 n);
int comparray_construct_set(comparray *da, i64 i, int x);
i64 comparray_construct_end(comparray *da, ushort L, int opt);
void comparray_maketbl(void);

i64 comparray_write(comparray *da, FILE *f);
void comparray_read(comparray *da, uchar **map);
int comparray_getbit(comparray *da, i64 i);
i64 comparray_rank(comparray *da, i64 i);
i64 comparray_rank0(comparray *da, i64 i);
i64 comparray_rank_and_bit(comparray *da, i64 i, int *c);

typedef struct {
  i64 n; // ベクトルの長さ
  i64 m; // 1の数
  i64 bufsize; // ベクトルのサイズ (bitvec_tの数)

  uchar opt;
  ushort L;
  Huffman *h;
  bitvec_t *buf; // ベクトル
// pointer to buf
  word *ps; // large block 内での buf のオフセット．ビット単位
//for rank
  word *rs;

} comparray_sb;


void comparray_sb_construct(comparray_sb *da, i64 n, bitvec_t *buf, ushort L, int opt);
int comparray_sb_construct_init(comparray_sb *da, i64 n);
int comparray_sb_construct_set(comparray_sb *da, i64 i, int x);
i64 comparray_sb_construct_end(comparray_sb *da, ushort L, int opt);

i64 comparray_sb_write(comparray_sb *da, FILE *f);
void comparray_sb_read(comparray_sb *da, uchar **map);
int comparray_sb_getbit(comparray_sb *da, i64 i);
i64 comparray_sb_rank(comparray_sb *da, i64 i);
i64 comparray_sb_rank0(comparray_sb *da, i64 i);


#endif // _COMPARRAY_H_
