/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _LF_WT_H_
#define _LF_WT_H_

#include "csa.h"
#include "densearray.h"
#include "sparsearray2.h"
#include "comparray.h"
#include "diskbuf.h"
#include "huffman.h"

#define ID_COMPPTR 0x80
#define ID_RR 0x40
#define ID_SUC 0x20


//typedef u64 bitvec_t;
#define DD (sizeof(bitvec_t)*8)

typedef struct lf_wt {
  i64 n;
  i64 last;
  int k; // number of bytes in an integer
//  comparray *da[SIGMA+2];
  void *da[SIGMA+2];
  Huffman *huffman;
  i64 id;
  ushort L;
  i64 psize;
  int opt;

  MMAP *map;

  int (*getbit)(void *da, i64 i);
  rank_t (*rank1)(void *da, rank_t i);
  rank_t (*rank0)(void *da, rank_t i);
  rank_t (*rank_bit)(void *da, rank_t i, int *c);
  i64 (*write)(void *da, FILE *f);
//  void (*read)(void *da, uchar **map);
  int (*init)(void *da, i64 n);
  int (*set)(void *da, i64 i, int x);
  i64 (*end)(void *da, ushort L, int opt);

} lf_wt;

i64 lf_wt_makeindex(CSA *csa, char *fname);
void lf_wt_read(CSA *sa, char *fname);
//static int lf_wt_BW_sub(lf_wt *lf,i64 i, int node);
//static int lf_wt_BW(CSA *csa,i64 i);
//static i64 lf_wt_rankc_sub(lf_wt *lf, i64 i, u64 x, int depth, int node);
//static i64 lf_wt_rankc(CSA *csa, i64 i, int c);
void lf_wt_options(CSA *csa, char *p);

int lf_wt_child_l(CSA *csa, i64 l, i64 r, uchar *head, rank_t *ll, rank_t *rr);


#endif // _LF_WT_H_
