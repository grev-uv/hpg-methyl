/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _LF_DNA_H_
#define _LF_DNA_H_

#include <stdbool.h>
#include "typedef.h"

#ifdef __SSE4_2__
#include <smmintrin.h>
#define POPCOUNT2(x) _mm_popcnt_u64(x)
#else
#define POPCOUNT2(x) popcount2(x)
#endif


//#ifndef _BITVEC_T_
//#define _BITVEC_T_
//typedef u64 bitvec_t;
//#endif
#define DD (sizeof(bitvec_t)*8)

#define logSIGMADNA 2
#define SIGMADNA (1<<logSIGMADNA)

//#define SB 32
//#define MB 256
#define logLB 16
#define LB (1<<logLB)


typedef struct {
  i64 n;
  i64 last;
  int k; // number of bytes in an integer
  int l; /* interval between two psi values stored explicitly */
  int lmb;
  bitvec_t *BW; /* pointer to the bit stream encoding psi */
  uchar *RL;
  u16 *RM;
  i64 *C;
  i64 psize;
  i64 id; // type of encoding
  MMAP *mapbwt,*mapidx;
} lf_dna;

i64 lf_dna_makeindex(CSA *csa, char *fname, bool coded);
void lf_dna_read(CSA *sa, char *fname);
//static int lf_dna_BW_sub(lf_dna *lf,i64 i);
//static int lf_dna_BW(CSA *csa,i64 i);
//static i64 lf_dna_rankc_sub(lf_dna *lf, i64 i, int c);
//static i64 lf_dna_rankc(CSA *csa, i64 i, int c);
void lf_dna_options(CSA *csa, char *p);



#endif // _LF_DNA_H_
