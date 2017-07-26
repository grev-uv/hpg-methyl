/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _LF_BIT_H_
#define _LF_BIT_H_

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
} lf_bit;

i64 lf_bit_makeindex(CSA *csa, char *fname);
void lf_bit_read(CSA *sa, char *fname);
void lf_bit_options(CSA *csa, char *p);



#endif // _LF_DNA_H_
