/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _PSI2_H_
#define _PSI2_H_

#include "csa.h"
#include "sparsearray2.h"
#include "diskbuf.h"


typedef struct {
  i64 n;
  i64 last;
  int k; // number of bytes in an integer
  i64 psize;
  i64 id; // type of encoding
  sparsearray4 *sa;
  MMAP *mappsi;
} psi2;

void psi2_options(CSA *csa, char *p);
i64 psi2_makeindex(CSA *csa, char *fname);
i64 psi2_read(CSA *csa, char *fname);

//static i64 psi2_psi(CSA *csa, i64 i);
//static i64 psi2_succ(CSA *csa, i64 pr, i64 l, i64 r);
//static i64 psi2_pred(CSA *csa, i64 pr, i64 l, i64 r);

//static int encodegamma(unsigned short *B,i64 p,i64 x);
//static int decodegamma(unsigned short *B,i64 p,i64 *ans);

#endif
