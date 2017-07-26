/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _PSI1_H_
#define _PSI1_H_

#include "csa.h"
#include "sparsearray.h"
#include "sparsearray2.h"
#include "diskbuf.h"

#define ID_COMPPTR 0x80

#define SPARSEARRAY sparsearray
#define SPARSEARRAY_construct_init sparsearray_construct_init
#define SPARSEARRAY_construct_set sparsearray_construct_set
#define SPARSEARRAY_construct_end sparsearray_construct_end
#define SPARSEARRAY_getbit sparsearray_getbit
#define SPARSEARRAY_select sparsearray_select
#define SPARSEARRAY_rank sparsearray_rank
#define SPARSEARRAY_write sparsearray_write
#define SPARSEARRAY_read sparsearray_read

#define ENCODENUM encodegamma
#define DECODENUM decodegamma
//#define ENCODENUM encodedelta
//#define DECODENUM decodedelta


typedef struct {
  i64 n;
  i64 last;
  uchar *R;
  unsigned short *B; /* pointer to the bit stream encoding psi */
  int k; // number of bytes in an integer
  int L; /* interval between two psi values stored explicitly */
  SPARSEARRAY *sx, *sb;
  i64 maxrun;
  i64 psize;
  i64 id; // type of encoding
  MMAP *mappsi,*mappsd;
} psi1;

typedef struct {
  psi1 *ps;
  i64 pos; // Ÿ‚Éo‚Ä‚­‚é—v‘f‚Ì“Y‚¦š
  i64 n; // —v‘f‚ÌŒÂ” (“Y‚¦š‚Í 0 ‚©‚ç n-1)
  i64 page; // Œ»İŠi”[‚µ‚Ä‚¢‚éƒy[ƒW
  i64 *buf; // ƒLƒƒƒbƒVƒ…
} psi1_iterator;

void psi1_options(CSA *csa, char *p);
i64 psi1_makeindex(CSA *csa, char *fname);
i64 psi1_read(CSA *csa, char *fname);
i64 psi12_makeindex(CSA *csa, char *fname);

//static i64 psi1_psi(CSA *csa, i64 i);
//static i64 psi1_succ(CSA *csa, i64 pr, i64 l, i64 r);
//static i64 psi1_pred(CSA *csa, i64 pr, i64 l, i64 r);

//static int encodegamma(unsigned short *B,i64 p,i64 x);
//static int decodegamma(unsigned short *B,i64 p,i64 *ans);


void psi1_release(psi1 *ps);
psi1_iterator *psi1_iterator_new(psi1 *ps, i64 start);
i64 psi1_iterator_next(psi1_iterator *pi);
i64 psi1_iterator_hasnext(psi1_iterator *pi);
i64 psi1_iterator_get(psi1_iterator *pi, i64 i);
void psi1_iterator_remove(psi1_iterator *pi);


#endif
