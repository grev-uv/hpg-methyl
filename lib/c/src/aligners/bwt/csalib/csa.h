/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
/* compressed suffix arrays 

 */

#ifndef _CSA_H_
#define _CSA_H_

#include <stdio.h>
#include "typedef.h"
#include "../bwt_commons.h"

#define VERSION 2009071800

#define ID_HEADER 0x00
#define ID_SA     0x01
#define ID_ISA    0x02
#define ID_PSI    0x03
#define ID_LF     0x04
#define ID_HUFFMAN     0x05

#define ID_ARRAY_FL 0x00
#define ID_DIFF_GAMMA 0x01
#define ID_DIFF_GAMMA_RL 0x02
#define ID_BWT_DNA 0x03
#define ID_BWT_WT 0x04
#define ID_BWT_WT_HUF 0x05
#define ID_SPARSE4 0x06
#define ID_BWT_WT_DENSE 0x07
#define ID_DIFF_GAMMA_SPARSE 0x08
#define ID_DIFF_GAMMA_RL_SPARSE 0x09
#define ID_BWT_WT_SPARSE4 0x0a
#define ID_DIFF_GAMMA_RR 0x0b
#define ID_BWT_WT_RR 0x0c
#define ID_BWT_HUF 0x0d
#define ID_BWT_BIT 0x0e

#ifndef min
#define min(x,y) ((x)<(y)?(x):(y))
#endif
#ifndef max
#define max(x,y) ((x)>(y)?(x):(y))
#endif

typedef i64 rank_t, pos_t;

#define logSIGMA 8
#define SIGMA (1<<logSIGMA)

typedef struct csa {
  i64 n; /* length of the text */
  int k; /* number of bytes in an integer */
  int m; /* the number of distinct characters in the text */
  int D; /* interval between two SA values stored explicitly */
  int D2; /* interval between two inverse SA values stored explicitly */
  i64 C[SIGMA+2]; /* frequency of characters */
  i64 K[SIGMA+2]; /* table of cumulative frequency */
  int CtoA[SIGMA];
  uchar AtoC[SIGMA]; /* table of character codes */
  uchar *SA,*ISA;
  int id; /* format of psi/bw */
  void *psi_struc;
  void *mapidx;
  i64 psize, isize;

  pos_t (*lookup)(struct csa *csa, rank_t i);
  rank_t (*psi)(struct csa *csa,rank_t i);
  rank_t (*psi_succ)(struct csa *csa, rank_t pr, rank_t l, rank_t r);
  rank_t (*psi_pred)(struct csa *csa, rank_t pl, rank_t l, rank_t r);
  rank_t (*LF)(struct csa *csa,rank_t i);
  rank_t (*rankc)(struct csa *csa,rank_t i, int c);
  rank_t (*selectc)(struct csa *csa,rank_t i, int c);
  rank_t (*inverse)(struct csa *csa, pos_t suf);
  void (*text)(uchar *p,struct csa *csa,pos_t i,pos_t j);
  i64 (*substring)(uchar *p,struct csa *csa,rank_t r,i64 len);
  i64 (*substring_lf)(uchar *p,struct csa *csa,rank_t r,i64 len);
  int (*T)(struct csa *csa,rank_t i);
  int (*BW)(struct csa *csa,rank_t i);
  int (*BW_rank)(struct csa *csa,i64 i, rank_t *r);
  int (*head)(struct csa *csa,rank_t i);
  i64 (*search)(uchar *key,i64 keylen,struct csa *csa,rank_t *li,rank_t *ri);
  i64 (*searchsub)(int c, struct csa *csa, rank_t *ll, rank_t *rr);
  int (*child_l)(struct csa *csa, rank_t l, rank_t r, uchar *tail, rank_t *ll, rank_t *rr);
  int (*child_r)(struct csa *csa, i64 len, rank_t l, rank_t r, uchar *tail, rank_t *ll, rank_t *rr);

} CSA;

/* calculate SA[i] */
pos_t csa_lookup(CSA *csa, rank_t i);

/* calculate Psi[i] */
rank_t csa_psi(CSA *csa,rank_t i);

/* calculate SA^{-1}[i] */
rank_t csa_inverse(CSA *csa, pos_t suf);

/* decode T[i..j] into p */
void csa_text(unsigned char *p,CSA *csa, pos_t i, pos_t j);

/* decode T[SA[pos]..SA[pos]+len-1] into p */
i64 csa_substring(unsigned char *p,CSA *csa,rank_t r,i64 len);
i64 csa_substring_lf(uchar *p,CSA *csa,rank_t r,i64 len);
i64 csa_substring_lf_naive(uchar *p,CSA *csa,rank_t r,i64 len);

void csa_new_from_bwt_gnu_bwt_wrapper(const char *directory, const char *name);
void csa_new_from_bwt_wrapper(int argc, char *argv[]);
int csa_read(CSA *SA, int argc, char *argv[]);

i64 csa_search(uchar *key,i64 keylen,CSA *csa,rank_t *li,rank_t *ri);
i64 csa_search_lf(uchar *key, i64 keylen, CSA *csa, rank_t *ll, rank_t *rr);
rank_t csa_searchsub(int c, CSA *csa, rank_t *ll, rank_t *rr);
rank_t csa_searchsub_lf(int c, CSA *csa, rank_t *ll, rank_t *rr);
i64 csa_search_r(i64 keylen,int c, CSA *csa,rank_t *li,rank_t *ri);
int csa_left_diverse(CSA *csa, rank_t l, rank_t r);
int csa_right_diverse(CSA *csa, rank_t l, rank_t r, i64 length);

int csa_child_l(CSA *csa, rank_t l, rank_t r, uchar *head, rank_t *ll, rank_t *rr);
int csa_child_r(CSA *csa, i64 len, rank_t l, rank_t r, uchar *tail, rank_t *ll, rank_t *rr);
int csa_child_r0(CSA *csa, i64 len, rank_t l, rank_t r, uchar *tail, rank_t *ll, rank_t *rr);
pos_t csa_lookup_lf(CSA *csa, rank_t i);
rank_t csa_inverse_lf(CSA *csa, pos_t suf);
void csa_text_lf(uchar *p,CSA *csa, pos_t s, pos_t t);
rank_t csa_LF(CSA *csa, rank_t i);
rank_t csa_psi_by_rankc_naive(CSA *csa, rank_t i);
rank_t csa_selectc(CSA *csa, i64 i, int c);
rank_t csa_psi_pred_naive(CSA *csa, rank_t pr, rank_t l, rank_t r);
rank_t csa_psi_succ_naive(CSA *csa, rank_t pr, rank_t l, rank_t r);
rank_t csa_BW_LF_by_psi(CSA *csa, rank_t i, int *cc);
rank_t csa_LF_by_psi(CSA *csa, rank_t i);
int csa_BW_rank(CSA *csa,i64 i, rank_t *r);
int csa_T(CSA *csa,rank_t i);
int csa_head_rank(CSA *csa,rank_t i);
//void csa_error(void);

void bw_to_psi(FILE *out, CSA *csa, char *fbw, char *flst, int *k);

#endif // _CSA_H_
