/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "typedef.h"
#include "mman.h"
#include "csa.h"
#include "lf_bit.h"

#ifdef __SSE4_2__
#include <smmintrin.h>
#define POPCOUNT(x) _mm_popcnt_u64(x)
#else
#define POPCOUNT(x) popcount(x)
#endif

#define DD (sizeof(bitvec_t)*8)

//#define SB 32
//#define MB 256
#define logLB 16
#define LB (1<<logLB)




static int blog(i64 x)
{
i64 l;
  l = -1;
  while (x>0) {
    x>>=1;
    l++;
  }
  return l;
}

static u64 getuint(uchar *s, i64 i, i64 w)
{
  u64 x;
  i64 j;
  s += i*w;
  x = 0;
  for (j=0; j<w; j++) {
    x += ((u64)(*s++)) << (j*8);
  }
  return x;
}

static void putuint(uchar *s, i64 i, i64 x, i64 w)
{
  i64 j;
  s += i*w;
  for (j=0; j<w; j++) {
    *s++ = x & 0xff;
    x >>= 8;
  }
}

static void writeint(int k,i64 x,FILE *f)
{
  int i;
  for (i=k-1; i>=0; i--) {
    fputc(x & 0xff,f);
    x >>= 8;
  }
}

static void writeuint(int k,u64 x,FILE *f)
{
  int i;
  for (i=k-1; i>=0; i--) {
    fputc(x & 0xff,f);
    x >>= 8;
  }
}



static int setbit(bitvec_t *B, i64 i,bitvec_t x)
{
  i64 j,l;
  j = i / DD;
  l = i % DD;
  if (x==0) B[j] &= (~(1L<<(DD-1-l)));
  else if (x==1) B[j] |= (1L<<(DD-1-l));
  else {
    printf("error setbit x=%ld\n",x);
    exit(1);
  }
  return x;
}

/*
static int setbits(bitvec_t *B, i64 i, bitvec_t x, int w)
{
  i64 j;
  
  for (j=0; j<w; j++) {
    setbit(B,i+j,(x>>(w-1-j)) & 1);
  }
  return 0;
}
*/

static int lf_bit_BW_sub(lf_bit *lf,i64 i)
{
  bitvec_t x;
  x = lf->BW[i/DD];
  x = (x >> (DD-1-(i % DD))) & 1;
  return x;
}

static bitvec_t popcount(bitvec_t x)
{
  bitvec_t r;

  r = x;
  r = ((r & 0xaaaaaaaaaaaaaaaa)>>1) + (r & 0x5555555555555555);
  r = ((r & 0xcccccccccccccccc)>>2) + (r & 0x3333333333333333);
  r = ((r>>4) + r) & 0x0f0f0f0f0f0f0f0f;

  r *= 0x0101010101010101;
  r >>= 64-8;

  return r;
}

static i64 lf_bit_rank1_sub(lf_bit *lf, i64 i)
{
  i64 j,r;
  i64 i2;
//  i64 i3;
  int lmb;
  bitvec_t x,*p;

  if (i >= lf->last) i--;
  if (i < 0) return 0;

  lmb = lf->lmb;
  r = getuint(lf->RL,i>>logLB,lf->k);
  r += lf->RM[i>>lmb];
  i2 = (i>>lmb) << lmb;
  p = &lf->BW[i2/DD];

  for (j=0; j+DD <= (i-i2); j+=DD) {
    x = *p++;
    r += POPCOUNT(x);
  }
  x = *p;
  x >>= (DD - ((i-i2) - j + 1));
  r += POPCOUNT(x);

  return r;
}


static int convertchar(uchar t)
{
  int c = 0;
  switch (t) {
  case '0':  c = 0;  break;
  case '1':  c = 1;  break;
  default:  printf("error char = %c [%02x]\n",t,t);
  }
  return c;
}

static void make_tbl(lf_bit *lf)
{
  i64 n;
  i64 i;
  i64 C;
  int c,mb;
  bitvec_t *BW;
  int k;

  n = lf->n;
  k = lf->k;
  BW = lf->BW;
  mb = 1 << lf->lmb;

  lf->RL = (uchar *) mymalloc(k*((n+LB-1)/LB+1));
  lf->RM = (u16 *) mymalloc(sizeof(lf->RM[0])*((n+mb-1)/mb+1));

  for (i=0;i<((n+LB-1)/LB);i++) {
    putuint(lf->RL,i,0,k);
  }
  for (i=0;i<((n+mb-1)/mb);i++) lf->RM[i] = 0;

  C = 0;

  fprintf(stderr,"making tables...\n");
  for (i = 0; i < n; i++) {
    if (i % 1000000 == 0) {
      fprintf(stderr,"%ld\r",i/1000);
      fflush(stderr);
    }
    if (i % LB == 0) {
      putuint(lf->RL,(i / LB),C,k);
    }
    if (i % mb == 0) {
      lf->RM[i / mb] = (u16)(C - getuint(lf->RL,i / LB,k));
    }
    c = lf_bit_BW_sub(lf,i);
    if (c == 1) C++;
  }
  
}

void lf_bit_options(CSA *csa, char *p)
{
  lf_bit *lf;
  csa->psi_struc = lf = (lf_bit *) mymalloc(sizeof(lf_bit));
  lf->l = 128;
  if (p[0] == 0) goto end;
  p++;
  if (p[0] != ':') {
    sscanf(p,"%d",&lf->l);
//    printf("L = %d\n",lf->l);
    while (p[0] != ':' && p[0] != 0) p++;
  }
  if (p[0] == 0) goto end;
  p++;
end:
  printf("BW_BIT L=%d\n",lf->l);
}

static int lf_bit_BW(CSA *csa, i64 i)
{
  lf_bit *lf;
  lf = (lf_bit *)csa->psi_struc;

  if (i == lf->last) return -1; // EOF
  if (i > lf->last) i--;

  return lf_bit_BW_sub(lf,i);
}

static i64 lf_bit_rankc(CSA *csa, i64 i, int c)
{
  lf_bit *lf;
  i64 r;

  lf = (lf_bit *)csa->psi_struc;
  r = lf_bit_rank1_sub(lf,i);
  if (c == 0) {
    if (i >= lf->last) i--;
    r = i+1-r;
  }
  return r;
}


i64 lf_bit_makeindex(CSA *csa, char *fname)
{
  i64 last;
  FILE *in,*out;
  bitvec_t *BW;
  char *fbw, *flst, *fbw2, *fidx;
  int k,L,c,c2;
  i64 size;
  i64 n,i;
  lf_bit *lf;
  
  lf = (lf_bit *)csa->psi_struc;

  size = 0;

  L = lf->l;
  lf->lmb = blog(L);
  L = 1 << lf->lmb;
  lf->l = L;

  k = strlen(fname);
  fbw = (char *) mymalloc(k+5);
  flst = (char *) mymalloc(k+5);
  fbw2 = (char *) mymalloc(k+5);
  fidx = (char *) mymalloc(k+5);
  sprintf(fbw,"%s.bw",fname);
  sprintf(flst,"%s.lst",fname);
  sprintf(fbw2,"%s.bw1",fname);
  sprintf(fidx,"%s.bwb",fname);

  in = fopen(flst,"r");
  if (in == NULL) {
    perror("lf_bit_makeindex:");  exit(1);
  }
  int res = fscanf(in,"%ld",&last);
  if (!res) {  }
  printf("last = %ld\n",last);
  lf->last = last;
  fclose(in);


  in = fopen(fbw,"rb");
  if (in == NULL) {
    printf("lf_bit_makeindex: cannot open %s\n",fbw);
    exit(1);
  }
  fseek(in,0,SEEK_END);
  n = ftell(in);
  fseek(in,0,SEEK_SET);
  printf("n=%ld\n",n);
  csa->n = lf->n = n;

  for (i=0; i<SIGMA; i++) csa->C[i] = 0;

  BW = (bitvec_t *) mymalloc(sizeof(*BW) * (n/DD+1));
  lf->BW = BW;

  fprintf(stderr,"packing...\n");
  for (i = 0; i < n; i++) {
    if (i % 1000000 == 0) {
      fprintf(stderr,"%ld\r",i/1000);
      fflush(stderr);
    }
    c = fgetc(in);
    csa->C[c]++;
    c2 = convertchar(c);
    setbit(BW,i,c2);
  }
  fclose(in);

  out = fopen(fbw2,"w");
  if (out == NULL) {
    printf("lf_bit_makeindex: cannot open %s\n",fbw2);
  }
  for (i=0; i<n/DD+1; i++) {
    writeuint(sizeof(*BW),BW[i],out);
    size += sizeof(*BW);
  }
  fclose(out);

  lf->k = k = (blog(n+1)+1+8-1)/8;

  make_tbl(lf);

  out = fopen(fidx,"w");
  if (out == NULL) {
    printf("lf_bit_makeindex: cannot open %s\n",fidx);
  }

  writeint(1,ID_LF,out);
  writeint(1,k,out); /* #bytes of integer */
  writeint(k,n,out);
  writeint(k,last,out);
  writeint(k,L,out);
  size += 1+1+3*k;

  writeint(1,ID_BWT_BIT,out);
  size += 1;

  for (i=0;i<((n+LB-1)/LB);i++) {
    writeint(k,getuint(lf->RL,i,k),out);
    size += k;
  }
  for (i=0;i<((n+L-1)/L);i++) {
    writeint(sizeof(lf->RM[0]),lf->RM[i],out);
    size += sizeof(lf->RM[0]);
  }
  fclose(out);

  csa->BW = lf_bit_BW;
  csa->rankc = lf_bit_rankc;

  csa->LF = csa_LF;

  free(fbw);
  free(flst);
  free(fbw2);
  free(fidx);

  return size;
}


void lf_bit_read(CSA *csa, char *fname)
{
  char *fbw, *fbwi, *fname2;
  int k,l,id;
  i64 psize1,psize2;
  i64 n;
  lf_bit *lf;
  uchar *p, *q;
  
  csa->psi_struc = lf = (lf_bit *)mymalloc(sizeof(lf_bit));

  k = strlen(fname);
  fname2 = (char *) mymalloc(k-4+1);
  strncpy(fname2,fname,k-4);
  fname2[k-4] = 0;
  k -= 5;

  fbw = (char *) mymalloc(k+5+1);
  fbwi = (char *) mymalloc(k+5+1);

  sprintf(fbwi,"%s.bwb",fname2);

//  printf("psi_read: read %s\n",fbwi);

  lf->mapidx = mymmap(fbwi);
  if (lf->mapidx->addr==NULL) {
    perror("psi1_read: mmap2\n");
    exit(1);
  }
  p = q = (uchar *)lf->mapidx->addr;
  psize1 = lf->mapidx->len;

  id = getuint(p,0,1);  p += 1;
  if (id != ID_LF) {
    printf("lf_bit_read: id = %d is not supported.\n",id);
    exit(1);
  }
  lf->k = k = getuint(p,0,1);  p += 1;
  lf->n = n = getuint(p,0,k);  p += k;
  lf->last = getuint(p,0,k);  p += k;
  lf->l = l = getuint(p,0,k);  p += k;
  lf->lmb = blog(l);
  if ((1 << lf->lmb) != l) {
    printf("L=%d must be a power of 2.\n",l);
    exit(1);
  }

  id = getuint(p,0,1);  p += 1;
  lf->id = id;

//  printf("lf_bit_read: psi_id = %d L = %d\n",id,l);
  switch (id) {
    case ID_BWT_BIT:
      printf("#lf format = BWT_BIT\n");
      sprintf(fbw,"%s.bw1",fname2);
      break;
    default:
      printf("lf_bit_read: ID %d is not supported.\n",id);
      break;
  }

  lf->RL = (uchar *)p;  p += ((n+LB-1)/LB) * k;
  lf->RM = (u16 *)p;  p += ((n+l-1)/l) * sizeof(lf->RM[0]);

//  printf("lf_bit_read: map %s\n",fbw);
  lf->mapbwt = mymmap(fbw);
  if (lf->mapbwt->addr==NULL) {
    perror("psi1_read: mmap1\n");
    exit(1);
  }
  lf->BW = (bitvec_t *)lf->mapbwt->addr;
  psize2 = lf->mapbwt->len;
//  printf("psize2 = %ld\n",psize2);

//  printf("lf_bit_read: psize1 = %ld psize2 = %ld\n",psize1,psize2);
  lf->psize = psize1 + psize2;

  free(fbw);
  free(fbwi);
  free(fname2);

// user-specific functions
  csa->BW = lf_bit_BW;
  csa->rankc = lf_bit_rankc;

// default functions
  csa->LF = csa_LF;
  csa->selectc = csa_selectc;
  csa->psi = csa_psi_by_rankc_naive;
  csa->psi_succ = csa_psi_succ_naive;
  csa->psi_pred = csa_psi_pred_naive;
  csa->lookup = csa_lookup_lf;
  csa->inverse = csa_inverse_lf;
  csa->text = csa_text_lf;
  csa->search = csa_search_lf;
  csa->searchsub = csa_searchsub_lf;

}

