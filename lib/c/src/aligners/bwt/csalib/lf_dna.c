/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mman.h"
#include "csa.h"
#include "lf_dna.h"

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

static int setbits(bitvec_t *B, i64 i, bitvec_t x, int w)
{
  i64 j;
  
  for (j=0; j<w; j++) {
    setbit(B,i+j,(x>>(w-1-j)) & 1);
  }
  return 0;
}

static int lf_dna_BW_sub(lf_dna *lf,i64 i)
{
  bitvec_t x;
  x = lf->BW[i/(DD/logSIGMADNA)];
  x = (x >> (DD-logSIGMADNA-(i % (DD/logSIGMADNA))*logSIGMADNA)) & (SIGMADNA-1);
  return x;
}

#if 0
static i64 lf_dna_rankc_sub_slow(lf_dna *lf, i64 i, int c)
{
  int j;
  i64 r;
  i64 i2;
  int lmb,mb;

//  printf("rankc(%ld)\n",i);
  if (i >= lf->last) i--;

  lmb = lf->lmb;
  mb = 1 << lmb;
  r = getuint(lf->RL,(i/LB)*SIGMADNA + c,lf->k);
  if (c < SIGMADNA-1) {
    r += lf->RM[(i>>lmb)*(SIGMADNA-1) + c];
  } else {
    r += ((i & (LB-1))>>lmb)<<lmb;
    r -= lf->RM[(i>>lmb)*(SIGMADNA-1) + 0];
    r -= lf->RM[(i>>lmb)*(SIGMADNA-1) + 1];
    r -= lf->RM[(i>>lmb)*(SIGMADNA-1) + 2];
  }
//  printf("true r=%d ",r);
  i2 = (i/mb)*mb;
  for (j=0; j <= (i % mb); j++) {
//    printf("W[%ld] = %d\n",i2+j,csa_W(lf,i2+j));
    if (lf_dna_BW_sub(lf,i2+j) == c) r++;
  }
//  printf("true r=%d\n",r);
  return r;
}
#endif

static bitvec_t popcount2(bitvec_t x)
{
  bitvec_t r;
  r = x;
  r = ((r >> 2) + r) & 0x3333333333333333L;
  r = ((r>>4) + r) & 0x0f0f0f0f0f0f0f0fL;
  r = (r>>8) + r;
  r = (r>>16) + r;
  r = (r>>32) + r;
  r = r & 127;
  return r;
}

static i64 lf_dna_rankc_sub(lf_dna *lf, i64 i, int c)
{
  i64 j,r;
  i64 i2;
//  i64 i3;
  int lmb;
  bitvec_t x,m,*p;
  static bitvec_t masktbl[4] = {0,0x5555555555555555L,0xAAAAAAAAAAAAAAAAL,0xFFFFFFFFFFFFFFFFL};

//    i3 = i;
  if (i >= lf->last) i--;
  if (i < 0) return 0;

  lmb = lf->lmb;
  r = getuint(lf->RL,(i>>logLB)*SIGMADNA + c,lf->k);
  if (c < SIGMADNA-1) {
    r += lf->RM[(i>>lmb)*(SIGMADNA-1) + c];
  } else {
    r += ((i & (LB-1))>>lmb)<<lmb;
    r -= lf->RM[(i>>lmb)*(SIGMADNA-1) + 0];
    r -= lf->RM[(i>>lmb)*(SIGMADNA-1) + 1];
    r -= lf->RM[(i>>lmb)*(SIGMADNA-1) + 2];
  }
  i2 = (i>>lmb) << lmb;
  p = &lf->BW[i2/(DD/logSIGMADNA)];

  for (j=0; j+(DD/logSIGMADNA)-1 <= (i-i2); j+=(DD/logSIGMADNA)) {
    x = (*p++) ^ masktbl[c];
    x = (x | (x>>1)) & 0x5555555555555555L;
    r += (DD/logSIGMADNA) - POPCOUNT2(x);
  }
  x = (*p) ^ masktbl[c];
  x = (x | (x>>1)) & 0x5555555555555555L;
  m = 0x5555555555555555L >> (((i-i2) - j + 1) * logSIGMADNA);
  x |= m;
  r += (DD/logSIGMADNA) - POPCOUNT2(x);

#if 0
  if (r != lf_dna_rankc_sub_slow(lf,i3,c)) {
	  printf("rankc: error i=%d c=%d rank=%d (%d)\n",i3,c,r,lf_dna_rankc_sub_slow(lf,i3,c));
  }
#endif
  return r;
}


static int convertchar(uchar t)
{
  int c = 0;
  switch (t) {
  case 'a':
  case 'A':  c = 0;  break;
  case 'c':
  case 'C':  c = 1;  break;
  case 'g':
  case 'G':  c = 2;  break;
  case 't':
  case 'T':  c = 3;  break;
  default:  printf("error char = %c [%02x]\n",t,t);
  }
  return c;
}

static void make_tbl(lf_dna *lf)
{
  i64 n;
  i64 i;
  i64 C[SIGMADNA];
  int c,mb;
  bitvec_t *BW;
  int k;

  n = lf->n;
  k = lf->k;
  BW = lf->BW;
  mb = 1 << lf->lmb;

  lf->RL = (uchar *) mymalloc(k*((n+LB-1)/LB+1)*SIGMADNA);
  lf->RM = (u16 *) mymalloc(sizeof(lf->RM[0])*((n+mb-1)/mb+1)*(SIGMADNA-1));

  for (i=0;i<((n+LB-1)/LB)*SIGMADNA;i++) {
    putuint(lf->RL,i,0,k);
  }
  for (i=0;i<((n+mb-1)/mb)*(SIGMADNA-1);i++) lf->RM[i] = 0;

  for (c = 0; c < SIGMADNA; c++) C[c] = 0;

  fprintf(stderr,"making tables...\n");
  for (i = 0; i < n; i++) {
    if (i % 1000000 == 0) {
      fprintf(stderr,"%ld\r",i/1000);
      fflush(stderr);
    }
    if (i % LB == 0) {
      for (c = 0; c < SIGMADNA; c++) {
        putuint(lf->RL,(i / LB)*SIGMADNA + c,C[c],k);
      }
    }
    if (i % mb == 0) {
      for (c = 0; c < SIGMADNA-1; c++) {
        lf->RM[(i / mb)*(SIGMADNA-1) + c]
            = (u16)(C[c] - getuint(lf->RL,(i / LB)*SIGMADNA + c,k));
      }
    }
    c = lf_dna_BW_sub(lf,i);
    C[c]++;
  }
  
}

void lf_dna_options(CSA *csa, char *p)
{
  lf_dna *lf;
  csa->psi_struc = lf = (lf_dna *) mymalloc(sizeof(lf_dna));
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
  printf("BW_DNA L=%d\n",lf->l);
}

static int lf_dna_BW(CSA *csa, i64 i)
{
  lf_dna *lf;
  lf = (lf_dna *)csa->psi_struc;

  if (i == lf->last) return -1; // EOF
  if (i > lf->last) i--;

  return lf_dna_BW_sub(lf,i);
}

static i64 lf_dna_rankc(CSA *csa, i64 i, int c)
{
  lf_dna *lf;

  lf = (lf_dna *)csa->psi_struc;
  return lf_dna_rankc_sub(lf,i,c);
}


i64 lf_dna_makeindex(CSA *csa, char *fname, bool coded)
{
  i64 last;
  FILE *in,*out;
  bitvec_t *BW;
  char *fbw, *flst, *fbw2, *fidx;
  int k,L,c,c2;
  i64 size;
  i64 n,i;
  lf_dna *lf;
  
  lf = (lf_dna *)csa->psi_struc;

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
  sprintf(fbw2,"%s.bw2",fname);
  sprintf(fidx,"%s.bwd",fname);

  in = fopen(flst,"r");
  if (in == NULL) {
    perror("lf_dna_makeindex:");  exit(1);
  }
  int res = fscanf(in,"%ld",&last);
  if (!res) {  }
  printf("last = %ld\n",last);
  lf->last = last;
  fclose(in);


  in = fopen(fbw,"rb");
  if (in == NULL) {
    printf("lf_dna_makeindex: cannot open %s\n",fbw);
    exit(1);
  }
  fseek(in,0,SEEK_END);
  n = ftell(in);
  fseek(in,0,SEEK_SET);
  printf("n=%ld\n",n);
  csa->n = lf->n = n;

  for (i=0; i<SIGMA; i++) csa->C[i] = 0;

  BW = (bitvec_t *) mymalloc(sizeof(*BW) * (n/(DD/logSIGMADNA)+1));
  lf->BW = BW;

  fprintf(stderr,"packing...\n");
  for (i = 0; i < n; i++) {
    if (i % 1000000 == 0) {
      fprintf(stderr,"%ld\r",i/1000);
      fflush(stderr);
    }
    c = fgetc(in);
    csa->C[c]++;
		if (coded==false) {
    	c2 = convertchar(c);
			setbits(BW,i*logSIGMADNA,c2,logSIGMADNA);
		} else {
    	setbits(BW,i*logSIGMADNA,c,logSIGMADNA);
		}
  }
  fclose(in);

  out = fopen(fbw2,"w");
  if (out == NULL) {
    printf("lf_dna_makeindex: cannot open %s\n",fbw2);
  }
  for (i=0; i<n/(DD/logSIGMADNA)+1; i++) {
    writeuint(sizeof(*BW),BW[i],out);
    size += sizeof(*BW);
  }
  fclose(out);

	lf->k = k = (blog(n+1)+1+8-1)/8;

	make_tbl(lf);

	out = fopen(fidx,"w");
  if (out == NULL) {
    printf("lf_dna_makeindex: cannot open %s\n",fidx);
  }

  writeint(1,ID_LF,out);
  writeint(1,k,out); /* #bytes of integer */
  writeint(k,n,out);
  writeint(k,last,out);
  writeint(k,L,out);
  size += 1+1+3*k;

  writeint(1,ID_BWT_DNA,out);
  size += 1;

  for (i=0;i<((n+LB-1)/LB)*SIGMADNA;i++) {
    writeint(k,getuint(lf->RL,i,k),out);
    size += k;
  }
  for (i=0;i<((n+L-1)/L)*(SIGMADNA-1);i++) {
    writeint(sizeof(lf->RM[0]),lf->RM[i],out);
    size += sizeof(lf->RM[0]);
  }
  fclose(out);

  csa->BW = lf_dna_BW;
  csa->rankc = lf_dna_rankc;

  csa->LF = csa_LF;

  free(fbw);
  free(flst);
  free(fbw2);
  free(fidx);

  return size;
}

void lf_dna_read(CSA *csa, char *fname)
{
  char *fbw, *fbwi, *fname2;
  int k,l,id;
  i64 psize1,psize2;
  i64 n;
  lf_dna *lf;
  uchar *p, *q;
  
  csa->psi_struc = lf = (lf_dna *)mymalloc(sizeof(lf_dna));

  k = strlen(fname);
  fname2 = (char *) mymalloc(k-4+1);
  strncpy(fname2,fname,k-4);
  fname2[k-4] = 0;
  k -= 5;

  fbw = (char *) mymalloc(k+5+1);
  fbwi = (char *) mymalloc(k+5+1);

  sprintf(fbwi,"%s.bwd",fname2);

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
    printf("lf_dna_read: id = %d is not supported.\n",id);
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

//  printf("lf_dna_read: psi_id = %d L = %d\n",id,l);
  switch (id) {
    case ID_BWT_DNA:
      printf("#lf format = BWT_DNA\n");
      sprintf(fbw,"%s.bw2",fname2);
      break;
    default:
      printf("lf_dna_read: ID %d is not supported.\n",id);
      break;
  }

  lf->RL = (uchar *)p;  p += ((n+LB-1)/LB)*SIGMADNA * k;
  lf->RM = (u16 *)p;  p += ((n+l-1)/l)*(SIGMADNA-1) * sizeof(lf->RM[0]);

//  printf("lf_dna_read: map %s\n",fbw);
  lf->mapbwt = mymmap(fbw);
  if (lf->mapbwt->addr==NULL) {
    perror("psi1_read: mmap1\n");
    exit(1);
  }
  lf->BW = (bitvec_t *)lf->mapbwt->addr;
  psize2 = lf->mapbwt->len;
//  printf("psize2 = %ld\n",psize2);

//  printf("lf_dna_read: psize1 = %ld psize2 = %ld\n",psize1,psize2);
  lf->psize = psize1 + psize2;

  free(fbw);
  free(fbwi);
  free(fname2);

// user-specific functions
  csa->BW = lf_dna_BW;
  csa->rankc = lf_dna_rankc;

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
