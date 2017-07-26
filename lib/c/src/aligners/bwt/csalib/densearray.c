/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/timeb.h>
#include "densearray.h"

#define CHECK
//#define RANDOM

#define logD 6

#define PBS (sizeof(bitvec_t)*8)
#define D (1<<logD)
#define logM 5
#define M (1<<logM)
#define logP 8
#define P (1<<logP)
#define logLL 16    // size of word
#define LL (1<<logLL)
//#define logLLL 7
//#define LLL 128
//#define LLL 32
#define logLLL 5
//#define logLLL 2
#define LLL (1<<logLLL)
//#define logL 10
//#define logL (logLL-3)
#define logL (logLL-1-5)
#define L (1<<logL)

#define logRR 16
#define RR (1<<logRR)
#define logRRR 9
#define RRR (1<<logRRR)

#define DA2_logSB 9
#define DA2_SB (1<<DA2_logSB)
#define DA2_logLB 11
#define DA2_LB (1<<DA2_logLB)
#define DA2_K 4

#ifndef min
 #define min(x,y) ((x)<(y)?(x):(y))
#endif

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


static void writeuint(int k,u64 x,FILE *f)
{
  int i;
  for (i=k-1; i>=0; i--) {
    fputc(x & 0xff,f);
    x >>= 8;
  }
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


static int setbit(bitvec_t *B, i64 i,int x)
{
  i64 j,l;

  j = i / D;
  l = i % D;
  if (x==0) B[j] &= (~(1L<<(D-1-l)));
  else if (x==1) B[j] |= (1L<<(D-1-l));
  else {
    printf("error setbit x=%d\n",x);
    exit(1);
  }
  return x;
}

static i64 setbits(bitvec_t *B, i64 i, int d, i64 x)
{
  int j;

  for (j=0; j<d; j++) {
    setbit(B,i+j,(x>>(d-j-1))&1);
  }
  return x;
}

static int getbit(bitvec_t *B, i64 i)
{
  i64 j,l;

  //j = i / D;
  //l = i % D;
  j = i >> logD;
  l = i & (D-1);
  return (B[j] >> (D-1-l)) & 1;
}

static const unsigned int popCount[] = {
0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};

static unsigned int selecttbl[8*256];

static unsigned int popcount(bitvec_t x)
{
  bitvec_t r;
//  __uint128_t rr;
#if 1
  r = x;
  r = ((r & 0xaaaaaaaaaaaaaaaa)>>1) + (r & 0x5555555555555555);
  r = ((r & 0xcccccccccccccccc)>>2) + (r & 0x3333333333333333);
  r = ((r>>4) + r) & 0x0f0f0f0f0f0f0f0f;
//  r = (r>>8) + r;
//  r = (r>>16) + r;
//  r = ((r>>32) + r) & 127;

  r *= 0x0101010101010101;
  r >>= 64-8;
//  printf("r1 %016lx\n",r);
//  rr = r;
//  rr *= 0x0101010101010101;
//  r = rr >> 56;
//  printf("r2 %016lx\n",r);
//  r &= 0xff;
#else
  r = popCount[x & 0xff];
  x >>= 8;
  r += popCount[x & 0xff];
  x >>= 8;
  r += popCount[x & 0xff];
  x >>= 8;
  r += popCount[x & 0xff];
#if 1
  x >>= 8;
  r += popCount[x & 0xff];
  x >>= 8;
  r += popCount[x & 0xff];
  x >>= 8;
  r += popCount[x & 0xff];
  x >>= 8;
  r += popCount[x & 0xff];
#endif
#endif
  return r;
}

#if 0
static int select_sub(bitvec_t x, int s)
{
  bitvec_t r, w;
  __uint128_t rr;
  int p;

  p = 0;

  r = x;
  r = ((r & 0xaaaaaaaaaaaaaaaa)>>1) + (r & 0x5555555555555555);
  r = ((r & 0xcccccccccccccccc)>>2) + (r & 0x3333333333333333);
  r = ((r>>4) + r) & 0x0f0f0f0f0f0f0f0f;

  rr = r;
  rr *= 0x0101010101010101;
  w = rr >> 56;
  w += 0x8080808080808080;
  w -= s * 0x0101010101010101;
  w &= 0x8080808080808080;
  if ((w & 0xffffffff00000000) == 0) {
    p += 32;
    x <<= 32;
  }
  if ((w & 0xffff000000000000) == 0) {
    p += 16;
    x <<= 16;
  }
  if ((w & 0xff00000000000000) == 0) {
    p += 8;
    x <<= 8;
  }
  return p;
}
#endif

void densearray_make_selecttbl(void)
{
  i64 i,x,r;
  bitvec_t buf[1];
  buf[0] = 0;
  for (x = 0; x < 256; x++) {
    setbits(buf,0,8,x);
    for (r=0; r<8; r++) selecttbl[(r<<8)+x] = -1;
    r = 0;
    for (i=0; i<8; i++) {
      if (getbit(buf,i)) {
        selecttbl[(r<<8)+x] = i;
        r++;
      }
    }
  }

}

int densearray_construct_init(densearray *da, i64 n)
{
  i64 i;

  da->n = n;
  da->buf = (bitvec_t *) mymalloc(((n+D-1)/D) * sizeof(bitvec_t));
  for (i=0; i<(n+D-1)/D; i++) da->buf[i] = 0;

  return 0;
}

int densearray_construct_set(densearray *da, i64 i, int x)
{
  if (x > 0) setbit(da->buf,i,x);
  return 0;
}

i64 densearray_construct_end(densearray *da, ushort l, int opt)
{
  densearray_construct(da, da->n, da->buf, opt);
//  return da->size;
  return da->size + sizeof(*da->buf) * ((da->n+D-1)/D);
}

void densearray_construct(densearray *da, i64 n, bitvec_t *buf, int opt)
{
  i64 i,j,m,k;
  i64 nl;
  i64 p,pp;
  i64 il,is,ml,ms;
  i64 r;
  i64 size;
//  dword *s;
  i64 p1,m1,p2,m2;
  i64 rrr;
  
#if 0
  i64 b;
  i64 freq[RRR+1];
#endif

  size = sizeof(densearray);
//  printf("densearray-size:0 %ld\n",size);

  densearray_make_selecttbl();

  da->rrr = rrr = opt >> 16;
  opt &= 0xffff;

  if (L/LLL == 0) {
    printf("ERROR: L=%d LLL=%d\n",L,LLL);
    exit(1);
  }

  m = 0;
  for (i=0; i<n; i++) m += getbit(buf,i);
  da->n = n;
  da->m = m;

  da->k = k = (blog(m+1)+1+8-1)/8;


//  printf("n=%d m=%d (%f)\n",n,m,(double)m/n);

  da->buf = buf;

//  mymalloc(s,m,0);
  //da->s = s;
  m = 0;
  for (i=0; i<n; i++) {
    if (getbit(buf,i)) {
//      s[m] = i;
      m++;
    }
  }

  if (opt & SDARRAY_SELECT1) {
    nl = (m-1) / L + 1;
    da->lp = (dword *) mymalloc((nl+1) * sizeof(dword));  size += sizeof(*(da->lp))*(nl+1);
    da->p  = (i64 *) mymalloc((nl+1) * sizeof(i64));  size += sizeof(*(da->p))*(nl+1);
//    printf("densearray-size:1 %ld\n",size);

    for (r = 0; r < 2; r++) {
      ml = ms = 0;
      p1 = p2 = -1;
      m1 = m2 = -1;
      for (il = 0; il < nl; il++) {
  //      pp = s[il*L];
        while (m1 < il*L) m1 += getbit(buf,++p1);
        if (p1 >= n) printf("???1 p1 = %ld\n",p1);
  //      if (pp != p1) printf("pp = %ld  p1 = %ld\n",pp,p1);
        pp = p1;
        da->lp[il] = pp;
        i = min((il+1)*L-1,m-1);
  //      p = s[i];
        while (m2 < i) m2 += getbit(buf,++p2);
        if (p2 >= n) printf("???2 p2 = %ld\n",p2);
  //      if (p != p2) printf("p = %ld  p2 = %ld\n",p,p2);
        p = p2;
        //printf("%d ",p-pp);
        if (p - pp >= LL) {
          if (r == 1) {
            for (is = 0; is < L; is++) {
              if (il*L+is >= m) break;
              while (m1 < il*L+is) m1 += getbit(buf,++p1);
              if (p1 >= n) printf("???3 p1 = %ld\n",p1);
  //          if (s[il*L+is] != p1) printf("1:s = %ld  p1 = %ld\n",s[il*L+is],p1);
//            da->sl[ml*L+is] = s[il*L+is];
              da->sl[ml*L+is] = p1;
            }
          }
          da->p[il] = -(ml+1);
          ml++;
        } else {
          if (r == 1) {
            for (is = 0; is < L/LLL; is++) {
              if (il*L+is*LLL >= m) break;
              while (m1 < il*L+is*LLL) m1 += getbit(buf,++p1);
              if (p1 >= n) printf("???4 p1 = %ld\n",p1);
//            if (s[il*L+is*LLL] != p1) 
//              printf("2:s = %ld  p1 = %ld\n",s[il*L+is*LLL],p1);
//            da->ss[ms*(L/LLL)+is] = s[il*L+is*LLL] - pp;
              da->ss[ms*(L/LLL)+is] = p1 - pp;
            }
          }
          da->p[il] = ms;
          ms++;
        }
      }
      if (r == 0) {
        da->ml = ml;  da->ms = ms;
        da->sl = (dword *) mymalloc((ml*L+1) * sizeof(dword));  size += sizeof(*(da->sl))*(ml*L+1);
        da->ss = (word *) mymalloc((ms*(L/LLL)+1) * sizeof(word));  
        size += sizeof(*(da->ss))*(ms*(L/LLL)+1);
//        printf("densearray-size:2 %ld\n",size);
      }
    }
  } else {
    da->lp = NULL;  da->p = NULL;
    da->sl = NULL;  da->ss = NULL;
  }

// rank index
  if (opt & SDARRAY_RANK1) {
//    mymalloc(da->rl,(n+RR-1)/RR,1);  
//      size += sizeof(*(da->rl))*((n+RR-1)/RR);
    da->rl = (uchar *) mymalloc(((n+RR-1)/RR*k) * sizeof(uchar));  
      size += sizeof(*(da->rl))*((n+RR-1)/RR)*k;
    da->rs = (word *) mymalloc(((n+rrr-1)/rrr) * sizeof(word));
      size += sizeof(*(da->rs))*((n+rrr-1)/rrr);
//  printf("densearray-size:3 %ld\n",size);
    r = 0;
    for (i=0; i<n; i+=RR) {
//      da->rl[i/RR] = r;
      putuint(da->rl,i/RR,r,k);
      m = 0;
      for (j=0; j<RR; j++) {
        if (j % rrr == 0 && i+j < n) {
          da->rs[(i+j)/rrr] = m;
        }
        if (i+j < n && getbit(buf,i+j)==1) m++;
      }
      r += m;
    }
  } else {
    da->rl = NULL;  da->rs = NULL;
  }
//
  da->opt = opt;
  da->size = size;
  return;
}

i64 densearray_write(densearray *da, FILE *f)
{
  i64 i,k;
  i64 nl;
  i64 rrr;

  i64 size;

  k = da->k;
  rrr = da->rrr;
  size = 0;
  writeuint(1, da->k, f);
  writeuint(sizeof(da->n), da->n, f);
  writeuint(sizeof(da->m), da->m, f);
  writeuint(1, da->opt, f);
  size += 1 + sizeof(da->n) + sizeof(da->m) + 1;
  if (da->opt & SDARRAY_RANK1) {
    writeuint(sizeof(int), RR, f);
    writeuint(sizeof(int), rrr, f);
    size += 2*sizeof(int);
  }
  if (da->opt & SDARRAY_SELECT1) {
    writeuint(sizeof(da->ml), da->ml, f);  size += sizeof(da->ml);
    writeuint(sizeof(da->ms), da->ms, f);  size += sizeof(da->ms);
  }

//  writeuint(sizeof(i64), size, f);

  if (da->opt & SDARRAY_RANK1) {
    for (i=0; i<(da->n+RR-1)/RR; i++) {
      writeuint(k, getuint(da->rl,i,k), f); size += k;
    }
    for (i=0; i<(da->n+rrr-1)/rrr; i++) {
      writeuint(sizeof(*da->rs), da->rs[i], f);
      size += sizeof(*da->rs);
    }
  }

  if (da->opt & SDARRAY_SELECT1) {
    nl = (da->m-1) / L + 1;
    for (i=0; i<nl+1; i++) {
      writeuint(sizeof(*da->lp), da->lp[i], f);  size += sizeof(*da->lp);
    }
    for (i=0; i<nl+1; i++) {
      writeuint(sizeof(*da->p), da->p[i], f);  size += sizeof(*da->p);
    }
    for (i=0; i<da->ml*L+1; i++) {
      writeuint(sizeof(*da->sl), da->sl[i], f);  size += sizeof(*da->sl);
    }
    for (i=0; i<da->ms*(L/LLL)+1; i++) {
      writeuint(sizeof(*da->ss), da->ss[i], f);  size += sizeof(*da->ss);
    }
  }

  if (!(da->opt & SDARRAY_NOBUF)) {
    for (i=0; i<(da->n+D-1)/D; i++) {
      writeuint(sizeof(*da->buf), da->buf[i], f);
      size += sizeof(*da->buf);
    }
  }
  return size;
}

void densearray_read(densearray *da, uchar **map)
{
  i64 nl;
  uchar *p;

  
  i64 rr = 0, rrr = 0;

//  make_selecttbl();

  p = *map;

  da->k = getuint(p,0,1);  p += 1;
  da->n = getuint(p,0,sizeof(da->n));  p += sizeof(da->n);
  da->m = getuint(p,0,sizeof(da->m));  p += sizeof(da->m);
  da->opt = getuint(p,0,1);  p += 1;
//  printf("densearray_read: n=%ld m=%ld opt=%d\n",da->n,da->m,da->opt);
  if (da->opt & SDARRAY_RANK1) {
    rr = getuint(p,0,sizeof(int));  p += sizeof(int);
    if (rr != RR) {
      printf("error2 RR=%ld must be %d\n",rr,RR);
    }
    rrr = getuint(p,0,sizeof(int));  p += sizeof(int);
    da->rrr = rrr;
//    printf("RRR = %ld\n",rrr);
#if 0
    if (rrr != RRR) {
      printf("error RRR=%ld must be %d\n",rrr,RRR);
    }
#endif
  }
  if (da->opt & SDARRAY_SELECT1) {
    da->ml = getuint(p,0,sizeof(da->ml));  p += sizeof(da->ml);
    da->ms = getuint(p,0,sizeof(da->ms));  p += sizeof(da->ms);
  }
//  size = getuint(p,0,sizeof(size));  p += sizeof(size);
//  printf("size %ld\n",size);

  if (da->opt & SDARRAY_RANK1) {
    da->rl = p;
    p += sizeof(*da->rl) * ((da->n+rr-1)/rr) * da->k;
    da->rs = (word *)p;
    p += sizeof(*da->rs) * ((da->n+rrr-1)/rrr);
  }

  if (da->opt & SDARRAY_SELECT1) {
    nl = (da->m-1) / L + 1;
    da->lp = (dword *)p;
    p += sizeof(*da->lp) * (nl+1);
    da->p = (i64 *)p;
    p += sizeof(*da->p) * (nl+1);
    da->sl = (dword *)p;
    p += sizeof(*da->sl) * (da->ml*L+1);
    da->ss = (word *)p;
    p += sizeof(*da->ss) * (da->ms*(L/LLL)+1);
  }
  
  if (!(da->opt & SDARRAY_NOBUF)) {
    da->buf = (bitvec_t *)p;
    p += sizeof(*da->buf) * ((da->n+D-1)/D);
  }
  *map = p;
}

int densearray_getbit(densearray *da, i64 i)
{
  return getbit(da->buf,i);
}

i64 densearray_rank(densearray *da, i64 i)
{
  i64 r,j;
  bitvec_t *p;
  i64 rrr;
  
  rrr = da->rrr;
//  r = da->rl[i>>logRR] + da->rs[i>>logRRR];
//  r = getuint(da->rl,i>>logRR,da->k) + da->rs[i>>logRRR];
  r = getuint(da->rl,i>>logRR,da->k) + da->rs[i/rrr];
//  p = da->buf + ((i>>logRRR)<<(logRRR-logD));
  p = da->buf + ((i & (~(rrr-1))) >> logD);
  j = i & (rrr-1);
//  if (j < D) r += popcount(*p >> (D-1-j));
//  else r += popcount(*p) + popcount(p[1] >> (D-1-(j-D)));
//  r += popcount(*p >> (D-1-j));
  while (j >= D) {
    r += POPCOUNT(*p++);
    j -= D;
  }
  r += POPCOUNT(*p >> (D-1-j));
  return r;
}

i64 densearray_rank_and_bit(densearray *da, i64 i, int *c)
{
  i64 r,j;
  bitvec_t *p,x;
  i64 rrr;

  rrr = da->rrr;
//  r = da->rl[i>>logRR] + da->rs[i>>logRRR];
//  r = getuint(da->rl,i>>logRR,da->k) + da->rs[i>>logRRR];
  r = getuint(da->rl,i>>logRR,da->k) + da->rs[i/rrr];
//  p = da->buf + ((i>>logRRR)<<(logRRR-logD));
  p = da->buf + ((i & (~(rrr-1))) >> logD);
  j = i & (rrr-1);
//  if (j < D) r += popcount(*p >> (D-1-j));
//  else r += popcount(*p) + popcount(p[1] >> (D-1-(j-D)));
//  r += popcount(*p >> (D-1-j));
  while (j >= D) {
    r += POPCOUNT(*p++);
    j -= D;
  }
  x = *p >> (D-1-j);
  r += POPCOUNT(x);
  *c = (int)(x & 1);

  return r;
}

i64 densearray_rank0(densearray *da, i64 i)
{
  return i+1 - densearray_rank(da,i);
}

i64 densearray_select(densearray *da, i64 i,int f)
{
  
  //  dword *s;
  i64 p,r;
  i64 il;
  i64 rr;
  bitvec_t x;
  bitvec_t *q;

  if (i == 0) return -1;

#if 0
  if (i > da->m) {
    printf("ERROR: m=%d i=%d\n",da->m,i);
    exit(1);
  }
#endif

  i--;

  il = da->p[i>>logL];
  if (il < 0) {
    il = -il-1;
    p = da->sl[(il<<logL)+(i & (L-1))];
  } else {
    p = da->lp[i>>logL];
    p += da->ss[(il<<(logL-logLLL))+(i & (L-1))/LLL];
    r = i - (i & (LLL-1));

    q = &(da->buf[p>>logD]);

    if (f == 1) {
      rr = p & (D-1);
      r -= POPCOUNT(*q >> (D-1-rr));
      p = p - rr;
      
      while (1) {
        rr = POPCOUNT(*q);
        if (r + rr >= i) break;
        r += rr;
        p += D;
        q++;
      }
      
      x = *q;
      while (1) {
        //rr = popcount(x >> (D-8));
//        rr = popCount[x >> (D-8)];
        rr = POPCOUNT(x >> (D-8));
        //rr = popcount8(x >> (D-8));
        if (r + rr >= i) break;
        r += rr;
        p += 8;
        x <<= 8;
      }
      p += selecttbl[((i-r-1)<<8)+(x>>(D-8))];
    } else {
      rr = p & (D-1);
      r -= POPCOUNT((~(*q))  >> (D-1-rr));
      p = p - rr;
      
      while (1) {
        rr = POPCOUNT(~(*q));
        if (r + rr >= i) break;
        r += rr;
        p += D;
        q++;
      }
      
      x = ~(*q);

      while (1) {
        //rr = popcount(x >> (D-8));
//        rr = popCount[x >> (D-8)];
        rr = POPCOUNT(x >> (D-8));
        //rr = popcount8(x >> (D-8));
        if (r + rr >= i) break;
        r += rr;
        p += 8;
        x <<= 8;
      }
      p += selecttbl[((i-r-1)<<8)+(x>>(D-8))];
    }
  }
  return p;
}

i64 densearray_select2(densearray *da, i64 i,int f)
{
  
  //  dword *s;
  i64 p,r;
  i64 il;
  i64 rr;
  bitvec_t x;
  bitvec_t *q;
  

  if (i == 0) return -1;

#if 0
  if (i > da->m) {
    printf("ERROR: m=%d i=%d\n",da->m,i);
    exit(1);
  }
#endif

  i--;

  il = da->p[i>>logL];
  if (il < 0) {
    il = -il-1;
    p = da->sl[(il<<logL)+(i & (L-1))];
  } else {
    p = da->lp[i>>logL];
    p += da->ss[(il<<(logL-logLLL))+(i & (L-1))/LLL];
    r = i - (i & (LLL-1));

    q = &(da->buf[p>>logD]);

    if (f == 1) {
      rr = p & (D-1);
      r -= POPCOUNT(*q >> (D-1-rr));
      p = p - rr;
      
      while (1) {
        rr = POPCOUNT(*q);
        if (r + rr >= i) break;
        r += rr;
        p += D;
        q++;
      }
      
      x = *q;
#if 1
      while (1) {
        //rr = popcount(x >> (D-8));
//        rr = popCount[x >> (D-8)];
        rr = POPCOUNT(x >> (D-8));
        //rr = popcount8(x >> (D-8));
        if (r + rr >= i) break;
        r += rr;
        p += 8;
        x <<= 8;
      }
      p += selecttbl[((i-r-1)<<8)+(x>>(D-8))];
#elif 0
     mask = 0xffffffff00000000L;
     rr = POPCOUNT(x & mask);
     if (r + rr < i) {
       r += rr;
       x <<= 32;
       p += 32;
     }
     mask = 0xffff000000000000L;
     rr = POPCOUNT(x & mask);
     if (r + rr < i) {
       r += rr;
       x <<= 16;
       p += 16;
     }
     mask = 0xff00000000000000L;
     rr = POPCOUNT(x & mask);
     if (r + rr < i) {
       r += rr;
       x <<= 8;
       p += 8;
     }
     p += selecttbl[((i-r-1)<<8)+(x>>(D-8))];
#else
     p += select_sub(x,i-r);
#endif
    } else {
      rr = p & (D-1);
      r -= POPCOUNT((~(*q))  >> (D-1-rr));
      p = p - rr;
      
      while (1) {
        rr = POPCOUNT(~(*q));
        if (r + rr >= i) break;
        r += rr;
        p += D;
        q++;
      }
      
      x = ~(*q);

      while (1) {
        //rr = popcount(x >> (D-8));
//        rr = popCount[x >> (D-8)];
        rr = POPCOUNT(x >> (D-8));
        //rr = popcount8(x >> (D-8));
        if (r + rr >= i) break;
        r += rr;
        p += 8;
        x <<= 8;
      }
      p += selecttbl[((i-r-1)<<8)+(x>>(D-8))];
    }
  }
  return p;
}



// rank(min{i | B[i]=1 and i >= x})
i64 densearray_succ_rank(densearray *da, i64 x)
{
  i64 r;
  if (x == 0) r = 1;
  else r = densearray_rank(da,x-1) + 1;
  
  return r;
}

// min{i | B[i]=1 and i >= x}
i64 densearray_succ(densearray *da, i64 x)
{
  i64 s0;
  i64 s;

// DEBUG
  if (x == 0) s0 = 1;
  else s0 = densearray_select(da,densearray_rank(da,x-1) + 1,1);
//

  s = s0;
  
  return s;
}

#define CHILD(i,k) (((i)+1)*DA2_K+(k))
#define PARENT(i) ((i)/DA2_K-1)

void densearray2_construct(densearray2 *da, i64 n, bitvec_t *buf, int opt)
{
  i64 i,j,m,d,k,c,r;
  i64 size;
  i64 ns, nl, nl2, ofs;

  size = sizeof(densearray2);
  printf("densearray2-size:0 %ld\n",size);

  densearray_make_selecttbl();

  m = 0;
  for (i=0; i<n; i++) m += getbit(buf,i);
  da->n = n;
  da->m = m;
  printf("n=%ld m=%ld\n",n,m);

  da->buf = buf;

  ns = (n+DA2_SB-1)/DA2_SB; // SBÇÃêî
  nl = (n+DA2_LB-1)/DA2_LB; // LBÇÃêî

  d = -1; // ñÿÇÃê[Ç≥
  k = 1; // Ç†ÇÈê[Ç≥ÇÃéqÇÃêî
  ofs = 0;
  while (k < nl) {
    d++;
    k *= DA2_K;
    ofs += k;
  }
  da->leaf_ofs = ofs;
  nl2 = nl + ofs;

  da->rs = (word *) mymalloc(ns * sizeof(word));  size += sizeof(*(da->rs))*ns;
  da->rl = (dword *) mymalloc(nl2 * sizeof(dword));  size += sizeof(*(da->rl))*nl2;
  printf("densearray2-size:1 %ld\n",size);

#if 0
  for (i=0; i<n; i+=DA2_SB) {
    r = 0;
    for (j=0; j<DA2_SB; j++) {
      if (i+j < n) r += getbit(buf,i+j);
    }
    da->rs[i/DA2_SB] = r;
  }

  k = DA2_LB / DA2_SB;  // 1Ç¬ÇÃLBÇ†ÇΩÇËÇÃSBÇÃêî
  for (i=0; i<ns; i+=k) {
    r = 0;
    for (j=0; j<k; j++) {
      if (i+j < ns) r += da->rs[i+j];
    }
    da->rl[ofs + i/k] = r;
  }

  for (i=ofs-1; i>=0; i--) {
    r = 0;
    for (j=0; j<DA2_K; j++) {
      c = CHILD(i,j);
      if (c < nl2) r += da->rl[c];
    }
    da->rl[i] = r;
  }
#else

  r = 0;
  for (i=0; i<n; i+=DA2_LB) {
    da->rl[ofs + i/DA2_LB] = r;
    m = 0;
    for (j=0; j<DA2_LB; j++) {
      if (j % DA2_SB == 0 && i+j < n) {
        da->rs[(i+j)/DA2_SB] = m;
      }
      if (i+j < n && getbit(buf,i+j)==1) m++;
    }
    r += m;
  }

  for (i=ofs; i<nl2; i++) {
//    printf("rl[%ld] = %ld\n",i,da->rl[i]);
  }
  for (i=ofs-1; i>=0; i--) {
    c = CHILD(i,0);
  //  printf("p=%ld c=%ld\n",i,c);
    if (c < nl2) da->rl[i] = da->rl[c];
    else da->rl[i] = n+1;
    //printf("rl[%ld] = %ld\n",i,da->rl[i]);
  }
  for (i=0; i<ofs; i++) {
//    printf("rl[%ld] = %ld\n",i,da->rl[i]);
  }
#endif
  da->opt = opt;
  da->size = size;
  return;
}

i64 densearray2_rank(densearray2 *da, i64 i)
{
  i64 r,j;
  
  bitvec_t *buf;

#if 0
  il = (i >> DA2_logLB) + da->leaf_ofs;
  r = 0;
  while (il >= 0) {
    p = PARENT(il);
    c = CHILD(p,0);
    while (c < il) {
      r += da->rl[c];
      c++;
    }
    il = p;
  }

  is = i >> DA2_logSB;
  c = (i & (~(DA2_LB-1))) >> DA2_logSB;
  while (c < is) {
    r += da->rs[c];
    c++;
  }
#else
  r = da->rl[da->leaf_ofs + (i>>DA2_logLB)] + da->rs[i>>DA2_logSB];
#endif

  buf = da->buf + ((i>>DA2_logSB)<<(DA2_logSB-logD));
  j = i & (DA2_SB-1);
#if DA2_SB > D
  while (j >= D) {
    r += POPCOUNT(*buf++);
    j -= D;
  }
#endif
  r += POPCOUNT(*buf >> (D-1-j));
  return r;
}

i64 densearray2_select(densearray2 *da, i64 i,int f)
{
  i64 j,k;
  //  dword *s;
  i64 p,r,c,ns;
  
  i64 rr;
  bitvec_t x;
  bitvec_t *q,*q0;
  
  i64 ofs;
  i64 i0;

  i0 = i;

  if (i == 0) return -1;

#if 0
  if (i > da->m) {
    printf("ERROR: m=%d i=%d\n",da->m,i);
    exit(1);
  }
#endif

//  i--;

  ofs = da->leaf_ofs;
  p = -1;
  while (p < ofs) {
    c = CHILD(p,0);
    for (j=0; j<DA2_K-1 && c+j+1 < ofs; j++) {
      if (da->rl[c+j+1] >= i) {
//        j--;
        break;
      }
    }
//    if (j == DA2_K) j--;
    p = c+j;
  }
  i -= da->rl[p];
  p -= ofs;

  p <<= (DA2_logLB - DA2_logSB);
  k = DA2_LB / DA2_SB;
  ns = (da->n+DA2_SB-1)/DA2_SB;
  for (j=0; j<k-1 && p+j+1 < ns; j++) {
    if (da->rs[p+j+1] >= i) {
//      j--;
      break;
    }
  }
//  if (j == k) j--;
  p += j;
  i -= da->rs[p];
  p <<= DA2_logSB;

    r = 0;

    q = &(da->buf[p>>logD]);
    q0 = q;

    if (f == 1) {
      while (1) {
        rr = POPCOUNT(*q);
        if (r + rr >= i) break;
        r += rr;
        p += D;
        q++;
      }
#if 0
      if (q - q0 > DA2_SB) {
        printf("too long i=%ld q0=%p q=%p\n",i,q0,q);
      }
#endif      
      x = *q;
      while (1) {
        rr = POPCOUNT(x >> (D-8));
        if (r + rr >= i) break;
        r += rr;
        p += 8;
        x <<= 8;
      }
      p += selecttbl[((i-r-1)<<8)+(x>>(D-8))];
    } else {
      while (1) {
        rr = POPCOUNT(~(*q));
        if (r + rr >= i) break;
        r += rr;
        p += D;
        q++;
      }
      
      x = ~(*q);

      while (1) {
        rr = POPCOUNT(x >> (D-8));
        if (r + rr >= i) break;
        r += rr;
        p += 8;
        x <<= 8;
      }
      p += selecttbl[((i-r-1)<<8)+(x>>(D-8))];
    }

  return p;
}





#ifdef DENSEMAIN
typedef struct timeb mytimestruct;

void mygettime(mytimestruct *t)
{
  ftime(t);
}

double mylaptime(mytimestruct *before,mytimestruct *after)
{
  double t;
  t = after->time - before->time;
  t += (double)(after->millitm - before->millitm)/1000;
  return t;
}


#define N 10485760

int main(int argc, char *argv[])
{
  densearray s1,s2;
//  densearray2 s2;

  i64 i,r,n,m,rr;
  u64 hoge,sum;
  FILE *infp = NULL;
  FILE *out;

  bitvec_t *B;
  byte *B2;
  dword *S,*R;
  MMAP *map;
  uchar *mapp;

  double t;
  mytimestruct before,after;

#ifdef __SSE4_2__
  printf("SSE4.2 is available.\n");
#else
  printf("no SSE4.2\n");
#endif

  srand(2);

  n = N; // length of bit vector
  r = 2; // ratio of ones
  m = 0; // number of ones
  if (argc >= 2){
    r = atoi(argv[1]);
    if (r == 0) {
      infp = fopen(argv[2],"rb");
      if (infp == NULL){
        printf("cannot open %s\n",argv[2]);
        return -1;
      }
      fseek(infp,0,SEEK_END);
      n = ftell(infp);
      rewind(infp);
      //printf("n: %d\n",n);
    } else if (argc >= 3) {
      n = atoi(argv[2]);
    }
  }
  rr = r;

  mymalloc(B,(n+PBS-1+P)/PBS,0);

  if (!infp) {
    m = 0;
    for (i = 0; i < n; i++) {
      if (rand() % 100 < r) {
        setbit(B,i,1);
        m++;
      } else {
        setbit(B,i,0);
      }
    }
  } else {
    m = 0;
    for (i = 0; i < n; i++){
      int c = fgetc(infp);
      if (c == EOF){
        printf("unexpected error at reading from %s\n",argv[2]);
        return -1;
      }
      if (c == '1' || c == '(') {
        setbit(B,i,1);
        m++;
      } else if (c == '0' || c == ')') {
        setbit(B,i,0);
      } else {
        printf("unexpected error (2) at reading from %s\n",argv[2]);
        return -1;
      }
    }
  }


  densearray_construct(&s2,n,B, (128<<16)+(SDARRAY_RANK1 | SDARRAY_SELECT1));
  printf("da: used memory: %d bytes (%lf bpc)\n",s2.size,(double)s2.size*8/n);

#if 1
  out = fopen("sparsearraytmp.dat","w");
  densearray_write(&s2, out);
  fclose(out);

  map = mymmap("sparsearraytmp.dat");
  if (map->addr==NULL) {
    perror("mmap2\n");
    exit(1);
  }
  mapp = (uchar *)map->addr;

  densearray_read(&s1, &mapp);
#endif


#if 0
  densearray2_construct(&s2,n,B, 0);
  printf("da2: used memory: %d bytes (%lf bpc)\n",s2.size,(double)s2.size*8/n);
#endif


#ifdef CHECK
  mymalloc(S,n+1,0);
  mymalloc(R,n+1,0);
  r = 0;
  S[r] = -1;
  for (i=0; i<n; i++) {
    if (getbit(B,i)) {
      r++;
      S[r] = i;
    }
    R[i] = r;
  }
#endif

#if 1
  srand(3);

  mygettime(&before);
  hoge = rand();
  sum = 0;
  for (i = 0; i < 100000000; i++) {
    i64 j;
    //j = (rand() % r)+1;
#ifdef RANDOM
    j = hoge % m + 1;
#else
    j = i % m + 1;
#endif
#ifdef CHECK
    if (densearray_select(&s1,j,1) != S[j]) {
      printf("ERROR: S[%d] = %d, s = %d\n",j,S[j],densearray_select(&s1,j,1));
    }
    sum += densearray_select(&s1,j,1);
#else
    sum += densearray_select(&s1,j,1);
#endif
    hoge = hoge * 1566083941U + 1;
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("%ld %f densearray select time sum=%ld\n",rr,t,sum);
#endif

#if 0
  srand(3);

  mygettime(&before);
  hoge = rand();
  sum = 0;
  for (i = 0; i < 100000000; i++) {
    i64 j;
    //j = (rand() % r)+1;
#ifdef RANDOM
    j = hoge % m + 1;
#else
    j = i % m + 1;
#endif
#ifdef CHECK
    if (densearray_select2(&s1,j,1) != S[j]) {
      printf("ERROR: S[%d] = %d, s = %d\n",j,S[j],densearray_select2(&s1,j,1));
    }
    sum += densearray_select2(&s1,j,1);
#else
    sum += densearray_select2(&s1,j,1);
#endif
    hoge = hoge * 1566083941U + 1;
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("%ld %f densearray select2 time sum=%ld\n",rr,t,sum);
#endif

#if 0
  srand(3);

  mygettime(&before);
  hoge = rand();
  sum = 0;
  for (i = 0; i < 100000000; i++) {
    i64 j;
    //j = (rand() % r)+1;
#ifdef RANDOM
    j = hoge % m + 1;
#else
    j = i % m + 1;
#endif
#ifdef CHECK
    if (densearray2_select(&s2,j,1) != S[j]) {
      printf("ERROR: S[%d] = %d, s = %d\n",j,S[j],densearray2_select(&s2,j,1));
    }
    sum += densearray2_select(&s2,j,1);
#else
    sum += densearray2_select(&s2,j,1);
#endif
    hoge = hoge * 1566083941U + 1;
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("%ld %f densearray2 select time sum=%ld\n",rr,t,sum);
#endif

#if 1
  srand(4);

  mygettime(&before);
  hoge = rand();
  sum = 0;
  for (i = 0; i < 100000000; i++) {
    int j;
    //j = (rand() % n);
#ifdef RANDOM
    j = hoge % n;
#else
    j = i % n;
#endif
#ifdef CHECK
    if (densearray_rank(&s1,j) != R[j]) {
      printf("ERROR: (%d) R[%d] = %d, r = %d\n", getbit(B,i),j, R[j],
                                                 densearray_rank(&s1,j));
    }
    sum += densearray_rank(&s1,j);
#else
    sum += densearray_rank(&s1,j);
#endif
    hoge = hoge * 1566083941U + 1;
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("%ld %f densearray rank time sum=%ld\n",rr,t,sum);
#endif

#if 0
  srand(4);

  mygettime(&before);
  hoge = rand();
  sum = 0;
  for (i = 0; i < 100000000; i++) {
    int j;
    //j = (rand() % n);
#ifdef RANDOM
    j = hoge % n;
#else
    j = i % n;
#endif
#ifdef CHECK
    if (densearray2_rank(&s2,j) != R[j]) {
      printf("ERROR: (%d) R[%d] = %d, r = %d\n", getbit(B,i),j, R[j],
                                                 densearray2_rank(&s2,j));
    }
    sum += densearray2_rank(&s2,j);
#else
    sum += densearray2_rank(&s2,j);
#endif
    hoge = hoge * 1566083941U + 1;
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("%ld %f densearray2 rank time sum=%ld\n",rr,t,sum);
#endif

  return 0;
}

#endif //  DENSEMAIN
