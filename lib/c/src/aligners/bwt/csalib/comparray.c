/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/timeb.h>
#include "comparray.h"

#define CHECK
//#define RANDOM

#define USE_HUFTBL 1

//#define USE_COMPRANK 1
//#define USE_COMPPTR 1


#define NCK (sizeof(bitvec_t)*8)
#define ENUMDECTBL (1<<20)
static bitvec_t *nCkDEC[NCK+1];

#define D (sizeof(bitvec_t)*8)
#define logD 6
//#define logRR 16
#define logRR 15
#define RR (1<<logRR)
#define logRRR 9
#define RRR (1<<logRRR)

#define MAXGRP 8

//static bitvec_t nCk[NCK+1][NCK+1];
static bitvec_t nCk[NCK*MAXGRP+1][NCK*MAXGRP+1];
static int nCk_len[NCK+1];

#define ENUMDECW 8
//#define ENUMDECS (1<<ENUMDECW)

typedef struct {
  int len; // ���������r�b�g��
  int b; // ���������r�b�g����1�̐�
  bitvec_t dec;  // ������������
} enumdectbl;

typedef enumdectbl *enumdectblp;



#define DD 16
#define TBLSIZE (1<<DD)
static int R4[TBLSIZE];



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


static void writeuint(int k,u64 x,FILE *f)
{
  int i;
  for (i=k-1; i>=0; i--) {
    fputc(x & 0xff,f);
    x >>= 8;
  }
}


static int blog(i64 x) // [0,n] �̐����i�[�����ɂ� blog(n)+1 �r�b�g�K�v
{
int l;
  l = -1;
  while (x>0) {
    x>>=1;
    l++;
  }
  return l;
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

  j = i / D;
  l = i % D;
  return (B[j] >> (D-1-l)) & 1;
}

/*
#if 1
static bitvec_t getbits0(bitvec_t *B, i64 i, i64 d)
{
  bitvec_t x;
  int j;
  x = 0;
  for (j=0; j<d; j++) {
    x <<= 1;
    x += getbit(B,i+j);
  }
  return x;
}
#endif


static bitvec_t getbits2(bitvec_t *B, i64 i, i64 d) // �x��
{
  bitvec_t x;

  B += (i >>logD);
  i &= (D-1);
  if (i <= D-d) {
    x = B[0] >> (D-d-i);
  } else {
    x = B[0] << (i+d-D);
    x += B[1] >> (2*D-d-i);
  }
  return x & ((1<<d)-1);
}
*/

static bitvec_t getbits(bitvec_t *B, i64 i, i64 d)
{
  bitvec_t x,z;



//  x0 = getbits0(B,i,d);

  B += (i >>logD);
  i &= (D-1);
  if (i+d <= D) {
    x = B[0];
    x <<= i;
    x >>= (D-d);  // D==64, d==0 ���Ɠ����Ȃ�
  } else {
    x = B[0] << i;
    x >>= D-d;
    z = B[1] >> (D-(i+d-D));
    x += z;
  }
#if 0
  if (x != x0) {
    printf("getbits: x=%lx x0=%lx\n",x,x0);
  }
  return x0;
#else
  return x;
#endif
}


bitvec_t getbits1(bitvec_t *B, i64 i, i64 d)
{
  bitvec_t x,z;

  B += (i >>logD);
  i &= (D-1);
  if (i+d <= D) {
    x = B[0];
    x <<= i;
    x >>= (D-d);
  } else {
    x = B[0] << i;
    x >>= D-d;
    z = B[1] >> (D-(i+d-D));
    x += z;
  }
  return x;
}

static int getbitDD(bitvec_t *B, i64 i)
{
  bitvec_t x,z;

  
//  x0 = getbits(B,i,DD);

  B += (i >>logD);
  i &= (D-1);
  if (i <= D-DD) {
    x = B[0] >> (D-DD-i);
  } else {
    x = B[0] << (i+DD-D);
    z = B[1] >> (2*D-DD-i);
    x += z;
  }
  x &= ((1<<DD)-1);
#if 0
  if (x != x0) {
    printf("getDD: x=%lx x0=%lx\n",x,x0);
  }
#endif
  return x;
}


static int encodegamma(bitvec_t *B,i64 p,i64 x) /* x >= 1 */
{
i64 j,w;
  if (x==0) {
    fprintf(stderr,"encodegamma %ld\n",x);  exit(1);
  }
  w = blog(x)+1;
  for (j=0;j<w-1;j++) setbit(B,p+j,0);
  for (j=w-1;j>=0;j--) setbit(B,p+(w-1)+(w-1)-j,(x >> j)&1);
  return 2*w-1;
}

static int initranktables(void)
{
  i64 i,j;
  for (i = 0; i < DD; i++) {
    for (j = (1<<i); j < (2<<i); j++) {
      R4[j] = DD-1-i;
    }
  }
  R4[0] = DD;
  
  return 0;
}


static int getzerorun(bitvec_t *B,i64 p)
{
i64 w,w2;
#if 0
  w = 0;
  while (getbit(B,p+w)==0) w++;
#else
  w = 0;
  while (1) {
    w2 = R4[getbitDD(B,p)];
    w += w2;
    if (w2 < DD) break;
    p += DD;
  }
#endif
  return w;
}

static int decodegamma(bitvec_t *B,i64 p,i64 *ans)
{
i64 w,w2;
i64 x;

  w = getzerorun(B,p);
#if 0
  x = 1;
  for (i=0;i<w;i++) {
    x <<= 1;
    x += getbit(B,p+w+1+i);
  }
#else
  p += w+1;
//  q = p;
#if 1
  x = 1;
  w2 = w;
  while (w2 > DD) {
    x <<= DD;
    x += getbitDD(B,p);
    p += DD;
    w2 -= DD;
  }
  x <<= w2;
  x += (getbitDD(B,p)>>(DD-w2));
#endif
#if 0
  y = 1L << w;
  if (w > 0) y += getbits(B,q,w);
  if (y != x) {
    printf("y=%lx x=%lx w=%ld\n",y,x,w);
  }
#endif
#endif
  *ans = x;
  return 2*w+1;
}


static unsigned int popcount(bitvec_t x)
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


static bitvec_t *mkenumdecodetable(int n, int m)
{
  bitvec_t i;
  int j,k;
  bitvec_t *p;
  bitvec_t x,y,w;

  p = (bitvec_t*) mymalloc(nCk[m][n] * sizeof(bitvec_t));

  for (i=0; i<nCk[m][n]; i++) {
    x = i;
    w = 0;
    j = m;
    y = x;
    for (k = 0; k < n; k++) {
      w <<= 1;
      if (x >= nCk[j][n-1-k]) {
        x -= nCk[j][n-1-k];
        w++;
        j--;
      }
    }
    p[i] = w;
//    printf("n=%d m=%d x=%ld w=%s\n",n,m,y,bins(w));
  }
  return p;
}

#if 0
static enumdectblp mkenumdecodetable2(int n, int m, int w)
{
  bitvec_t i;
  int j,k,l;
  enumdectblp p;
  bitvec_t x1,x2;
  bitvec_t y1,y2;
  int b,len;

//  printf("mkenumdecodetable2 n=%d m=%d w=%d\n",n,m,w);
  mymalloc(p,1<<w);

  l = blog(nCk[m][n]-1)+1;
//  printf("nck=%ld l=%d\n",nCk[m][n],l);
//  if (l < 20) return NULL;

  for (i=0; i<(1<<w); i++) {
    x1 = i << (l-w);
    x2 = (i << (l-w)) + (1<<(l-w))-1;
    if (x1 > nCk[m][n]-1) {
//      printf("x1 = %ld nck = %ld\n",x1,nCk[m][n]);
      break;
    }
    if (x2 > nCk[m][n]-1) x2 = nCk[m][n]-1;
//    printf("i=%ld x1=%lx x2=%lx\n",i,x1,x2);
    y1 = 0;  j = m;
    for (k = 0; k < n; k++) {
      y1 <<= 1;
      if (x1 >= nCk[j][n-1-k]) {
        x1 -= nCk[j][n-1-k];
        j--;
        y1++;
      }
    }
    y2 = 0;  j = m;
    for (k = 0; k < n; k++) {
      y2 <<= 1;
      if (x2 >= nCk[j][n-1-k]) {
        x2 -= nCk[j][n-1-k];
        j--;
        y2++;
      }
    }
//    printf("y1=%lx y2=%lx\n",i,y1,y2);
    len = 0;  b = 0;
    for (j=l-1; j>=0; j--) {
      if ((y1 & (1<<j)) != (y2 & (1<<j))) break;
      len++;
      if (y1 & (1<<j)) b++;
    }
    if (len == 0) {
//      printf("len == 0\n");
    }
    y1 >>= l-len;
    y1 <<= l-len;
    p[i].len = len;
    p[i].b = b;
    p[i].dec = y1;
//    printf("len=%d b=%d dec=%lx\n",len,b,y1);
  }
  return p;
}
#endif

int encode_enum(bitvec_t *out, i64 i, int n, bitvec_t *in, i64 j)
{
  int k,d,w,m;
  bitvec_t x;
  int len;
  
  len = 0;

  m = 0;
  for (k=0; k<n; k++) m += getbit(in,j+k);

  if (m == D || m == 0) {
    setbits(out, i, blog(D-1)+1, 0);
    len += blog(D-1)+1;
    setbits(out, i+len, 1, m/D);
    len++;
  } else {
    setbits(out, i, blog(D-1)+1, m);
    len += blog(D-1)+1;

    d = 0;
    x = 0;
    for (k = 0; k < n; k++) {
      if (getbit(in,j+k)==1) {
        x += nCk[m][n-1-k];
        d++;
        m--;
      }
    }

    w = blog(nCk[d][n]-1)+1;
    setbits(out,i+len,w,x);
    len += w;
  }
  return len;
}

int encode_enum2(comparray *da, i64 i, int n, bitvec_t *in, i64 j)
{
  int k,d,w,m;
  bitvec_t x;
  int len;
  int clen;
  bitvec_t *out;
  bitvec_t code;

  out = da->buf;
  
  len = 0;

  m = 0;
  for (k=0; k<n; k++) m += getbit(in,j+k);

  clen = da->h->clen[m];
  code = da->h->code[m];
  code >>= sizeof(code)*8 - clen;
  setbits(out, i, clen, code);
  len += clen;

  if (m == D || m == 0) {
    // nothing
  } else {
//    setbits(out, i, blog(D-1)+1, m);
//    len += blog(D-1)+1;

    d = 0;
    x = 0;
    for (k = 0; k < n; k++) {
      if (getbit(in,j+k)==1) {
        x += nCk[m][n-1-k];
        d++;
        m--;
      }
    }

    w = blog(nCk[d][n]-1)+1;
    setbits(out,i+len,w,x);
    len += w;
  }
  return len;
}

int encode_enum3(bitvec_t *out, i64 i, int n, bitvec_t *in, i64 j)
{
  int k,d,w,m;
  bitvec_t x;

  
  m = 0;
  for (k=0; k<n; k++) m += getbit(in,j+k);

  if (m == D || m == 0) {
    w = 0;
  } else {
    d = 0;
    x = 0;
    for (k = 0; k < n; k++) {
      if (getbit(in,j+k)==1) {
        x += nCk[m][n-1-k];
        d++;
        m--;
      }
    }

    w = blog(nCk[d][n]-1)+1;
//    printf("enc i=%ld w=%ld x=%ld\n",i,w,x);
    setbits(out,i,w,x);
  }
  return w;
}



int decode_enum(bitvec_t *in, i64 i, int n, bitvec_t *out, int *r)
{
  int k,d,m,m2,n2;
  bitvec_t x,y,x2;
  int len;


  m = getbits(in,i,logD);
  len = logD;
  if (m == 0) {
    if (getbit(in,i+len)==1) {
      *out = (bitvec_t)-1;
      *r = D;
    } else {
      *out = 0;
      *r = 0;
    }
    len++;
    return len;
  }
  *r = m;

//  d = blog(nCk[m][n]-1)+1;
  d = nCk_len[m];
  x = getbits(in,i+len,d);
  len += d;

  x2 = x;  m2 = m;  n2 = n;
  if (nCkDEC[m] != NULL) {
    y = nCkDEC[m][x];
  } else {
#if 0
    y2 = 0;
    while (n > 0) {
      enumdectblp p;
      p = nCkDEC2[m2][n2];
      n2 -= p[x2 >> (d - ENUMDECW)].len;
      m2 -= p[x2 >> (d - ENUMDECW)].b;
      y2 += p[x2 >> (d - ENUMDECW)].dec;
      x2 -= p[x2 >> (d - ENUMDECW)].dec;
    }
#endif
#if 0
    x2 = x;  m2 = m;
    y = 0;
    for (k = 0; k < n; k++) {
      y <<= 1;
      if (x2 >= nCk[m2][n-1-k]) {
        x2 -= nCk[m2][n-1-k];
        y++;
        m2--;
      }
    }
#endif
#if 1
    x2 = x;  m2 = m;
    y = 0;
    for (k = 0; k < n; k++) {
      int f;
      y <<= 1;
      f = (x2 >= nCk[m2][n-1-k]);
      x2 -= f*nCk[m2][n-1-k];
      y += f;
      m2 -= f;
    }
#endif
#if 0
    y2 = 0;
//    printf("x2 = %lu\n",x2);
    k = 0;
    while (m2 > 0) {
      ll = k;  rr = n-1;
      while (ll <= rr) {
        mm = (ll+rr) / 2;
  //      printf("mm=%d %lu\n",mm,nCk[m2][n-1-mm]);
        if (nCk[m2][n-1-mm] <= x2) rr = mm-1; else ll = mm+1;
      }
      k = rr+1;
    //  printf("hit %d  nck=%lu\n",k,nCk[m2][n-1-k]);
      x2 -= nCk[m2][n-1-k];
      y2 += 1UL << (n-1-k);
      m2--;
      k++;
    }
    y = y2;
//    if (y != y2) {
//      printf("y = %ld y2 = %ld\n",y,y2);
//    }
#endif
  }
  *out = y;

  return len;
}

int decode_enum1(bitvec_t *in, i64 i, int n, bitvec_t *out, int *r)
{
  int k,d,m,m2;
  bitvec_t x,y,x2;
  int len;


  m = getbits(in,i,logD+1);
  len = logD;
  if ((m>>1) == 0) {
    if (m & 1) {
      *out = (bitvec_t)-1;
      *r = D;
    } else {
      *out = 0;
      *r = 0;
    }
    len++;
    return len;
  }
  m >>= 1;
  *r = m;

  d = nCk_len[m];
  x = getbits(in,i+len,d);
  len += d;

  x2 = x;  m2 = m; 
  if (nCkDEC[m] != NULL) {
    y = nCkDEC[m][x];
  } else {
    x2 = x;  m2 = m;
    y = 0;

#if 1
    for (k = 0; k < n; k++) {
      int f;
      y <<= 1;
      f = (x2 >= nCk[m2][n-1-k]);
      x2 -= f*nCk[m2][n-1-k];
      y += f;
      m2 -= f;
    }
#else
    for (k = 0; k < n; k++) {
      y <<= 1;
      if (x2 >= nCk[m2][n-1-k]) {
        x2 -= nCk[m2][n-1-k];
        y++;
        m2--;
      }
    }
#endif
  }
  *out = y;

  return len;
}

int decode_enum2(comparray *da, i64 i, int n, bitvec_t *out, int *r)
{
  int k,d,m,m2;
  bitvec_t x,y,x2;
  int len;

  
  bitvec_t *in;
  
  in = da->buf;

  x = getbits(in,i,32)<<32;
#if USE_HUFTBL
//  m = DecodeHuffman_tbl(da->h, x);
  m = da->h->tbl[x>>(sizeof(u64)*8-HUFTBLWIDTH)];
  if (m < 0) m = DecodeHuffman(da->h, x);
#else
  m = DecodeHuffman(da->h, x);
#endif
//  if (m != m2) printf("decode_enum2: m=%d m2=%d\n",m,m2);
  len = da->h->clen[m];
  
  if (m == 0) {
    *out = 0;
    *r = 0;
    return len;
  }
  if (m == D) {
    *out = (bitvec_t)-1;
    *r = D;
    return len;
  }
  *r = m;

//  d = blog(nCk[m][n]-1)+1;
  d = nCk_len[m];
  x = getbits(in,i+len,d);
  len += d;

  x2 = x;  m2 = m;
  if (nCkDEC[m] != NULL) {
    y = nCkDEC[m][x];
  } else {
#if 1
    y = 0;
#if 1
    for (k = 0; k < n; k++) {
      int f;
      y <<= 1;
      f = (x2 >= nCk[m2][n-1-k]);
      x2 -= f*nCk[m2][n-1-k];
      y += f;
      m2 -= f;
    }
#else
    for (k = 0; k < n; k++) {
      y <<= 1;
      if (x >= nCk[m][n-1-k]) {
        x -= nCk[m][n-1-k];
        y++;
        m--;
      }
    }
#endif
  //  y <<= D-n;
#else
    y2 = 0;
//    printf("x2 = %lu\n",x2);
    k = 0;
    while (m2 > 0) {
      ll = k;  rr = n-1;
      while (ll <= rr) {
        mm = (ll+rr) / 2;
  //      printf("mm=%d %lu\n",mm,nCk[m2][n-1-mm]);
        if (nCk[m2][n-1-mm] <= x2) rr = mm-1; else ll = mm+1;
      }
      k = rr+1;
    //  printf("hit %d  nck=%lu\n",k,nCk[m2][n-1-k]);
      x2 -= nCk[m2][n-1-k];
      y2 += 1UL << (n-1-k);
      m2--;
      k++;
    }
    y = y2;
//    if (y != y2) {
//      printf("y = %ld y2 = %ld\n",y,y2);
//    }
#endif
  }
  *out = y;

  return len;
}

int decode_enum3(bitvec_t *in, i64 i, int n, bitvec_t *out, int m)
{
  int k,d,m2;
  bitvec_t x,y,x2;
  int len;


  if (m == D) {
    *out = (bitvec_t)-1;
    len=0;
    return len;
  }
  if (m == 0) {
    *out = 0;
    len=0;
    return len;
  }
  len = 0;

//  d = blog(nCk[m][n]-1)+1;
  d = nCk_len[m];
  x = getbits(in,i+len,d);
  len += d;

  x2 = x;  m2 = m; 
  if (nCkDEC[m] != NULL) {
    y = nCkDEC[m][x];
  } else {
    x2 = x;  m2 = m;
    y = 0;
    for (k = 0; k < n; k++) {
      y <<= 1;
      if (x2 >= nCk[m2][n-1-k]) {
        x2 -= nCk[m2][n-1-k];
        y++;
        m2--;
      }
    }
  }
  *out = y;

  return len;
}


int skip_enum(bitvec_t *in, i64 i, int n, int *r)
{
  int d,m;
  
  int len;

  m = getbits(in,i,logD);
  len = logD;
  if (m == 0) {
    if (getbit(in,i+len)==1) {
      *r = D;
    } else {
      *r = 0;
    }
    len++;
    return len;
  }
  *r = m;

//  d = blog(nCk[m][n]-1)+1;
  d = nCk_len[m];
  len += d;

  return len;
}

int skip_enum1(bitvec_t *in, i64 i, int n, int *r)
{
  int d,m;
  
  int len;

  m = getbits(in,i,logD+1);
  len = logD;
  if ((m>>1) == 0) {
    if (m & 1) {
      *r = D;
    } else {
      *r = 0;
    }
    len++;
    return len;
  }
  m >>= 1;
  *r = m;

//  d = blog(nCk[m][n]-1)+1;
  d = nCk_len[m];
  len += d;

  return len;
}


int skip_enum2(comparray *da, i64 i, int n, int *r)
{
  int d,m;
  bitvec_t x;
  int len;
  bitvec_t *in;
  
  in = da->buf;

  x = getbits(in,i,32)<<32;
#if USE_HUFTBL
//  m = DecodeHuffman_tbl(da->h, x);
  m = da->h->tbl[x>>(sizeof(u64)*8-HUFTBLWIDTH)];
  if (m < 0) m = DecodeHuffman(da->h, x);
#else
  m = DecodeHuffman(da->h, x);
#endif
//  if (m != m2) printf("skip_enum2: m=%d m2=%d\n",m,m2);
  len = da->h->clen[m];

  if (m == 0) {
    *r = 0;
    return len;
  }
  if (m == D) {
    *r = D;
    return len;
  }
  *r = m;

//  d = blog(nCk[m][n]-1)+1;
  d = nCk_len[m];
  len += d;

  return len;
}

int *enumtbl_w[MAXGRP+1][NCK*MAXGRP];
int *enumtbl_p[MAXGRP+1][NCK*MAXGRP];
int enumtbl_kwp[MAXGRP+1][NCK*MAXGRP];

static void make_tables(void)
{
i64 i,j,k;
i64 n,w;
i64 x,y;
i64 p;
int l, l_min, l_max;
int f;

  if (enumtbl_w[1][1] != NULL) return;

  for (w = 0; w <= NCK*MAXGRP; w++) { // w: ���̍��v�l
    for (k = 1; k <= MAXGRP; k++) { // k: �v�f��
      if (w+k > NCK*MAXGRP) continue;
      if (k*sizeof(bitvec_t)*8 < w) continue;
      n = nCk[k-1][w+k-1]; // n: ���vw��k�g�̐�
      if (n > 0 && n <= (1<<16)) {
//        printf("[k=%ld w=%ld n=%ld] ",k,w,n);
        if (enumtbl_w[k][w]==NULL) enumtbl_w[k][w] = (int *) mymalloc((n*k) * sizeof(int)); // (k,w) �̂Ƃ��̊e�d��(1�̐�)(�擪����i�̘a (i=0,...,k-1)
        if (enumtbl_p[k][w]==NULL) enumtbl_p[k][w] = (int *) mymalloc((n*k) * sizeof(int)); // (k,w) �̂Ƃ��̊e�����̃I�t�Z�b�g
        l_max = -1;
        l_min = 1 << 30;
        for (x = 0; x < n; x++) { // x: �p�^��
//          printf("[k=%ld w=%ld x=%ld] ",k,w,x);
          l = blog(n-1)+1;
          y = x;
          j = k-1;
          p = w+k-1-1;
          f = 1;
          for (i = w+k-1-1; i>=0; i--) {
            if (y >= nCk[j][i]) {
//              printf("%ld ",p-i);
              if (p-i > NCK) {
                f = 0;
                break;
              }
//              enumtbl_w[k][w][x*k+(k-1-j)] = (w+k-1-1)-p - (k-1-j);
              enumtbl_w[k][w][x*k+(k-1-j)] = p-i;
              if (j < k-1) enumtbl_w[k][w][x*k+(k-1-j)] += enumtbl_w[k][w][x*k+(k-1-j)-1];
              enumtbl_p[k][w][x*k+(k-1-j)] = l;
//              printf("len=%d ",blog(nCk[p-i][NCK]-1)+1);
              l += blog(nCk[p-i][NCK]-1)+1;
              p = i-1;
              y -= nCk[j][i];
              j--;
            }
          }
//          printf("%ld ",p-i);
          if (p-i > NCK) {
            f = 0;
          } else {
//            enumtbl_w[k][w][x*k+(k-1-j)] = (w+k-1-1)-p - (k-1-j);
            enumtbl_w[k][w][x*k+(k-1-j)] = p-i;
            if (j < k-1) enumtbl_w[k][w][x*k+(k-1-j)] += enumtbl_w[k][w][x*k+(k-1-j)-1];
            enumtbl_p[k][w][x*k+(k-1-j)] = l;
//            printf("len=%d ",blog(nCk[p-i][NCK]-1)+1);
            l += blog(nCk[p-i][NCK]-1)+1;
          }
          if (f) {
//            printf(" len = %d ",l);
//            printf("w = ");
//            for (i=0; i<k; i++) printf("%d ",enumtbl_w[k][w][x*k+i]);
//            printf("\n");
//            printf("p = ");
//            for (i=0; i<k; i++) printf("%d ",enumtbl_p[k][w][x*k+i]);
//            printf("l=%d\n",l);
            if (l > l_max) l_max = l;
            if (l < l_min) l_min = l;
          }
//          printf("\n");
        }
        enumtbl_kwp[k][w] = l_max;
//        printf("k=%ld w=%ld l_max = %d l_min = %d %f bpc\n",
//                k,w,l_max,l_min,(double)l_max / (k*NCK));
      } else {
        enumtbl_w[k][w] = enumtbl_p[k][w] = NULL;
      }
    }
  }
  
//  fflush(stdout);
}



void comparray_maketbl(void)
{
  int i,j;

  for (i=0; i<=MAXGRP*NCK; i++) for (j=0; j<=MAXGRP*NCK; j++) nCk[i][j] = 0;
  for (i=0; i<=MAXGRP*NCK; i++) nCk[0][i] = 1;
  for (i=0; i<=MAXGRP*NCK; i++) nCk[i][i] = 1;
  for (j=2; j<=MAXGRP*NCK; j++) {
    for (i=1; i<j; i++) {
      if (nCk[i-1][j-1]>0 && nCk[i][j-1]>0) {
        nCk[i][j] = nCk[i-1][j-1]+nCk[i][j-1];
        if (nCk[i][j] > (1L<<(8*sizeof(bitvec_t)-2))) nCk[i][j] = 0;
        if (nCk[j-i][j] == 0) nCk[j-i][j] = nCk[i][j];
      }
    }
  }
  for (i=0; i<=NCK; i++) nCk_len[i] = blog(nCk[i][NCK]-1)+1;
#if 0
  for (j=0; j<=NCK*MAXGRP; j++) {
    for (i=0; i<=j; i++) {
      if (nCk[i][j]>0) printf("(%d,%d) %ld\n",j,i,nCk[i][j]);
    }
    printf("\n");
  }
  for (i=0; i<=NCK; i++) {
    printf("%d ",nCk_len[i]);
  }
  printf("\n");
#endif

#if 1
  for (i=0; i<=NCK; i++) {
//    if (nCkDEC[i] == NULL && nCk[i][NCK] <= ENUMDECTBL) {
    if (nCk[i][NCK] <= ENUMDECTBL) {
      nCkDEC[i] = mkenumdecodetable(NCK,i);
//      printf("DEC[%d] %ld %p\n",i,nCk[i][NCK], nCkDEC[i]);
    } else {
      nCkDEC[i] = NULL;
    }
  }
#endif

#if 0
  for (j=0; j<=NCK; j++) {
    for (i=0; i<=j; i++) {
      nCkDEC2[i][j] = mkenumdecodetable2(j,i,ENUMDECW);
    }
  }
#endif

  make_tables();
  initranktables();

}




void comparray_construct(comparray *da, i64 n, bitvec_t *buf, ushort L, int opt)
{
  i64 i,j,k,k1,k2,m,r;
  
  i64 size;
  i64 pass;
  i64 bs, bss;
  i64 f;
  i64 freq[D+1];
  double ff[D+1];
  bitvec_t buftmp[1],buftmp2[2];
  int runlen, d;
  int *num_ones = NULL, *num_blk = NULL, *w_blk = NULL, blk = 0;
  int w2;


  comparray_maketbl();

  size = sizeof(comparray);
//  printf("densearray-size:0 %ld\n",size);

  m = 0;
  for (i=0; i<n; i++) m += getbit(buf,i);
  da->n = n;
  da->m = m;
  da->L = L;
  da->opt = opt;

//  printf("n=%d m=%d (%f)\n",n,m,(double)m/n);

  da->k2 = ((blog(m)+1)+7)/8;

// rank index
  if (opt & SDARRAY_COMPRANK) {
    da->rsa = (sparsearray * ) mymalloc((1) * sizeof(sparsearray));
    sparsearray_construct_init(da->rsa, n, (n+L-1)/L);
  } else {
    da->rl = (uchar *) mymalloc((da->k2*((n+RR-1)/RR)) * sizeof(uchar));
            size += da->k2 * ((n+RR-1)/RR);
    da->rs = (word *) mymalloc(((n+L-1)/L) * sizeof(word));
            size += sizeof(*(da->rs))*((n+L-1)/L);
  }
  r = 0;
  for (i=0; i<n; i+=RR) {
    if (!(opt & SDARRAY_COMPRANK)) {
//    da->rl[i/RR] = r;
    putuint(da->rl, i/RR, r, da->k2);
    }
    m = 0;
    for (j=0; j<RR; j++) {
      if (j % L == 0 && i+j < n) {
        if (opt & SDARRAY_COMPRANK) {
          sparsearray_construct_set(da->rsa,(i+j)/L, r+m);
        } else {
          da->rs[(i+j)/L] = m;
        }
      }
      if (i+j < n && getbit(buf,i+j)==1) m++;
    }
    r += m;
  }
  if (opt & SDARRAY_COMPRANK) {
    sparsearray_construct_end(da->rsa, SDARRAY_SELECT1);
  }

// pointer
//  printf("densearray-size:3 %ld\n",size);

  for (i=0; i<D+1; i++) freq[i] = 0;

  if (opt & SDARRAY_SUC) {
    num_ones = (int *) mymalloc((L/D+1) * sizeof(int));
    num_blk = (int *) mymalloc((L/D+1) * sizeof(int));
    w_blk = (int *) mymalloc((L/D+1) * sizeof(int));
  }

  for (pass = 1; pass <= 2; pass++) {
    bs = 0;

    for (i=0; i<n; i+=RR) {
      if (pass==2) {
        if (!(opt & SDARRAY_COMPPTR)) {
//        da->pl[i/RR] = bs; // large block�̐擪
          putuint(da->pl, i/RR, bs, da->k1);
        }
      }
      bss = 0;
      for (j=0; j<RR && i+j<n; j+=L) {
        if (pass==2) {
          if (opt & SDARRAY_COMPPTR) {
            sparsearray_construct_set(da->psa, (i+j)/L, bs+bss);
          } else {
            da->ps[(i+j)/L] = bss; // small block�̐擪
          }
        }
        if (opt & SDARRAY_RR) {
          d = getbit(buf,i+j+0);
//          printf("B[%ld] = %d\n",i+j,d);
          if (pass==2) setbit(da->buf,bs+bss, d);
          bss++;
          runlen = 1;
          for (k=1; k<L && i+j+k<n; k++) {
  //          printf("B[%ld] = %d\n",i+j+k,getbit(buf,i+j+k));
            if (getbit(buf,i+j+k) == d) {
              runlen++;
//              printf("runlen = %d\n",runlen);
            } else {
              if (pass==1) {
    //            printf("encode runlen = %d\n",runlen);
                bss += encodegamma(buftmp,0,runlen);
              } else {
      //          printf("encode runlen = %d\n",runlen);
                bss += encodegamma(da->buf,bs+bss,runlen);
              }
              d = 1-d;
              runlen = 1;
            }
          }
          if (runlen > 0) {
            if (pass==1) {
              bss += encodegamma(buftmp,0,runlen);
            } else {
              bss += encodegamma(da->buf,bs+bss,runlen);
            }
          }
        } else if (opt & SDARRAY_SUC) {
          i64 bss2;
          int b;
// �e SB ��1�̐��𐔂���
          for (k=0; k<L && i+j+k<n; k+=D) {
            buftmp[0] = 0;
            f = 0;
            for (k2=0; k2<D && i+j+k+k2<n; k2++) {
              if (getbit(buf,i+j+k+k2)==1) {
                setbit(buftmp,k2,1);
                f++;
              }
            }
            num_ones[k/D] = f;
//              printf("num_ones[%d] = %d\n",k/D,num_ones[k/D]);
          }
// SB�̃O���[�v������
//          printf("grouping i=%ld j=%ld\n",i,j);
          blk = 0;
          k2 = 0;  w2 = 0;
          k1 = 0;
          for (k=0;k<L && i+j+k<n ;k+=D) {
            // [k2,k] ��1�O���[�v�ɂł��邩����
            f = num_ones[k/D];
            if ((k1-k2+1 > MAXGRP) || (w2+f > NCK*MAXGRP)
                || (enumtbl_w[k1-k2+1][w2+f] == NULL)) { // �ł��Ȃ�
              //if (pass==1) bss += enumtbl_kwp[k1-k2][w2];
              num_blk[blk] = k1-k2;
              w_blk[blk] = w2;
//              printf("group %d #blk=%d #ones=%d\n",
//                      blk,num_blk[blk],w_blk[blk]);
              blk++;
              k2 = k1; // ���̃O���[�v�̐擪
              w2 = 0;
            }
            w2 += f;
            k1++;
          }
          if (k1-k2 > 0) { // �Ō���1�O���[�v��
            if ((k1-k2 > MAXGRP) || (w2 > NCK*MAXGRP)
                || (enumtbl_w[k1-k2][w2] == NULL)) { // �ł��Ȃ�
              //printf("??? k1=%d k2=%d w2=%d\n",k1,k2,w2);
              exit(1);
            }
            //if (pass==1) bss += enumtbl_kwp[k1-k2][w2];
            num_blk[blk] = k1-k2;
            w_blk[blk] = w2;
//            printf("group %d #blk=%d #ones=%d\n",
//                    blk,num_blk[blk],w_blk[blk]);
            blk++;
          }
          if (pass == 1) {
            // �w�b�_�[�̃T�C�Y
//            printf("group blk=%d [",blk);
            bss += encodegamma(buftmp2,0,blk);
            for (k=0; k<blk; k++) {
//              printf("(k=%d w=%d) ",num_blk[k],w_blk[k]);
              bss += encodegamma(buftmp2,0,num_blk[k]);
              bss += encodegamma(buftmp2,0,w_blk[k]+1);
            // �{�̂̃T�C�Y
//              printf("body size=%ld ",enumtbl_kwp[num_blk[k]][w_blk[k]]);
              bss += enumtbl_kwp[num_blk[k]][w_blk[k]];
//              printf("block bss=%ld\n",bss);
            }
//            printf("]\n");
          } else {
            // group��encode
//            printf("encode num-group %d bss=%d \n",blk,bss);
            bss += encodegamma(da->buf,bs+bss,blk);
            for (k=0; k<blk; k++) {
//              printf("encode (k,w)=(%d,%d) bss=%d \n",
//                     num_blk[k],w_blk[k],bss);
              bss += encodegamma(da->buf,bs+bss,num_blk[k]);
              bss += encodegamma(da->buf,bs+bss,w_blk[k]+1);
            }
            // group�w�b�_�[��encode
            w2 = 0;
            for (b=0; b<blk; b++) {
              bitvec_t x;
              int ww, kk, p, nn;
              bss2 = bss;

              kk = num_blk[b];
              ww = w_blk[b];
              nn = nCk[kk-1][ww+kk-1];
              x = 0;
              p = ww+kk-1;
//              printf("#ones ");
              for (k=0; k<num_blk[b]; k++) {
//                printf("%d ",num_ones[w2+k]);
                p -= num_ones[w2+k]+1;
                kk--;
//                printf("x += %ld\n",nCk[kk][p]);
                x += nCk[kk][p];
              }
#if 0
{
  i64 y,p,l,i,j,w,k;
          k = num_blk[b];
          w = w_blk[b];
          l = blog(n-1)+1;
          y = x;
          j = k-1;
          p = w+k-1-1;
          printf("decode ");
          for (i = w+k-1-1; i>=0; i--) {
//            printf("y=%ld try (%ld %ld) %ld\n",y,j,i,nCk[j][i]);
            if (y >= nCk[j][i]) {
              printf("%ld ",p-i);
              if (p-i > NCK) {
                break;
              }
              l += blog(nCk[p-i][NCK]-1)+1;
              p = i-1;
              y -= nCk[j][i];
              j--;
            }
          }
          printf("%ld \n",p-i);
  
}
#endif
//              printf("encode header bss=%d len=%d x=%ld\n",bss,blog(n-1)+1,x);
              setbits(da->buf,bs+bss,blog(nn-1)+1,x);
              bss += blog(nn-1)+1;
            // SB��encode
              f = 0;
              for (k=0; k<num_blk[b]; k++) {
                buftmp[0] = 0;
                for (k2=0; k2<D && i+j+w2*D+k2<n; k2++) {
                  if (getbit(buf,i+j+w2*D+k2)==1) {
                    setbit(buftmp,k2,1);
                    f++;
                  }
                }
//                printf("encode enum bss=%d\n",bss);
                bss += encode_enum3(da->buf,bs+bss, D, buftmp, 0);
                w2++;
              }
//              printf("#blk=%ld #ones=%d l_max=%ld  actual=%ld\n",
//                     num_blk[b],w_blk[b],enumtbl_kwp[num_blk[b]][w_blk[b]],bss-bss2);

              bss = bss2 + enumtbl_kwp[num_blk[b]][w_blk[b]];
//              printf("block bss=%ld\n",bss);
            }
          }
        } else { // �ʏ��̕�����
          for (k=0; k<L && i+j+k<n; k+=D) {
            buftmp[0] = 0;
            f = 0;
            for (k2=0; k2<D && i+j+k+k2<n; k2++) {
//            printf("i=%ld j=%ld k=%ld k2=%ld\n",i,j,k,k2);
              if (getbit(buf,i+j+k+k2)==1) {
                if (k2 >= D) {
                  printf("??? k2 = %ld\n",k2);
                }
                setbit(buftmp,k2,1);
                f++;
              }
            }
            if (pass == 1) {
              freq[f]++;
              if (opt & SDARRAY_LENHUF) {
                bss += blog(nCk[f][D]-1)+1;
              } else {
                bss += encode_enum(buftmp2,0, D, buftmp, 0);
              }
            } else {
              if (opt & SDARRAY_LENHUF) {
                bss += encode_enum2(da,bs+bss, D, buftmp, 0);
              } else {
                bss += encode_enum(da->buf,bs+bss, D, buftmp, 0);
              }
            }
          }
        }
      }
      bs += bss;
    }
    if (pass == 1) {
      if (opt & SDARRAY_LENHUF) {
        for (i=0; i<D+1; i++) ff[i] = (double)(freq[i]+1) / n;
        da->h = MakeHuffmanTree(D+1, ff);

        for (i=0; i<D+1; i++) {
          bs += da->h->clen[i] * freq[i];
       if (da->h->clen[i]>32) printf("%lu: len = %i code = %lx\n",i,da->h->clen[i],da->h->code[i]);
        }
      }
      bs = (bs+D-1)/D;
      da->buf = (bitvec_t *) mymalloc(bs * sizeof(bitvec_t));
      da->bufsize = bs;
      size += sizeof(*(da->buf)) * bs;
      if (opt & SDARRAY_COMPPTR) {
        da->psa = (sparsearray *) mymalloc(1 * sizeof(sparsearray));
        sparsearray_construct_init(da->psa, bs*D, (n+L-1)/L);
      } else {
        da->k1 = ((blog(sizeof(*(da->buf)) * bs*8)+1)+7)/8;
        da->pl = (uchar *) mymalloc((da->k1 * ((n+RR-1)/RR)) * sizeof(uchar));
        size += da->k1 * ((n+RR-1)/RR);
        da->ps = (word *) mymalloc(((n+L-1)/L) * sizeof(word));
        size += sizeof(*(da->ps))*((n+L-1)/L);
      }
    }
  }
  if (opt & SDARRAY_COMPPTR) {
    sparsearray_construct_end(da->psa, SDARRAY_SELECT1);
  }

  if (opt & SDARRAY_SUC) {
    free(num_ones);
    free(num_blk);
    free(w_blk);
  }

//
//  printf("k1 = %d k2 = %d\n",da->k1,da->k2);
  da->size = size;
  return;
}


i64 comparray_write(comparray *da, FILE *f)
{
  i64 i;
  ushort L;
  i64 size,size2;
  
  L = da->L;

  size = 0;  size2 = 0;
  writeuint(sizeof(da->n), da->n, f);  size2 += sizeof(da->n);
  writeuint(sizeof(da->m), da->m, f);  size2 += sizeof(da->m);
  writeuint(sizeof(da->opt), da->opt, f);  size2 += sizeof(da->opt);
  writeuint(sizeof(int), RR, f);    size2 += sizeof(int);
  writeuint(sizeof(da->L), da->L, f);  size2 += sizeof(da->L);
//  writeuint(sizeof(i64), size, f);
  writeuint(sizeof(i64), da->bufsize, f);  size2 += sizeof(i64);

  if (da->opt & SDARRAY_LENHUF) Huffman_write(da->h, f);

//  printf("k1 = %d k2 = %d\n",da->k1,da->k2);

  if (!(da->opt & SDARRAY_COMPRANK)) {
    size += da->k2 * ((da->n+RR-1)/RR)
         + sizeof(*da->rs) * ((da->n+L-1)/L);
  }
  if (!(da->opt & SDARRAY_COMPPTR)) {
    size += da->k1 * ((da->n+RR-1)/RR)
         + sizeof(*da->ps) * ((da->n+L-1)/L);
  }
  size += sizeof(*da->buf) * da->bufsize;

  if (da->opt & SDARRAY_COMPRANK) {
    size += sparsearray_write(da->rsa, f);
  } else {
    for (i=0; i<(da->n+RR-1)/RR; i++) {
//    writeuint(sizeof(*da->rl), da->rl[i], f);
//    printf("RL[%ld] = %ld\n",i,getuint(da->rl,i,da->k2));
      writeuint(da->k2, getuint(da->rl,i,da->k2), f);
    }
    for (i=0; i<(da->n+L-1)/L; i++) {
      writeuint(sizeof(*da->rs), da->rs[i], f);
    }
  }
  if (da->opt & SDARRAY_COMPPTR) {
    size += sparsearray_write(da->psa, f);
  } else {
    for (i=0; i<(da->n+RR-1)/RR; i++) {
//    writeuint(sizeof(*da->pl), da->pl[i], f);
      writeuint(da->k1, getuint(da->pl,i,da->k1), f);
    }
    for (i=0; i<(da->n+L-1)/L; i++) {
      writeuint(sizeof(*da->ps), da->ps[i], f);
    }
  }
  for (i=0; i<da->bufsize; i++) {
    writeuint(sizeof(*da->buf), da->buf[i], f);
  }
  return size + size2;
}

void comparray_read(comparray *da, uchar **map)
{
  i64 i;
  uchar *p,*q;
  ushort L;
  i64 rr;

//  comparray_maketbl();
  p = q = *map;

  da->n = getuint(p,0,sizeof(da->n));  p += sizeof(da->n);
  da->m = getuint(p,0,sizeof(da->m));  p += sizeof(da->m);
//  printf("n = %ld m = %ld\n",da->n,da->m);
  da->k2 = ((blog(da->m)+1)+7)/8;

  da->opt = getuint(p,0,sizeof(da->opt));  p += sizeof(da->opt);
//  printf("densearray_read: n=%ld m=%ld\n",da->n,da->m);
//  printf("n = %ld m = %ld opt = %d\n",da->n,da->m,da->opt);

  if (da->opt & SDARRAY_RR || 1) {
//    initranktables();
  }

  rr = getuint(p,0,sizeof(int));  p += sizeof(int);
  if (rr != RR) {
    printf("error1 RR=%ld must be %d\n",rr,RR);
  }
  da->L = L = getuint(p,0,sizeof(da->L));  p += sizeof(da->L);
//  size = getuint(p,0,sizeof(size));  p += sizeof(size);
//  printf("size %ld\n",size);
  da->bufsize = getuint(p,0,sizeof(da->bufsize));  p += sizeof(da->bufsize);
  da->k1 = ((blog(sizeof(*(da->buf)) * da->bufsize*8)+1)+7)/8;

//  printf("k1 = %d k2 = %d\n",da->k1,da->k2);

  if (da->opt & SDARRAY_LENHUF) {
    *map = p;
    da->h = Huffman_read2(map);
    p = *map;
    for (i=0; i<D+1; i++) {
//    printf("%d: len = %d code = %lx\n",i,da->h->clen[i],da->h->code[i]);
      if (da->h->clen[i] > 32) {
        printf("%ld: len = %d code = %lx\n",i,da->h->clen[i],da->h->code[i]);
      }
    }
  }

  if (da->opt & SDARRAY_COMPRANK) {
    da->rsa = (sparsearray *) mymalloc(1 * sizeof(sparsearray));
    *map = p;
    sparsearray_read(da->rsa, map);
    p = *map;
  } else {
    da->rl = (uchar *)p;
    p += da->k2 * ((da->n+rr-1)/rr);
    da->rs = (word *)p;
    p += sizeof(*da->rs) * ((da->n+L-1)/L);
  }
  if (da->opt & SDARRAY_COMPPTR) {
    da->psa = (sparsearray *) mymalloc(1 * sizeof(sparsearray));
    *map = p;
    sparsearray_read(da->psa, map);
    p = *map;
  } else {
    da->pl = (uchar *)p;
    p += da->k1 * ((da->n+rr-1)/rr);
    da->ps = (word *)p;
    p += sizeof(*da->ps) * ((da->n+L-1)/L);
  }
  da->buf = (bitvec_t *)p;
  p += sizeof(*da->buf) * da->bufsize;
  *map = p;
  
  da->size = p - q;

}

int comparray_getbit(comparray *da, i64 i)
{
  bitvec_t x, *buf;
  
  i64 j,k;
  i64 ofs;
  int r,c;
  i64 d,runlen;

  buf = da->buf;
//  ofs = da->pl[i>>logRR];
  if (da->opt & SDARRAY_COMPPTR) {
    ofs = sparsearray_select(da->psa, i/da->L +1);
  } else {
    ofs = getuint(da->pl,i>>logRR, da->k1);
    ofs += da->ps[i/da->L];
  }
  i &= da->L -1;

  if (da->opt & SDARRAY_RR) {
    d = getbit(buf,ofs);
    ofs++;
    if (d == 1) runlen = 1;  else runlen = 0;
//    r += d;
    k = 0;
    c = 0;
    while (k <= i) {
      if (runlen == 1) {
        ofs += decodegamma(buf,ofs,&d);
        if (k+d-1 > i) {
          c = 1;
          break;
        }
        k += d;
      }
      if (k > i) {
        c = 1;
        break;
      }

      ofs += decodegamma(buf,ofs,&d);
      k += d;
      runlen = 1;
    }
    return c;
  } else if (da->opt & SDARRAY_SUC) {
    i64 blk,b;
    i64 k,w,k2 = 0,w2 = 0,j2,nn,l,m;
    i64 ofs2;

    ofs += decodegamma(buf,ofs,&blk);
    ofs2 = 0;  j2 = 0;  j = -1;
    for (b=0; b<blk; b++) {
      ofs += decodegamma(buf,ofs,&k);
      ofs += decodegamma(buf,ofs,&w);  w -= 1;
      if ((j < 0) && (j2+k*D > i)) {
        j = j2;
        k2 = k;  w2 = w;
      }
      if (j < 0) ofs2 += enumtbl_kwp[k][w];
      j2 += k*D;
    }
    ofs += ofs2;
    nn = nCk[k2-1][w2+k2-1];
    l = blog(nn-1)+1;
    x = 0;
    if (l > 0) x = getbits(buf,ofs,l);
//    ofs += l;
    ofs += enumtbl_p[k2][w2][x*k2+(i-j)/D];
    m = enumtbl_w[k2][w2][x*k2+(i-j)/D];
    if ((i-j)/D>0) m -= enumtbl_w[k2][w2][x*k2+(i-j)/D-1];
    decode_enum3(buf, ofs, D, &x, m);
    j = i & (D-1);
    return (x >> (D-1-j)) & 1;
  } else {
    if (da->opt & SDARRAY_LENHUF) {
      for (j=D; j<=i; j+=D) {
        ofs += skip_enum2(da, ofs, D, &r);
      }
      decode_enum2(da, ofs, D, &x, &r);
    } else {
//    buf += (ofs >>logD);
//    ofs &= (D-1);
      for (j=D; j<=i; j+=D) {
        ofs += skip_enum1(buf, ofs, D, &r);
      }
      decode_enum1(buf, ofs, D, &x, &r);
    }
    j = i & (D-1);
    return (x >> (D-1-j)) & 1;
  }
}

i64 comparray_rank(comparray *da, i64 i)
{
  bitvec_t x, *buf;
  
  i64 j,k;
  i64 ofs;
  i64 r;
  int rr;
  i64 d,runlen;

  if (da->opt & SDARRAY_COMPRANK) {
    r = sparsearray_select(da->rsa, i/da->L+1);
  } else {
//  r = da->rl[i>>logRR] + da->rs[i>>logRRR];
    r = getuint(da->rl,i>>logRR,da->k2);
    r += da->rs[i/da->L];
  }
  rr = 0;
  
  buf = da->buf;
//  ofs = da->pl[i>>logRR];
  if (da->opt & SDARRAY_COMPPTR) {
    ofs = sparsearray_select(da->psa, i/da->L +1);
  } else {
    ofs = getuint(da->pl,i>>logRR, da->k1);
    ofs += da->ps[i/da->L];
  }
  i &= (da->L - 1);

  if (da->opt & SDARRAY_RR) {
    d = getbit(buf,ofs);
    ofs++;
    if (d == 1) runlen = 1;  else runlen = 0;
//    r += d;
    k = 0;
    while (k <= i) {
      if (runlen == 1) {
        ofs += decodegamma(buf,ofs,&d);
        if (k+d-1 > i) {
          r += i-k+1;
          break;
        }
        r += d;
        k += d;
      }
      if (k > i) break;

      ofs += decodegamma(buf,ofs,&d);
      k += d;
      runlen = 1;
    }
  } else if (da->opt & SDARRAY_SUC) {
    i64 blk,b;
    i64 k,w,k2 = 0, w2 = 0,j2,nn,l,m;
    i64 ofs2;

    ofs += decodegamma(buf,ofs,&blk);
    ofs2 = 0;  j2 = 0;  j = -1;
    for (b=0; b<blk; b++) {
      ofs += decodegamma(buf,ofs,&k);
      ofs += decodegamma(buf,ofs,&w);  w -= 1;
      if ((j < 0) && (j2+k*D > i)) {
        j = j2;
        k2 = k;  w2 = w;
      }
      if (j < 0) {
        ofs2 += enumtbl_kwp[k][w];
        r += w;
      }
      j2 += k*D;
    }
    ofs += ofs2;
    nn = nCk[k2-1][w2+k2-1];
    l = blog(nn-1)+1;
    x = 0;
    if (l>0) x = getbits(buf,ofs,l);
//    ofs += l;
    ofs += enumtbl_p[k2][w2][x*k2+(i-j)/D];
    m = enumtbl_w[k2][w2][x*k2+(i-j)/D];
    if ((i-j)/D>0) m -= enumtbl_w[k2][w2][x*k2+(i-j)/D-1];
    if ((i-j)/D>0) r += enumtbl_w[k2][w2][x*k2+(i-j)/D-1];
    decode_enum3(buf, ofs, D, &x, m);
    j = i & (D-1);
    r += POPCOUNT(x >> (D-1-j));
  } else {
    if (da->opt & SDARRAY_LENHUF) {
      for (j=D; j<=i; j+=D) {
        ofs += skip_enum2(da, ofs, D, &rr);
        r += rr;
      }
      decode_enum2(da, ofs, D, &x, &rr);
    } else {
      for (j=D; j<=i; j+=D) {
        ofs += skip_enum1(buf, ofs, D, &rr);
        r += rr;
      }
      decode_enum1(buf, ofs, D, &x, &rr);
    }
    j = i & (D-1);
    r += POPCOUNT(x >> (D-1-j));
  }
  return r;
}

i64 comparray_rank_and_bit(comparray *da, i64 i, int *c)
{
  bitvec_t x, *buf;
  
  i64 j,k;
  i64 ofs;
  i64 r;
  int rr;
  i64 d,runlen;

  if (da->opt & SDARRAY_COMPRANK) {
    r = sparsearray_select(da->rsa, i/da->L+1);
  } else {
//  r = da->rl[i>>logRR] + da->rs[i>>logRRR];
    r = getuint(da->rl,i>>logRR,da->k2);
    r += da->rs[i/da->L];
  }
  rr = 0;
  
  buf = da->buf;
//  ofs = da->pl[i>>logRR];
  if (da->opt & SDARRAY_COMPPTR) {
    ofs = sparsearray_select(da->psa, i/da->L +1);
  } else {
    ofs = getuint(da->pl,i>>logRR, da->k1);
    ofs += da->ps[i/da->L];
  }
  i &= da->L -1;
  
  if (da->opt & SDARRAY_RR) {
    d = getbit(buf,ofs);
    ofs++;
    if (d == 1) runlen = 1;  else runlen = 0;
//    r += d;
    k = 0;
    *c = 0;
    while (k <= i) {
      if (runlen == 1) {
        ofs += decodegamma(buf,ofs,&d);
        if (k+d-1 > i) {
          r += i-k+1;
          *c = 1;
          break;
        }
        r += d;
        k += d;
      }
      if (k > i) {
        *c = 1;
        break;
      }

      ofs += decodegamma(buf,ofs,&d);
      k += d;
      runlen = 1;
    }
  } else if (da->opt & SDARRAY_SUC) {
    i64 blk,b;
    i64 k,w,k2 = 0, w2 = 0,j2,nn,l,m;
    i64 ofs2;

    ofs += decodegamma(buf,ofs,&blk);
    ofs2 = 0;  j2 = 0;  j = -1;
    for (b=0; b<blk; b++) {
      ofs += decodegamma(buf,ofs,&k);
      ofs += decodegamma(buf,ofs,&w);  w -= 1;
      if ((j < 0) && (j2+k*D > i)) {
        j = j2;
        k2 = k;  w2 = w;
      }
      if (j < 0) {
        ofs2 += enumtbl_kwp[k][w];
        r += w;
      }
      j2 += k*D;
    }
    ofs += ofs2;
    nn = nCk[k2-1][w2+k2-1];
    l = blog(nn-1)+1;
    x = 0;
    if (l>0) x = getbits(buf,ofs,l);
//    ofs += l;
    ofs += enumtbl_p[k2][w2][x*k2+(i-j)/D];
    m = enumtbl_w[k2][w2][x*k2+(i-j)/D];
    if ((i-j)/D>0) m -= enumtbl_w[k2][w2][x*k2+(i-j)/D-1];
    if ((i-j)/D>0) r += enumtbl_w[k2][w2][x*k2+(i-j)/D-1];
    decode_enum3(buf, ofs, D, &x, m);
    j = i & (D-1);
    x >>= (D-1-j);
    r += POPCOUNT(x);
    *c = (int)(x & 1);
  } else {
    if (da->opt & SDARRAY_LENHUF) {
      for (j=D; j<=i; j+=D) {
        ofs += skip_enum2(da, ofs, D, &rr);
        r += rr;
      }
      decode_enum2(da, ofs, D, &x, &rr);
    } else {
      buf += (ofs >>logD);
      ofs &= (D-1);
      for (j=D; j<=i; j+=D) {
        ofs += skip_enum1(buf, ofs, D, &rr);
        r += rr;
      }
      decode_enum1(buf, ofs, D, &x, &rr);
    }
    j = i & (D-1);
    x >>= (D-1-j);
    r += POPCOUNT(x);
    *c = (int)(x & 1);
  }

  return r;
}



i64 comparray_rank0(comparray *da, i64 i)
{
  return i+1 - comparray_rank(da,i);
}

int comparray_construct_init(comparray *da, i64 n)
{
  i64 i;

  da->n = n;
  da->buf = (bitvec_t *) mymalloc(((n+D-1)/D) * sizeof(bitvec_t));
  for (i=0; i<(n+D-1)/D; i++) da->buf[i] = 0;

  return 0;
}

int comparray_construct_set(comparray *da, i64 i, int x)
{
  if (x > 0) setbit(da->buf,i,x);
  return 0;
}

i64 comparray_construct_end(comparray *da, ushort L, int opt)
{
  bitvec_t *tmp;
  tmp = da->buf;
  comparray_construct(da, da->n, tmp, L, opt);
  free(tmp);
  return da->size;
}


#if 0
void comparray_sb_construct(comparray_sb *da, i64 n, bitvec_t *buf, ushort L, int opt)
{
  i64 i,j,k,k2,m,r;
  
  i64 size;
  i64 pass;
  i64 bs, bss;
  i64 f;
  i64 freq[D+1];
  double ff[D+1];
  bitvec_t buftmp[1],buftmp2[2];

  da->opt = opt;
  da->L = L;

  comparray_maketbl();

  size = sizeof(comparray_sb);

  m = 0;
  for (i=0; i<n; i++) m += getbit(buf,i);

// rank index
  mymalloc(da->rs,(n+L-1)/L);
          size += sizeof(*(da->rs))*((n+L-1)/L);

  for (j=0; j<n; j++) {
    if (j % L == 0) {
      da->rs[j/L] = m;
    }
    if (getbit(buf,j)==1) m++;
  }

// pointer
  for (i=0; i<D+1; i++) freq[i] = 0;
  for (pass = 1; pass <= 2; pass++) {
    bs = 0;

    bss = 0;
    for (j=0; j<n; j+=L) {
      if (pass==2) da->ps[(i+j)/L] = bss; // small block�̐擪
      for (k=0; k<L && j+k<n; k+=D) {
        buftmp[0] = 0;
        f = 0;
        for (k2=0; k2<D && j+k+k2<n; k2++) {
          if (getbit(buf,j+k+k2)==1) {
            setbit(buftmp,k2,1);
           f++;
          }
        }
        if (pass == 1) {
          freq[f]++;
          if (opt & SDARRAY_LENHUF) {
            bss += blog(nCk[f][D]-1)+1;
          } else {
            bss += encode_enum(buftmp2,0, D, buftmp, 0);
          }
        } else {
          if (opt & SDARRAY_LENHUF) {
            bss += encode_enum2(da,bs+bss, D, buftmp, 0);
          } else {
            bss += encode_enum(da->buf,bs+bss, D, buftmp, 0);
          }
        }
      }
    }
    bs += bss;

    if (pass == 1) {
      if (opt & SDARRAY_LENHUF) {
        for (i=0; i<D+1; i++) ff[i] = (double)(freq[i]+1) / n;
        da->h = MakeHuffmanTree(D+1, ff);

        for (i=0; i<D+1; i++) {
          bs += da->h->clen[i] * freq[i];
        }
      }
      bs = (bs+D-1)/D;
      mymalloc(da->buf,bs);
      da->bufsize = bs;
      size += sizeof(*(da->buf)) * bs;
      mymalloc(da->ps,(n+L-1)/L);
      size += sizeof(*(da->ps))*((n+L-1)/L);
    }
  }

}

i64 comparray_sb_write(comparray_sb *da, FILE *f)
{
  i64 i;
  ushort L;
  i64 size,size2;
  
  L = da->L;

  size2 = 0;

  if (da->opt & SDARRAY_LENHUF) Huffman_write(da->h, f);

  size = sizeof(*da->rs) * ((da->n+L-1)/L)
       + sizeof(*da->ps) * ((da->n+L-1)/L)
       + sizeof(*da->buf) * da->bufsize;

  for (i=0; i<(da->n+L-1)/L; i++) {
    writeuint(sizeof(*da->rs), da->rs[i], f);
  }
  for (i=0; i<(da->n+L-1)/L; i++) {
    writeuint(sizeof(*da->ps), da->ps[i], f);
  }
  for (i=0; i<da->bufsize; i++) {
    writeuint(sizeof(*da->buf), da->buf[i], f);
  }
  return size + size2;
}

void comparray_sb_read(comparray_sb *da, uchar **map)
{
  i64 i;
  uchar *p,*q;
  ushort L;
  i64 rr, rrr;

  comparray_maketbl();

  p = q = *map;

  da->n = getuint(p,0,sizeof(da->n));  p += sizeof(da->n);
  da->m = getuint(p,0,sizeof(da->m));  p += sizeof(da->m);
//  printf("n = %ld m = %ld\n",da->n,da->m);

  da->opt = getuint(p,0,sizeof(da->opt));  p += sizeof(da->opt);
//  printf("densearray_read: n=%ld m=%ld\n",da->n,da->m);
  da->L = L = getuint(p,0,sizeof(da->L));  p += sizeof(da->L);
//  size = getuint(p,0,sizeof(size));  p += sizeof(size);
//  printf("size %ld\n",size);
  da->bufsize = getuint(p,0,sizeof(da->bufsize));  p += sizeof(da->bufsize);

//  printf("k1 = %d k2 = %d\n",da->k1,da->k2);

  if (da->opt & SDARRAY_LENHUF) {
    *map = p;
    da->h = Huffman_read2(map);
    p = *map;
    for (i=0; i<D+1; i++) {
//    printf("%d: len = %d code = %lx\n",i,da->h->clen[i],da->h->code[i]);
      if (da->h->clen[i] > 32) {
        printf("%ld: len = %d code = %lx\n",i,da->h->clen[i],da->h->code[i]);
      }
    }
  }

  da->rs = (word *)p;
  p += sizeof(*da->rs) * ((da->n+L-1)/L);
  da->ps = (word *)p;
  p += sizeof(*da->ps) * ((da->n+L-1)/L);
  da->buf = (bitvec_t *)p;
  p += sizeof(*da->buf) * da->bufsize;
  *map = p;
  
//  da->size = p - q;

}



int comparray_sb_getbit(comparray_sb *da, i64 i)
{
  bitvec_t x, *buf;
  
  i64 j;
  i64 ofs;
  int r;

  buf = da->buf;
//  ofs = da->pl[i>>logRR];
  ofs = da->ps[i/da->L];
  i &= da->L -1;
  
  for (j=D; j<=i; j+=D) {
    if (da->opt & SDARRAY_LENHUF) {
      ofs += skip_enum2(da, ofs, D, &r);
    } else {
      ofs += skip_enum(buf, ofs, D, &r);
    }
  }
  if (da->opt & SDARRAY_LENHUF) {
    decode_enum2(da, ofs, D, &x, &r);
  } else {
    decode_enum1(buf, ofs, D, &x, &r);
  }
  j = i & (D-1);
  return (x >> (D-1-j)) & 1;
}

i64 comparray_sb_rank(comparray_sb *da, i64 i)
{
  bitvec_t x, *buf;
  
  i64 j;
  i64 ofs;
  i64 r;
  int rr;

  r = da->rs[i/da->L];
  rr = 0;
  
  buf = da->buf;
//  ofs = da->pl[i>>logRR];
  ofs = da->ps[i/da->L];
  i &= da->L -1;
  
  for (j=D; j<=i; j+=D) {
    if (da->opt & SDARRAY_LENHUF) {
      ofs += skip_enum2(da, ofs, D, &rr);
    } else {
      ofs += skip_enum(buf, ofs, D, &rr);
    }
    r += rr;
  }
  if (da->opt & SDARRAY_LENHUF) {
    decode_enum2(da, ofs, D, &x, &rr);
  } else {
    decode_enum1(buf, ofs, D, &x, &rr);
  }
  j = i & (D-1);
  r += POPCOUNT(x >> (D-1-j));
  return r;
}

i64 comparray_sb_rank0(comparray_sb *da, i64 i)
{
  return i+1 - comparray_sb_rank(da,i);
}

int comparray_sb_construct_init(comparray_sb *da, i64 n)
{
  i64 i;

  da->n = n;
  mymalloc(da->buf, (n+D-1)/D);
  for (i=0; i<(n+D-1)/D; i++) da->buf[i] = 0;

  return 0;
}

int comparray_sb_construct_set(comparray_sb *da, i64 i, int x)
{
  if (x > 0) setbit(da->buf,i,x);
  return 0;
}

i64 comparray_sb_construct_end(comparray_sb *da, ushort L, int opt)
{
  bitvec_t *tmp;
  tmp = da->buf;
  comparray_sb_construct(da, da->n, tmp, L, opt);
  free(tmp);
  return 0;
}

#endif

#ifdef COMPARRAYMAIN
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
  comparray s1,s2;

  i64 i,r,n,m,rr;
  u64 hoge,sum;
  FILE *infp = NULL;
  FILE *out;

  bitvec_t *B;
  byte *B2;
  dword *S,*R;
  MMAP *map;
  uchar *mapp;
  
  int L;

  double t;
  mytimestruct before,after;

#ifdef __SSE4_2__
  printf("SSE4.2 is available.\n");
#else
  printf("no SSE4.2\n");
#endif

  srand(2);

  n = N; // length of bit vector
//  n = 1<<24;
  L = 2048;
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
//      n = atoi(argv[2]);
      L = atoi(argv[2]);
    }
  }
  rr = r;

  mymalloc(B,(n+D-1)/D);

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

  comparray_maketbl();
  initranktables();

// comparray_construct(&s2,n,B, 512, SDARRAY_RANK1 | SDARRAY_LENHUF);
//  comparray_construct(&s2,n,B, L, SDARRAY_RANK1);
//  comparray_construct(&s2,n,B, L, SDARRAY_RANK1 | SDARRAY_RR);
  comparray_construct(&s2,n,B, L, SDARRAY_RANK1 | SDARRAY_SUC);



//  comp
#if 1
  out = fopen("comparraytmp.dat","w");
  comparray_write(&s2, out);
  fclose(out);

  map = mymmap("comparraytmp.dat");
  if (map->addr==NULL) {
    perror("mmap2\n");
    exit(1);
  }
  mapp = (uchar *)map->addr;

  comparray_read(&s1, &mapp);
#endif
  printf("da: used memory: %d bytes (%lf bpc)\n",s1.size,(double)s1.size*8/n);

#ifdef CHECK
  mymalloc(S,n+1);
  mymalloc(R,n+1);
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
  for (i = 0; i < n; i++) {
    if (getbit(B,i) != comparray_getbit(&s1,i)) {
      printf("B[%ld]=%d getbit=%d\n",i,getbit(B,i),comparray_getbit(&s1,i));
    }
  }
#endif

#if 1
  srand(4);

  mygettime(&before);
  hoge = rand();
  sum = 0;
  for (i = 0; i < 100000000; i++) {
//  for (i = 0; i < n; i++) {
    int j;
    //if (i % 1000 == 0) {printf("%ld\r",i);  fflush(stdout);}
    //j = (rand() % n);
#ifdef RANDOM
    j = hoge % n;
#else
    j = i % n;
#endif
#ifdef CHECK
    if (comparray_rank(&s1,j) != R[j] || 0) {
      printf("ERROR: (%d) R[%d] = %d, r = %d\n", getbit(B,i),j, R[j],
                                                 comparray_rank(&s1,j));
    }
    sum += comparray_rank(&s1,j);
#else
    sum += comparray_rank(&s1,j);
#endif
    hoge = hoge * 1566083941U + 1;
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
//  printf("%ld %f comparray rank time sum=%ld\n",rr,t,sum);
  printf("%f %f r=%d L=%d comparray rank time sum=%ld\n",(double)s1.size*8/n,t,rr,L,sum);
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
    if (comparray_select(&s1,j,1) != S[j]) {
      printf("ERROR: S[%d] = %d, s = %d\n",j,S[j],comparray_select(&s1,j,1));
    }
    sum += comparray_select(&s1,j,1);
#else
    sum += comparray_select(&s1,j,1);
#endif
    hoge = hoge * 1566083941U + 1;
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("%ld %f comparray select time sum=%ld\n",rr,t,sum);
#endif



  return 0;
}

#endif //  DENSEMAIN
