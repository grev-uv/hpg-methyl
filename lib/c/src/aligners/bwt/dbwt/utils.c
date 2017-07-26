#include "utils.h"

int blog(ulong x)
{
int l;
  l = -1;
  while (x>0) {
    x>>=1;
    l++;
  }
  return l;
}

int setbits0(pb *B, ulong i, int d, ulong x)
{
  int j;

  for (j=0; j<d; j++) {
    setbit(B,i+j,(x>>(d-j-1))&1);
  }
  return x;
}

int getbit(pb *B, ulong i)
{
  ulong j,l;

  //j = i / D;
  //l = i % D;
  j = i >> _logD_;
  l = i & (_D_-1);
  return (B[j] >> (_D_-1-l)) & 1;
}

#if 1
dword getbits(pb *B, ulong i, int d)
{
  dword x,z;

  B += (i >> _logD_);
  i &= (_D_-1);
  if (i+d <= 2*_D_) {
    x = (((dword)B[0]) << _D_) + B[1];
    x <<= i;
    x >>= (_D_*2-1-d);
    x >>= 1;
  } else {
    x = (((dword)B[0])<<_D_)+B[1];
    z = (x<<_D_)+B[2];
    x <<= i;
    x &= ((1L<<_D_)-1)<<_D_;
    z <<= i;
    z >>= _D_;
    x += z;
    x >>= (2*_D_-d);
  }

	return x;
}
#else
dword getbits(pb *B, ulong i, int d)
{
  dword j,x;

  x = 0;
  for (j=0; j<d; j++) {
    x <<= 1;
    x += getbit(B,i+j);
  }
  return x;
}
#endif

int setbits(pb *B, ulong i, int d, ulong x)
{
  //ulong j;
  ulong y,m;
  int d2;
  //ulong ii,xx;
  //int dd;
  //pb *BB;

  //  BB = B;  ii = i;  dd = d;  xx = x;

  B += (i>>_logD_);
  i %= _D_;

  while (i+d > _D_) {
    d2 = _D_-i; // x の上位 d2 ビットを格納
    y = x >> (d-d2);
    m = (1<<d2)-1;
    *B = (*B & (~m)) | y;
    B++;  i=0;
    d -= d2;
    x &= (1<<d)-1; // x の上位ビットを消去
  }

  m = (1<<d)-1;
  y = x << (_D_-i-d);
  m <<= (_D_-i-d);
  *B = (*B & (~m)) | y;

#if 0
  if (getbits(BB,ii,dd) != xx) {
    printf("setbits ??? x=%ld %ld\n",xx,getbits(BB,ii,dd));
  }
#endif
  return x;
}

pb *allocate_vector(ulong n)
{
  ulong i,x;
  pb *b;

  x = (n+PBS-1)/PBS;
  b = (pb *) mymalloc(x*sizeof(pb));
  for (i=0; i<x; i++) b[i] = 0;
  return b;
}

packed_array *allocate_packed_array(ulong n, int w)
{
  ulong i,x;
  packed_array *p;

  if (w >= 32) {
    printf("warning: w=%d\n",w);
  }

	p = (packed_array *) mymalloc(sizeof(packed_array));
  p->n = n;  p->w = w;
  x = (n / PBS)*w + ((n % PBS)*w + PBS-1) / PBS;
  p->b = (pb *) mymalloc((x+1)*sizeof(pb));
  for (i=0; i<x+1; i++) p->b[i] = 0;

  return p;
}

void free_packed_array(packed_array *p)
{
  ulong x;
  x = (p->n/PBS+1) * p->w;
  myfree(p->b,x*sizeof(pb));
  myfree(p,sizeof(packed_array));
}

ulong pa_get(packed_array *p, ulong i)
{
  int w;
  pb *b;

  w = p->w;
  b = p->b + (i>>_logD_)*w;
  i = (i % _D_)*w;

  return (ulong) getbits(b,i,p->w);
}

void pa_set(packed_array *p, ulong i, long x)
{
  int w;
  pb *b;
#if 0
  if (x < 0 || x > (1<<p->w)) {
    printf("pa_set: x=%ld w=%d\n",x,p->w);
  }
  if (i < 0 || i >= p->n) {
    printf("pa_set: i=%ld n=%d\n",i,p->n);
  }
#endif
  w = p->w;
  b = p->b + (i>>_logD_)*w;
  i = (i % _D_)*w;

  setbits(b,i,p->w,x);
}

int setbit(pb *B, ulong i,int x)
{
  ulong j,l;

  j = i / _D_;
  l = i % _D_;
  if (x==0) B[j] &= (~(1<<(_D_-1-l)));
  else if (x==1) B[j] |= (1<<(_D_-1-l));
  else {
    printf("error setbit x=%d\n",x);
    exit(1);
  }
  return x;
}
/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
