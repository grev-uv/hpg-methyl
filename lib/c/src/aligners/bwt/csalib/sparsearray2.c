/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>
#include <math.h>
#include "sparsearray2.h"

#define CHECK
//#define RANDOM

#define BLOCK (1<<12)  

static int blog(i64 x) // [0,n] の数を格納するには blog(n)+1 ビット必要
{
int l;
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

static void putuint(uchar *s, i64 i, i64 x, i64 w)
{
  i64 j;
  s += i*w;
  for (j=0; j<w; j++) {
    *s++ = x & 0xff;
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


static int setbit(bitvec_t *B, i64 i,int x)
{
  i64 j,l;

  j = i / PBS;
  l = i % PBS;
  if (x==0) B[j] &= (~(1L<<(PBS-1-l)));
  else if (x==1) B[j] |= (1L<<(PBS-1-l));
  else {
    printf("error setbit x=%d\n",x);
    exit(1);
  }
  return x;
}

static int setbits(bitvec_t *B, i64 i, i64 d, i64 x)
{
  i64 j;

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
  l = i & (PBS-1);
  return (B[j] >> (PBS-1-l)) & 1;
}

static bitvec_t getbits(bitvec_t *B, i64 i, i64 d)
{
  bitvec_t x,z;

  B += (i >> logD);
  i &= (PBS-1);
  if (i+d <= PBS) {
    x = B[0];
    x <<= i;
    x >>= (PBS-1-d);
    x >>= 1;
  } else {
    x = B[0] << i;
    x >>= PBS-d;
    z = B[1] >> (PBS-(i+d-PBS));
    x += z;
  }
  
  return x;
}

i64 sparsearray_sb_construct(sparsearray_sb *sa, i64 n, bitvec_t *buf, int opt)
{
  i64 i,m,r,size;

  m = 0;
  for (i=0; i<n; i++) {
    if (getbit(buf,i)==1) m++;
  }

  sparsearray_sb_construct_init(sa, n, m);
  r = 0;
  for (i=0; i<n; i++) {
    if (getbit(buf,i)==1) {
      sparsearray_sb_construct_set(sa, r, i);
      r++;
    }
  }
  size = sparsearray_sb_construct_end(sa,opt);
  return size;
}

int sparsearray_sb_construct_init(sparsearray_sb *sa, i64 n, i64 m)
{
  i64 i;
  
  
  
  i64 d,mm;
  

  sa->n = n;
  sa->m = m;

  if (m >= n/4) {
    d = 0;
    sa->hi = (bitvec_t *) mymalloc(((n+PBS-1)/PBS+1) * sizeof(bitvec_t));
//    printf("hi = %p size = %ld\n",sa->hi,sizeof(*sa->hi)*((n+PBS-1)/PBS+1));
    sa->low = NULL;
//    for (i=0; i<n; i++) setbit(sa->hi,i,0);
    for (i=0; i<(n+PBS-1)/PBS+1; i++) sa->hi[i] = 0;
  } else {
    mm = m;
    d = 0;
    while (mm < n) {
      mm <<= 1;
      d++;
    }
    sa->hi = (bitvec_t*) mymalloc(((2*m+PBS-1)/PBS+1) * sizeof(bitvec_t));
    sa->low = (bitvec_t*) mymalloc(((d*m+PBS-1)/PBS+1) * sizeof(bitvec_t));
//    for (i=0; i<m*2; i++) setbit(sa->hi,i,0);
//    for (i=0; i<m; i++) setbits(sa->low,i*d,d,0);
    for (i=0; i<(2*m+PBS-1)/PBS+1; i++) sa->hi[i] = 0;
    for (i=0; i<(d*m+PBS-1)/PBS+1; i++) sa->low[i] = 0;
  }
//  printf("sparsearray_sb_construct_init: n=%ld m=%ld d=%ld\n",n,m,d);
  sa->d = d;

  return 0;
}

int sparsearray_sb_construct_set(sparsearray_sb *sa, i64 i, i64 x)
{
  i64 d;
  d = sa->d;

//  printf("setbit(%ld)\n",(x>>d)+i);
#if 0
  if (i >= 0 && i <= 10) {
    printf("3:set i=%ld x=%ld\n",i,x);
  }
#endif
  if (d > 0) {
    setbit(sa->hi,(x>>d)+i,1);
    setbits(sa->low,i*d,d,x & ((1L<<d)-1));
  } else {
//    printf("setbit %ld\n",x);
    setbit(sa->hi,x,1);
  }
  return 0;
}

i64 sparsearray_sb_construct_end(sparsearray_sb *sa, int opt)
{
  bitvec_t *buf2;
  i64 i,n,m,d,n2;
  i64 size;
  int opt_da1, opt_da0;

  n = sa->n;
  m = sa->m;
  d = sa->d;
  size = sizeof(sparsearray_sb);
//  printf("sparsearray-size:0 %ld\n",size);
  opt_da1 = opt_da0 = 0;
  if (d > 0) {
    n2 = 2*m;
    size += sizeof(*(sa->low)) * ((d*m+PBS-1)/PBS);
    //printf("sparsearray-size:1 %ld\n",size);
    if (opt & SDARRAY_SELECT1) opt_da1 |= SDARRAY_SELECT1;
    if (opt & SDARRAY_RANK1) opt_da0 |= SDARRAY_SELECT1 | SDARRAY_NOBUF;
  } else {
    n2 = n;
    if (opt & SDARRAY_SELECT1) opt_da1 |= SDARRAY_SELECT1;
    if (opt & SDARRAY_RANK1) opt_da1 |= SDARRAY_RANK1;
  }
  size += sizeof(*(sa->hi)) * ((n2+PBS-1)/PBS);
  //printf("sparsearray-size:2 %ld\n",size);

  buf2 = sa->hi;

  sa->sd1 = NULL;
  sa->sd1 = (densearray_sb *) mymalloc(1 * sizeof(densearray_sb));  size += sizeof(*(sa->sd1));
  size += densearray_sb_construct(sa->sd1,n2,buf2,opt_da1);

  sa->sd0 = NULL;
  if ((opt & SDARRAY_RANK1) && (d>0)) {
    for (i=0; i<n2; i++) setbit(buf2,i,1-getbit(buf2,i));
    sa->sd0 = (densearray_sb *) mymalloc(1 * sizeof(densearray_sb));  size += sizeof(*(sa->sd0));
    size += densearray_sb_construct(sa->sd0,n2,buf2,opt_da0);

    for (i=0; i<n2; i++) setbit(buf2,i,1-getbit(buf2,i));
  }
  //printf("sparsearray-size:4 %ld\n",size);
  return size;
}

i64 sparsearray_sb_write(sparsearray_sb *sa, int opt, FILE *f)
{

  i64 i,n,m,d,n2;
  i64 size;
  int opt_da1, opt_da0;

  size = 0;
  writeuint(sizeof(sa->n), sa->n, f);  size += sizeof(sa->n);
  writeuint(sizeof(sa->m), sa->m, f);  size += sizeof(sa->m);
  writeuint(sizeof(sa->d), sa->d, f);  size += sizeof(sa->d);
  n = sa->n;  m = sa->m;  d = sa->d;
//  printf("sparsearray_sb_write: n=%ld m=%ld d=%ld opt=%d\n",n,m,d,opt);

  opt_da1 = opt_da0 = 0;
  if (d > 0) {
    n2 = 2*m;
    for (i=0; i<(d*m+PBS-1)/PBS; i++) {
      writeuint(sizeof(*(sa->low)), sa->low[i], f);
    }
    size += sizeof(*(sa->low)) * ((d*m+PBS-1)/PBS);
    if (opt & SDARRAY_SELECT1) opt_da1 |= SDARRAY_SELECT1;
    if (opt & SDARRAY_RANK1) opt_da0 |= SDARRAY_SELECT1 | SDARRAY_NOBUF;
  } else {
    n2 = n;
    if (opt & SDARRAY_SELECT1) opt_da1 |= SDARRAY_SELECT1;
    if (opt & SDARRAY_RANK1) opt_da1 |= SDARRAY_RANK1;
  }
  
  size += densearray_sb_write(sa->sd1,n2,opt_da1, f);
  if ((opt & SDARRAY_RANK1) && (d>0)) {
    size += densearray_sb_write(sa->sd0,n2,opt_da0, f);
  }
  return size;
}

i64 sparsearray_sb_read(sparsearray_sb *sa, int opt, uchar **map)
{
  
  i64 n,m,d,n2;
  uchar *p;
  int opt_da1, opt_da0;

  p = *map;

  sa->n = n = getuint(p,0,sizeof(sa->n));  p += sizeof(sa->n);
  sa->m = m = getuint(p,0,sizeof(sa->m));  p += sizeof(sa->m);
  sa->d = d = getuint(p,0,sizeof(sa->d));  p += sizeof(sa->d);
//  printf("sparsearray_sb_read: n=%ld m=%ld d=%ld opt=%d\n",n,m,d,opt);

  opt_da1 = opt_da0 = 0;
  if (d > 0) {
    n2 = 2*m;
    sa->low = (bitvec_t *)p;
    p += sizeof(*(sa->low)) * ((d*m+PBS-1)/PBS);
    if (opt & SDARRAY_SELECT1) opt_da1 |= SDARRAY_SELECT1;
    if (opt & SDARRAY_RANK1) opt_da0 |= SDARRAY_SELECT1 | SDARRAY_NOBUF;
  } else {
    n2 = n;
    if (opt & SDARRAY_SELECT1) opt_da1 |= SDARRAY_SELECT1;
    if (opt & SDARRAY_RANK1) opt_da1 |= SDARRAY_RANK1;
  }

  *map = p;

  sa->sd1 = (densearray_sb *) mymalloc(1 * sizeof(densearray_sb));
  densearray_sb_read(sa->sd1,n2,opt_da1, map);
  sa->hi = sa->sd1->b;

  if ((opt & SDARRAY_RANK1) && (d>0)) {
    sa->sd0 = (densearray_sb *) mymalloc(1 * sizeof(densearray_sb));
    densearray_sb_read(sa->sd0,n2,opt_da0, map);
    sa->sd0->b = sa->sd1->b;
  }

  return 0;
}



i64 sparsearray_sb_select(sparsearray_sb *sa, i64 i)
{
  i64 d,x;


#if 0
  if ((sa->opt & DENSEARRAY_SELECT1) == 0) {
    printf("sparsearray_select: select1 is not supported.\n");
  }
  if (i > sa->m) {
    printf("ERROR: m=%d i=%d\n",sa->m,i);
    exit(1);
  }
#endif

  if (i == 0) return -1;

  d = sa->d;
  if (d > 0) {
    x = densearray_sb_select(sa->sd1,i,1) - (i-1);
    x <<= d;
    x += getbits(sa->low,(i-1)*d,d);
  } else {
    x = densearray_sb_select(sa->sd1,i,1);
  }
  return x;

}

i64 sparsearray_sb_rank(sparsearray_sb *sa, i64 i)
{
  i64 d,x,w,y;
  i64 j;

#if 0
  if ((sa->opt & DENSEARRAY_RANK1) == 0) {
    printf("sparsearray_select: rank1 is not supported.\n");
  }
  if (i > sa->n) {
    printf("ERROR: n=%d i=%d\n",sa->n,i);
    exit(1);
  }
#endif

  d = sa->d;
  if (d > 0) {
    y = densearray_sb_select(sa->sd0,i>>d,0)+1;
    x = y - (i>>d);
    j = i - ((i>>d)<<d);

    while (1) {
      if (getbit(sa->hi,y)==0) break;
      w = getbits(sa->low,x*d,d);
      if (w >= j) {
        if (w == j) x++;
        break;
      }
      x++;
      y++;
    }
  } else {
    x = densearray_sb_rank(sa->sd1,i);
  }
  return x;
}

i64 sparsearray_sb_rank0(sparsearray_sb *sa, i64 i)
{
  return i+1 - sparsearray_sb_rank(sa,i);
}

i64 sparsearray_sb_rank_bit(sparsearray_sb *sa, i64 i, int *c)
{
  i64 d,x,w,y;
  i64 j;

#if 0
  if ((sa->opt & DENSEARRAY_RANK1) == 0) {
    printf("sparsearray_select: rank1 is not supported.\n");
  }
  if (i > sa->n) {
    printf("ERROR: n=%d i=%d\n",sa->n,i);
    exit(1);
  }
#endif

  d = sa->d;
  if (d > 0) {
    y = densearray_sb_select(sa->sd0,i>>d,0)+1;
    x = y - (i>>d);
    j = i - ((i>>d)<<d);

    *c  = 0;
    while (1) {
      if (getbit(sa->hi,y)==0) break;
      w = getbits(sa->low,x*d,d);
      if (w >= j) {
        if (w == j) {x++; *c = 1;}
        break;
      }
      x++;
      y++;
    }
  } else {
    x = densearray_sb_rank_bit(sa->sd1,i,c);
  }
  return x;
}


int sparsearray_sb_getbit(sparsearray_sb *sa, i64 i)
{
  i64 d,x,w,y;
  i64 j;
  int b;

#if 0
  if ((sa->opt & DENSEARRAY_RANK1) == 0) {
    printf("sparsearray_select: rank1 is not supported.\n");
  }
  if (i > sa->n) {
    printf("ERROR: n=%d i=%d\n",sa->n,i);
    exit(1);
  }
#endif

  d = sa->d;
  if (d > 0) {
    y = densearray_sb_select(sa->sd0,i>>d,0)+1;
    x = y - (i>>d);
    j = i - ((i>>d)<<d);
    b = 0;

    while (1) {
      if (getbit(sa->hi,y)==0) break;
      w = getbits(sa->low,x*d,d);
      if (w >= j) {
        if (w == j) b = 1;
        break;
      }
      x++;
      y++;
    }
  } else {
    b = densearray_sb_getbit(sa->sd1,i);
  }
  return b;
}


i64 sparsearray4_construct(sparsearray4 *sa, i64 n, bitvec_t *buf, int opt)
{
  i64 i,m,r,size;

  m = 0;
  for (i=0; i<n; i++) {
    if (getbit(buf,i)==1) m++;
  }

  sparsearray4_construct_init(sa, n, m);
  r = 0;
  for (i=0; i<n; i++) {
    if (getbit(buf,i)==1) {
      sparsearray4_construct_set(sa, r, i);
      r++;
    }
  }
  size = sparsearray4_construct_end(sa,0, opt);
  return size;
}

int sparsearray4_construct_init(sparsearray4 *sa, i64 n, i64 m)
{
  i64 i;
  
  
  
  i64 d,mm;
  

  sa->n = n;
  sa->m = m;
  sa->k = (blog(n+1)+1+8-1)/8;

  mm = m;
  d = 0;
  while (mm < n) {
    mm <<= 1;
    d++;
  }
//  printf("sparsearray4_construct_init: n=%ld m=%ld d=%ld\n",n,m,d);
  sa->d = d;

  sa->hi = (bitvec_t*) mymalloc(((2*m+PBS-1)/PBS+1) * sizeof(bitvec_t));
  sa->low = (bitvec_t*) mymalloc(((d*m+PBS-1)/PBS+1) * sizeof(bitvec_t));

//  for (i=0; i<m*2; i++) setbit(sa->hi,i,0);
//  for (i=0; i<m; i++) setbits(sa->low,i*d,d,0);
  for (i=0; i<(2*m+PBS-1)/PBS+1; i++) sa->hi[i] = 0;
  for (i=0; i<(d*m+PBS-1)/PBS+1; i++) sa->low[i] = 0;

  return 0;
}

int sparsearray4_construct_set(sparsearray4 *sa, i64 i, i64 x)
{
  i64 d;
  d = sa->d;

//  printf("setbit(%ld)\n",(x>>d)+i);
#if 0
  if (i >= 0 && i <= 10) {
    printf("4:set i=%ld x=%ld hi=%ld low=%ld \n",i,x,(x>>d)+i,x & ((1L<<d)-1));
  }
#endif
  setbit(sa->hi,(x>>d)+i,1);
  setbits(sa->low,i*d,d,x & ((1L<<d)-1));

  return 0;
}

i64 sparsearray4_construct_end(sparsearray4 *sa, ushort L, int opt)
{
  
  i64 i,n,m,d;
  i64 size;
  i64 r,r2;
  i64 b,bi,s,e,x,xs,xe,r0;
  i64 k;

  size = 0;

  sa->opt = opt;

  n = sa->n;
  m = sa->m;
  d = sa->d;
  k = sa->k;

//  printf("sparsearray4_end: n=%ld m=%ld d=%ld opt=%d\n",n,m,d,opt);
  r0 = 0;
  for (i=0; i<2*m; i++)  r0 += getbit(sa->hi,i);
//  r0 = m;
//  printf("r0=%ld\n",r0);


  b = (m+BLOCK-1)/BLOCK;
  sa->sa = (sparsearray_sb *) mymalloc(b * sizeof(sparsearray_sb));  size += sizeof(*sa->sa) * b;
  sa->s  = (uchar *) mymalloc((b*k) * sizeof(uchar));  size += sizeof(*sa->s) * b*k;

  r = 0;  s = -1;  xs = 0;
  for (bi = 0; bi < b; bi++) {
    if (bi % 100 == 0) {printf("%ld/%ld \r",bi,b);  fflush(stdout);}
    s++; // bi 番目のブロックの先頭のビット位置
//    sa->s[bi] = xs;
    putuint(sa->s,bi,xs,sa->k);
    //printf("s[%ld] = %ld\n",bi,xs);
    e = s;  r2 = r;
    while (r2 < (bi+1)*BLOCK && r2 < r0 && e < 2*m) {
      //printf("B[%ld] = %ld\n",e,getbit(sa->hi,e));
      r2 += getbit(sa->hi,e++);
    }
    e--;
    xe = ((e - (r2-1)) << d) + getbits(sa->low,(r2-1)*d,d);
    //printf("[%ld,%ld] n=%ld m=%ld\n",xs,xe,xe-xs+1,r2-r);
    // e は bi番目のブロックの最後の位置
    // bi番目のブロック [s,e] 1の数 r2-r
    sparsearray_sb_construct_init(&sa->sa[bi], xe-xs+1, r2-r);
  //  printf("[%ld,%ld] d = %ld\n",s,e,sa->sa[bi].d);

    e = s; r2 = r;
    for (i = 1; i <= BLOCK; i++) {
      while (r2 < bi*BLOCK + i && r2 < r0 && e < 2*m) {
        //printf("r2=%ld B[%ld] = %ld\n", r2, e,getbit(sa->hi,e));
        r2 += getbit(sa->hi,e++);
      }
      e--;
      if (r2 < bi*BLOCK + i) break; // ベクトルの最後まで行った
      x = ((e - (r2-1)) << d) + getbits(sa->low,(r2-1)*d,d);
      if (x > xe) {
        printf("??? x = %ld xe = %ld\n",x,xe);
        exit(1);
      }
      sparsearray_sb_construct_set(&sa->sa[bi],i-1, x-xs);
      e++;
    }
    e--;
    size += sparsearray_sb_construct_end(&sa->sa[bi],opt);
    if (bi == 0) {
//      for (i=0; i<10; i++) printf("b[%ld] = %ld\n",i,sa->sa[bi].sd1->b[i]);
    }
    s = e;
    xs = xe+1;
    r = r2;
  }


  free(sa->hi);
  free(sa->low);
  return size;
}

i64 sparsearray4_write(sparsearray4 *sa, FILE *f)
{
  i64 n,m,d,k;
  i64 size;
  i64 b,bi;

  size = 0;

  k = sa->k;
  n = sa->n;
  m = sa->m;
  d = sa->d;
  writeuint(1,k,f);  size += 1;
  writeuint(k, sa->n, f);  size += k;
  writeuint(k, sa->m, f);  size += k;
  writeuint(k, sa->d, f);  size += k;
  writeuint(1, sa->opt, f);  size += 1;

  printf("sparsearray4_write: n=%ld m=%ld d=%ld opt=%d\n",n,m,d,sa->opt);


  b = (m+BLOCK-1)/BLOCK;

  for (bi = 0; bi < b; bi++) {
//    writeuint(sizeof(sa->s[0]), sa->s[bi], f);  size += sizeof(sa->s[0]);
    writeuint(k, getuint(sa->s,bi,k), f);  size += k;
  }

  for (bi = 0; bi < b; bi++) {
    size += sparsearray_sb_write(&sa->sa[bi],sa->opt, f);
  }

  return size;
}

i64 sparsearray4_read(sparsearray4 *sa, uchar **map)
{
  
  i64 i,n,m,b,d,k;
  
  uchar *p;
  
  p = *map;

  sa->k = k = getuint(p,0,1);  p += 1;
  sa->n = n = getuint(p,0,k);  p += k;
  sa->m = m = getuint(p,0,k);  p += k;
  sa->d = d = getuint(p,0,k);  p += k;
  sa->opt = getuint(p,0,1);  p += 1;
//  printf("sparsearray4_read: n=%ld m=%ld d=%ld opt=%d\n",n,m,d,sa->opt);

  b = (m+BLOCK-1)/BLOCK;
  sa->sa = (sparsearray_sb *) mymalloc(b * sizeof(sparsearray_sb));

  sa->s = p;
  p += sizeof(sa->s[0]) * b * k;
  *map = p;

  for (i = 0; i < b; i++) {
    sparsearray_sb_read(&sa->sa[i],sa->opt, map);
  }
//  for (i=0; i<10; i++) printf("b[%ld] = %ld\n",i,sa->sa[0].sd1->b[i]);

  return 0;
}


i64 sparsearray4_select(sparsearray4 *sa, i64 i)
{
  i64 x,a,r;

#if 0
  if ((sa->opt & DENSEARRAY_SELECT1) == 0) {
    printf("sparsearray_select: select1 is not supported.\n");
  }
  if (i > sa->m) {
    printf("ERROR: m=%d i=%d\n",sa->m,i);
    exit(1);
  }
#endif

  if (i == 0) return -1;

  a = (i+BLOCK-1) / BLOCK - 1;
  r = (i-1) % BLOCK + 1;
//  if (r == 0) r = BLOCK;
//  x = sa->s[a] + sparsearray_sb_select(&sa->sa[a],r);
  x = getuint(sa->s,a,sa->k) + sparsearray_sb_select(&sa->sa[a],r);
  
  return x;

}

i64 sparsearray4_rank(sparsearray4 *sa, i64 i)
{
  i64 b,bi,r, ll,rr,mm;

  b = (sa->m+BLOCK-1)/BLOCK;
  ll = 0;  rr = b-1;
  while (ll <= rr) {
    mm = (ll + rr) / 2;
    if (getuint(sa->s, mm, sa->k) <= i) ll = mm + 1; else rr = mm - 1;
  }
  bi = ll - 1;
  r = getuint(sa->s, bi, sa->k);
  i -= r;
  return bi*BLOCK + sparsearray_sb_rank(&sa->sa[bi],i);
}

i64 sparsearray4_rank0(sparsearray4 *da, i64 i)
{
  return i+1 - sparsearray4_rank(da,i);
}

int sparsearray4_getbit(sparsearray4 *sa, i64 i)
{
  i64 b,bi,r, ll,rr,mm;

  b = (sa->m+BLOCK-1)/BLOCK;
  ll = 0;  rr = b-1;
  while (ll <= rr) {
    mm = (ll + rr) / 2;
    if (getuint(sa->s, mm, sa->k) <= i) ll = mm + 1; else rr = mm - 1;
  }
  bi = ll - 1;
  r = getuint(sa->s, bi, sa->k);
  i -= r;
  return sparsearray_sb_getbit(&sa->sa[bi],i);
}

i64 sparsearray4_rank_and_bit(sparsearray4 *sa, i64 i, int *c)
{
  i64 b,bi,r, ll,rr,mm;

  b = (sa->m+BLOCK-1)/BLOCK;
  ll = 0;  rr = b-1;
  while (ll <= rr) {
    mm = (ll + rr) / 2;
    if (getuint(sa->s, mm, sa->k) <= i) ll = mm + 1; else rr = mm - 1;
  }
  bi = ll - 1;
  r = getuint(sa->s, bi, sa->k);
  i -= r;
//  *c = sparsearray_sb_getbit(&sa->sa[bi],i);
//  return bi*BLOCK + sparsearray_sb_rank(&sa->sa[bi],i);
  return bi*BLOCK + sparsearray_sb_rank_bit(&sa->sa[bi],i,c);
}



#ifdef SPARSEMAIN
//#if 1
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
//#define N (1<<20)

int main(int argc, char *argv[])
{
  sparsearray4 s1,s2;
  i64 i,r,n,rr,m;
  u64 hoge,sum;
  i64 size;
  FILE *infp = NULL;
  FILE *out;
  MMAP *map;
  uchar *mapp;

  bitvec_t *B;
  dword *S,*R;

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

  mymalloc(B,(n+PBS-1)/PBS,0);

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

  sparsearray4_construct_init(&s2, n, m);
  r = 0;
  for (i=0; i<n; i++) {
    if (getbit(B,i)==1) {
//      printf("B[%ld] = 1\n",i);
      sparsearray4_construct_set(&s2, r, i);
      r++;
    }
  }
  size = sparsearray4_construct_end(&s2, 0, SDARRAY_SELECT1 | SDARRAY_RANK1);
  printf("used memory: %d bytes (%lf bpc)\n",size,(double)size*8/n);

#if 1
  out = fopen("sparsearraytmp.dat","w");
  sparsearray4_write(&s2, out);
  fclose(out);

  map = mymmap("sparsearraytmp.dat");
  if (map->addr==NULL) {
    perror("mmap2\n");
    exit(1);
  }
  mapp = (uchar *)map->addr;

  sparsearray4_read(&s1, &mapp);
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

#if 0
  for (i = 0; i < n; i++) {
    if (getbit(B,i) != sparsearray4_getbit(&s1,i)) {
     printf("B[%ld]=%d getbit=%d\n",i,getbit(B,i),sparsearray4_getbit(&s1,i));
    }
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
    if (sparsearray4_select(&s1,j) != S[j]) {
      printf("ERROR: S[%d] = %d, s = %d\n",j,S[j],sparsearray4_select(&s1,j));
    }
    sum += sparsearray4_select(&s1,j);
#else
    sum += sparsearray4_select(&s1,j);
#endif
    hoge = hoge * 1566083941U + 1;
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("%ld %f sparsearray4 select time sum=%ld\n",rr,t,sum);
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
    if (sparsearray4_rank(&s1,j) != R[j]) {
      printf("ERROR: (%d) R[%d] = %d, r = %d\n", getbit(B,i),j, R[j],
                                                 sparsearray4_rank(&s1,j));
    }
    sum += sparsearray4_rank(&s1,j);
#else
    sum += sparsearray4_rank(&s1,j);
#endif
    hoge = hoge * 1566083941U + 1;
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("%ld %f sparsearray4 rank time sum=%ld\n",rr,t,sum);
#endif

  return 0;
}

#endif //  SPARSEMAIN
