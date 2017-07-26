/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/timeb.h>
#include "densearray2.h"

//#define CHECK
#define RANDOM

#define logD 6

#define PBS (sizeof(bitvec_t)*8)
#define D (1<<logD)

#define logL 8
#define L (1<<logL)
#define logR 7
#define R (1<<logR)


#ifndef min
 #define min(x,y) ((x)<(y)?(x):(y))
#endif


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

static unsigned int selecttbl[8*256];

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

static void make_selecttbl(void)
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

i64 densearray_sb_construct(densearray_sb *da, i64 n, bitvec_t *buf, int opt)
{
  i64 i,j,m;
  i64 nl;
  
  

  i64 size;
  i64 p1,m1;

  if (n > (1<<16)) {
    printf("densearray_sb_construct: n = %ld\n",n);
    exit(1);
  }

  size = sizeof(densearray_sb);
//  printf("densearray_sb-size:0 %ld\n",size);

  make_selecttbl();

  m = 0;
  for (i=0; i<n; i++) m += getbit(buf,i);

  da->b = buf;

// select index
  if (opt & SDARRAY_SELECT1) {
    nl = (m-1) / L + 1;
    da->s = (word *) mymalloc((nl+1) * sizeof(word));  size += sizeof(*(da->s))*(nl+1);
//    printf("densearray_sb-size:2 %ld\n",size);

    p1 = -1;
    m1 = -1;
    for (i = 0; i < nl; i++) {
      while (m1 < i*L && p1<n-1) m1 += getbit(buf,++p1);
      if (m1 < i*L) printf("???1 m1 = %ld i=%ld\n",m1,i);
      da->s[i] = p1;
    }
  } else {
    da->s = NULL;
  }

// rank index
  if (opt & SDARRAY_RANK1) {
    da->r = (word *) mymalloc((n/R+1) * sizeof(word));  size += sizeof(*(da->r))*(n/R+1);
//    printf("densearray_sb-size:3 %ld\n",size);
    m = 0;
    for (j=0; j<n; j+=R) {
      da->r[j/R] = m;
      for (i=0; i<R; i++) {
        if (i+j < n && getbit(buf,i+j)==1) m++;
      }
    }
    opt |= SDARRAY_RANK1;
  } else {
    da->r = NULL;
  }
//
  return size;
}

i64 densearray_sb_write(densearray_sb *da, i64 n, int opt, FILE *f)
{
  i64 i,j,m;
  i64 nl;
  
  
  
  i64 size;
  

  size = 0;
//  printf("densearray_sb-size:0 %ld\n",size);

  m = 0;
  for (i=0; i<n; i++) m += getbit(da->b,i);
  writeuint(sizeof(i64), m, f);  size += sizeof(i64);
  //printf("densearray_sb_write: n=%ld m=%ld\n",n,m);

// select index
  if (opt & SDARRAY_SELECT1) {
    nl = (m-1) / L + 1;
    for (i = 0; i < nl; i++) {
      writeuint(sizeof(da->s[0]), da->s[i], f);  size += sizeof(da->s[0]);
    }
    //printf("densearray_sb_write: writing select index (%ld)\n",nl);
  }

// rank index
  if (opt & SDARRAY_RANK1) {
    m = 0;
    for (j=0; j<n; j+=R) {
      writeuint(sizeof(da->r[0]), da->r[j/R], f);  size += sizeof(da->s[0]);
    }
  //  printf("densearray_sb_write: writing rank index (%ld)\n",n/R);
  }
//
  if (!(opt & SDARRAY_NOBUF)) {
    m = (n+sizeof(da->b)*8-1) / (sizeof(da->b)*8);
//    printf("write buf m=%ld\n",m);
    for (i=0; i<m; i++) {
      writeuint(sizeof(*da->b), da->b[i], f);  size += sizeof(*da->b);
    }
  }
  return size;
}

i64 densearray_sb_read(densearray_sb *da, i64 n, int opt, uchar **map)
{
  i64 m;
  i64 nl;
  
  
  uchar *p;

  make_selecttbl();

  p = *map;

  m = getuint(p,0,sizeof(i64));  p += sizeof(i64);
//  printf("densearray_sb_read: n=%ld m=%ld opt=%d\n",n,m,opt);

// select index
  if (opt & SDARRAY_SELECT1) {
    nl = (m-1) / L + 1;
    da->s = (word *)p;
    p += sizeof(*da->s) * nl;
//    printf("densearray_sb_read: reading select index (%ld)\n",nl);
  }

// rank index
  if (opt & SDARRAY_RANK1) {
    da->r = (word *)p;
    p += sizeof(*da->r) * ((n+R-1)/R);
  //  printf("densearray_sb_read: reading rank index (%ld)\n",n/R);
  }
//
  if (!(opt & SDARRAY_NOBUF)) {
    m = (n+sizeof(da->b)*8-1) / (sizeof(da->b)*8);
    //printf("read buf m=%ld\n",m);
    da->b = (bitvec_t *)p;
    p += sizeof(*da->b) * m;
  }
  *map = p;
  return 0;
}


i64 densearray_sb_rank(densearray_sb *da, i64 i)
{
  i64 r,j;
  bitvec_t *p;
  r = da->r[i>>logR];
  p = da->b + ((i>>logR)<<(logR-logD));
  j = i & (R-1);
#if R > D
  while (j >= D) {
    r += POPCOUNT(*p++);
    j -= D;
  }
#endif
  r += POPCOUNT(*p >> (D-1-j));
  return r;
}

i64 densearray_sb_rank0(densearray_sb *da, i64 i)
{
  return i+1 - densearray_sb_rank(da,i);
}

i64 densearray_sb_rank_bit(densearray_sb *da, i64 i, int *c)
{
  i64 r,j;
  bitvec_t *p,x;
  r = da->r[i>>logR];
  p = da->b + ((i>>logR)<<(logR-logD));
  j = i & (R-1);
#if R > D
  while (j >= D) {
    r += POPCOUNT(*p++);
    j -= D;
  }
#endif
  x = *p >> (D-1-j);
  r += POPCOUNT(x);
  *c = (int)(x & 1);
  return r;
}


int densearray_sb_getbit(densearray_sb *da, i64 i)
{
  return getbit(da->b, i);
}


i64 densearray_sb_select(densearray_sb *da, i64 i,int f)
{
  
  //  dword *s;
  i64 p,r;
  
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

  p = da->s[i>>logL];
  r = i - (i & (L-1));

  q = &(da->b[p>>logD]);

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
      rr = POPCOUNT(x >> (D-8));
      if (r + rr >= i) break;
      r += rr;
      p += 8;
      x <<= 8;
    }
    p += selecttbl[((i-r-1)<<8)+(x>>(D-8))];
  } else {
    rr = p & (D-1);
    r -= POPCOUNT((~*q) >> (D-1-rr));
    p = p - rr;
    
    while (1) {
      rr = POPCOUNT(~*q);
      if (r + rr >= i) break;
      r += rr;
      p += D;
      q++;
    }
    
    x = ~*q;
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
//#define N (1<<14)

int main(int argc, char *argv[])
{
  densearray_sb s1;

  i64 i,r,n,m,rr,size;
  u64 hoge,sum;
  FILE *infp = NULL;

  bitvec_t *B;
  dword *sel,*rank;

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

  size = densearray_sb_construct(&s1, n,B, SDARRAY_RANK1 | SDARRAY_SELECT1);
  printf("da: used memory: %d bytes (%lf bpc)\n",size,(double)size*8/n);

#ifdef CHECK
  mymalloc(sel,n+1,0);
  mymalloc(rank,n+1,0);
  r = 0;
  sel[r] = -1;
  for (i=0; i<n; i++) {
    if (getbit(B,i)) {
      r++;
      sel[r] = i;
    }
    rank[i] = r;
  }
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
    if (densearray_sb_rank(&s1,j) != rank[j]) {
      printf("ERROR: (%d) R[%d] = %d, r = %d\n", getbit(B,i),j, rank[j],
                                                 densearray_sb_rank(&s1,j));
    }
    sum += densearray_sb_rank(&s1,j);
#else
    sum += densearray_sb_rank(&s1,j);
#endif
    hoge = hoge * 1566083941U + 1;
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("%ld %f densearray_sb rank time sum=%ld\n",rr,t,sum);
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
    if (densearray_sb_select(&s1,j,1) != sel[j]) {
      printf("ERROR: S[%d] = %d, s = %d\n",j,sel[j],densearray_sb_select(&s1,j,1));
    }
    sum += densearray_sb_select(&s1,j,1);
#else
    sum += densearray_sb_select(&s1,j,1);
#endif
    hoge = hoge * 1566083941U + 1;
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("%ld %f densearray_sb select time sum=%ld\n",rr,t,sum);
#endif

  return 0;
}

#endif //  DENSEMAIN
