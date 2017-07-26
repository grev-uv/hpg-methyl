/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>
#include <math.h>
#include "sparsearray.h"

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
#define logLLL 5
//#define LLL 128
//#define LLL 32
#define LLL (1<<logLLL)
//#define logL 10
//#define logL (logLL-3)
#define logL (logLL-1-5)
#define L (1<<logL)


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

static void writeuint(int k,u64 x,FILE *f)
{
  int i;
  for (i=k-1; i>=0; i--) {
    fputc(x & 0xff,f);
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
  l = i & (D-1);
  return (B[j] >> (D-1-l)) & 1;
}

static bitvec_t getbits(bitvec_t *B, i64 i, i64 d)
{
  bitvec_t x,z;

  B += (i >> logD);
  i &= (D-1);
  if (i+d <= D) {
    x = B[0];
    x <<= i;
    x >>= (D-1-d);
    x >>= 1;
  } else {
    x = B[0] << i;
    x >>= D-d;
    z = B[1] >> (D-(i+d-D));
    x += z;
  }
  
  return x;
}

int sparsearray_construct(sparsearray *sa, i64 n, bitvec_t *buf, int opt)
{
  i64 i,m,r;

  m = 0;
  for (i=0; i<n; i++) {
    if (getbit(buf,i)==1) m++;
  }

  sparsearray_construct_init(sa, n, m);
  r = 0;
  for (i=0; i<n; i++) {
    if (getbit(buf,i)==1) {
      sparsearray_construct_set(sa, r, i);
      r++;
    }
  }
  sparsearray_construct_end(sa, opt);
  return 0;
}

int sparsearray_construct_init(sparsearray *sa, i64 n, i64 m)
{
  i64 i;
  
  
  
  i64 d,mm;
#if LOWCHAR
  byte *low;
#else
  bitvec_t *low;
#endif
  bitvec_t *buf2;
  
  
  
  i64 size;

//  sa->opt = opt;

  size = sizeof(sparsearray);

  sa->n = n;
  sa->m = m;

#if LOWCHAR
  d = 8;
#else
  mm = m;
  d = 0;
  while (mm < n) {
    mm <<= 1;
    d++;
  }
#endif
//  printf("n=%ld m=%ld d=%ld\n",n,m,d);
  sa->d = d;

  buf2 = (bitvec_t*) mymalloc(((2*m+PBS-1)/PBS+1) * sizeof(bitvec_t));  size += sizeof(*buf2) * ((2*m+PBS-1)/PBS+1);
#if LOWCHAR
  low = (byte *) mymalloc(m * sizeof(byte));  size += sizeof(*low) * m;
#else
  low = (bitvec_t *) mymalloc(((d*m+PBS-1)/PBS+1) * sizeof(bitvec_t));  size += sizeof(*low) * ((d*m+PBS-1)/PBS+1);
#endif
//  printf("sparsearray-size:0 %ld\n",size);
  sa->hi = buf2;
  sa->low = low;

//  for (i=0; i<m*2; i++) setbit(buf2,i,0);
//  for (i=0; i<m; i++) setbits(low,i*d,d,0);
  for (i=0; i<(2*m+PBS-1)/PBS+1; i++) buf2[i] = 0;
#if LOWCHAR
  for (i=0; i<m; i++) low[i] = 0;
#else
  for (i=0; i<(d*m+PBS-1)/PBS+1; i++) low[i] = 0;
#endif

  sa->size = size;
  return 0;
}

int sparsearray_construct_set(sparsearray *sa, i64 i, i64 x)
{
  i64 d;
  d = sa->d;

//  printf("setbit(%ld)\n",(x>>d)+i);
  setbit(sa->hi,(x>>d)+i,1);
#if LOWCHAR
  sa->low[i] = x & 0xff;
#else
  setbits(sa->low,i*d,d,x & ((1<<d)-1));
#endif
  return 0;
}

int sparsearray_construct_end(sparsearray *sa, int opt)
{
  bitvec_t *buf2;
  i64 i,m;

  sa->opt = opt;

  buf2 = sa->hi;
  m = sa->m;

  if (sa->opt & SDARRAY_SELECT1) {
    sa->sd1 = (densearray *) mymalloc(1 * sizeof(densearray));
    densearray_construct(sa->sd1,m*2,buf2,SDARRAY_SELECT1);  
    sa->size += sa->sd1->size;
  } else {
    sa->sd1 = NULL;
  }

  if (sa->opt & SDARRAY_RANK1) {
    for (i=0; i<m*2; i++) setbit(buf2,i,1-getbit(buf2,i));
    sa->sd0 = (densearray *) mymalloc(1 * sizeof(densearray));
    densearray_construct(sa->sd0,m*2,buf2,SDARRAY_SELECT1 | SDARRAY_NOBUF);
    sa->size += sa->sd0->size;

    for (i=0; i<m*2; i++) setbit(buf2,i,1-getbit(buf2,i));
  } else {
    sa->sd0 = NULL;
  }
//  printf("sparsearray-size:1 %ld\n",sa->size);
  return 0;
}

i64 sparsearray_write(sparsearray *sa, FILE *f)
{
  i64 n,m,d,k;
  i64 size;
  i64 i;

  size = 0;

//  k = sa->k;
  n = sa->n;
  m = sa->m;
  d = sa->d;

  k = (blog(n+1)+1+8-1)/8;

  writeuint(1,k,f);  size += 1;
  writeuint(k, sa->n, f);  size += k;
  writeuint(k, sa->m, f);  size += k;
  writeuint(k, sa->d, f);  size += k;
  writeuint(1, sa->opt, f);  size += 1;

  for (i=0; i<(d*m+PBS-1)/PBS+1; i++) {
    writeuint(sizeof(sa->low[0]), sa->low[i], f);  size += sizeof(sa->low[0]);
  }

  size += densearray_write(sa->sd1, f);

  if (sa->opt & SDARRAY_RANK1) {
    size += densearray_write(sa->sd0, f);
  }

  return size;
}

void sparsearray_read(sparsearray *sa, uchar **map)
{
  uchar *p;
  i64 k;

  p = *map;

  k = getuint(p,0,1);  p += 1;
  sa->n = getuint(p,0,k);  p += k;
  sa->m = getuint(p,0,k);  p += k;
  sa->d = getuint(p,0,k);  p += k;
  sa->opt = getuint(p,0,1);  p += 1;
//  printf("sparsearray_read k=%ld n=%ld m=%ld d=%ld\n",k,sa->n, sa->m, sa->d);
#if LOWCHAR 
  sa->low = (byte *) p;
#else
  sa->low = (bitvec_t *) p;
#endif

  p += sizeof(*sa->low) * ((sa->d*sa->m+PBS-1)/PBS+1);

  *map = p;
  
  sa->sd1 = (densearray *) mymalloc(1 * sizeof(densearray));
  densearray_read(sa->sd1, map);
  sa->hi = sa->sd1->buf;

  if (sa->opt & SDARRAY_RANK1) {
    sa->sd0 = (densearray *) mymalloc(1 * sizeof(densearray));
    densearray_read(sa->sd0, map);
    sa->sd0->buf = sa->sd1->buf;
  }

}


i64 sparsearray_select(sparsearray *sa, i64 i)
{
  i64 d,x;

#if 0
  if ((sa->opt & SDARRAY_SELECT1) == 0) {
    printf("sparsearray_select: select1 is not supported.\n");
  }
  if (i > sa->m) {
    printf("ERROR: m=%d i=%d\n",sa->m,i);
    exit(1);
  }
#endif

  if (i == 0) return -1;

  d = sa->d;

  x = densearray_select(sa->sd1,i,1) - (i-1);
  x <<= d;
#if LOWCHAR
  x += sa->low[i-1];
#else
  x += getbits(sa->low,(i-1)*d,d);
#endif
  return x;

}

i64 sparsearray_rank(sparsearray *sa, i64 i)
{
  i64 d,x,w,y;
  i64 j;

#if 0
  if ((sa->opt & SDARRAY_RANK1) == 0) {
    printf("sparsearray_select: rank1 is not supported.\n");
  }
  if (i > sa->n) {
    printf("ERROR: n=%d i=%d\n",sa->n,i);
    exit(1);
  }
#endif

  d = sa->d;

  y = densearray_select(sa->sd0,i>>d,0)+1;
  x = y - (i>>d);

  j = i - ((i>>d)<<d);

  while (1) {
    if (getbit(sa->hi,y)==0) break;
#if LOWCHAR
    w = sa->low[x];
#else
    w = getbits(sa->low,x*d,d);
#endif
    if (w >= j) {
      if (w == j) x++;
      break;
    }
    x++;
    y++;
  }
  return x;
}

#ifdef SPARSEMAIN
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
  sparsearray s1,s2;
  i64 i,r,n,rr,m;
  u64 hoge,sum;
  FILE *infp = NULL;
  FILE *out;

  bitvec_t *B;
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

//  sparsearray_construct(&s1,n,B, SDARRAY_SELECT1 | SDARRAY_RANK1);
  sparsearray_construct_init(&s2, n, m);
  r = 0;
  for (i=0; i<n; i++) {
    if (getbit(B,i)==1) {
      sparsearray_construct_set(&s2, r, i);
      r++;
    }
  }
  sparsearray_construct_end(&s2, SDARRAY_SELECT1 | SDARRAY_RANK1);
//  sparsearray_construct_end(&s2, SDARRAY_SELECT1);
  printf("used memory: %d bytes (%lf bpc)\n",s2.size,(double)s2.size*8/n);

#if 1
  out = fopen("sparsearraytmp.dat","w");
  sparsearray_write(&s2, out);
  fclose(out);

  map = mymmap("sparsearraytmp.dat");
  if (map->addr==NULL) {
    perror("mmap2\n");
    exit(1);
  }
  mapp = (uchar *)map->addr;

  sparsearray_read(&s1, &mapp);
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
    if (sparsearray_select(&s1,j) != S[j]) {
      printf("ERROR: S[%d] = %d, s = %d\n",j,S[j],sparsearray_select(&s1,j));
    }
    sum += sparsearray_select(&s1,j);
#else
    sum += sparsearray_select(&s1,j);
#endif
    hoge = hoge * 1566083941U + 1;
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("%ld %f sparsearray select time sum=%ld\n",rr,t,sum);
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
    if (sparsearray_rank(&s1,j) != R[j]) {
      printf("ERROR: (%d) R[%d] = %d, r = %d\n", getbit(B,i),j, R[j],
                                                 sparsearray_rank(&s1,j));
    }
    sum += sparsearray_rank(&s1,j);
#else
    sum += sparsearray_rank(&s1,j);
#endif
    hoge = hoge * 1566083941U + 1;
  }
  mygettime(&after);
  t = mylaptime(&before,&after);
  printf("%ld %f sparsearray rank time sum=%ld\n",rr,t,sum);
#endif

  return 0;
}

#endif //  SPARSEMAIN
