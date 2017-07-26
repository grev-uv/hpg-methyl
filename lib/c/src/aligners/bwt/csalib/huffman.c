/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include "huffman.h"

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

static Huffman *newHuffman(int n)
{
  Huffman *p;
  p = (Huffman *) mymalloc(1 * sizeof(Huffman));
  p->n = n;
  p->v = (int *) mymalloc(n * sizeof(int));
  p->left = (int *) mymalloc((2*n+3) * sizeof(int));
  p->right = (int *) mymalloc((2*n+3) * sizeof(int));
  p->clen = (byte *) mymalloc((n*2) * sizeof(byte));
  p->code = (u64 *) mymalloc((n*2) * sizeof(u64));
  return p;
}

void freeHuffman(Huffman *p)
{
  free(p->v);
  free(p->left);
  free(p->right);
  free(p->clen);
  free(p->code);
  free(p);
}

static void maketree(int n, int r, int d, Huffman *h, u64 c)
{
  h->clen[r] = d;
  if (h->left[r] >= 0) {
    maketree(n,h->left[r],d+1,h, c);
  } else {
    //h->v[r] = r;
    h->code[r] = c;
  }
  if (h->right[r] >= 0) {
    maketree(n,h->right[r],d+1,h, c + (1L<<(sizeof(c)*8-1-d)));
  } else {
    //h->v[r] = r;
    h->code[r] = c;
  }
}

void mkhufdecodetable(Huffman *h)
{
  u64 i;
  u64 v,w;

  for (i = 0; i < HUFTBLSIZ; i++) {
    w = i << (sizeof(u64)*8-HUFTBLWIDTH); // bit stream
    v = DecodeHuffman(h, w);
    if (h->clen[v] <= HUFTBLWIDTH) {
      h->tbl[i] = v;
    } else {
      h->tbl[i] = -1;
    }
  }
}

int DecodeHuffman(Huffman *h, u64 x)
{
  unsigned int i = sizeof(x)*8 - 1;
  int p = h->n * 2 - 2;
  while (1) {
    if ((x >> i) & 1) p = h->right[p]; else p = h->left[p];
    i--;
    //if (p < h->n) return h->v[p];
    if (p < h->n) return p;
  }
}

int DecodeHuffman_tbl(Huffman *h, u64 x)
{
  u64 w;
  
  w = x >> (sizeof(u64)*8-HUFTBLWIDTH);
  if (h->tbl[w] >= 0) return h->tbl[w];
  return DecodeHuffman(h,x);
}

Huffman *MakeHuffmanTree(int n, double *freq)
{
  int i,j;
  Huffman *h;
  double *tmpf;
  int l,r;
  int m1,m2;

  h = newHuffman(n);
  tmpf = (double *) mymalloc((2*n+3) * sizeof(double));

  for (i=0; i<n; i++) tmpf[i] = freq[i];
  for (i=0; i<n*2+3; i++) h->left[i] = h->right[i] = -1;
  for (i=n; i<n*2-1+2; i++) tmpf[i] = 0.0;
//  tmpf[n*2+1] = tmpf[n*2+2] = 10000.0;
  tmpf[n*2-1] = tmpf[n*2] = 1000000.0;

  l = 0; r = n-1;
  for (j=0; j<n-1; j++) {
//    m1 = n*2+1;  m2 = n*2+2;
    m1 = n*2-1;  m2 = n*2;
    for (i=l; i<=r; i++) {
      if ((tmpf[i] > 0) && (tmpf[i] < tmpf[m2])) m2 = i;
      if ((tmpf[i] > 0) && (tmpf[i] < tmpf[m1])) {m2 = m1; m1 = i;}
    }
    h->left[r+1] = m1;  h->right[r+1] = m2;
    tmpf[r+1] = tmpf[m1] + tmpf[m2];
    tmpf[m1] = tmpf[m2] = 0.0;
    r++;
  }
  maketree(n, r, 0, h, 0);

  free(tmpf);

  mkhufdecodetable(h);
  return h;
}

void Huffman_write(Huffman *h, FILE *out)
{
  int k,n,i,len,d;
  u64 x;
//  printf("Huffman_write\n");
  writeuint(1,ID_HUFFMAN,out);
  n = h->n;
  k = (blog(2*n+4)+1+8-1)/8;
  writeuint(1,k,out);
  writeuint(k,n,out);

  for (i=n; i<=2*n-2; i++) {
//    printf("i=%d left=%d right=%d\n",i,h->left[i],h->right[i]);
    writeuint(k,h->left[i],out);
    writeuint(k,h->right[i],out);
  }
  for (i=0; i<n; i++) {
    len = h->clen[i];
    writeuint(1,len,out);
    d = (len+7)/8;
    x = h->code[i];
    x >>= (sizeof(x)-d)*8;
//    printf("i=%d len=%d d=%d code=%lx\n",i,len,d,x);
    writeuint(d,x,out);
  }
}

#if 0
static u64 readuint(int k,FILE *f)
{
  u64 x;
  int i;
  x = 0;
   for (i=0; i<k; i++) {
    x += ((u64)fgetc(f)<<(8*i));
  }
 return x;
}

Huffman *Huffman_read(FILE *in)
{
  int id;
  int k,n,i,len,d;
  u64 x;
  Huffman *h;
//  printf("Huffman_read\n");
  if ((id = readuint(1,in)) != ID_HUFFMAN) {
//    printf("Huffman_read: id = %d\n",id);
//    exit(1);
  }
  k = readuint(1,in);
//  printf("Huffman_read k=%d\n",k);
  n = readuint(k,in);
//  printf("Huffman_read n=%d\n",n);

  h = newHuffman(n);

  for (i=n; i<=2*n-2; i++) {
    h->left[i] = readuint(k,in);
    h->right[i] = readuint(k,in);
//    printf("i=%d left=%d right=%d\n",i,h->left[i],h->right[i]);
  }
  for (i=0; i<n; i++) {
    len = h->clen[i] = readuint(1,in);
    d = (len+7)/8;
    x = readuint(d,in);
//    printf("i=%d len=%d d=%d code=%lx\n",i,len,d,x);
    x <<= (sizeof(x)-d)*8;
    h->code[i] = x;
  }
  mkhufdecodetable(h);
  return h;
}
#endif

Huffman *Huffman_read2(uchar **map)
{
  int id;
  int k,n,i,len,d;
  u64 x;
  Huffman *h;
  uchar *p;
  
  p = *map;

  if ((id = getuint(p,0,1)) != ID_HUFFMAN) {
//    printf("Huffman_read: id = %d\n",id);
//    exit(1);
  }
  p += 1;
  k = getuint(p,0,1);  p += 1;
//  printf("Huffman_read k=%d\n",k);
  n = getuint(p,0,k);  p += k;
//  printf("Huffman_read n=%d\n",n);

  h = newHuffman(n);

  for (i=n; i<=2*n-2; i++) {
    h->left[i] = getuint(p,0,k);  p += k;
    h->right[i] = getuint(p,0,k);  p += k;
//    printf("i=%d left=%d right=%d\n",i,h->left[i],h->right[i]);
  }
  for (i=0; i<n; i++) {
    len = h->clen[i] = getuint(p,0,1); p += 1;
    d = (len+7)/8;
    x = getuint(p,0,d);  p += d;
//    printf("i=%d len=%d d=%d code=%lx\n",i,len,d,x);
    x <<= (sizeof(x)-d)*8;
    h->code[i] = x;
  }
  mkhufdecodetable(h);
  *map = p;
  return h;
}


