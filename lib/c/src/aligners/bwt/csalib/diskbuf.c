/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include "diskbuf.h"

static void *mymalloc(size_t n)
{
  void *p;

//  printf("allocating %ld bytes (%ld bytes available)\n",n,available_memory());

  p = malloc(n);
  if (p == NULL) {
    printf("malloc failed.\n");
    exit(1);
  }
  if (n == 0) {
    printf("warning: 0 bytes allocated p=%p\n",p);
  }

  return p;
}


FILE *create_tmp(int idx)
{
  FILE *f;
  char fname[20];
  sprintf(fname,"tmpfile.%08x",idx);
//  printf("create %s\n",fname);
  f = fopen(fname,"w+");
  if (f == NULL) {
    perror("create_tmp: creating failed.\n");
    exit(1);
  }

  return f;
}

FILE *open_tmp(int idx)
{
  FILE *f;
  char fname[20];
  sprintf(fname,"tmpfile.%08x",idx);
//  printf("open %s\n",fname);
  f = fopen(fname,"r+");
  if (f == NULL) {
    perror("open_tmp: open failed.\n");
    printf("filename %s\n",fname);
    exit(1);
  }

  return f;
}

void *remove_tmp(int idx)
{

  int r;
  char fname[20];
  sprintf(fname,"tmpfile.%08x",idx);
//  printf("remove %s\n",fname);
  r = unlink(fname);
  if (r == -1) {
    perror("remove_tmp");
    printf("remove %s\n",fname);
//    exit(1);
  }
  return NULL;
}



diskbuf *open_diskbuf(FILE *f, int width)
{
  diskbuf *b;

//  width = (width+7)/8;

  b = (diskbuf *)mymalloc(sizeof(diskbuf));
  b->fp = f;
  b->width = width;
  b->m = PAGESIZE / width;
  b->ptr = 0;
  fseek(b->fp,b->ptr*PAGESIZE,SEEK_SET);
  int res = fread(b->buf,1,PAGESIZE,b->fp);
  if (!res) {  }
  b->dirty = 0;
  b->minidx = PAGESIZE;  b->maxidx = -1;

  return b;
}

long getint_diskbuf(diskbuf *b, long idx)
{
  long i,m,w;
  long x;
  long q,r;
  m = b->m;
  w = b->width;
  q = idx / m;  r = idx % m;
  if (q != b->ptr) {
    if (b->dirty) {
      fseek(b->fp,b->ptr*PAGESIZE+b->minidx*w,SEEK_SET);
//      fwrite(b->buf,1,PAGESIZE,b->fp);
      fwrite(b->buf+b->minidx*w,1,(b->maxidx - b->minidx+1)*w,b->fp);
    }
    b->ptr = q;
    fseek(b->fp,b->ptr*PAGESIZE,SEEK_SET);
    int res = fread(b->buf,1,PAGESIZE,b->fp);
    if (!res) {  }
    b->dirty = 0;
    b->minidx = PAGESIZE;  b->maxidx = -1;
  }
  x = 0;
  for (i=0; i<w; i++) {
    x <<= 8;
    x += b->buf[r*w+i];
  }
  return x;
}

void getstr_diskbuf(unsigned char *buf, diskbuf *b, long idx)
{
  long i,m,w;
  long q,r;
  m = b->m;
  w = b->width;
  q = idx / m;  r = idx % m;
  if (q != b->ptr) {
    if (b->dirty) {
      fseek(b->fp,b->ptr*PAGESIZE+b->minidx*w,SEEK_SET);
//      fwrite(b->buf,1,PAGESIZE,b->fp);
      fwrite(b->buf+b->minidx*w,1,(b->maxidx - b->minidx+1)*w,b->fp);
    }
    b->ptr = q;
    fseek(b->fp,b->ptr*PAGESIZE,SEEK_SET);
    int res = fread(b->buf,1,PAGESIZE,b->fp);
    if (!res) {  }
    b->dirty = 0;
    b->minidx = PAGESIZE;  b->maxidx = -1;
  }
  for (i=0; i<w; i++) {
    buf[i] = b->buf[r*w+i];
  }
}


int setint_diskbuf(diskbuf *b, long idx, long x)
{
  long i,w,m;
  long q,r;
  m = b->m;
  w = b->width;
  q = idx / m;  r = idx % m;
  if (q != b->ptr) {
    if (b->dirty) {
      fseek(b->fp,b->ptr*PAGESIZE+b->minidx*w,SEEK_SET);
//      fwrite(b->buf,1,PAGESIZE,b->fp);
      fwrite(b->buf+b->minidx*w,1,(b->maxidx - b->minidx+1)*w,b->fp);
    }
    b->ptr = q;
    fseek(b->fp,b->ptr*PAGESIZE,SEEK_SET);
    int res = fread(b->buf,1,PAGESIZE,b->fp);
    if (!res) {  }
    b->dirty = 0;
    b->minidx = PAGESIZE;  b->maxidx = -1;
  }
  for (i=w-1; i>=0; i--) {
    b->buf[r*w+i] = x & 0xff;
    x >>= 8;
  }
  b->dirty = 1;
  if (r > b->maxidx) b->maxidx = r;
  if (r < b->minidx) b->minidx = r;
  return x;
}

int setstr_diskbuf(unsigned char *buf, diskbuf *b, long idx)
{
  long i,w,m;
  long q,r;
  m = b->m;
  w = b->width;
  q = idx / m;  r = idx % m;
  if (q != b->ptr) {
    if (b->dirty) {
      fseek(b->fp,b->ptr*PAGESIZE+b->minidx*w,SEEK_SET);
//      fwrite(b->buf,1,PAGESIZE,b->fp);
      fwrite(b->buf+b->minidx*w,1,(b->maxidx - b->minidx+1)*w,b->fp);
    }
    b->ptr = q;
    fseek(b->fp,b->ptr*PAGESIZE,SEEK_SET);
    int res = fread(b->buf,1,PAGESIZE,b->fp);
    if (!res) {  }
    b->dirty = 0;
  }
  for (i=0; i<w; i++) {
    b->buf[r*w+i] = buf[i];
  }
  b->dirty = 1;
  if (r > b->maxidx) b->maxidx = r;
  if (r < b->minidx) b->minidx = r;
  return 0;
}

int close_diskbuf(diskbuf *b)
{
  long w;
  w = b->width;
  if (b->dirty) {
    fseek(b->fp,b->ptr*PAGESIZE+b->minidx*w,SEEK_SET);
//    fwrite(b->buf,1,PAGESIZE,b->fp);
    fwrite(b->buf+b->minidx*w,1,(b->maxidx - b->minidx+1)*w,b->fp);
  }
  //fclose(b->fp);
  free(b);
  return 0;
}
