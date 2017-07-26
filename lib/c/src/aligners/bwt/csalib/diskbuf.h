/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _DISKBUF_H_
#define _DISKBUF_H_

#include <unistd.h>
#include <stdio.h>

#define PAGESIZE (1<<16)

typedef struct {
  long n; // 要素数
  int dirty; // 現在のバッファ内が更新されたか
  int width; // 格納する整数の桁数
  int m; // ページ内の整数の数
  long ptr; // 現在のバッファのページ番号 (単位:PAGESIZE)
  int minidx, maxidx; // ページ内で更新したアドレスの最大と最小
  FILE *fp;
  unsigned char buf[PAGESIZE];
} diskbuf;

diskbuf *open_diskbuf(FILE *f, int width);
long getint_diskbuf(diskbuf *b, long idx);
int setint_diskbuf(diskbuf *b, long idx, long x);
int close_diskbuf(diskbuf *b);

void getstr_diskbuf(unsigned char *buf, diskbuf *b, long idx);
int setstr_diskbuf(unsigned char *buf, diskbuf *b, long idx);
FILE *create_tmp(int idx);
FILE *open_tmp(int idx);
void *remove_tmp(int idx);
//void *mymalloc(size_t n);

#endif
