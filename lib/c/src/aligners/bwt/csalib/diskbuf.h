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
  long n; // �v�f��
  int dirty; // ���݂̃o�b�t�@�����X�V���ꂽ��
  int width; // �i�[���鐮���̌���
  int m; // �y�[�W���̐����̐�
  long ptr; // ���݂̃o�b�t�@�̃y�[�W�ԍ� (�P��:PAGESIZE)
  int minidx, maxidx; // �y�[�W���ōX�V�����A�h���X�̍ő�ƍŏ�
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
