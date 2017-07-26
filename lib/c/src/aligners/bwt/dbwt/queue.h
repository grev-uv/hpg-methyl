#ifndef _DBWT_QUEUE_
#define _DBWT_QUEUE_

#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "queue.h"

#define QSIZ 1024

typedef struct qblock {
  struct qblock *prev, *next;
  packed_array *b;
} qblock;

typedef struct {
  long n; // 要素数
  int w; // 格納する値のビット数
  qblock *sb, *eb;
  long s_ofs, e_ofs; // 0 <= s_ofs, e_ofs < QSIZ
} queue;

queue *init_queue(int w);
void enqueue(queue *que, long x);
void enqueue_l(queue *que, long x);
long dequeue(queue *que);
int emptyqueue(queue *que);
void free_queue(queue *que);
void printqueue(queue *que);

#endif

/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
