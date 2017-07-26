#include "queue.h"

queue *init_queue(int w)
{
  queue *que;

  que = (queue *)mymalloc(sizeof(queue));
  que->n = 0;
  que->w = w;
  que->sb = que->eb = NULL;
  que->s_ofs = 0;
  que->e_ofs = QSIZ-1;
  return que;
}

void enqueue(queue *que, long x)
{
  qblock *qb;
  if (que->e_ofs == QSIZ-1) { // 現在のブロックが一杯
    qb = (qblock *)mymalloc(sizeof(qblock));
    qb->b = allocate_packed_array(QSIZ,que->w);
    if (que->eb == NULL) { // ブロックが1つもない
      que->sb = que->eb = qb;
      que->s_ofs = 0;
      qb->prev = qb->next = NULL;
    } else {
      que->eb->next = qb; // 現在の最後のブロックの次に追加
      qb->prev = que->eb;
      qb->next = NULL;
      que->eb = qb;
    }
    que->e_ofs = -1;
  }
  pa_set(que->eb->b, ++que->e_ofs, x);
  que->n++;
}

void enqueue_l(queue *que, long x) // 先頭に追加
{
  qblock *qb;
  if (que->s_ofs == 0) { // 現在のブロックが一杯
    qb = (qblock *)mymalloc(sizeof(qblock));
    qb->b = allocate_packed_array(QSIZ,que->w);
    if (que->sb == NULL) { // ブロックが1つもない
      que->sb = que->eb = qb;
      que->e_ofs = QSIZ-1;
      qb->prev = qb->next = NULL;
    } else {
      que->sb->prev = qb; // 現在の最初のブロックの前に追加
      qb->next = que->sb;
      qb->prev = NULL;
      que->sb = qb;
    }
    que->s_ofs = QSIZ;
  }
  pa_set(que->sb->b, --que->s_ofs, x);
  que->n++;
}

long dequeue(queue *que)
{
  long x;
  qblock *qb;
  x = pa_get(que->sb->b, que->s_ofs++);
  if (que->s_ofs == QSIZ) { // 現在のブロックが空に
    qb = que->sb;           // 現在の先頭ブロック
    que->sb = qb->next;
    free_packed_array(qb->b);
    myfree(qb,sizeof(qblock));
    if (que->sb == NULL) { // ブロックが無くなった
      que->eb = NULL;
      que->e_ofs = QSIZ-1;
      que->s_ofs = 0;
    } else {
      que->sb->prev = NULL;
      que->s_ofs = 0;
    }
  }
  que->n--;
  return x;
}

long dequeue_r(queue *que) // 最後から削除
{
  long x;
  qblock *qb;
  x = pa_get(que->eb->b, que->e_ofs--);
  if (que->e_ofs < 0) { // 現在のブロックが空に
    qb = que->eb;           // 現在の最終ブロック
    que->eb = qb->prev;
    free_packed_array(qb->b);
    myfree(qb,sizeof(qblock));
    if (que->eb == NULL) { // ブロックが無くなった
      que->sb = NULL;
      que->e_ofs = QSIZ-1;
      que->s_ofs = 0;
    } else {
      que->eb->next = NULL;
      que->e_ofs = QSIZ-1;
    }
  }
  que->n--;
  return x;
}

int emptyqueue(queue *que)
{
  return (que->n == 0);
}

void free_queue(queue *que)
{
  qblock *qb,*q;
  qb = que->sb;
  while (qb != NULL) {
    q = qb->next;
    free_packed_array(qb->b);
    myfree(qb,sizeof(qblock));
    qb = q;
  }
  myfree(que,sizeof(que));
}

void printqueue(queue *que)
{
  long i,s,t;
  qblock *qb;
  printf("printqueue que=%p\n",que);
  if (que == NULL) {
    printf("printqueue: que == NULL\n");
    return;
  }
  qb = que->sb;
  if (qb == NULL) {
    printf("printqueue: empty\n");
    return;
  }
  while (qb != NULL) {
    printf("qblock %p prev=%p next=%p (",qb,qb->prev,qb->next);
    s = 0;  t = QSIZ-1;
    if (qb == que->sb) s = que->s_ofs;
    if (qb == que->eb) t = que->e_ofs;
    for (i=s; i<=t; i++) {
      printf("%ld ",pa_get(qb->b,i));
    }
    printf(")\n");
    qb = qb->next;
  }
  printf("end\n");
}

#if 0
int main(void)
{
  int i,x;
  queue *q;
  q = init_queue(10);
  for (i=0; i<100; i++) {
    x = rand() % 1000;
    printf("%d: %d\n",i,x);
    enqueue(q,x);
  }
  for (i=0; i<100; i++) {
    printf("%d: %d\n",i,dequeue(q));
  }
  getchar();
  for (i=0; i<100; i++) {
    x = rand() % 1000;
    printf("%d: %d\n",i,x);
    enqueue(q,x);
  }
  for (i=0; i<100; i++) {
    printf("%d: %d\n",i,dequeue(q));
  }
  getchar();
}
#endif

#if 0
test()
{
  queue *Q[2][SIGMA+2];

#if 00
  bw2 = (uchar *)mymalloc(n+1);
  for (i=0; i<=SIGMA; i++) {
    Q[TYPE_L][i] = Q[TYPE_S][i] = NULL;
    //printf("i=%d N=%d C=%d\n",i,N[TYPE_L][i],C[i]);
    if (N[TYPE_L][i] > 0) Q[TYPE_L][i] = init_queue(N[TYPE_L][i]);
    if (C[i] > 0) Q[TYPE_S][i] = init_queue(C[i]);
  }
  for (i=0; i<n1; i++) {
    int l;
    uchar *q;
    p = sa[i];
    q = S[T2[p-1]];
    l = getlen(&q);
    q += 1; // for sentinel
    q += l-1;
    c = q[0];
    if (i == 0) {
      c = -1;
    }
    enqueue(Q[TYPE_S][c+1],q);
  }
#endif

#if 00
  printf("induced-sorting 1a\n");
  cc = 0;
  for (i=0; i<=SIGMA; i++) {
    M2[i] = cc; // M2[i] is the first index for i
    cc += M[i];
  }

  for (i = 0; i <= SIGMA; i++) {
    uchar *q;
    int c1,c2;
    C3[i] = 0;
    while (!emptyqueue(Q[TYPE_L][i])) {
      q = quehead(Q[TYPE_L][i]);
      if (q == lastptr) {
        last = i;
	dequeue(Q[TYPE_L][i]);
        continue;
      }
      c1 = q[-1];  c2 = q[0];
      printf("%d: %d(%c) %d(%c)\n",i,c1,c1,c2,c2);
      if (c1 >= c2) { // TYPE_L
        printf("copy %d to %d\n",i,M2[c1+1]);
	q = dequeue(Q[TYPE_L][i]);
        enqueue(Q[TYPE_L][c1+1],q-1);
        bw2[M2[c2+1]++] = c1;
      }
    }
    while (!emptyqueue(Q[TYPE_S][i])) {
      q = quehead(Q[TYPE_S][i]);
      if (q == lastptr) {
	last = i;
	dequeue(Q[TYPE_S][i]);
	continue;
      }
      c1 = q[-1];  c2 = q[0];
      if (i == 0) c2 = -1;
      printf("%d: %d(%c) %d(%c)\n",i,c1,c1,c2,c2);
      if (c1 >= c2) { // TYPE_L
        printf("copy %d to %d\n",i,M2[c1+1]);
	q = dequeue(Q[TYPE_S][i]);
        enqueue(Q[TYPE_L][c1+1],q-1);
        LMS_BW[c2+1][C3[c2+1]++] = c1;
      }
    }
  }
  for (i=0; i<=SIGMA; i++) {
    free_queue(Q[TYPE_S][i]);
  }
#endif

#if 00
  printf("induced-sorting 2a\n");
  for (i=1; i<=SIGMA; i++) M2[i] = 0;
  for (i=1; i<=SIGMA; i++) {
    M2[i] = M2[i-1] + M[i]; // M2[i] is the last index for i
  }
  for (i=0; i<=SIGMA; i++) {
    Q[TYPE_S][i] = NULL;
    if (N[TYPE_S][i] > 0) Q[TYPE_S][i] = init_queue(N[TYPE_S][i]);
  }
  for (i = SIGMA; i>0; i--) {
    uchar *q;
    int c1,c2;

    while (!emptyqueue(Q[TYPE_S][i])) {
      q = quetail(Q[TYPE_S][i]);
      if (q == lastptr) {
        last = i;
	dequeue_r(Q[TYPE_S][i]);
        continue;
      }
      c1 = q[-1];  c2 = q[0];
      printf("%d: %d(%c) %d(%c)\n",i,c1,c1,c2,c2);
      if (c1 <= c2) { // TYPE_S
	//printf("copy %d to %d\n",i,M2[c1+1]);
	q = dequeue_r(Q[TYPE_S][i]);
	enqueue_l(Q[TYPE_S][c1+1],q-1);
	bw2[M2[c2+1]--] = c1;
      } else {
	bw2[i] = LMS_BW[c2+1][--C3[c2+1]];
      }
    }
    while (!emptyqueue(Q[TYPE_L][i])) {
      q = quetail(Q[TYPE_L][i]);
      if (q == lastptr) {
        last = i;
	dequeue_r(Q[TYPE_S][i]);
        continue;
      }
      c1 = q[-1];  c2 = q[0];
      printf("%d: %d(%c) %d(%c)\n",i,c1,c1,c2,c2);
      if (c1 <= c2) { // TYPE_S
	//printf("copy %d to %d\n",i,M2[c1+1]);
	q = dequeue_r(Q[TYPE_L][i]);
	enqueue_l(Q[TYPE_S][c1+1],q-1);
	bw2[M2[c2+1]--] = c1;
      } else {
	bw2[i] = LMS_BW[c2+1][--C3[c2+1]];
      }
    }
  }
#endif



}
#endif

/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
