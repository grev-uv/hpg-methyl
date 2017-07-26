/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include "csa.h"
#include "approx.h"

#define INDEX(x,y) ((x)*(1+keylen)+(y))

#define GAP_PENALTY (1)

int nnodes;

#define BLOCKSIZE (1<<20)
approx_list *approx_list_init(void)
{
  approx_list *list;
  list = (approx_list *) malloc(sizeof(approx_list));
  if (list == NULL) {
    printf("approx_list_init: malloc failed.\n");
    exit(1);
  }
  list->p = (approx_set *) malloc(sizeof(approx_set)*BLOCKSIZE);
  if (list->p == NULL) {
    printf("approx_list_init: malloc failed.\n");
    exit(1);
  }
  list->len = BLOCKSIZE;
  list->num = 0;
  return list;
}

void approx_list_append(approx_list *list, int error, int len, i64 l, i64 r)
{
  if (list->num == list->len) {
    list->len += BLOCKSIZE;
//    printf("list->p = %p ",list->p);
    list->p = (approx_set *) realloc(list->p, sizeof(approx_set)*list->len);
//    printf("realloc list->p = %p\n",list->p);
  }
  list->p[list->num].error = error;
  list->p[list->num].len = len;
  list->p[list->num].l = l;
  list->p[list->num].r = r;
  list->num++;
}

void approx_list_print(CSA *csa, approx_list *list)
{
  int i;
  uchar *text;
  for (i=0; i<list->num; i++) {
    printf("[%ld,%ld] error=%d ",list->p[i].l,list->p[i].r,list->p[i].error);
    text = (uchar *) mymalloc(list->p[i].len+1);
    csa->substring(text, csa, list->p[i].l, list->p[i].len);
    text[list->p[i].len] = 0;
    printf("%s\n",text);
    free(text);
  }
}

void approx_list_print2(CSA *csa, approx_list *list)
{
  int i;
  uchar *text;
  for (i=0; i<list->num; i++) {
    printf("[%ld,%ld] error=%d ",list->p[i].l,list->p[i].r,list->p[i].error);
    text = (uchar *) mymalloc(list->p[i].len+1);
    csa->text(text, csa, list->p[i].l, list->p[i].l+list->p[i].len-1);
    text[list->p[i].len] = 0;
    printf("%s\n",text);
    free(text);
  }
}

void approx_list_clear(approx_list *list)
{
  list->num = 0;
}

void approx_list_free(approx_list *list)
{
  free(list->p);
  free(list);
}

static void csa_approxsearch_sub(uchar *key,int keylen,int k, int *score, i64 l, i64 r, int depth, int dofs, uchar *buf, int lr, CSA *csa, approx_list *list)
{
  i64 i,j,c,cc;
  i64 lll, rrr;
  uchar set[SIGMA];
  rank_t ll[SIGMA], rr[SIGMA];
  int mm;
  int num;

  depth++;
  score[INDEX(depth,0)] = score[INDEX(depth-1,0)]+1;

  if (lr == 0) {
    num = csa->child_l(csa, l, r, set, ll, rr);
  } else {
    num = csa->child_r(csa, dofs+depth-1, l, r, set, ll, rr);
  }
  for (i=0; i<num; i++) {
    nnodes++;
//    if (nnodes % 10000 == 0) printf("%d\n",nnodes);
    c = set[i];
    lll = ll[i];  rrr = rr[i];
    buf[depth-1] = c;
    mm = score[INDEX(depth,0)];
    for (j=1; j<=keylen; j++) {
//        printf("key[%d] = %c\n",keylen-j,key[keylen-j]);
      if (lr == 0) cc = key[keylen-j]; else cc = key[j-1];
      score[INDEX(depth,j)]
       = min(score[INDEX(depth-1,j-1)] + (c != cc),
         min(score[INDEX(depth-1,j)] +GAP_PENALTY,
             score[INDEX(depth,j-1)] +GAP_PENALTY));
        mm = min(mm, score[INDEX(depth,j)]);
    }
    if (score[INDEX(depth,keylen)] <= k) {
#if 0
      printf("match k=%d ",score[INDEX(depth,keylen)]);
      for (j=depth-1; j>=0; j--) printf("%c",buf[j]);
      printf(" %ld [%ld,%ld]\n",rrr-lll+1,lll,rrr);
#endif
      approx_list_append(list, score[INDEX(depth,keylen)], dofs+depth, lll, rrr);
    }
    if (mm <= k) {
      csa_approxsearch_sub(key,keylen,k,score,lll,rrr,depth,dofs,buf,lr,csa, list);
    }
  }
}

approx_list *csa_approxsearch(unsigned char *key,int keylen,int k,CSA *csa)
{

  
  i64 i,l,r;
  int *score;
  uchar *buf;
  int buflen;
  approx_list *list;

  score = (int *) malloc(sizeof(int)*(1+keylen)*(1+keylen+k+1));
  buflen = keylen+k+1;
  buf = (uchar *) malloc(buflen);
//  printf("buf %d\n",buflen);
  for (i=0; i<=keylen; i++) score[INDEX(0,i)] = i * GAP_PENALTY;
  for (i=0; i<=keylen+k; i++) score[INDEX(i,0)] = i * GAP_PENALTY;
  for (i=0; i<=keylen+k; i++) buf[i] = 0;

  l = 1;  r = csa->n;

  list = approx_list_init();

  nnodes = 0;
  csa_approxsearch_sub(key, keylen, k, score, l, r, 0, 0, buf, 0, csa, list);
  printf("#visited nodes = %d\n",nnodes);
//  approx_list_print(csa, list);
//  approx_list_free(list);

  free(score);
  free(buf);

  return list;
}

approx_list *csa_approxsearch2(unsigned char *key,int keylen,int k,CSA *csa)
{


  i64 i,j,l,r;
  int *score;
  uchar *buf;
  int buflen;
  approx_list *list, *list1, *list2;
  
  int l1, l2, e1, e2;

  e1 = k/2;
  e2 = k-e1;
  l1 = keylen/2;
  l2 = keylen-l1;

  score = (int *) malloc(sizeof(int)*(1+keylen)*(1+keylen+k+1));
  buflen = keylen+k+1;
  buf = (uchar *) malloc(buflen);
//  printf("buf %d\n",buflen);
  for (i=0; i<=keylen; i++) score[INDEX(0,i)] = i;
  for (i=0; i<=keylen+k; i++) score[INDEX(i,0)] = i;
  for (i=0; i<=keylen+k; i++) buf[i] = 0;

  l = 1;  r = csa->n;

  list = approx_list_init();

  nnodes = 0;

  list1 = approx_list_init();
  for (j=0; j<=l1; j++) score[INDEX(0,j)] = j;
  csa_approxsearch_sub(key, l1, e1, score, l, r, 0, 0, buf, 0, csa, list1);
//  printf("#visited nodes = %d\n",nnodes);

  list2 = approx_list_init();
  for (i=0; i<list1->num; i++) {
    int e, dofs;
    i64 ll, rr;
    e = list1->p[i].error;
    dofs = list1->p[i].len;
    ll = list1->p[i].l;
    rr = list1->p[i].r;
    for (j=0; j<=l2; j++) score[INDEX(0,j)] = e+j;
    csa_approxsearch_sub(key+l1, l2, k, score, ll, rr, 0, dofs, buf, 1, csa, list2);
    approx_list_clear(list2);
  }
  approx_list_clear(list1);
//  printf("#visited nodes = %d\n",nnodes);

  for (j=0; j<=l2; j++) score[INDEX(0,j)] = j;
  csa_approxsearch_sub(key+l1, l2, e2, score, l, r, 0, 0, buf, 0, csa, list1);
//  printf("#visited nodes = %d\n",nnodes);

  for (i=0; i<list1->num; i++) {
    int e, dofs;
    i64 ll,rr;
    e = list1->p[i].error;
    dofs = list1->p[i].len;
    ll = list1->p[i].l;
    rr = list1->p[i].r;
    for (j=0; j<=l1; j++) score[INDEX(0,j)] = e+j;
    csa_approxsearch_sub(key, l1, k, score, ll, rr, 0, dofs, buf, 0, csa, list2);
    approx_list_clear(list2);
  }

//  approx_list_free(list);
  approx_list_free(list1);
  approx_list_free(list2);


  printf("#visited nodes = %d\n",nnodes);

  free(score);
  free(buf);
  
  return list;
}
