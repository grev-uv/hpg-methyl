/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
typedef struct approxmatch_ {
  int err;
  int l,r;
} approxmatch;

typedef struct {
  int error;
  int len;
  i64 l,r;
} approx_set;

typedef struct {
  i64 len;
  i64 num;
  approx_set *p;
} approx_list;

void approx_list_print(CSA *csa, approx_list *list);
void approx_list_print2(CSA *csa, approx_list *list);
approx_list *approx_list_init(void);
void approx_list_append(approx_list *list, int error, int len, i64 l, i64 r);
void approx_list_clear(approx_list *list);
void approx_list_free(approx_list *list);

approx_list *csa_approxsearch(unsigned char *key,int keylen,int k,CSA *csa);
approx_list *csa_approxsearch2(unsigned char *key,int keylen,int k,CSA *csa);
