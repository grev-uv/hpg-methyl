/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "psi2.h"

#define display_progressbar(str,i,n) if (i % 10000000 == 0) {fprintf(stderr,"%s %ld/%ld                       \r",str,i/10000000,n/10000000);  fflush(stderr); }

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

static int blog(i64 x) // [0,n] の数を格納するには blog(n)+1 ビット必要
{
int l;
  l = -1;
  while (x>0) {
    x>>=1;
    l++;
  }
  return l;
}

static void writeint(int k,i64 x,FILE *f)
{
  int i;
  for (i=k-1; i>=0; i--) {
    fputc(x & 0xff,f);
    x >>= 8;
  }
}

void psi2_options(CSA *csa, char *p)
{
  psi2 *ps;

  csa->psi_struc = ps = (psi2 *) mymalloc(sizeof(psi2));
  ps->id = csa->id;
  if (p[0] == 0) goto end;
  p++;
#if 0
  if (p[0] != ':') {
    sscanf(p,"%d",&ps->L);
    printf("L = %d\n",ps->L);
    while (p[0] != ':' && p[0] != 0) p++;
  }
#endif
  if (p[0] == 0) goto end;
  p++;
end:
  printf("psi_id = %ld\n",ps->id);
}

static i64 psi2_psi(CSA *csa, i64 i)
{
  i64 x;
  psi2 *ps;

  ps = (psi2 *)csa->psi_struc;
  x = sparsearray4_select(ps->sa, i+1);
  x %= csa->n + 1;

  return x;
}

i64 psi2_makeindex(CSA *csa, char *fname)
{
i64 psize;
i64 i,j,x;
 i64 y,d;
int k;
FILE *f2;
char *fpsi;
//psi1_iterator *pi;
i64 n,m;
i64 nn, mm;
psi2 *ps;
int id;
FILE *out;
diskbuf *psi;
char *fbw, *flst;
sparsearray4 sa;

  ps = (psi2 *)csa->psi_struc;
  id = csa->id;

  k = strlen(fname);
  fbw = (char *) mymalloc(k+5);
  flst = (char *) mymalloc(k+5);
  fpsi = (char *) mymalloc(k+5);
  sprintf(fbw,"%s.bw",fname);
  sprintf(flst,"%s.lst",fname);

  sprintf(fpsi,"%s.psa",fname);

  out = create_tmp(0);
  bw_to_psi(out,csa,fbw,flst,&k);

  m = 0;
  for (i=0; i<SIGMA; i++) {
    if (csa->C[i] > 0) m++;
  }

  psi = open_diskbuf(out,k);
  ps->last = getint_diskbuf(psi,0);
  printf("last = %ld\n",ps->last);

  n = csa->n;

  mm = n+1;  nn = (m+1) * (n+1);
  sparsearray4_construct_init(&sa, nn, mm);

  y = 0;  d = 0;
  for (j=0; j<=n; j++) {
    display_progressbar("compressing psi ",j,n);
    x = getint_diskbuf(psi,j);
    if (x <= y) d++;
    sparsearray4_construct_set(&sa, j, d*(n+1)+x);
    y = x;
  }
  sparsearray4_construct_end(&sa,0, SDARRAY_SELECT1);

  f2 = fopen(fpsi,"wb");
  psize = 0;

  ps->k = k = (blog(n+1)+1+8-1)/8;

  writeint(1,ID_PSI,f2);
  writeint(1,k,f2); /* #bytes of integer */
  writeint(k,n,f2);
  psize += 1+1+k;

  writeint(1,id,f2);
  psize += 1;

  psize += sparsearray4_write(&sa, f2);

  printf("size %ld (%1.3f bpc)\n",psize,(double)psize*8 / n);

  fclose(f2);

//  psi1_iterator_remove(pi);

  close_diskbuf(psi);
  fclose(out);
  remove_tmp(0);

  psi2_read(csa, fpsi);

  free(fpsi);
  free(fbw);
  free(flst);

  return psize;
}



i64 psi2_read(CSA *csa, char *fname)
{
  i64 psize;
  i64 n;
  int k,id;
  psi2 *ps;
  uchar *p, *q;
 
  csa->psi_struc = ps = (psi2 *) mymalloc(sizeof(psi2));

  printf("psi_read: map %s\n",fname);
  ps->mappsi = mymmap(fname);
  if (ps->mappsi->addr==NULL) {
    perror("psi2_read: mmap2\n");
    exit(1);
  }
  p = q = (uchar *)ps->mappsi->addr;

  id = getuint(p,0,1);  p += 1;
  if (id != ID_PSI) {
    printf("read_psi: id = %d is not supported.\n",id);
    exit(1);
  }
  ps->k = k = getuint(p,0,1);  p += 1;
  ps->n = n = getuint(p,0,k);  p += k;

  id = getuint(p,0,1);  p += 1;
//  printf("read_psi: psi_id = %d\n",id);
  csa->id = ps->id = id;
  switch (id) {
  case ID_SPARSE4:
    printf("#psi format = SPARSE4\n");
    break;
  default:
    printf("read_csa: ID %d is not supported.\n",id);
    break;
  }


  ps->sa = (sparsearray4 *) mymalloc(sizeof(*ps->sa));
  sparsearray4_read(ps->sa, &p);
  psize = p - q;

//  printf("psi2_read: psize = %ld\n",psize);
  ps->psize = psize;


// user-specific functions
  csa->psi = psi2_psi;
  csa->psi_succ = csa_psi_succ_naive;
  csa->psi_pred = csa_psi_pred_naive;

// default functions
  csa->lookup = csa_lookup;
  csa->inverse = csa_inverse;
  csa->text = csa_text;
  csa->substring = csa_substring;
  csa->T = csa_T;
  csa->head = csa_head_rank;
  csa->search = csa_search;


  return psize;
}

