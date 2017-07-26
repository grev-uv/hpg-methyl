/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lf_wt.h"

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

static void writeint(int k,i64 x,FILE *f)
{
  int i;
  for (i=k-1; i>=0; i--) {
    fputc(x & 0xff,f);
    x >>= 8;
  }
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

void lf_wt_options(CSA *csa, char *p)
{
  lf_wt *lf;
  int L,opt;

  csa->psi_struc = lf = (lf_wt *) mymalloc(sizeof(lf_wt));
  lf->L = 512;
  lf->opt = opt = 0;

  if (p[0] == 0) goto end;
  p++;
  if (p[0] != ':') {
    if (sscanf(p,"%d",&L)==1) lf->L = L;
    if (L < sizeof(u64)*8) {
      printf("L must be >= %lu\n",sizeof(u64)*8);
      exit(1);
    }
//    printf("L = %d\n",lf->l);
    while (p[0] != ':' && p[0] != 0) p++;
  }

  if (p[0] == 0) goto end;
  p++;
  if (p[0] != ':') {
    if (sscanf(p,"%d",&opt)==1) {
      if (opt & 1) {
        lf->opt |= ID_COMPPTR;
        printf("COMPPTR\n");
      }
      if (opt & 2) {
        lf->opt |= ID_RR;
        printf("RR\n");
      }
      if (opt & 4) {
        lf->opt |= ID_SUC;
        printf("SUC\n");
      }
    }
    while (p[0] != ':' && p[0] != 0) p++;
  }
end:
  ;
//  printf("BW_WT L=%d\n",lf->L);
}

static void make_huffman_tree(CSA *csa, lf_wt *lf)
{
  int i,m;
  double *freq;

  for (m=0,i=0; i<SIGMA; i++) {
    if (csa->C[i]>0) {
//      printf("C[%d] = %d\n",i,csa->C[i]);
      csa->AtoC[m] = i;
      csa->CtoA[i] = m;
      m++;
    }
  }
  csa->m = m;

  freq = (double *) mymalloc(m*sizeof(freq[0]));
  for (i=0; i<m; i++) {
    freq[i] = (double)(csa->C[csa->AtoC[i]]+1) / csa->n;
//    printf("freq[%d] = %lf\n",i,freq[i]);
  }
  lf->huffman = MakeHuffmanTree(m, freq);

#if 0
  for (i=0; i<m; i++) {
    printf("%d: len = %d code = %lx\n",
           i,lf->huffman->clen[i],lf->huffman->code[i]);
  }
#endif
  free(freq);
}

static void make_wavelet_sub(lf_wt *lf, i64 n, int depth, int node)
{
  i64 i,nl,nr = 0;
  int c,d,m,w;
  FILE *in,*fl,*fr;
  Huffman *h;
  u64 x;
  i64 stmp;
  int opt;

  printf("sub node=%d depth=%d n=%ld\n",node,depth,n);
  fflush(stdout);
  h = lf->huffman;
  m = h->n;

  if (node < m) {
    remove_tmp(node);
    return;
  }

  w = sizeof(x)*8;

//  lf->init(lf->da[node-m], n);
  switch (lf->id) {
  case ID_BWT_WT:
  case ID_BWT_WT_HUF:
  case ID_BWT_WT_RR:
    lf->da[node-m] = (comparray *) mymalloc(sizeof(comparray));
    comparray_construct_init( (comparray *) lf->da[node-m], n);
    break;
  case ID_BWT_WT_DENSE:
    lf->da[node-m] = (densearray *) mymalloc(sizeof(densearray));
    densearray_construct_init((densearray *) lf->da[node-m], n);
    break;
  case ID_BWT_WT_SPARSE4:
    lf->da[node-m] = (sparsearray4 *) mymalloc(sizeof(sparsearray4));
    in = open_tmp(node);
    for (i=0; i<n; i++) {
      c = fgetc(in);
      x = h->code[c];
      d = (x >> (w-1-depth)) & 1;
      if (d == 1) nr++;
    }
    fclose(in);
    sparsearray4_construct_init((sparsearray4 *) lf->da[node-m], n, nr);
    break;
  }

  in = open_tmp(node);

  fl = create_tmp(h->left[node]);
  fr = create_tmp(h->right[node]);

  nl = nr = 0;
  for (i=0; i<n; i++) {
    c = fgetc(in);
//    printf("BW[%ld] = %d\n",i,c);
    x = h->code[c];
    d = (x >> (w-1-depth)) & 1;
//    comparray_construct_set(lf->da[node-m], i, d);
    if (lf->id != ID_BWT_WT_SPARSE4) {
      lf->set(lf->da[node-m], i, d);
    } else {
      if (d==1) sparsearray4_construct_set((sparsearray4 *) lf->da[node-m], nr, i);
    }
    if (d == 0) {
      fputc(c, fl);
      nl++;
    } else {
      fputc(c, fr);
      nr++;
    }
  }
  fclose(in);
  fclose(fl);  fclose(fr);
  remove_tmp(node);

  opt = lf->opt;
  opt |= SDARRAY_RANK1;
  if (lf->id == ID_BWT_WT_HUF) {
    opt |= SDARRAY_LENHUF;
  }
  if (lf->id == ID_BWT_WT_RR) {
    opt |= SDARRAY_RR;
  }
  if (lf->opt & ID_COMPPTR) {
    opt |= SDARRAY_COMPRANK | SDARRAY_COMPPTR;
  }
  if (lf->opt & ID_RR) {
    opt |= SDARRAY_RR;
  }
  if (lf->opt & ID_SUC) {
    opt |= SDARRAY_SUC;
  }
  if (lf->id == ID_BWT_WT_DENSE) {
    opt |= (lf->L << 16);
  }
  stmp = lf->end(lf->da[node-m], lf->L, opt);
  printf("comparray size=%ld (%1.3f bpc)\n",stmp,(double)stmp*8/n);
  fflush(stdout);
  lf->psize += stmp;
//  lf->psize += (n+DD-1)/8;
  make_wavelet_sub(lf, nl, depth+1, h->left[node]);
  make_wavelet_sub(lf, nr, depth+1, h->right[node]);
}

static void make_wavelet(CSA *csa, lf_wt *lf, char *fbw)
{
  FILE *in,*out;
  i64 i,n;
  int c,root;

  lf->id = csa->id;

  in = fopen(fbw,"rb");
  if (in == NULL) {
    printf("lf_wt:make_wavelet: cannot open %s\n",fbw);
    exit(1);
  }
  fseek(in,0,SEEK_END);
  n = ftell(in);
  fseek(in,0,SEEK_SET);
  printf("n=%ld\n",n);
  csa->n = lf->n = n;

  for (i=0; i<SIGMA; i++) csa->C[i] = 0;

  fprintf(stderr,"counting...\n");
  for (i = 0; i < n; i++) {
    if (i % 1000000 == 0) {
      fprintf(stderr,"%ld\r",i/1000);
      fflush(stderr);
    }
    c = fgetc(in);
    csa->C[c]++;
  }
  fseek(in,0,SEEK_SET);

  make_huffman_tree(csa, lf);
  
  root = 2 * lf->huffman->n - 2;
  out = create_tmp(root);
  fprintf(stderr,"packing...\n");
  for (i = 0; i < n; i++) {
    if (i % 1000000 == 0) {
      fprintf(stderr,"%ld\r",i/1000);
      fflush(stderr);
    }
    c = fgetc(in);
    fputc(csa->CtoA[c],out);
  }
  fclose(out);
  
  for (i=0; i<SIGMA; i++) lf->da[i] = NULL;

  lf->psize = 0;
  make_wavelet_sub(lf, n, 0, root);
  printf("wavelet: size = %ld (%1.3f bpc) \n",lf->psize,
          (double)lf->psize*8/n);
  fclose(in);
}

#if 0
static int lf_wt_BW_sub_old(lf_wt *lf,i64 i, int node)
{
  int c,m;
  int c2;
  i64 r;
  i64 r2;
  m = lf->huffman->n;

  if (node < m) return node;

  c = lf->getbit(lf->da[node-m], i);
  if (c == 0) {
    r = lf->rank0(lf->da[node-m], i) - 1;
    return lf_wt_BW_sub_old(lf, r, lf->huffman->left[node]);
  } else {
    r = lf->rank1(lf->da[node-m], i) - 1;
    return lf_wt_BW_sub_old(lf, r, lf->huffman->right[node]);
  }
  return node;
}

static int lf_wt_BW_old(CSA *csa,i64 i)
{
  lf_wt *lf;
  int root;

  lf = (lf_wt *)csa->psi_struc;
  if (i == lf->last) return -1;
  if (i >= lf->last) i--;
  root = 2 * lf->huffman->n - 2;
  return lf_wt_BW_sub_old(lf, i, root);
}
#endif

static int lf_wt_BW(CSA *csa, i64 i)
{
  lf_wt *lf;
  int node;
  int c,m;
  i64 r;

  lf = (lf_wt *)csa->psi_struc;
  if (i == lf->last) return -1;
  if (i >= lf->last) i--;

  m = lf->huffman->n;
  node = 2 * m - 2;

  while (node >= m && i >= 0) {
    r = lf->rank_bit(lf->da[node-m], i, &c);
    if (c == 0) {
      i = (i+1) - r - 1;
      node = lf->huffman->left[node];
    } else {
      i = r - 1;
      node = lf->huffman->right[node];
    }
//    if (i < 0) {
//      printf("BW: i=%ld\n",i);
//   }
  }
  return node;
}

static int lf_wt_BW_rank(CSA *csa, i64 i, rank_t *rank)
{
  lf_wt *lf;
  int node;
  int c,m;
  i64 r;

  lf = (lf_wt *)csa->psi_struc;
  if (i == lf->last) return -1;
  if (i >= lf->last) i--;
  m = lf->huffman->n;
  node = 2 * m - 2;

  while (node >= m && i >= 0) {
    r = lf->rank_bit(lf->da[node-m], i, &c);
    if (c == 0) {
      i = (i+1) - r - 1;
      node = lf->huffman->left[node];
    } else {
      i = r - 1;
      node = lf->huffman->right[node];
    }
//    if (i < 0) {
//      printf("BW_rank: i=%ld\n",i);
//    }
  }
  *rank = i+1;
  return node;
}

static rank_t lf_wt_LF(CSA *csa, rank_t i)  /* 0 <= i <= n */
{
  int c;
  i64 rank;

  c = lf_wt_BW_rank(csa, i, &rank);
  if (c == -1) return 0; // last
  return csa->K[c+1]-1 + rank;
}

/********************
static i64 lf_wt_rankc_sub_old(lf_wt *lf, i64 i, u64 x, int depth, int node)
{
  int d,m;
  i64 r;

  if (i < 0) {
//    printf("rankc_sub: i=%ld\n",i);
    return 0;
  }

  m = lf->huffman->n;
  if (node < m) return i+1;

  d = (x >> (sizeof(x)*8-1-depth)) & 1;
  if (d == 0) {
//    r = comparray_rank0(lf->da[node-m], i) - 1;
    r = lf->rank0(lf->da[node-m], i) - 1;
    return lf_wt_rankc_sub_old(lf, r, x, depth+1, lf->huffman->left[node]);
  } else {
//    r = comparray_rank(lf->da[node-m], i) - 1;
    r = lf->rank1(lf->da[node-m], i) - 1;
    return lf_wt_rankc_sub_old(lf, r, x, depth+1, lf->huffman->right[node]);
  }
}

static i64 lf_wt_rankc_old(CSA *csa, i64 i, int c)
{
  lf_wt *lf;
  int root;
  u64 x;

  lf = (lf_wt *)csa->psi_struc;
  if (i >= lf->last) i--;
  root = 2 * lf->huffman->n - 2;
  x = lf->huffman->code[c];
  return lf_wt_rankc_sub_old(lf,i, x, 0, root);
}
*********************/

static i64 lf_wt_rankc(CSA *csa, i64 i, int c)
{
  lf_wt *lf;
  int node;
  u64 x;
  int d,m;


  lf = (lf_wt *)csa->psi_struc;
  if (i >= lf->last) i--;
  if (i < 0) return 0;
  m = lf->huffman->n;
  node = 2 * m - 2;
  x = lf->huffman->code[c];

  while (node >= m && i >= 0) {
    d = (x >> (sizeof(x)*8-1)) & 1;
    if (d == 0) {
      i = lf->rank0(lf->da[node-m], i) - 1;
      node = lf->huffman->left[node];
    } else {
      i = lf->rank1(lf->da[node-m], i) - 1;
      node = lf->huffman->right[node];
    }
    x <<= 1;
//    if (i < 0) {
//      printf("rankc: i=%ld\n",i);
//    }
  }
  return i+1;
}

static void lf_wt_text(uchar *p,CSA *csa, pos_t s, pos_t t)
{
  i64 i;
  i64 j,c;
  i64 rank = 0;

#if DEBUG
  if (s < 1 || t > csa->n || s > t) {
	  printf("csa_text: out of range [%d,%d] n=%d\n",s,t,csa->n);
  }
#endif

  i = csa->inverse(csa,(t+1));
  for (j=t-s; j>=0; j--) {
    c = lf_wt_BW_rank(csa, i, &rank);
#if 0
    if (c != lf_wt_BW_old(csa, i)) {
      printf("error c = %d old = %d\n",c,lf_wt_BW_old(csa,i));
    }
    if (rank != lf_wt_rankc_old(csa, i, c)) {
      printf("error rank = %ld old = %ld\n",rank,lf_wt_rankc_old(csa, i, c));
    }
#if 1
    if (c == -1) {
      printf("lf_wt_text: ??? j = %ld i = %ld c = %ld  s=%ld t=%ld \n",
              j,i,c,s,t);
      exit(1);
    }
#endif
    printf("j=%ld rank=%ld c=%d\n",j,rank,c);
#endif
    p[j] = csa->AtoC[c];
    i = csa->K[c+1]-1 + rank;
  }
}

void lf_wt_child_l_sub(CSA *csa, i64 l, i64 r, int node, int *num, uchar *head, rank_t *ll, rank_t *rr)
{
  lf_wt *lf;
  int m;
  i64 rank, rank1, rank2;


  lf = (lf_wt *)csa->psi_struc;
  m = lf->huffman->n;

#if 0
  if (node >= m) {
    for (i=l; i<=r; i++) {
      printf("b[%ld] = %d \n",i,lf->getbit(lf->da[node-m], i));
    }
  }
#endif

  if (node < m) {
//    printf("found node = %d  [%ld,%ld] num=%d\n",node, l, r, *num);
    head[*num] = csa->AtoC[node];
    ll[*num] = csa->K[node+1]+l;
    rr[*num] = csa->K[node+1]+r;
    *num = *num+1;
    return;
  }

  rank1 = lf->rank1(lf->da[node-m], r);
  if (l > 0) rank2 = lf->rank1(lf->da[node-m], l-1); else rank2 = 0;
  rank = rank1 - rank2;

  if (rank < r-l+1) { // there exists at least one 0
    lf_wt_child_l_sub(csa, l-rank2, r-rank1, lf->huffman->left[node], num, head, ll, rr);
  }
  if (rank > 0) { // there exists at least one 1
    lf_wt_child_l_sub(csa, rank2, rank1-1, lf->huffman->right[node], num, head, ll, rr);
  }


}

int lf_wt_child_l(CSA *csa, i64 l, i64 r, uchar *head, rank_t *ll, rank_t *rr)
{
  lf_wt *lf;
  int node;
  int i,c,m, num;
  rank_t *ltmp, *rtmp;
  uchar *headtmp;
  int *hc;

  lf = (lf_wt *)csa->psi_struc;
  m = lf->huffman->n;
  if (l == lf->last) return 0;
  if (l >= lf->last) l--;
  if (r >= lf->last) r--;
  node = 2 * m - 2;

  headtmp = (uchar*) mymalloc(m * sizeof(*headtmp));
  ltmp = (rank_t *) mymalloc(m * sizeof(*ltmp));
  rtmp = (rank_t *) mymalloc(m * sizeof(*ltmp));
  hc = (int *) mymalloc(SIGMA * sizeof(*hc));

  num = 0;
  lf_wt_child_l_sub(csa, l, r, node, &num, headtmp, ltmp, rtmp);

  for (c=0; c<SIGMA; c++) hc[c] = -1;
  for (i=0; i<num; i++) hc[headtmp[i]] = i;
  i = 0;
  for (c=0; c<SIGMA; c++) {
    if (hc[c] >= 0) {
      head[i] = headtmp[hc[c]];
      ll[i] = ltmp[hc[c]];
      rr[i] = rtmp[hc[c]];
      i++;
    }
  }

  free(hc);  free(rtmp);  free(ltmp);  free(headtmp);

  return num;
}



i64 lf_wt_makeindex(CSA *csa, char *fname)
{
  i64 last;
  FILE *in,*out;
  
  char *fbw, *flst, *fbw2, *fidx;
  int k;
  i64 size;
  i64 n,i;
  lf_wt *lf;
  int nn;
  
  comparray_maketbl();

  
  lf = (lf_wt *)csa->psi_struc;

  k = strlen(fname);
  fbw = (char *) mymalloc(k+5);
  flst = (char *) mymalloc(k+5);
  fbw2 = (char *) mymalloc(k+5);
  fidx = (char *) mymalloc(k+5);
  sprintf(fbw,"%s.bw",fname);
  sprintf(flst,"%s.lst",fname);
  switch (csa->id) {
  case ID_BWT_WT:
    sprintf(fidx,"%s.wtd",fname);
    printf("BW_WT L=%d\n",lf->L);
    break;
  case ID_BWT_WT_HUF:
    sprintf(fidx,"%s.whd",fname);
    printf("BW_WT_HUF L=%d\n",lf->L);
    break;
  case ID_BWT_WT_DENSE:
    sprintf(fidx,"%s.wda",fname);
    printf("BW_WT_DENSE\n");
    break;
  case ID_BWT_WT_SPARSE4:
    sprintf(fidx,"%s.wsa",fname);
    printf("BW_WT_SPARSE4\n");
    break;
  case ID_BWT_WT_RR:
    sprintf(fidx,"%s.wxd",fname);
    printf("BW_WT_RR\n");
    break;
  }

  switch (csa->id) {
  case ID_BWT_WT:
  case ID_BWT_WT_HUF:
  case ID_BWT_WT_RR:
    lf->init = (int (*)(void *, i64))comparray_construct_init;
    lf->set = (int (*)(void *, i64, int))comparray_construct_set;
    lf->end = (i64 (*)(void *, ushort, int))comparray_construct_end;
    lf->write = (i64 (*)(void *, FILE *))comparray_write;
    lf->getbit = (int (*)(void *, i64))comparray_getbit;
    lf->rank0 = (rank_t (*)(void *, rank_t))comparray_rank0;
    lf->rank1 = (rank_t (*)(void *, rank_t))comparray_rank;
    lf->rank_bit = (rank_t (*)(void *, rank_t, int *))comparray_rank_and_bit;
    break;
  case ID_BWT_WT_DENSE:
    lf->init = (int (*)(void *, i64))densearray_construct_init;
    lf->set = (int (*)(void *, i64, int))densearray_construct_set;
    lf->end = (i64 (*)(void *, ushort, int))densearray_construct_end;
    lf->write = (i64 (*)(void *, FILE *))densearray_write;
    lf->getbit = (int (*)(void *, i64))densearray_getbit;
    lf->rank0 = (rank_t (*)(void *, rank_t))densearray_rank0;
    lf->rank1 = (rank_t (*)(void *, rank_t))densearray_rank;
    lf->rank_bit = (rank_t (*)(void *, rank_t, int *))densearray_rank_and_bit;
    break;
  case ID_BWT_WT_SPARSE4:
    lf->init = (int (*)(void *, i64))sparsearray4_construct_init;
    lf->set = (int (*)(void *, i64, int))sparsearray4_construct_set;
    lf->end = (i64 (*)(void *, ushort, int))sparsearray4_construct_end;
    lf->write = (i64 (*)(void *, FILE *))sparsearray4_write;
    lf->getbit = (int (*)(void *, i64))sparsearray4_getbit;
    lf->rank0 = (rank_t (*)(void *, rank_t))sparsearray4_rank0;
    lf->rank1 = (rank_t (*)(void *, rank_t))sparsearray4_rank;
    lf->rank_bit = (rank_t (*)(void *, rank_t, int *))sparsearray4_rank_and_bit;
    break;
  }

  in = fopen(flst,"r");
  if (in == NULL) {
    perror("lf_wt_makeindex:");  exit(1);
  }
  int res = fscanf(in,"%ld",&last);
  if (!res) {}
  printf("last = %ld\n",last);
  lf->last = last;
  fclose(in);


  make_wavelet(csa, lf, fbw);
  n = lf->n;
  lf->k = k = (blog(n+1)+1+8-1)/8;

  out = fopen(fidx,"w");
  if (out == NULL) {
    printf("lf_dna_makeindex: cannot open %s\n",fidx);
  }

  size = 0;

  writeint(1,ID_LF,out);
  writeint(1,k,out); /* #bytes of integer */
  writeint(k,n,out);
  writeint(k,last,out);
  writeint(k,1,out); // dummy L
  size += 1+1+3*k;

  writeint(1,lf->id,out);
  size += 1;

  Huffman_write(lf->huffman, out);

  nn = lf->huffman->n;
  for (i=nn; i<=2*nn-2; i++) {
//    printf("write node %d\n",i);
    //size += densearray_write(lf->da[i-nn], out);
//    size += comparray_write(lf->da[i-nn], out);
    size += lf->write(lf->da[i-nn], out);
  }

  fclose(out);

  csa->BW = lf_wt_BW;
  csa->rankc = lf_wt_rankc;

//  csa->LF = csa_LF;
  csa->LF = lf_wt_LF;

  free(fbw);
  free(flst);
  free(fidx);

  return size;
}

void lf_wt_read(CSA *csa, char *fname)
{
  
  int k;
  i64 bsize;
  i64 n,i,id;
  lf_wt *lf;
  int nn;
  char *fidx;
  uchar *p, *q;

  comparray_maketbl();
  
  csa->psi_struc = lf = (lf_wt *)mymalloc(sizeof(lf_wt));

  fidx = fname;
//  printf("psi_read: read %s\n",fidx);

  lf->map = mymmap(fidx);
  if (lf->map->addr==NULL) {
    perror("lf_wt_read: mmap2\n");
    exit(1);
  }
  p = (uchar *)lf->map->addr;
  q = p;

  id = getuint(p,0,1);  p += 1;
  if (id != ID_LF) {
    printf("lf_dna_read: id = %ld is not supported.\n",id);
    exit(1);
  }
  lf->k = k = getuint(p,0,1);  p += 1;
  lf->n = n = getuint(p,0,k);  p += k;
  lf->last = getuint(p,0,k);  p += k;
  p += k; // dymmy L

  id = getuint(p,0,1);  p += 1;
  lf->id = id;

//  printf("lf_wt_read: lf_id = %ld\n",id);
  switch (id) {
    case ID_BWT_WT:
      printf("#lf format = BWT_WT\n");
      break;
    case ID_BWT_WT_HUF:
      printf("#lf format = BWT_WT_HUF\n");
      break;
    case ID_BWT_WT_DENSE:
      printf("#lf format = BWT_WT_DENSE\n");
      break;
    case ID_BWT_WT_SPARSE4:
      printf("#lf format = BWT_WT_SPARSE4\n");
      break;
    case ID_BWT_WT_RR:
      printf("#lf format = BWT_WT_RR\n");
      break;
    default:
      printf("lf_wt_read: ID %ld is not supported.\n",id);
      break;
  }

  switch (id) {
  case ID_BWT_WT:
  case ID_BWT_WT_HUF:
  case ID_BWT_WT_RR:
//    lf->init = (int (*)(void *, i64))comparray_construct_init;
//    lf->set = (int (*)(void *, i64, int))comparray_construct_set;
//    lf->end = (i64 (*)(void *, ushort, int))comparray_construct_end;
//    lf->write = (i64 (*)(void *, FILE *))comparray_write;
    lf->getbit = (int (*)(void *, i64))comparray_getbit;
    lf->rank0 = (rank_t (*)(void *, rank_t))comparray_rank0;
    lf->rank1 = (rank_t (*)(void *, rank_t))comparray_rank;
    lf->rank_bit = (rank_t (*)(void *, rank_t, int *))comparray_rank_and_bit;
    break;
  case ID_BWT_WT_DENSE:
//    lf->init = (int (*)(void *, i64))densearray_construct_init;
//    lf->set = (int (*)(void *, i64, int))densearray_construct_set;
//    lf->end = (i64 (*)(void *, ushort, int))densearray_construct_end;
//    lf->write = (i64 (*)(void *, FILE *))densearray_write;
    lf->getbit = (int (*)(void *, i64))densearray_getbit;
    lf->rank0 = (rank_t (*)(void *, rank_t))densearray_rank0;
    lf->rank1 = (rank_t (*)(void *, rank_t))densearray_rank;
    lf->rank_bit = (rank_t (*)(void *, rank_t, int *))densearray_rank_and_bit;
    break;
  case ID_BWT_WT_SPARSE4:
//    lf->init = (int (*)(void *, i64))sparsearray4_construct_init;
//    lf->set = (int (*)(void *, i64, int))sparsearray4_construct_set;
//    lf->end = (i64 (*)(void *, ushort, int))sparsearray4_construct_end;
//    lf->write = (i64 (*)(void *, FILE *))sparsearray4_write;
    lf->getbit = (int (*)(void *, i64))sparsearray4_getbit;
    lf->rank0 = (rank_t (*)(void *, rank_t))sparsearray4_rank0;
    lf->rank1 = (rank_t (*)(void *, rank_t))sparsearray4_rank;
    lf->rank_bit = (rank_t (*)(void *, rank_t, int *))sparsearray4_rank_and_bit;
    break;
  }

  for (i=0; i<SIGMA; i++) lf->da[i] = NULL;

  lf->huffman = Huffman_read2(&p);

  switch (lf->id) {
  case ID_BWT_WT:
  case ID_BWT_WT_HUF:
  case ID_BWT_WT_RR:
    comparray_maketbl();
    break;
  case ID_BWT_WT_DENSE:
    densearray_make_selecttbl();
    break;
  case ID_BWT_WT_SPARSE4:
    break;
  }

  bsize = 0;
  nn = lf->huffman->n;
  for (i=nn; i<=2*nn-2; i++) {
//    printf("read node %d\n",i);

    switch (lf->id) {

    case ID_BWT_WT:
    case ID_BWT_WT_HUF:
    case ID_BWT_WT_RR:
      lf->da[i-nn] = mymalloc(sizeof(comparray));
      comparray_read((comparray *)lf->da[i-nn], &p);
      bsize += ((comparray *)lf->da[i-nn])->bufsize;
      break;
    case ID_BWT_WT_DENSE:
      lf->da[i-nn] = mymalloc(sizeof(densearray));
      densearray_read((densearray *)lf->da[i-nn], &p);
      break;
    case ID_BWT_WT_SPARSE4:
      lf->da[i-nn] = mymalloc(sizeof(sparsearray4));
      sparsearray4_read((sparsearray4 *)lf->da[i-nn], &p);
      break;
    }

//    printf("len %ld\n",lf->da[i-nn]->n);
  }
  lf->psize = p - q;
//  printf("lf_wt: size = %ld\n",lf->psize);
//  printf("bufsize %ld (%1.3f bpc)\n",bsize*sizeof(bitvec_t),
//               (double)bsize*sizeof(bitvec_t)*8/n);


  csa->BW = lf_wt_BW;
  csa->rankc = lf_wt_rankc;
//  csa->rankc = lf_wt_rankc_old;
  csa->LF = lf_wt_LF;
  csa->text = lf_wt_text;
//  csa->text = csa_text_lf;
  csa->BW_rank = lf_wt_BW_rank;
  csa->child_l = lf_wt_child_l;
//  csa->child_l = csa_child_l;

// default functions
//  csa->LF = csa_LF;
  csa->substring_lf = csa_substring_lf;
  csa->selectc = csa_selectc;
  csa->psi = csa_psi_by_rankc_naive;
  csa->psi_succ = csa_psi_succ_naive;
  csa->psi_pred = csa_psi_pred_naive;
  csa->lookup = csa_lookup_lf;
  csa->inverse = csa_inverse_lf;
  csa->search = csa_search_lf;
  csa->searchsub = csa_searchsub_lf;
  csa->child_r = csa_child_r;

}

