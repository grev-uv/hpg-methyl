/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "psi1.h"

#define display_progressbar(str,i,n) if (i % 10000000 == 0) {fprintf(stderr,"%s %ld/%ld                       \r",str,i/10000000,n/10000000);  fflush(stderr); }

#define DD (8*sizeof(short))
#define TBLSIZE (1<<DD)


static int R4[TBLSIZE];
static int R5n[TBLSIZE],R5b[TBLSIZE],R5x[TBLSIZE];

static int getbitD(unsigned short *B, i64 i)
{
  i64 j,l,x;
  i--;
  j = i / DD;
  l = i & (DD-1);
  x = (B[j]<<DD)+B[j+1];
  return (x >> (DD-l)) & 0xffff;
}

static int getbit(unsigned short *B, i64 i)
{
  i64 j,l;
  i--;
  j = i / DD;
  l = i & (DD-1);
  return (B[j] >> (DD-1-l)) & 1;
}

static int setbit(unsigned short *B, i64 i,int x)
{
  i64 j,l;
  i--;
  j = i / DD;
  l = i % DD;
  if (x==0) B[j] &= (~(1<<(DD-1-l)));
  else if (x==1) B[j] |= (1<<(DD-1-l));
  else {
    printf("error setbit x=%d\n",x);
    exit(1);
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

static void writeint(int k,i64 x,FILE *f)
{
  int i;
  for (i=k-1; i>=0; i--) {
    fputc(x & 0xff,f);
    x >>= 8;
  }
}

static int encodegamma(unsigned short *B,i64 p,i64 x) /* x >= 1 */
{
i64 j,w;
  if (x==0) {
    fprintf(stderr,"encodegamma %ld\n",x);  exit(1);
  }
  w = blog(x)+1;
  for (j=0;j<w-1;j++) setbit(B,1+p+j,0);
  for (j=w-1;j>=0;j--) setbit(B,1+p+(w-1)+(w-1)-j,(x >> j)&1);
  return 2*w-1;
}

/*
static int encodedelta(unsigned short *B,i64 p,i64 x) // x >= 1 
{
i64 j,l,w;
  if (x==0) {
    fprintf(stderr,"encodedelta %ld\n",x);  exit(1);
  }
  w = blog(x)+1;
  l = encodegamma(B,p,w);
  for (j=w-2;j>=0;j--) setbit(B,1+p+l+(w-2)-j,((x >> j)&1));
  return (w-1)+l;
}
*/

static int initranktables(void)
{
  i64 i,j;
  for (i = 0; i < DD; i++) {
    for (j = (1<<i); j < (2<<i); j++) {
      R4[j] = DD-1-i;
    }
  }
  R4[0] = DD;
  
  return 0;
}


static int getzerorun(unsigned short *B,i64 p)
{
i64 w,w2;
#if 0
  w = 0;
  while (getbit(B,1+p+w)==0) w++;
#else
  w = 0;
  while (1) {
    w2 = R4[getbitD(B,1+p)];
    w += w2;
    if (w2 < DD) break;
    p += DD;
  }
#endif
  return w;
}

static int decodegamma(unsigned short *B,i64 p,i64 *ans)
{
i64 w,w2;
i64 x;
  w = getzerorun(B,p);
#if 0
  x = 1;
  for (i=0;i<w;i++) {
    x <<= 1;
    x += getbit(B,1+p+w+1+i);
  }
#else
  p += w+1;
  x = 1;
  w2 = w;
  while (w2 > DD) {
    x <<= DD;
    x += getbitD(B,1+p);
    p += DD;
    w2 -= DD;
  }
  x <<= w2;
  x += (getbitD(B,1+p)>>(DD-w2));
#endif
  *ans = x;
  return 2*w+1;
}

/*
static int decodedelta(unsigned short *B,i64 p,i64 *ans)
{
i64 l,w,w2;
i64 x;
  l = decodegamma(B,p,&w);
#if 0
  x = 1;
  for (i=1;i<w;i++) {
    x <<= 1;
    x += getbit(B,1+p+l+i-1);
  }
#else
  p += l;
  x = 1;
  w2 = w-1;
  while (w2 > DD) {
    x <<= DD;
    x += getbitD(B,1+p);
    p += DD;
    w2 -= DD;
  }
  x <<= w2;
  x += (getbitD(B,1+p)>>(DD-w2));
#endif
  *ans = x;
  return l+w-1;
}
*/

static void mkdecodetable(void)
{
  unsigned short B[256];
  i64 i,j,b,b2,d,x;

  for (i=0; i<256; i++) B[i] = 0xffff;
  for (i = 0; i < TBLSIZE; i++) {
    B[0] = i;
    b = 0;  j = 0;  x = 0;
    while (1) {
      b2 = DECODENUM(B,b,&d);
      if (b+b2 > DD) break;
      b += b2;
      x += d;
      j++;
    }
    R5n[i] = j;  R5b[i] = b;  R5x[i] = x;
  }
}

void psi1_options(CSA *csa, char *p)
{
  psi1 *ps;
  int opt;

  csa->psi_struc = ps = (psi1 *) mymalloc(sizeof(psi1));
  ps->id = csa->id;
  ps->L = 128;
  if (p[0] == 0) goto end;
  p++;
  if (p[0] != ':') {
    sscanf(p,"%d",&ps->L);
    printf("L = %d\n",ps->L);
    while (p[0] != ':' && p[0] != 0) p++;
  }
  if (p[0] == 0) goto end;
  p++;
  if (p[0] != ':') {
    if (sscanf(p,"%d",&opt)==1) {
      if (opt & 1) {
        csa->id |= ID_COMPPTR;
        ps->id |= ID_COMPPTR;
        printf("COMPPTR\n");
      }
    }
    while (p[0] != ':' && p[0] != 0) p++;
  }
  if (p[0] == 0) goto end;
  p++;
end:
  printf("psi_id = %ld L = %d\n",ps->id,ps->L);
}

static i64 psi1_psi(CSA *csa, i64 i)
{
  i64 j,k;
  i64 x;
  i64 k2,p,n;
  i64 L;
  i64 b,d,sp;
  i64 maxrun;
  unsigned short *B;
  psi1 *ps;

  ps = (psi1 *)csa->psi_struc;

//  printf("psi1_psi[%ld] (no run)\n",i);

#ifdef DEBUG
  if (i > csa->n || i < 1) {
    printf("error csa2_psi i=%u n=%u\n",i,csa->n);
    exit(1);
  }
#endif

  L = ps->L;
  if (ps->id & ID_COMPPTR) {
    x = SPARSEARRAY_select(ps->sx, (i/L)+1) % (csa->n+1);
    sp = SPARSEARRAY_select(ps->sb, (i/L)+1);
  } else {
    x = getuint(ps->R,(i / L)*2,ps->k);
    sp = getuint(ps->R,(i / L)*2+1,ps->k);
  }
  maxrun = L;

  b = 0;
  j = i % L;

  n = ps->n;
  B = ps->B;

  switch (ps->id & 0x3f) {
    case ID_DIFF_GAMMA:
    case ID_DIFF_GAMMA_SPARSE:
      k = 0;
      while (k < j) {
        p = getbitD(B+sp,1+b);
        k2 = R5n[p];
        if (k2 == 0) {
          b += DECODENUM(B+sp,b,&d);
          x += d;
          x %= (n+1);
          k++;
        } else {
          if (k+k2 > j) break;
          k += k2;
          b += R5b[p];
          x += R5x[p];
          x %= (n+1);
        }
      }

      for (; k<j; k++) {
        b += DECODENUM(B+sp,b,&d);
        x += d;
        x %= (n+1);
      }
      break;
    case ID_DIFF_GAMMA_RL:
    case ID_DIFF_GAMMA_RL_SPARSE:
//      psi1_decbuf[0] = x;

      b = 0;
      for (k=1; k<=j; k++) {
        b += DECODENUM(B+sp,b,&d);
        if (d <= maxrun*2) {
          if (d % 2 == 0) {
#if 0
            for (i=0; i<d/2; i++) {
              x += 1;
              x %= (n+1);
              if (k+i >= L) {
                printf("psi1_psi: error k=%d i=%d l=%d\n",k,i,L);
              }
//              psi1_decbuf[k+i] = x;
                if (k+i == j) break;
            }
#else
            if (k+d/2-1 >= j) {
              x += j-k+1;
              x %= (n+1);
              break;
            }
            x += d/2;
#endif
            k += (d/2)-1;
          } else {
            x += (d+3)/2;
            x %= (n+1);
//            psi1_decbuf[k] = x;
          }
        } else {
          x += d-maxrun+1;
          x %= (n+1);
//          psi1_decbuf[k] = x;
        }
      }
//      x = psi1_decbuf[j];
      break;
  }

#ifdef DEBUG
  if (x < 0 || x > csa->n) {
    printf("error csa2_psi(%u) %u\n",i,x);
  }
#endif
//  printf("psi1_psi: psi[%ld] = %ld\n",i,x);
  return x;
}

static i64 psi12_psi(CSA *csa, i64 i)
{
  i64 j,k;
  i64 x;
  i64 n;
  i64 L;
  i64 b,d,sp;
  unsigned short *B;
  psi1 *ps;
  int runlen;

  ps = (psi1 *)csa->psi_struc;

  L = ps->L;
  if (ps->id & ID_COMPPTR) {
    x = SPARSEARRAY_select(ps->sx, (i/L)+1) % (csa->n+1);
    sp = SPARSEARRAY_select(ps->sb, (i/L)+1);
  } else {
    x = getuint(ps->R,(i / L)*2,ps->k);
    sp = getuint(ps->R,(i / L)*2+1,ps->k);
  }

  b = 0;
  j = i % L;

  n = ps->n;
  B = ps->B;

  b = 0;
  d = getbit(B+sp,b+1);
  if (d == 1) runlen = 1;  else runlen = 0;
  b++;
  k = 1;
  while (k<=j) {
    if (runlen == 1) {
      b += DECODENUM(B+sp,b,&d);
//      printf("1. k = %d j = %d x = %ld d = %d runlen = %d\n",k,j,x,d,runlen);
      if (k+d-1 > j) {
        x += j-k+1;
        x %= (n+1);
        break;
      }
      x += d-1;
      k += d-1;
//      printf("11. k = %d j = %d x = %ld d = %d runlen = %d\n",k,j,x,d,runlen);
    }
    if (k > j) break;

    b += DECODENUM(B+sp,b,&d);
//    printf("2. k = %d j = %d x = %ld d = %d runlen = %d\n",k,j,x,d,runlen);
    x += d+1;
    x %= (n+1);
    k++;
//    printf("22. k = %d j = %d x = %ld d = %d runlen = %d\n",k,j,x,d,runlen);
    runlen = 1;
  }

  return x;
}


static i64 psi1_pred(CSA *csa, i64 pr, i64 l, i64 r)
{
  i64 m;
  i64 x;
  i64 sp,L,d,n;
  i64 lb, rb, b;
  uchar *R;
  unsigned short *B, *buf, offset;
  i64 maxrun, run;
  int k;
  psi1 *ps;
  i64 ans;

  ps = (psi1 *)csa->psi_struc;


  R = ps->R;  B = ps->B;  L = ps->L;  k = ps->k;
  n = ps->n;
  maxrun = L;

  lb = (l+L-1) / L;
  rb = r / L;

  b = lb-1;
  while (lb <= rb) {
    m = (lb+rb) / 2;
    if (ps->id & ID_COMPPTR) {
      x = SPARSEARRAY_select(ps->sx, m+1) % (csa->n+1);
    } else {
      x = getuint(R,m*2,k);
    }
    if (x == pr) {
      b = m;
      break;
    }
    if (x < pr) {
      b = m;
      lb = m+1;
    } else {
      rb = m-1;
    }
  }
  
// ブロック b ([b*L..(b+1)*L-1]) の中で [l..r] の範囲を探す
  if (ps->id & ID_COMPPTR) {
    x = SPARSEARRAY_select(ps->sx, b+1) % (csa->n+1);
    sp = SPARSEARRAY_select(ps->sb, b+1);
  } else {
    x = getuint(R,b*2,k); // ブロックの先頭の値
    sp = getuint(R,b*2+1,k);
  }
  ans = l-1;

  if (r > (b+1)*L-1) r = (b+1)*L-1; // ブロックの右端


  m = b*L;
  buf = B+sp;
  offset = 0;
  switch (ps->id & 0x3f) {
    case ID_DIFF_GAMMA:
      while (m < l) { // l より左はスキップ
        offset += DECODENUM(buf,offset,&d);
        x += d;
        x %= (n+1);
        m++;
      }
      while (1) {
        if (x > pr) break;
        ans = m;
        // 次のPsiを計算
        m++;
        if (m > r) break;
        offset += DECODENUM(buf,offset,&d);
        x += d;
        x %= (n+1);
      }
      break;
    case ID_DIFF_GAMMA_RL:
      run = 0;
      while (1) {
        if (m >= l && x > pr) break;
        ans = m;
        // 次のPsiを計算
        m++;
        if (m > r) break;
        if (run > 0) {
          x += 1;
          x %= (n+1);
          run--;
        } else {
          offset += DECODENUM(buf,offset,&d);
          if (d <= maxrun*2) {
            if (d % 2 == 0) {
              run = d/2-1;
              x += 1;
              x %= (n+1);
            } else {
              x += (d+3)/2;
              x %= (n+1);
            }
          } else {
            x += d-maxrun+1;
            x %= (n+1);
          }
        }
      }
      break;
    case ID_DIFF_GAMMA_RR:
      d = getbit(B+sp,offset+1);
      offset++;
      if (d == 1) {
//        printf("runlen = 1\n");
        offset += DECODENUM(B+sp,offset,&d);
//        printf("d = %ld\n",d);
        run = d-1;
      } else {
//        printf("runlen = 0\n");
        run = 0;
      }
      while (1) {
        if (m >= l && x > pr) break;
        ans = m;
        // 次のPsiを計算
        m++;
        if (m > r) break;
        if (run > 0) {
          x += 1;
          x %= (n+1);
          run--;
        } else {
          offset += DECODENUM(B+sp,offset,&d); // 0の数
//          printf("zero = %ld\n",d);
          x += d+1;
          x %= (n+1);
          offset += DECODENUM(B+sp,offset,&d);
//          printf("run = %ld\n",d);
          run = d-1;
        }
      }
      break;
    default:
      printf("psi id %lu is not supported.\n", ps->id & 0x3f);
      exit(1);
  }
  return ans;
}

static i64 psi1_succ(CSA *csa, i64 pl, i64 l, i64 r)
{
  i64 m;
  i64 x;
  i64 sp,L,d,n;
  uchar *R;
  unsigned short *B, *buf, offset;
  i64 maxrun, run;
  int k;
  psi1 *ps;
  i64 lb, rb, b;
  i64 ans;

  ps = (psi1 *)csa->psi_struc;

  R = ps->R;  B = ps->B;  L = ps->L;  k = ps->k;
  n = ps->n;
  maxrun = L;

  lb = (l+L-1) / L;
  rb = r / L;

  b = rb+1;
  while (lb <= rb) {
    m = (lb+rb) / 2;
    if (ps->id & ID_COMPPTR) {
      x = SPARSEARRAY_select(ps->sx, m+1) % (csa->n+1);
    } else {
      x = getuint(R,m*2,k);
    }
    if (x == pl) {
      b = m;
      break;
    }
    if (x > pl) {
      b = m;
      rb = m-1;
    } else {
      lb = m+1;
    }
  }


// ブロック b-1 ([(b-1)*L..b*L-1]) の中で [l..r] の範囲を探す
  if (ps->id & ID_COMPPTR) {
    x = SPARSEARRAY_select(ps->sx, (b-1)+1) % (csa->n+1);
    sp = SPARSEARRAY_select(ps->sb, (b-1)+1);
  } else {
    x = getuint(R,(b-1)*2,k); // ブロックの先頭の値
    sp = getuint(R,(b-1)*2+1,k);
  }

  if (r > b*L-1) r = b*L-1; // ブロックの右端
  ans = r+1; // 探しても見つからないときはブロックの右端+1


  m = (b-1)*L;


  buf = B + sp;
  offset = 0;

  switch (ps->id & 0x3f) {
    case ID_DIFF_GAMMA:
      while (m < l) { // l より左はスキップ
        offset += DECODENUM(buf,offset,&d);
        x += d;
        x %= (n+1);
        m++;
      }
      while (1) {
        if (x >= pl) {
          ans = m;
          break;
        }
        // 次のPsiを計算
        m++;
        if (m > r) break;
        offset += DECODENUM(buf,offset,&d);
        x += d;
        x %= (n+1);
      }
      break;
    case ID_DIFF_GAMMA_RL:
      run = 0;
      while (1) {
        if (m >= l && x >= pl) {
          ans = m;
          break;
        }
        // 次のPsiを計算
        m++;
        if (m > r) break;
        if (run > 0) {
          x += 1;
          x %= (n+1);
          run--;
        } else {
          offset += DECODENUM(buf,offset,&d);
//          printf("d = %ld\n",d);
          if (d <= maxrun*2) {
            if (d % 2 == 0) {
              run = d/2-1;
              x += 1;
              x %= (n+1);
            } else {
              x += (d+3)/2;
              x %= (n+1);
            }
          } else {
            x += d-maxrun+1;
            x %= (n+1);
          }
        }
      }
      break;
    case ID_DIFF_GAMMA_RR:
      d = getbit(B+sp,offset+1);
      offset++;
      if (d == 1) {
//        printf("runlen = 1\n");
        offset += DECODENUM(B+sp,offset,&d);
//        printf("d = %ld\n",d);
        run = d-1;
      } else {
//        printf("runlen = 0\n");
        run = 0;
      }
      while (1) {
        if (m >= l && x >= pl) {
          ans = m;
          break;
        }
        // 次のPsiを計算
        m++;
        if (m > r) break;
        if (run > 0) {
          x += 1;
          x %= (n+1);
          run--;
        } else {
          offset += DECODENUM(B+sp,offset,&d); // 0の数
//          printf("zero = %ld\n",d);
          x += d+1;
          x %= (n+1);
          offset += DECODENUM(B+sp,offset,&d);
//          printf("run = %ld\n",d);
          run = d-1;
        }
      }
      break;
    default:
      printf("psi id %lu is not supported.\n", ps->id & 0x3f);
      exit(1);
      break;
  }    
  return ans;
}

#if 0
static i64 psi1_succ_tmp(CSA *csa, i64 pr, i64 l, i64 r)
{
  i64 p1,p2;
  p2 = csa_psi_succ_naive(csa, pr, l, r);
  p1 = psi1_succ(csa, pr, l, r);
  if (p1 != p2) {
    printf("pr = %ld [%ld,%ld] psi_succ_naive = %ld\n", pr, l, r, p2);
    printf("pr = %ld [%ld,%ld] psi1_succ = %ld\n", pr, l, r, p1);
//    csa_psi_succ_naive(csa, pr, l, r);
//    psi1_succ(csa, pr, l, r);
  }
  return p2;
}

static i64 psi1_pred_tmp(CSA *csa, i64 pr, i64 l, i64 r)
{
  i64 p1,p2;
  p2 = csa_psi_pred_naive(csa, pr, l, r);
  p1 = psi1_pred(csa, pr, l, r);
  if (p1 != p2) {
    printf("pr = %ld [%ld,%ld] psi_pred_naive = %ld\n", pr, l, r, p2);
    printf("pr = %ld [%ld,%ld] psi1_pred = %ld\n", pr, l, r, p1);
  }
  return p2;
}
#endif

static unsigned short Btmp[10240];
i64 psi1_makeindex(CSA *csa, char *fname)
{
i64 psize,psize1,psize2;
i64 b, b2;
i64 i,j,x,xx;
i64 y,d,w;
int k;
FILE *f1,*f2;
char *fpsi, *fpsd;
//psi1_iterator *pi;
i64 runlen;
i64 maxrun;
i64 n,L;
psi1 *ps;
int id,id2;
FILE *out;
diskbuf *psi;
char *fbw, *flst;
SPARSEARRAY sx, sb;
int mm;

  ps = (psi1 *)csa->psi_struc;
  id = ps->id;
  id2 = id & 0x3f;

  k = strlen(fname);
  fbw = (char *) mymalloc(k+5);
  flst = (char *) mymalloc(k+5);
  fpsi = (char *) mymalloc(k+5);
  fpsd = (char *) mymalloc(k+5);
  sprintf(fbw,"%s.bw",fname);
  sprintf(flst,"%s.lst",fname);

  switch (id2) {
  case ID_DIFF_GAMMA:
    sprintf(fpsi,"%s.psi",fname);
    sprintf(fpsd,"%s.psd",fname);
    break;
  case ID_DIFF_GAMMA_RL:
    sprintf(fpsi,"%s.pri",fname);
    sprintf(fpsd,"%s.prd",fname);
    break;
  case ID_DIFF_GAMMA_SPARSE:
    sprintf(fpsi,"%s.psi",fname);
    sprintf(fpsd,"%s.pss",fname);
    break;
  case ID_DIFF_GAMMA_RL_SPARSE:
    sprintf(fpsi,"%s.pri",fname);
    sprintf(fpsd,"%s.prs",fname);
    break;
  }

  out = create_tmp(0);
  bw_to_psi(out,csa,fbw,flst,&k);
  psi = open_diskbuf(out,k);
  ps->last = getint_diskbuf(psi,0);
  printf("last = %ld\n",ps->last);

  n = csa->n;
  L = ps->L;
  if (L >= n) {
    printf("L=%ld >= n=%ld\n",L,n);
    exit(0);
  }


  maxrun = L;


  mkdecodetable();

  f1 = fopen(fpsi,"wb");
  psize1 = 0;

  f2 = fopen(fpsd,"wb");
  psize2 = 0;

  ps->k = k = (blog(n+1)+1+8-1)/8;

//  pi = psi1_iterator_new(ps,0);


  writeint(1,ID_PSI,f2);
  writeint(1,k,f2); /* #bytes of integer */
  writeint(k,n,f2);
  writeint(k,L,f2);
  psize2 += 1+1+2*k;

//  writeint(1,ID_DIFF_GAMMA_RL,f2);
  writeint(1,id,f2);
  psize2 += 1;

  if (id & ID_COMPPTR) {
    mm = 0;
    for (i=0; i<SIGMA; i++) {
      if (csa->C[i] > 0) mm++;
    }
    SPARSEARRAY_construct_init(&sx, (mm+1)*(n+1), n/L+1);
    SPARSEARRAY_construct_init(&sb, n, n/L+1);
  }

  b = b2 = 0;
  mm = 0;  xx = 0;
  for (j=0; j<=n/L; j++) {
//    display_progressbar("writing psi ",j,n/L);
    if (j % 100000 == 0) {
      printf("%ld %1.3f bpc\r",j,(double)psize2*8/(j+1)/L);  fflush(stdout);
    }
//    x = psi1_iterator_get(pi,j*L);
    y = getint_diskbuf(psi,j*L);
    if (id & ID_COMPPTR) {
      if (y <= xx) {
        mm++;
      }
      SPARSEARRAY_construct_set(&sx, j, mm*(n+1) + y);
      SPARSEARRAY_construct_set(&sb, j, b);
      xx = y;
    } else {
//      printf("%ld   x=%ld   sp=%ld\n",j,y,b);
      writeint(k,y,f2);
      writeint(k,b,f2);
      psize2 += 2*k;
    }
    x = y;
    runlen = 0;
    b2 = 0;
    for (i=j*L+1; i<(j+1)*L && i <= n; i++) { /* psi[j*L] are not encoded */
      //y = psi1_iterator_get(pi,i);
      y = getint_diskbuf(psi,i);
      d = y - x;
      if (d <= 0) {
        d += n+1;
      }
      if (id2 == ID_DIFF_GAMMA_RL || id2 == ID_DIFF_GAMMA_RL_SPARSE) {
        if (d > 1) {
          if (runlen>0) {
            w = ENCODENUM(Btmp,b2,runlen*2);
            b2 += w;
            runlen = 0;
          }
          if (d >= maxrun+2) {
            w = ENCODENUM(Btmp,b2,d+maxrun-1);
          } else {
            w = ENCODENUM(Btmp,b2,d*2-3);
          }
          b2 += w;
        } else {
          runlen++;
          if (runlen == maxrun) {
            w = ENCODENUM(Btmp,b2,runlen*2);
            b2 += w;
            runlen = 0;
          }
        }
      } else {
        w = ENCODENUM(Btmp,b2,d);
        b2 += w;
      }
      x = y;
    }
    if (id2 == ID_DIFF_GAMMA_RL || id2== ID_DIFF_GAMMA_RL_SPARSE) {
      if (runlen>0) {
        w = ENCODENUM(Btmp,b2,runlen*2);
        b2 += w;
        runlen = 0;
      }
    }
    fwrite(Btmp,(b2+15) / 16,sizeof(short),f1);
    psize1 += (b2+15)/16*sizeof(short);
    b += (b2+15) / 16;
    b2 = 0;
  }
  if (b2 > 0) {
    fwrite(Btmp,(b2+15) / 16,sizeof(short),f1);
    psize1 += (b2+15)/16*sizeof(short);
  }
  fwrite(Btmp,1,sizeof(short),f1); // getbitDで１ワード余計に読むため
  psize1 += 1*sizeof(short);

  if (id & ID_COMPPTR) {
    SPARSEARRAY_construct_end(&sx, SDARRAY_SELECT1);
    SPARSEARRAY_construct_end(&sb, SDARRAY_SELECT1);
    SPARSEARRAY_write(&sx, f2);
    SPARSEARRAY_write(&sb, f2);
  }

  psize = psize1 + psize2;
  printf("size %ld (%1.3f bpc)\n",psize,(double)psize*8 / n);

  fclose(f1);
  fclose(f2);

//  psi1_iterator_remove(pi);

  close_diskbuf(psi);
  fclose(out);
  remove_tmp(0);

  psi1_read(csa, fpsd);

  free(fpsi);
  free(fpsd);
  free(fbw);
  free(flst);

  return psize;
}

i64 psi12_makeindex(CSA *csa, char *fname)
{
i64 psize,psize1,psize2;
i64 b, b2;
i64 i,j,x,xx;
i64 y,d,w;
int k;
FILE *f1,*f2;
char *fpsi, *fpsd;
//psi1_iterator *pi;
i64 runlen;
i64 maxrun;
i64 n,L;
psi1 *ps;
int id,id2;
FILE *out;
diskbuf *psi;
char *fbw, *flst;
SPARSEARRAY sx, sb;
int mm;

  ps = (psi1 *)csa->psi_struc;
  id = ps->id;
  id2 = id & 0x3f;

  k = strlen(fname);
  fbw = (char *) mymalloc(k+5);
  flst = (char *) mymalloc(k+5);
  fpsi = (char *) mymalloc(k+5);
  fpsd = (char *) mymalloc(k+5);
  sprintf(fbw,"%s.bw",fname);
  sprintf(flst,"%s.lst",fname);

  switch (id2) {
  case ID_DIFF_GAMMA_RR:
    sprintf(fpsi,"%s.pxi",fname);
    sprintf(fpsd,"%s.pxd",fname);
    break;
  }

  out = create_tmp(0);
  bw_to_psi(out,csa,fbw,flst,&k);
  psi = open_diskbuf(out,k);
  ps->last = getint_diskbuf(psi,0);
  printf("last = %ld\n",ps->last);

  n = csa->n;
  L = ps->L;
  if (L >= n) {
    printf("L=%ld >= n=%ld\n",L,n);
    exit(0);
  }


  maxrun = L;


  mkdecodetable();

  f1 = fopen(fpsi,"wb");
  psize1 = 0;

  f2 = fopen(fpsd,"wb");
  psize2 = 0;

  ps->k = k = (blog(n+1)+1+8-1)/8;

//  pi = psi1_iterator_new(ps,0);


  writeint(1,ID_PSI,f2);
  writeint(1,k,f2); /* #bytes of integer */
  writeint(k,n,f2);
  writeint(k,L,f2);
  psize2 += 1+1+2*k;

  writeint(1,id,f2);
  psize2 += 1;

  if (id & ID_COMPPTR) {
    mm = 0;
    for (i=0; i<SIGMA; i++) {
      if (csa->C[i] > 0) mm++;
    }
    SPARSEARRAY_construct_init(&sx, (mm+1)*(n+1), n/L+1);
    SPARSEARRAY_construct_init(&sb, n, n/L+1);
  }



  b = b2 = 0;
  mm = 0;  xx = 0;
  for (j=0; j<=n/L; j++) {
//    display_progressbar("writing psi ",j,n/L);
    if (j % 100000 == 0) {
      printf("%ld %1.3f bpc\r",j,(double)psize2*8/(j+1)/L);  fflush(stdout);
    }
    y = getint_diskbuf(psi,j*L);

    if (id & ID_COMPPTR) {
      if (y <= xx) {
        mm++;
      }
      SPARSEARRAY_construct_set(&sx, j, mm*(n+1) + y);
      SPARSEARRAY_construct_set(&sb, j, b);
      xx = y;
    } else {
//      printf("%ld   x=%ld   sp=%ld\n",j,y,b);
      writeint(k,y,f2);
      writeint(k,b,f2);
      psize2 += 2*k;
    }

    x = y;
    runlen = 0;
    b2 = 0;
    for (i=j*L+1; i<(j+1)*L && i <= n; i++) { /* psi[j*L] are not encoded */
      y = getint_diskbuf(psi,i);
      d = y - x;
      if (d <= 0) {
        d += n+1;
      }
      if (i == j*L+1) {
        if (d == 1) {setbit(Btmp,b2+1,1);  runlen = 1;}
        else {setbit(Btmp,b2+1,0);  runlen = 0;}
        b2++;
      }
      if (d > 1) {
        if (runlen>0) {
//          printf("d=%ld encode runlen=%ld\n",d,runlen);
          w = ENCODENUM(Btmp,b2,runlen);
          b2 += w;
          runlen = 0;
        }
//        printf("encode d-1=%ld\n",d-1);
        w = ENCODENUM(Btmp,b2,d-1);
        b2 += w;
        runlen = 1;
      } else {
        runlen++;
      }
      x = y;
    }
    if (runlen>0) {
//      printf("encode runlen=%ld\n",runlen);
      w = ENCODENUM(Btmp,b2,runlen);
      b2 += w;
      runlen = 0;
    }
    fwrite(Btmp,(b2+15) / 16,sizeof(short),f1);
    psize1 += (b2+15)/16*sizeof(short);
    b += (b2+15) / 16;
    b2 = 0;
  }
  if (b2 > 0) {
    fwrite(Btmp,(b2+15) / 16,sizeof(short),f1);
    psize1 += (b2+15)/16*sizeof(short);
  }
  fwrite(Btmp,1,sizeof(short),f1); // getbitDで１ワード余計に読むため
  psize1 += 1*sizeof(short);

  if (id & ID_COMPPTR) {
    SPARSEARRAY_construct_end(&sx, SDARRAY_SELECT1);
    SPARSEARRAY_construct_end(&sb, SDARRAY_SELECT1);
    SPARSEARRAY_write(&sx, f2);
    SPARSEARRAY_write(&sb, f2);
  }

  psize = psize1 + psize2;
  printf("size %ld (%1.3f bpc)\n",psize,(double)psize*8 / n);

  fclose(f1);
  fclose(f2);

//  psi1_iterator_remove(pi);

  close_diskbuf(psi);
  fclose(out);
  remove_tmp(0);

  psi1_read(csa, fpsd);

  free(fpsi);
  free(fpsd);
  free(fbw);
  free(flst);

  return psize;
}



i64 psi1_read(CSA *csa, char *fname)
{

  i64 psize1,psize2;
  i64 n;
  int k,l,id,id2;
  char *fpsi, *fpsd, *fname2;
  psi1 *ps;
  uchar *p,*q;
  
  csa->psi_struc = ps = (psi1 *) mymalloc(sizeof(psi1));

  k = strlen(fname);
  fname2 = (char *) mymalloc(k-4+1);
  strncpy(fname2,fname,k-4);
  fname2[k-4] = 0;
  k -= 5;

  initranktables();
  mkdecodetable();

  fpsi = (char *) mymalloc(k+5+1);
//  fpsd = mymalloc(k+5);

//  sprintf(fpsd,"%s.psd",fname2);
  fpsd = fname;
//  printf("psi_read: read %s\n",fpsd);

  ps->mappsd = mymmap(fpsd);
  if (ps->mappsd->addr==NULL) {
    perror("psi1_read: mmap2\n");
    exit(1);
  }
  p = q = (uchar *)ps->mappsd->addr;
  psize1 = ps->mappsd->len;

  id = getuint(p,0,1);  p += 1;
  if (id != ID_PSI) {
    printf("read_psi: id = %d is not supported.\n",id);
    exit(1);
  }
  ps->k = k = getuint(p,0,1);  p += 1;
  ps->n = n = getuint(p,0,k);  p += k;
  ps->L = l = getuint(p,0,k);  p += k;

  id = getuint(p,0,1);  p += 1;
//  printf("read_psi: psi_id = %d L = %d\n",id,l);
  csa->id = ps->id = id;
  id2 = id & 0x3f;
  switch (id2) {
    case ID_DIFF_GAMMA:
      printf("#psi format = GAMMA L=%d C=%d\n",l,(id>>7));
      sprintf(fpsi,"%s.psi",fname2);
      break;
    case ID_DIFF_GAMMA_RL:
      printf("#psi format = GAMMA_RL L=%d C=%d\n",l,(id>>7));
      sprintf(fpsi,"%s.pri",fname2);
      break;
    case ID_DIFF_GAMMA_SPARSE:
      printf("#psi format = GAMMA_SPARSE L=%d C=%d\n",l,(id>>7));
      sprintf(fpsi,"%s.psi",fname2);
      break;
    case ID_DIFF_GAMMA_RL_SPARSE:
      printf("#psi format = GAMMA_RL_SPARSE L=%d C=%d\n",l,(id>>7));
      sprintf(fpsi,"%s.pri",fname2);
      break;
    case ID_DIFF_GAMMA_RR:
      printf("#psi format = GAMMA_RR L=%d C=%d\n",l,(id>>7));
      sprintf(fpsi,"%s.pxi",fname2);
      break;
    default:
      printf("read_csa: ID %d is not supported.\n",id);
      break;
  }

  if (id & ID_COMPPTR) {
    printf("COMPPTR\n");
    ps->sx = (sparsearray *) mymalloc(sizeof(SPARSEARRAY));
    ps->sb = (sparsearray *) mymalloc(sizeof(SPARSEARRAY));
    SPARSEARRAY_read(ps->sx, &p);
    SPARSEARRAY_read(ps->sb, &p);
  } else {
    ps->R = p;
  }
//  printf("psize = %ld\n",psize);

////   read psi
//  printf("psi_read: map %s\n",fpsi);
  ps->mappsi = mymmap(fpsi);
  if (ps->mappsi->addr==NULL) {
    perror("psi1_read: mmap1\n");
    exit(1);
  }
  ps->B = (unsigned short *)ps->mappsi->addr;
  psize2 = ps->mappsi->len;
//  printf("psize2 = %ld\n",psize2);

//  printf("psi1_read: psize1 = %ld psize2 = %ld\n",psize1,psize2);
  ps->psize = psize1 + psize2;

  free(fpsi);
//  free(fpsd);
  free(fname2);

// user-specific functions
  csa->psi = psi1_psi;
  if (id2 == ID_DIFF_GAMMA_RR) {
    csa->psi = psi12_psi;
//    csa->psi_pred = csa_psi_pred_naive;
//    csa->psi_succ = csa_psi_succ_naive;
    csa->psi_pred = psi1_pred;
    csa->psi_succ = psi1_succ;
  } else {
    csa->psi_succ = psi1_succ;
    csa->psi_pred = psi1_pred;
//      csa->psi_succ = psi1_succ_tmp;
//      csa->psi_pred = psi1_pred_tmp;
//    csa->psi_succ = csa_psi_succ_naive;
//    csa->psi_pred = csa_psi_pred_naive;
  }

// default functions
  csa->LF = csa_LF_by_psi;
  csa->lookup = csa_lookup;
  csa->inverse = csa_inverse;
  csa->text = csa_text;
  csa->substring = csa_substring;
  csa->T = csa_T;
  csa->head = csa_head_rank;
  csa->search = csa_search;
  csa->searchsub = csa_searchsub;
//  csa->child_l = csa_child_l;
//  csa->child_r = csa_child_r;


  return psize1 + psize2;
}

void psi1_release(psi1 *ps)
{
  mymunmap(ps->mappsi);
  mymunmap(ps->mappsd);
}


psi1_iterator *psi1_iterator_new(psi1 *ps, i64 start)
{
  i64 L;
  psi1_iterator *pi;

  L = ps->L;
  pi = (psi1_iterator *) mymalloc(sizeof(*pi));
  pi->buf = (i64 *) mymalloc(L*sizeof(i64));

  pi->ps = ps;
  pi->n = ps->n + 1;
  pi->page = -1;
  pi->pos = start;

  return pi;
}

static void psi1_iterator_readpage(psi1_iterator *pi, i64 page)
{
  i64 L,i,j,k,id;
  i64 x,sp,n,b,d;
  unsigned short *B;
  i64 maxrun;
  psi1 *ps;

  n = pi->n;
  ps = pi->ps;
  L = ps->L;
  id = ps->id;
  maxrun = L;

  B = ps->B;
  j = L;
  if (page*L + j > n) j = n - page*L;
  x = getuint(ps->R,page*2,ps->k);
  sp = getuint(ps->R,page*2+1,ps->k);

  pi->buf[0] = x;

  b = 0;
  for (k=1; k<j; k++) {
    b += DECODENUM(B+sp,b,&d);
    if (id == ID_DIFF_GAMMA) {
      x += d;
      x %= n;
      pi->buf[k] = x;
    } else if (id == ID_DIFF_GAMMA_RL) {
      if (d <= maxrun*2) {
        if (d % 2 == 0) {
          for (i=0; i<d/2; i++) {
            x += 1;
            x %= n;
            if (k+i >= L) {
              printf("readpage: error k=%ld i=%ld l=%ld\n",k,i,L);
            }
            pi->buf[k+i] = x;
          }
          k += (d/2)-1;
        } else {
          x += (d+3)/2;
          x %= n;
          pi->buf[k] = x;
        }
      } else {
        x += d-maxrun+1;
        x %= n;
        pi->buf[k] = x;
      }
    } else {
      printf("??? id = %ld\n",id);
    }
  }
  pi->page = page;
}

i64 psi1_iterator_next(psi1_iterator *pi)
{
  i64 page,r,L,i,x,n;
  psi1 *ps;

  n = pi->n;
  i = pi->pos;
  if (i >= n) {
    printf("psi1_iterator_next: i=%ld n=%ld\n",i,n);
    exit(1);
  }

  ps = pi->ps;

  L = ps->L;
  page = i / L;
  r = i % L;

  if (page != pi->page) {
    psi1_iterator_readpage(pi,page);
  }

  x = pi->buf[r];
  pi->pos++;

  return x;
}

i64 psi1_iterator_hasnext(psi1_iterator *pi)
{
  if (pi->pos < pi->n) return 1; else return 0;
}

i64 psi1_iterator_get(psi1_iterator *pi, i64 i)
{
  i64 page,r,L,x,n;
  psi1 *ps;

  n = pi->n;
  if (i >= n) {
    printf("psi1_iterator_get: i=%ld n=%ld\n",i,n);
    exit(1);
  }

  ps = pi->ps;

  L = ps->L;
  page = i / L;
  r = i % L;

  if (page != pi->page) {
    psi1_iterator_readpage(pi,page);
  }

  x = pi->buf[r];
  pi->pos++;

  return x;
}

void psi1_iterator_remove(psi1_iterator *pi)
{
  free(pi->buf);
  free(pi);
}
