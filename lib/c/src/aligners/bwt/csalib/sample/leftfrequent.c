/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdlib.h>
#include <stdio.h>
#include "csa.h"

int utf_head(int c)
{
  return (0 <= c && c <= 0x7f);
}

int utf8_len(int c)
{
  int k;
  if (0 <= c && c <= 0x7f) k = 1;
  else if (0xc2 <= c && c <= 0xdf) k = 2;
  else if (0xe0 <= c && c <= 0xef) k = 3;
  else if (0xf0 <= c && c <= 0xf7) k = 4;
  else if (0xf8 <= c && c <= 0xfb) k = 5;
  else if (0xfc <= c && c <= 0xfd) k = 6;
  else k = -1;
  return k;
}

int valid_utf8(uchar *p, int len)
{
  int c, k, i;

  while (len > 0) {
    c = *p++;
    k = utf8_len(c);
    if (k <= 0) return 0;
    if (k > len) return 0;

    for (i=1; i<k; i++) {
      c = *p++;
      if (!(0x80 <= c && c <= 0xbf)) return 0;
    }
    len -= k;
  }
  return 1;
}

////////////////////////////////////////////////////////////////
// ���ɑ�ȕp�o�p�^����񋓂���
// �p�^���̒����� maxlen �ȉ�
// �e�L�X�g��Unicode (utf8) �Ɖ��肵�C�����̐؂�ڂ����������̂������߂�
////////////////////////////////////////////////////////////////
void frequent_leftmaximal(CSA *csa, i64 l, i64 r, i64 th, uchar *key, int len, int maxlen)
{
  int i,c,k;
  int maximal;
  i64 ll, rr, num;
  uchar head[SIGMA];
  i64 L[SIGMA], R[SIGMA];

  len++;  if (len > maxlen) return; // �p�^�����̍ő�l�𒴂�����I��

  k = csa->child_l(csa,l,r,head, L, R); // �p�^�� P �̍����ɏo�����镶���W�������߂�

  maximal = 1;
  for (i=0; i<k; i++) {
    c = head[i];
    ll = L[i];  rr = R[i];
    if (rr - ll + 1 >= th) { // �p�^�� cP �̕p�x��臒l�ȏ�Ȃ�ΒT���𑱂���
      *(key-1) = c;
      frequent_leftmaximal(csa, ll, rr, th, key-1, len, maxlen);
      maximal = 0; // cP ���p�o�Ȃ̂ŁCP �͋ɑ�ł͂Ȃ�
    }
  }

  if (maximal && valid_utf8(key,len)) { // P ���ɑ�Ő�����Unicode������Ȃ�Ε\��
    printf("[%ld,%ld] key = %s\n",l,r,key);
  }
}

#define MAXLEN 100
int main(int argc, char *argv[])
{
  i64 n;
  CSA SA;
  double t;
  i64 th;
  uchar key[MAXLEN+1];

  if (argc<3) {
    fprintf(stderr, "syntax: %s threshold {indexfiles}\n", argv[0]);
    return 1;
  }

  t = atof(argv[1]);

  csa_read(&SA,argc-2, argv+2);
  n = SA.n;
  if (t < 1.0) th = n * t; else th = t;
  printf("n = %ld threshold = %ld\n", n, th);

  key[MAXLEN] = 0;
  frequent_leftmaximal(&SA, 0, n, th, key+MAXLEN, 0, MAXLEN);

  return 0;
}
