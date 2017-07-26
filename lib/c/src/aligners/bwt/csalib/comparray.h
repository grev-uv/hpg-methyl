/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _COMPARRAY_H_
#define _COMPARRAY_H_

#include "typedef.h"
#include "sparsearray.h"
#include "huffman.h"

#ifdef __SSE4_2__
#include <smmintrin.h>
#define POPCOUNT(x) _mm_popcnt_u64(x)
#else
#define POPCOUNT(x) popcount(x)
#endif

//#define SDARRAY_RANK1 (1<<0)
//#define SDARRAY_RANK0 (1<<1)
//#define SDARRAY_SELECT1 (1<<2)
//#define SDARRAY_SELECT0 (1<<3)
//#define SDARRAY_LENHUF (1<<4)

#ifndef _BITVEC_T_
#define _BITVEC_T_
typedef u64 bitvec_t;
#endif


typedef struct {
  i64 n; // �x�N�g���̒���
  i64 m; // 1�̐�
  i64 size; // �����T�C�Y (�x�N�g���܂�)
  bitvec_t *buf; // �x�N�g��
  i64 bufsize; // �x�N�g���̃T�C�Y (bitvec_t�̐�)
  uchar opt; // �I�v�V���� (1�̐����n�t�}�������ŕ��������邩)
  ushort L;

  uchar k1; // pl �̃T�C�Y
  uchar k2; // rl �̃T�C�Y

  Huffman *h;

// pointer to buf
  sparsearray *psa;
  uchar *pl; // buf�̃I�t�Z�b�g�D�r�b�g�P�ʁD
  word *ps; // large block ���ł� buf �̃I�t�Z�b�g�D�r�b�g�P��

//for rank
  sparsearray *rsa;
  uchar *rl;
  word *rs;

} comparray;


void comparray_construct(comparray *da, i64 n, bitvec_t *buf, ushort L, int opt);
int comparray_construct_init(comparray *da, i64 n);
int comparray_construct_set(comparray *da, i64 i, int x);
i64 comparray_construct_end(comparray *da, ushort L, int opt);
void comparray_maketbl(void);

i64 comparray_write(comparray *da, FILE *f);
void comparray_read(comparray *da, uchar **map);
int comparray_getbit(comparray *da, i64 i);
i64 comparray_rank(comparray *da, i64 i);
i64 comparray_rank0(comparray *da, i64 i);
i64 comparray_rank_and_bit(comparray *da, i64 i, int *c);

typedef struct {
  i64 n; // �x�N�g���̒���
  i64 m; // 1�̐�
  i64 bufsize; // �x�N�g���̃T�C�Y (bitvec_t�̐�)

  uchar opt;
  ushort L;
  Huffman *h;
  bitvec_t *buf; // �x�N�g��
// pointer to buf
  word *ps; // large block ���ł� buf �̃I�t�Z�b�g�D�r�b�g�P��
//for rank
  word *rs;

} comparray_sb;


void comparray_sb_construct(comparray_sb *da, i64 n, bitvec_t *buf, ushort L, int opt);
int comparray_sb_construct_init(comparray_sb *da, i64 n);
int comparray_sb_construct_set(comparray_sb *da, i64 i, int x);
i64 comparray_sb_construct_end(comparray_sb *da, ushort L, int opt);

i64 comparray_sb_write(comparray_sb *da, FILE *f);
void comparray_sb_read(comparray_sb *da, uchar **map);
int comparray_sb_getbit(comparray_sb *da, i64 i);
i64 comparray_sb_rank(comparray_sb *da, i64 i);
i64 comparray_sb_rank0(comparray_sb *da, i64 i);


#endif // _COMPARRAY_H_
