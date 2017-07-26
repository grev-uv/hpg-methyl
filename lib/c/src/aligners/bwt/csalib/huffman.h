/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _HUFFMAN_H_
#define _HUFFMAN_H_

#include "csa.h"

#define HUFTBLWIDTH 8
#define HUFTBLSIZ (1<<HUFTBLWIDTH)

typedef struct {
  int n;  // alphabet size
  int *v;
  int *left, *right;
  byte *clen;
  u64 *code;
  short tbl[HUFTBLSIZ];
} Huffman;

void freeHuffman(Huffman *p);
int DecodeHuffman(Huffman *h, u64 x);
int DecodeHuffman_tbl(Huffman *h, u64 x);
Huffman *MakeHuffmanTree(int n, double *freq);
void Huffman_write(Huffman *h, FILE *out);
Huffman *Huffman_read(FILE *in);
Huffman *Huffman_read2(uchar **map);

#endif // _HUFFMAN_H_
