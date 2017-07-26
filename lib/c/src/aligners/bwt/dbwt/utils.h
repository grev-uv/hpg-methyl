#ifndef _DBWT_UTILS_
#define _DBWT_UTILS_

#include <stdio.h>
#include <stdlib.h>

#include "../csalib/typedef.h"
#include "../bwt_commons.h"

typedef word pb;
#define _logD_ 4

#define PBS (sizeof(pb)*8)
#define _D_ (1<<_logD_)

unsigned int getbits(unsigned short *B, unsigned long i, int d);

typedef struct {
  ulong n;
  int w;
  pb *b;
} packed_array;

int blog(ulong x);
packed_array *allocate_packed_array(ulong n, int w);
void free_packed_array(packed_array *p);
ulong pa_get(packed_array *p, ulong i);
void pa_set(packed_array *p, ulong i, long x);

pb *allocate_vector(ulong n);
int getbit(pb *B, ulong i);
int setbit(pb *B, ulong i,int x);
dword getbits(pb *B, ulong i, int d);
int setbits(pb *B, ulong i, int d, ulong x);

#endif

/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.
    
   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
