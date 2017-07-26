/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

/***************************************************************
To compile this, you need to obtain interface.h and run_queries.c from
Pizza&Chili Corpus:  http://pizzachili.dcc.uchile.cl/ or http://pizzachili.di.unipi.it/
***************************************************************/

#include "../csa.h"
#include "interface.h"

int load_index(char *filename, void **index)
{
  CSA *csa;
  char *fname1;
  int fnamelen;
  char **argv;

  fnamelen = strlen(filename);
  fname1 = malloc(fnamelen);
  if (fname1 == NULL) {
    printf("load_index:1. not enough mem.\n");
    exit(1);
  }
  strncpy(fname1,filename,fnamelen-4);
  strcpy(fname1+fnamelen-4,".idx");

  argv = malloc(sizeof(char *)*2);
  if (argv == NULL) {
    printf("load_index:3. not enough mem.\n");
    exit(1);
  }
  argv[0] = fname1;
  argv[1] = filename;

  csa = malloc(sizeof(CSA));
  if (csa == NULL) {
    printf("load_index:2. not enough mem.\n");
    exit(1);
  }
//  printf("csa_read %s %s\n",argv[1],argv[2]);
  csa_read(csa,2,argv);
  *index = (void *)csa;

  free(argv);
  free(fname1);
  return 0;
}



/*
 * Frees the memory occupied by index.
 */
int free_index(void *index)
{
  CSA *csa;

  csa = (CSA *)index;

  return 0;
}

int count(void *index, uchar *pattern, ulong length, ulong *numocc)
{
  CSA *csa;
  i64 l,r;
  i64 mlen;

  csa = (CSA *)index;
  mlen = csa->search(pattern, length, csa, &l, &r);
  if (mlen < length) {
    *numocc = 0;
    return 0;
  }
  *numocc = (ulong)(r - l + 1);

  return 0;
}

int locate(void *index, uchar *pattern, ulong length, ulong **occ, 
		   ulong *numocc)
{
  CSA *csa;
  i64 l,r;
  ulong *buf,*buf2;
  ulong i,oc;
  i64 mlen;

  csa = (CSA *)index;
  mlen = csa->search(pattern, length, csa, &l, &r);
  if (mlen < length) {
    *numocc = 0;
    return 0;
  }
  oc = (ulong)(r - l + 1);

  buf = malloc((*numocc) * sizeof(ulong));
  if (buf == NULL) {
    printf("locate: not enough mem.\n");
    exit(1);
  }

  for (i=0; i<oc; i++) {
    buf[i] = csa->lookup(csa,l + i);
  }

  *numocc = oc;
  *occ = buf;
  return 0;
}

int extract(void *index, ulong from, ulong to, uchar **snippet,
			ulong *snippet_length)
{
  CSA *csa;
  uchar *text;
  i64 i,len;

  csa = (CSA *)index;

  from++;  to++;
  if (to > csa->n) to = csa->n;

  len = to - from + 1;
  text = malloc(len);

  csa->text(text,csa,from,to);

  *snippet = text;
  *snippet_length = (ulong)len;
  return 0;
}

int display(void *index, uchar *pattern, ulong length, ulong numc, 
			ulong *numocc, 	uchar **snippet_text, ulong **snippet_len)
{
  printf("display: not supported.\n");
  return 0;
}

int get_length(void *index, ulong *length)
{
  *length = ((CSA *)(index))->n;
  return 0;
}

int index_size(void *index, ulong *size)
{
  *size = (ulong)(((CSA *)(index))->psize + ((CSA *)(index))->isize);
  return 0;
}

int index_size_count(void *index, ulong *size)
{
  *size = (ulong)(((CSA *)(index))->psize);
  return 0;
}

char* error_index(int e)
{
  printf("error_index: %d\n",e);
  return "hoge";
}

int build_index(uchar *text, ulong length, char *build_options, void **index)
{
  printf("build_index: not supported.\n");
  return 0;
}

int save_index(void *index, char *filename)
{
  printf("save_index: not supported.\n");
  return 0;
}

int save_index_mem(void *index, uchar *compress)
{
  printf("save_index_mem: not supported.\n");
  return 99;
}

int load_index_mem (void ** index, uchar *compress, ulong size)
{
  printf("load_index_mem: not supported.\n");
  return 99;
}

