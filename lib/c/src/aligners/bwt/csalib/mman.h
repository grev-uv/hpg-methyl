/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#ifndef _MYMMAP_H_
#define _MYMMAP_H_


#ifdef WIN32
#include <windows.h>
#include <inttypes.h>
#else
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <unistd.h>
#include <inttypes.h>
#endif

#ifdef WIN32
#define PAGE_READONLY          0x02     
#define SECTION_MAP_READ    0x0004
#define FILE_MAP_READ       SECTION_MAP_READ
#endif

typedef struct {
  void *addr;
  int64_t len;
#ifdef WIN32
  HANDLE h1,h2;
#else
  int fd;
#endif
} MMAP;

MMAP *mymmap (char *fname);
int mymunmap (MMAP *m);

#endif
