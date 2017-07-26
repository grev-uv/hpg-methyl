/*
   Copyright 2010, Kunihiko Sadakane, all rights reserved.

   This software may be used freely for any purpose.
   No warranty is given regarding the quality of this software.
*/
#include "typedef.h"
#include "mman.h"

#ifdef WIN32
MMAP *mymmap (char *fname)
{
void *base;
HANDLE fd,h;
i64 len;
MMAP *m;
  m = malloc(sizeof(*m));
  if (m==NULL) {perror("mymmap malloc");  exit(1);}
  fd = CreateFile(fname,GENERIC_READ,FILE_SHARE_READ,NULL,OPEN_EXISTING,
         FILE_ATTRIBUTE_NORMAL,NULL);
  if (fd==INVALID_HANDLE_VALUE) {
    printf("createfile\n");
    exit(1);
  }
  m->h1 = fd;
  len = GetFileSize(fd,0);
  m->len = len;
  //h = CreateFileMapping (fd, NULL, PAGE_READONLY, 0, len, NULL);
  h = CreateFileMapping (fd, NULL, PAGE_READONLY, 0, 0, NULL);
  if (h==NULL) {
    printf("createfilemapping\n");
    exit(1);
  }
  m->h2 = h;
  //base = MapViewOfFile (h, FILE_MAP_READ, 0, 0, len);
  base = MapViewOfFile (h, FILE_MAP_READ, 0, 0, 0);
  if (base==NULL) {
    printf("mapviewoffile\n");
    return NULL;
  }
  m->addr = base;
  return m;
}

int mymunmap (MMAP *m)
{
 UnmapViewOfFile (m->addr);
 CloseHandle(m->h2);
 CloseHandle(m->h1);
 return 0;
}                

#else

MMAP *mymmap (char *fname)
{
int fd;
i64 len;
MMAP *m;
struct stat statbuf;
void *base;
  m = (MMAP *) malloc(sizeof(*m));
  if (m==NULL) {perror("mymmap malloc");  exit(1);}

  stat(fname,&statbuf);
  len = statbuf.st_size;
  fd = open(fname,O_RDONLY);
  base = mmap(0,len,PROT_READ,MAP_SHARED,fd,0);
  if (base==(void *)-1) {
    perror("mmap1\n");
    exit(1);
  }
//  printf("mmap: fname=%s fd=%d base=%p len=%ld\n",fname,fd,base,len);
  m->addr = base;
  m->fd = fd;
  m->len = len;
  return m;
}

int mymunmap (MMAP *m)
{
  if (munmap(m->addr,m->len)==-1) {
    perror("munmap 1:");
  }
  close(m->fd);
  return 0;
}                
#endif
