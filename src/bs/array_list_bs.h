#ifndef ARRAY_LIST_BS_H
#define ARRAY_LIST_BS_H

#include <stdio.h>
#include <math.h>
#include <limits.h>

#include "commons/string_utils.h"
#include "commons/log.h"

#include "containers/containers.h"
#include "containers/cprops/hashtable.h"

//=====================================================
// structures
//=====================================================

typedef struct metil_data {
  char*  query_name;
  char   status;
  int    chromosome;
  size_t start;
  char   context;
  int    strand;
  int    zone;
} metil_data_t;


typedef struct array_list_bs {
  size_t capacity;
  size_t size;
  float realloc_factor;
  int mode;
  int flag;

  int (*compare_fn)(const void *, const void *);

  pthread_mutex_t lock;
  pthread_cond_t condition;

  metil_data_t *items;
} array_list_bs_t;


#endif /* array_list_bs_H */
