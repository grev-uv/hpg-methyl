#ifndef STATIC_QUEUE_H
#define STATIC_QUEUE_H

#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>
#include <omp.h>

#include "containers.h"

typedef struct static_queue {
  size_t last;
  size_t first;
  size_t capacity;
  size_t count;

  void** data;
} static_queue_t;

static_queue_t* static_queue_new(size_t initial_capacity);

void static_queue_free(static_queue_t* queue);

void static_queue_clear(static_queue_t* queue);

void static_queue_push(void* data, static_queue_t* queue);

void* static_queue_pop(static_queue_t* queue);

int8_t static_queue_is_full(static_queue_t* queue);

int8_t static_queue_is_empty(static_queue_t* queue);

size_t static_queue_size(static_queue_t* queue);

#endif // STATIC_QUEUE_H