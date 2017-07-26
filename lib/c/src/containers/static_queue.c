#include "static_queue.h"


#define NEXT(i,sz)        (i + 1) % sz

//-----------------------------------------------------

static_queue_t* static_queue_new(size_t initial_capacity) {
  static_queue_t* queue = calloc(1, sizeof(static_queue_t));

  queue->first = 0;
  queue->last = 0;
  queue->capacity = initial_capacity;
  queue->data = calloc(queue->capacity, sizeof(void*));
  queue->count = 0;

  return queue;
}

//-----------------------------------------------------

void static_queue_free(static_queue_t* queue) {
  if (queue) {
    if (queue->data) {
      free(queue->data);
    }

    free(queue);
  }
}

//-----------------------------------------------------

void static_queue_clear(static_queue_t* queue) {
  queue->first = 0;
  queue->last = 0;
  queue->count = 0;
}

//-----------------------------------------------------

void static_queue_push(void* data, static_queue_t* queue) {
  assert(NEXT(queue->last, queue->capacity) != queue->first);
  
  queue->data[queue->last] = data;
  queue->last = NEXT(queue->last, queue->capacity);
  queue->count++;
}

//-----------------------------------------------------

void* static_queue_pop(static_queue_t* queue) {
  void* result = NULL;

  if (queue->first != queue->last) {
    result = queue->data[queue->first];
    queue->first = NEXT(queue->first, queue->capacity);
    queue->count--;
  }

  return result;
}

//-----------------------------------------------------

int8_t static_queue_is_full(static_queue_t* queue) {
  return NEXT(queue->last, queue->capacity) == queue->first;
}

//-----------------------------------------------------

int8_t static_queue_is_empty(static_queue_t* queue) {
  return queue->first == queue->last;
}

//-----------------------------------------------------

size_t static_queue_size(static_queue_t* queue) {
  return queue->count;
}
