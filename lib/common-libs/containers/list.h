#ifndef LIST_H
#define LIST_H

#include <stdio.h>
#include <stdio.h>
#include <pthread.h>

#include <commons/log.h>

//=====================================================
// structures
//=====================================================

typedef struct list_item {
  int id;
  int type;
  void* data_p;

  struct list_item* prev_p;
  struct list_item* next_p;
} list_item_t;

//-----------------------------------------------------

typedef struct list {
  char* name;
  size_t length;
  size_t max_length;

  int writers;
  int inserting;
  int removing;

  pthread_mutex_t lock;
  pthread_cond_t condition;

  list_item_t* first_p;
  list_item_t* last_p;  
} list_t;

//=====================================================
// functions
//=====================================================


list_item_t* list_item_new(int id, int type, void* data_p);
void list_item_free(list_item_t* item_p);

//-----------------------------------------------------
void list_init(char* name, int writers, size_t max_length, list_t* list_p);
void list_free_deep(list_t* list_p, void* (*data_callback) (void* data));

int list_insert_item(list_item_t* item_p, list_t* list_p);
list_item_t* list_remove_item(list_t* list_p);

int list_insert_item_async(list_item_t* item_p, list_t* list_p);
list_item_t* list_remove_item_async(list_t* list_p);

int list_get_length(list_t* list_p);
int list_get_max_length(list_t* list_p);

int list_set_writers(int writers, list_t* list_p);
int list_get_writers(list_t* list_p);
int list_incr_writers(list_t* list_p);
int list_decr_writers(list_t* list_p);

void list_print(list_t * list_p);

void **list_to_array(list_t *list_p);

#endif /* LIST_H */
