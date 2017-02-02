#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "list.h"

//=====================================================
// functions to manage list item
//=====================================================

list_item_t* list_item_new(int id, int type, void* data_p) {
  list_item_t* item_p;
  
  item_p = (list_item_t*) calloc(1, sizeof(list_item_t));

  item_p->id = id;
  item_p->type = type;
  item_p->data_p = data_p;

  return item_p;
}

//-----------------------------------------------------

void list_item_free(list_item_t* item_p) {
  free(item_p);
}

//=====================================================
// functions to manage list
//=====================================================

//-----------------------------------------------------
// list_init
//
// init to zero the object and initialize the lock
//-----------------------------------------------------

void list_init(char* name, int writers, size_t max_length, list_t* list_p) {
   memset(list_p, 0, sizeof(list_t));

   pthread_mutex_init(&(list_p->lock), NULL);
   //pthread_cond_init(&(list_p->condition), NULL);
   
   list_p->name = name;
   list_p->writers = writers;
   list_p->max_length = max_length;

   list_p->inserting = 0;
   list_p->removing = 0;
}

void list_free_deep(list_t* list_p, void* (*data_callback) (void* data)) {
    list_item_t *freed = NULL;
    while ((freed = list_remove_item_async(list_p)) != NULL) {
        if (data_callback) {
            data_callback(freed->data_p);
        }
        list_item_free(freed);
    }
    free(list_p);
}

//-----------------------------------------------------
// list_insert_item asychronous
//
// insert a object in the end of the list,
// according to fifo order,
//-----------------------------------------------------

int list_insert_item_async(list_item_t* item_p, list_t* list_p) {

  if (list_p==NULL) return 0;

  pthread_mutex_lock(&list_p->lock);

  if (list_p->length >= list_p->max_length) {
    pthread_mutex_unlock(&list_p->lock);
    return 0;
  }

  if (list_p->first_p==NULL) {
    item_p->prev_p = NULL;
    item_p->next_p = NULL;
    list_p->first_p = item_p;
    list_p->last_p = item_p;
  } else {
    list_p->last_p->next_p = item_p;
    item_p->prev_p = list_p->last_p;
    item_p->next_p = NULL;
    list_p->last_p = item_p;
  }
  list_p->length++;

  pthread_mutex_unlock(&list_p->lock);

  return 1;
}


/**
 * list_remove_item asynchronously
 *
 * remove the first object in the begining 
 * of the list,
 * according to fifo order,
 */

list_item_t* list_remove_item_async(list_t* list_p) {
  
  if (list_p==NULL) return NULL;

  pthread_mutex_lock(&list_p->lock);

  // just get the first element, and if is not null 
  // update the first-element pointer
  //
  list_item_t* item_p = list_p->first_p;

  if (item_p!=NULL) {
    list_p->first_p = item_p->next_p;
    list_p->length--;
  }
  
  pthread_mutex_unlock(&list_p->lock);

  return item_p;
}


//-----------------------------------------------------
// syncrhonous insert item (blocking call)
//
// insert a object in the end of the list,
// according to fifo order,
//-----------------------------------------------------

int list_insert_item(list_item_t* item_p, list_t* list_p) {

  if (list_p==NULL) return 0;

  pthread_mutex_lock(&list_p->lock);
  list_p->inserting++;

  while (list_p->length >= list_p->max_length) {
//    LOG_DEBUG("-->Insert cond wait\n");
    pthread_cond_wait(&list_p->condition, &list_p->lock);
  }

//  LOG_DEBUG_F("--> inserting %s: length = %i, max_length = %i\n", list_p->name, list_p->length, list_p->max_length);

  if (list_p->first_p==NULL) {
    item_p->prev_p = NULL;
    item_p->next_p = NULL;
    list_p->first_p = item_p;
    list_p->last_p = item_p;
  } else {
    list_p->last_p->next_p = item_p;
    item_p->prev_p = list_p->last_p;
    item_p->next_p = NULL;
    list_p->last_p = item_p;
  }
  list_p->length++;
//  LOG_DEBUG_F("<-- inserting %s: length = %i, max_length = %i\n", list_p->name, list_p->length, list_p->max_length);

  list_p->inserting--;

  if (list_p->removing > 0) {
    pthread_cond_broadcast(&list_p->condition);
  }
  pthread_mutex_unlock(&list_p->lock);

  return 1;
}
//-----------------------------------------------------
// synchronous remove item (blocking call)
//
// remove the first object in the begining 
// of the list,
// according to fifo order,
//-----------------------------------------------------

list_item_t* list_remove_item(list_t* list_p) {
  
  if (list_p==NULL) return NULL;

  pthread_mutex_lock(&list_p->lock);
  list_p->removing++;

  // just get the first element, and if is not null 
  // update the first-element pointer
  //
  list_item_t* item_p;

  while ((item_p = list_p->first_p) == NULL) {
    if (list_p->writers == 0)
      break;
    pthread_cond_wait(&list_p->condition, &list_p->lock);
  }

  if (item_p!=NULL) {
    list_p->first_p = item_p->next_p;
    list_p->length--;
  }

  list_p->removing--;
  
  if (list_p->inserting > 0) {
    pthread_cond_broadcast(&list_p->condition);
  }

  pthread_mutex_unlock(&list_p->lock);

  return item_p;
}

//-----------------------------------------------------
// list_get_length
//
// returns list length
//-----------------------------------------------------

int list_get_length(list_t* list_p) {

  int length = 0;

  if (list_p==NULL) return length;

  pthread_mutex_lock(&list_p->lock);
  length = list_p->length;
  pthread_mutex_unlock(&list_p->lock);

  return length;
}

//-----------------------------------------------------
// list_get_max_length
//
// returns list maximum length
//-----------------------------------------------------

int list_get_max_length(list_t* list_p) {

  int length = 0;

  if (list_p==NULL) return length;

  pthread_mutex_lock(&list_p->lock);
  length = list_p->max_length;
  pthread_mutex_unlock(&list_p->lock);

  return length;
}

//-----------------------------------------------------
// list_set_producers
//
// set the number of writers of the input list
//-----------------------------------------------------

int list_set_writers(int writers, list_t* list_p) {

  if (list_p==NULL) return 0;

  pthread_mutex_lock(&list_p->lock);
  list_p->writers = writers;
  pthread_mutex_unlock(&list_p->lock);

  return writers;  
}

//-----------------------------------------------------
// list_get_producers
//
// get the number of writers of the input list
//-----------------------------------------------------

int list_get_writers(list_t* list_p) {

  int writers = 0;

  if (list_p==NULL) return writers;

  pthread_mutex_lock(&list_p->lock);
  writers = list_p->writers;
  pthread_mutex_unlock(&list_p->lock);

  return writers;  
}

//-----------------------------------------------------
// list_incr_writers
//
// increment the number of writers in the list
//-----------------------------------------------------

int list_incr_writers(list_t* list_p) {

  int writers = 0;

  if (list_p==NULL) return 0;

  pthread_mutex_lock(&list_p->lock);
  list_p->writers++;
  writers = list_p->writers;

  pthread_mutex_unlock(&list_p->lock);

  return list_p->writers;
}

//-----------------------------------------------------
// list_decr_writers
//
// decrement the number of writers in the list
//-----------------------------------------------------

int list_decr_writers(list_t* list_p) {

  int writers = 0;

  if (list_p==NULL) return 0;

  pthread_mutex_lock(&list_p->lock);
  list_p->writers--;
  writers = list_p->writers;

  if(list_p->removing > 0);
     pthread_cond_broadcast(&list_p->condition);
  pthread_mutex_unlock(&list_p->lock);

  return writers;
}

//-----------------------------------------------------
// print batch list
//
// 
//-----------------------------------------------------

void list_print(list_t* list_p) {

  if (list_p==NULL) return;
  
  pthread_mutex_lock(&list_p->lock);

/*
>>>>>>> c6a2ae1bccac4045850bc3dd6ed3f2a081272e73
  printf("Number of items: %i\n", list_p->length);
<<<<<<< HEAD
    
  bam_data_batch_list_item_t* item_p = list_p->first_p;
 
  while(item_p!=NULL) {
    printf("batch id: %i", item_p->id);
    printf("\talignments: %i\n", item_p->num_alignments);

    for(int i=0; i<item_p->num_alignments; i++) {
	//printf("i: %d, strand: %i\n", i, item_p->batch_p[i].strand);
	//printf("i: %d, sequence: %i\n", i, item_p->sequence[i]);
	//printf("i: %d, cigar: %i\n", i, item_p->cigar[i]);
    }    
    
    item_p = item_p->next_p;
  }

=======
  list_item_t *item_p;
  int read;
  for(item_p=list_p->first_p; item_p != NULL; item_p=item_p->next_p)
  {
     printf("Item nº %d\n", item_p->id);
     fastq_batch_t *batch;
     batch=(fastq_batch_t *)item_p->data_p;
     
     for(read=0; read< batch->num_reads;  read++)
     {
        printf("read %d", read);
	printf("  %s\n", &(batch->data[batch->data_indices[read]]) );
     }
     printf("TOTAL READS %d\n", read);
  }
>>>>>>> c6a2ae1bccac4045850bc3dd6ed3f2a081272e73
*/  
  pthread_mutex_unlock(&list_p->lock);
}

//-----------------------------------------------------
// bam_data_batch_list_items_free
//
// free list itmes
//
//-----------------------------------------------------
/*
void bam_data_batch_list_items_free(bam_data_batch_list_t* list_p) {
  bam_data_batch_list_item_t* item_p;
  
  while ((item_p = bam_data_batch_list_remove(list_p)) != NULL) {
    //printf("liberating bam data item %i\n", item_p->id);
    bam_data_batch_list_item_free(item_p, true);
  }
}
*/

void** list_to_array(list_t* list_p) {
    void **sample_data = calloc (list_p->length + 1, sizeof(void*));
    list_item_t *item = list_p->first_p;
    
    for (int j = 0; j < list_p->length && item != NULL; j++) {
        sample_data[j] = item->data_p;
        item = item->next_p;
    }
    
    return sample_data;
}
