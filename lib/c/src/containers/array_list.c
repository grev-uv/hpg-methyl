#include "array_list.h"


/**
 * array_list items functions
 */
//void* array_list_item_new(void *type, void *item_p) {
//	void *item = (void*) malloc(sizeof(void));
//	item->type = type;
//	item->data = item_p;
//
//	return item;
//}
//
//void array_list_item_free(void *item_p) {
//	if(item_p != NULL) {
//		if(item_p->type != NULL) {
//			free(item_p->type);
//		}
//		if(item_p->data != NULL) {
//			free(item_p->data);
//		}
//		free(item_p);
//	}
//}


/**
 * array_list functions
 */

int compare_items(const void *item1, const void *item2) {
	return item1 != item2;
}

array_list_t* array_list_new(size_t initial_capacity, float realloc_factor, int SYNC_MODE) {
	array_list_t *array_list_p = (array_list_t*) malloc(sizeof(array_list_t));
	array_list_p->capacity = initial_capacity;
	array_list_p->size = 0;
	array_list_p->realloc_factor = realloc_factor;
	array_list_p->mode = SYNC_MODE;
	array_list_p->compare_fn = compare_items;
	array_list_p->flag = 0;
	array_list_p->items = (void**) malloc(initial_capacity * sizeof(void*));

	if (array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
	  pthread_mutex_init(&(array_list_p->lock), NULL);
	}

	return array_list_p;
}

array_list_t* array_list_dup(array_list_t *array_list_p) {
        array_list_t* new_list = NULL;

        if(array_list_p != NULL) {
	  if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
	      pthread_mutex_lock(&array_list_p->lock);
	  }

	  new_list = array_list_new(array_list_p->size, array_list_p->realloc_factor, array_list_p->mode);
	  for (size_t i = 0; i < array_list_p->size; i++) {
	    array_list_insert(array_list_get(i, array_list_p), new_list);
	  }

	  if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
	    pthread_mutex_unlock(&array_list_p->lock);
	  }

	}

	return new_list;
}

static array_list_t *reallocate(array_list_t * array_list_p, size_t inc_size) {
	// Capacity is increased in factor.
	size_t new_capacity;// = array_list_t->capacity;
	if(!inc_size) {
		new_capacity = (int)ceil((float)array_list_p->capacity * array_list_p->realloc_factor);
		//new_capacity = (int)ceil((float)array_list_p->capacity * array_list_p->realloc_factor);
	}else {
		new_capacity = array_list_p->capacity + inc_size;
	}
	// Realloc items with the new capacity. Size remains equals.
	void **items_aux = (void**) realloc(array_list_p->items, new_capacity * sizeof(void*));
	if(items_aux != NULL) {
		array_list_p->items = items_aux;
		array_list_p->capacity = new_capacity;
	}else {
		LOG_ERROR("Error in reallocate");
	}
	return array_list_p;
}
 
//void array_list_init(size_t initial_capacity, float realloc_factor, int SYNC_MODE, array_list_t *array_list_p) {
//	array_list_p->capacity = initial_capacity;
//	array_list_p->size = 0;
//	array_list_p->realloc_factor = realloc_factor;
//	array_list_p->compare_fn = compare_items;
//	if (array_list_p->items) {
//		array_list_p->items = (void**) realloc(initial_capacity * sizeof(void*));
//	} else {
//		array_list_p->items = (void**) malloc(initial_capacity * sizeof(void*));
//	}
//}

int array_list_clear(array_list_t *array_list_p,  void (*data_callback) (void* data)) {
	if(array_list_p != NULL) {
		// Free c
	        if (data_callback != NULL) {
		        for(size_t i=0; i < array_list_p->size; i++) {
			        if(array_list_p->items != NULL && array_list_p->items[i] != NULL) {
				        data_callback(array_list_p->items[i]);
			        }
		         }
		}
		// Set default parameters
		array_list_p->size = 0;
//		array_list_init(10, 1.5, COLLECTION_MODE_SYNCHRONIZED, array_list_p);
		return 1;
	}
	return 0;
}

void array_list_free(array_list_t *array_list_p, void (*data_callback) (void* data)) {
	if(array_list_p != NULL) {
		array_list_clear(array_list_p, data_callback);
		free(array_list_p->items);
		free(array_list_p);
	}
}



size_t array_list_capacity(array_list_t *array_list_p) {
	if(array_list_p != NULL) {
		return array_list_p->capacity;
	}
	return ULONG_MAX;
}

size_t array_list_size(array_list_t *array_list_p) {
	size_t length;
	
	if(array_list_p != NULL) {
	  if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
	      pthread_mutex_lock(&array_list_p->lock);
	  }
	  
	  length = array_list_p->size;
	
	  if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
	      pthread_mutex_unlock(&array_list_p->lock);
	  }
	}else{
	  length = ULONG_MAX;
	}
	
	return length;
}

size_t array_list_index_of(void *item_p, array_list_t* array_list_p) {
	if(array_list_p != NULL) {
		if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&array_list_p->lock);
		}

		for(size_t i=0; i<array_list_p->size; i++) {
			if(array_list_p->compare_fn(array_list_p->items[i], item_p) == 0) {
				return i;
			}
		}

		if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&array_list_p->lock);
		}
	}
	return ULONG_MAX;
}

int array_list_contains(void *item_p, array_list_t *array_list_p) {
	if(array_list_p != NULL) {
		if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&array_list_p->lock);
		}

		for(size_t i=0; i<array_list_p->size; i++) {
			if(array_list_p->compare_fn(array_list_p->items[i], item_p) == 0) {
				return 1;
			}
		}

		if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&array_list_p->lock);
		}
	}
	return 0;
}



int array_list_insert(void *item_p, array_list_t *array_list_p) {
	if(array_list_p != NULL) {
		if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&array_list_p->lock);
		}

		// Check if array_list is full
		if(array_list_p->size >= array_list_p->capacity) {
			array_list_p = reallocate(array_list_p, 0);
		}
			
		array_list_p->items[array_list_p->size] = item_p;
		array_list_p->size++;
		if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&array_list_p->lock);
		}
		return 1;
	}
	return 0;
}

int array_list_insert_at(size_t index, void *item_p, array_list_t *array_list_p) {
	if(array_list_p != NULL) {
		if(index < array_list_p->size) {
			if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
				pthread_mutex_lock(&array_list_p->lock);
			}

			// Check if array_list is full
			if(array_list_p->size >= array_list_p->capacity) {
				array_list_p = reallocate(array_list_p, 0);
			}
			// move the items to the end
			for(size_t i=array_list_p->size; i>index; i--) {
				array_list_p->items[i] = array_list_p->items[i-1];
			}
			array_list_p->items[index] = item_p;
			array_list_p->size++;

			if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
				pthread_mutex_unlock(&array_list_p->lock);
			}
			return 1;
		}
	}
	return 0;
}


int array_list_insert_all(void **items_p, size_t num_items, array_list_t *array_list_p) {
	if(array_list_p != NULL) {
		if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&array_list_p->lock);
		}

		// First check if we need to realloc items.
		// This avoids to check size for each insert
		if(array_list_p->size + num_items > array_list_p->capacity) {
			array_list_p = reallocate(array_list_p, num_items+1);
		}
		for(size_t i=0; i<num_items; i++) {
			array_list_p->items[array_list_p->size++] = items_p[i];
		}

		if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&array_list_p->lock);
		}
		return 1;
	}
	return 0;
}

int array_list_insert_all_at(size_t index, void **items_p, size_t num_items, array_list_t *array_list_p) {
	if(array_list_p != NULL) {
		if(index < array_list_p->size) {
			if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
				pthread_mutex_lock(&array_list_p->lock);
			}

			// First check if we need to realloc items.
			// This avoids to check size for each insert
			if(array_list_p->size + num_items > array_list_p->capacity) {
				array_list_p = reallocate(array_list_p, num_items+1);
			}
			// move the items to the end
			for(size_t i=array_list_p->size; i>index; i--) {
				array_list_p->items[i] = array_list_p->items[i-1];
			}
			for(size_t i=0; i<num_items; i++) {
				array_list_p->items[array_list_p->size++] = items_p[i];
			}

			if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
				pthread_mutex_unlock(&array_list_p->lock);
			}
			return 1;
		}
	}
	return 0;
}



void* array_list_remove(void *item_p, array_list_t *array_list_p) {
	if(array_list_p != NULL) {
		if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&array_list_p->lock);
		}
		
		size_t index = array_list_index_of(item_p, array_list_p);
		
		if(index != ULONG_MAX) {
			return array_list_remove_at(index, array_list_p);
		}

		if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&array_list_p->lock);
		}
		return NULL;
	}
	return NULL;
}

void* array_list_remove_at(size_t index, array_list_t *array_list_p) {
	if(array_list_p != NULL && index >= 0 && index < array_list_p->size) {
		if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&array_list_p->lock);
		}

		void *aux = array_list_p->items[index];
		for(size_t i=index; i<array_list_p->size - 1; i++) {
			array_list_p->items[i] = array_list_p->items[i+1];
		}

		if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&array_list_p->lock);
		}
		
        array_list_p->size--;
		return aux;
	}
	return NULL;
}
//TODO
void** array_list_remove_range(size_t start, size_t end, array_list_t *array_list_p) {
	return NULL;
}



void* array_list_get(size_t index, array_list_t *array_list_p) {
	void *item_p;
	if(array_list_p != NULL && index >= 0 && index < array_list_p->size) {
	    if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
		    pthread_mutex_lock(&array_list_p->lock);
	    }
	    
		  item_p = array_list_p->items[index];

	    if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
		    pthread_mutex_unlock(&array_list_p->lock);
	    }

	}else{
	    item_p = NULL;
	}
	
	return item_p;
}

//TODO
array_list_t* array_list_sublist(size_t start, size_t end, array_list_t *array_list_p, array_list_t *sublist) {
    return NULL;
}

void* array_list_set(size_t index, void *item_p, array_list_t *array_list_p) {
	if(array_list_p != NULL && index >= 0 && index < array_list_p->size) {
		array_list_p->items[index] = item_p;
		return array_list_p->items[index];
	}
	return NULL;
}

int array_list_qsort(array_list_t *array_list_p, int (*compare_fn)(const void*,const void*))
{
	if(array_list_p != NULL && compare_fn != NULL)
	{
		if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_lock(&array_list_p->lock);
		}

		//Quicksort
		qsort(array_list_p->items, array_list_p->size, sizeof(void *), compare_fn);

		if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
			pthread_mutex_unlock(&array_list_p->lock);
		}
	}
	return 0;
}


void array_list_print(array_list_t *array_list_p) {
	printf("[");
	for(size_t i=0; i < array_list_p->size-1; i++) {
		printf("%s,", (char*)array_list_p->items[i]);
	}
	printf("%s", (char*)array_list_p->items[array_list_p->size-1]);
	printf("]");
}

// void **list_to_array(array_list_t *array_list_p) {
// 	return NULL;
// }


int array_list_swap(const int pos1, const int pos2, array_list_t *array_list_p){
	
	if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
		pthread_mutex_lock(&array_list_p->lock);
	}
	
	size_t size = array_list_p->size;
	
	if(((pos1 < 0) || (pos1 > size)) || ((pos2 < 0) || (pos2 > size)) ){
	  return 0;
	}else{
		void *tmp = array_list_p->items[pos1];
		array_list_p->items[pos1] = array_list_p->items[pos2];
		array_list_p->items[pos2] = tmp;
	}

	if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
		pthread_mutex_unlock(&array_list_p->lock);
	}
	
	return 1;
}

void array_list_set_flag(int flag, array_list_t *array_list_p) {
  if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
    pthread_mutex_lock(&array_list_p->lock);
  }

  array_list_p->flag = flag;

  if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
    pthread_mutex_unlock(&array_list_p->lock);
  }
}

int array_list_get_flag(array_list_t *array_list_p) {
  int flag;
  if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
    pthread_mutex_lock(&array_list_p->lock);
  }

  flag = array_list_p->flag;

  if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
    pthread_mutex_unlock(&array_list_p->lock);
  }

  return flag;
}


//int list_set_writers(int writers, list_t* list_p);
//int list_get_writers(list_t* list_p);
//int list_incr_writers(list_t* list_p);
//int list_decr_writers(list_t* list_p);

array_list_t* array_list_unique(array_list_t *orig, int (*compare)(const void *a, const void *b), array_list_t *dest) {
    char* item_description = NULL;
    size_t orig_list_length = array_list_size(orig);
    cp_hashtable* visited = cp_hashtable_create(orig_list_length, cp_hash_string, cp_hash_compare_string);
    int value; // variable for inserting a pointer as a value for key-value pair
    
    for (size_t i = 0; i < orig_list_length; i++) {
        item_description = (char*) array_list_get(i, orig);

        if (cp_hashtable_contains(visited, (void*) item_description) == 0) {
            array_list_insert(item_description, dest);
            cp_hashtable_put(visited, item_description, &value);  // inserting NULL does not allow future searches
        }
    }
    
    // free visited hashtable
    cp_hashtable_destroy(visited);
    
    return dest;
}

array_list_t* array_list_intersect(array_list_t *al1, array_list_t *al2, int (*compare)(const void *a, const void *b), array_list_t *dest) {
    char* item_description = NULL;
    array_list_t *al_little = NULL, *al_big = NULL;
    size_t al1_length = array_list_size(al1);
    size_t al2_length = array_list_size(al2);
    size_t al_little_length, al_big_length;
    int value; // variable for inserting a pointer as a value for key-value pair
    
    if (al1_length <= al2_length) {
        al_little = al1;
        al_big = al2;
        al_little_length = al1_length;
        al_big_length = al2_length;
    } else {
        al_little = al2;
        al_big = al1;
        al_little_length = al2_length;
        al_big_length = al1_length;
    }    
    
    // put into the hash the largest array list for time efficiency (i.e.: 10*log(70) < 70*log(10)
    cp_hashtable* visited = cp_hashtable_create(al_big_length, cp_hash_string, cp_hash_compare_string);

    // fill the hash with the biggest array list
    for (size_t i = 0; i < al_big_length; i++) {
        item_description = (char*) array_list_get(i, al_big);

        if (cp_hashtable_contains(visited, item_description) == 0) {
            cp_hashtable_put(visited, item_description, &value);
        }
    }
    
    // if item is contained in the hash/big array list then put it in the destination array list
    for (size_t i = 0; i < al_little_length; i++) {
        item_description = (char*) array_list_get(i, al_little);

        if (cp_hashtable_contains(visited, item_description) != 0) {
            array_list_insert(item_description, dest);        
        }
    }
    
    // free visited hashtable
    cp_hashtable_destroy(visited);
    
    return dest;
}

array_list_t* array_list_complement(array_list_t *al1, array_list_t *al2, int (*compare)(const void *a, const void *b), array_list_t *dest) {
    char* item_description = NULL;
    size_t al1_length = array_list_size(al1);
    size_t al2_length = array_list_size(al2);
    int value; // variable for inserting a pointer as a value for key-value pair

    // put into the hash the array list 1 (it will be the search space)
    cp_hashtable* visited = cp_hashtable_create(al1_length, cp_hash_string, cp_hash_compare_string);

    // fill the hash with array list 1
    for (size_t i = 0; i < al1_length; i++) {
        item_description = (char*) array_list_get(i, al1);

        // OPTIMIZATION: this condition can be eliminated if array list have unique values
        if (cp_hashtable_contains(visited, item_description) == 0) {
            cp_hashtable_put(visited, item_description, &value);
        }
    }
    
    // if item from array list 2 is NOT contained in the array list 1 then put it in the destination array list
    for (size_t i = 0; i < al2_length; i++) {
        item_description = (char*) array_list_get(i, al2);

        if (cp_hashtable_contains(visited, item_description) == 0) {
            array_list_insert(item_description, dest); 
        }
    }
    
    // free visited hashtable
    cp_hashtable_destroy(visited);
    
    return dest;
}

int compare(const void *a, const void *b) {
    return strcmp((char*)a, (char*)b);
}


int array_list_replace_at(size_t index, void *item_p, array_list_t *array_list_p) {
  if(array_list_p != NULL) {
    if(index < array_list_p->size) {
      if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
	pthread_mutex_lock(&array_list_p->lock);
      }
      
      array_list_p->items[index] = item_p;

      if(array_list_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
	pthread_mutex_unlock(&array_list_p->lock);
      }
      return 1;
    }
  }

  return 0;

}
