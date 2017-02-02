/*
 * linked_list.h
 *
 *  Created on: Nov 7, 2012
 *  Last modified: Mar 26, 2013
 *      Author: imedina, hmartinez, cgonzalez
 */

#ifndef LINKED_LIST_H
#define LINKED_LIST_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include "containers.h"
#include <commons/string_utils.h>
#include <commons/log.h>

typedef struct linked_list_item {
	void *item;

	struct linked_list_item *prev;
	struct linked_list_item *next;
} linked_list_item_t;

typedef struct linked_list {
	size_t size;
	int mode;
	int flag;

	// linked_list_compare_fn compare_fn;
	int (*compare_fn)(const void *, const void *);

	pthread_mutex_t lock;
	pthread_cond_t condition;

	linked_list_item_t *first;
	linked_list_item_t *last;
} linked_list_t;

typedef struct linked_list_iterator {
  linked_list_t *linked_list_p;
  linked_list_item_t *curr_pos;
} linked_list_iterator_t;



void print_item(void *item);

/*********************************************************************/
/***             LinkedList Iterator Implementation               ***/
/*******************************************************************/

linked_list_iterator_t* linked_list_iterator_new(linked_list_t *linked_list_p);

linked_list_iterator_t* linked_list_iterator_new_by_item(linked_list_t *linked_list_p, linked_list_item_t *item);

linked_list_iterator_t* linked_list_iterator_init(linked_list_t *linked_list_p, linked_list_iterator_t *iterator_p);

linked_list_iterator_t* linked_list_iterator_init_by_item(linked_list_t *linked_list_p, 
							  linked_list_item_t *item,
							  linked_list_iterator_t *iterator_p);

void linked_list_iterator_free(linked_list_iterator_t *iterator_p);

void* linked_list_iterator_curr(linked_list_iterator_t *iterator_p);

linked_list_item_t* linked_list_iterator_list_item_curr(linked_list_iterator_t *iterator_p);

void* linked_list_iterator_next(linked_list_iterator_t *iterator_p);

void* linked_list_iterator_prev(linked_list_iterator_t *iterator_p);

void* linked_list_iterator_last(linked_list_iterator_t *iterator_p);

void* linked_list_iterator_first(linked_list_iterator_t *iterator_p);

int linked_list_iterator_insert(void *item, linked_list_iterator_t *iterator_p);

void* linked_list_iterator_remove(linked_list_iterator_t *iterator_p);

linked_list_item_t* linked_list_iterator_remove_2(linked_list_iterator_t *iterator_p);

/*****************************************************************/

linked_list_t* linked_list_new(int SYNC_MODE);

void linked_list_free(linked_list_t* linked_list_p, void (*data_callback) (void* data));

linked_list_item_t* linked_list_item_new(void *item);

void linked_list_item_free(linked_list_item_t *linked_list_item, void (*data_callback) (void* data));


size_t linked_list_size(linked_list_t *linked_list_p);

size_t linked_list_index_of(void *item, linked_list_t *linked_list_p);

int linked_list_contains(void* item, linked_list_t *linked_list_p);

int linked_list_clear(linked_list_t *linked_list_p, void (*data_callback) (void* data));


int linked_list_insert(void* item_p, linked_list_t *linked_list_p);

int linked_list_insert_first(void* item_p, linked_list_t *linked_list_p);

int linked_list_insert_last(void* item_p, linked_list_t *linked_list_p);

int linked_list_insert_at(size_t index, void* item_p, linked_list_t *linked_list_p);

int linked_list_insert_all(void** item_p, size_t num_items, linked_list_t *linked_list_p);

int linked_list_insert_all_at(size_t index, void** item_p, size_t num_items, linked_list_t* linked_list_p);



void* linked_list_remove(void *item, linked_list_t *linked_list_p);

void* linked_list_remove_first(linked_list_t *linked_list_p);

void* linked_list_remove_last(linked_list_t *linked_list_p);

void* linked_list_remove_at(size_t index, linked_list_t *linked_list_p);

void** linked_list_remove_range(size_t start, size_t end, linked_list_t* linked_list_p);



void* linked_list_get(size_t index, linked_list_t *linked_list_p);

void* linked_list_get_first(linked_list_t *linked_list_p);

void* linked_list_get_last(linked_list_t *linked_list_p);

linked_list_t* linked_list_sublist(size_t start, size_t end, linked_list_t *linked_list_p, linked_list_t *sublist);

void* linked_list_set(size_t index, void* new_item, linked_list_t *linked_list_p);



void linked_list_print(linked_list_t *linked_list_p, void (*data_callback) (void* data));

// void **linked_list_to_array(linked_list_t *linked_list_p);

static int compare_items(const void *item1, const void *item2);


int linked_list_swap(const int pos1, const int pos2, linked_list_t *linked_list_p);


void linked_list_set_flag(int flag, linked_list_t *linked_list_p);

int linked_list_get_flag(linked_list_t *linked_list_p);


/**
 *  @brief Compare function for strings
 *  @param a pointer to string
 *  @param b pointer to string
 *  @return 0: if equal, <>0: if not equal
 *
 *  Compare function for strings
 */
int compare(const void *a, const void *b);


#endif /* LINKED_LIST_H_ */
