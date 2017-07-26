
#ifndef FASTQ_BATCH_LIST_H
#define FASTQ_BATCH_LIST_H

#include <stdio.h>
#include <pthread.h>

#include "fastq_batch.h"
#include "fastq_file.h"

/* **************************************
 *  		Structures		*
 * *************************************/

/**
* @brief Item containing fastq_batch structure
* 
* List item that contains a fastq batch
*/
typedef struct fastq_batch_list_item {
    int id;					/**< Batch id. */
    fastq_batch_t *batch_p;			/**< Pointer to the fastq_batch element. */
    struct fastq_batch_list_item *prev_p;	/**< Pointer to the previous fastq_batch_list_item. */
    struct fastq_batch_list_item *next_p;	/**< Pointer to the next fastq_batch_list_item. */
} fastq_batch_list_item_t;

/**
* @brief List of fastq_batch_list_item elements
* 
* List containing and linking fastq_batch_list_item elements
*/
typedef struct fastq_batch_list {
    int length;					/**< Total length of the list. */
    int length_by_source_id[MAX_NUM_PRODUCERS];	/**< Length of list by source id (pair 1 and pair2). */
    int producers;				/**< Number of threads inserting in the list. */
    fastq_batch_list_item_t *first_p;		/**< Pointer to the first item in the list. */
    fastq_batch_list_item_t *last_p;		/**< Pointer to the last item in the list. */
    pthread_mutex_t lock;			/**< Lock. */
} fastq_batch_list_t;

/* **************************************
 *  		Functions		*
 * *************************************/

/**
*  @brief Frees a fastq batch list item
*  @param[in,out] fastq_batch_list_item_p pointer to the fastq batch list
*  @param all flag to indicate if batch is freed
*  @return void
*  
*  Free a fastq batch list item (all = 0 -> fastq batch not freed, all = 1 -> fastq batch freed)
*/
void fastq_batch_list_item_free(fastq_batch_list_item_t* fastq_batch_list_item_p, int all);

/**
*  @brief Init a fastq batch list
*  @param[in,out] fastq_batch_list_p pointer to the list to initialize
*  @param producers number of threads feeding the list
*  @return void
*  
*  Performs fastq batch list initialization
*/
void fastq_batch_list_init(fastq_batch_list_t* fastq_batch_list_p, int producers);

/**
*  @brief Inserts a fastq batch item in the given list
*  @param fastq_batch_list_item_p pointer to the fastq batch item to insert in the list
*  @param[in,out] fastq_batch_list_p pointer to the list in where item will be inserted
*  @return void
*  
*  Inserts a fastq batch list item in a fastq batch list
*/
void fastq_batch_list_insert(fastq_batch_list_item_t* fastq_batch_list_item_p, fastq_batch_list_t* fastq_batch_list_p);

/**
*  @brief Removes a fastq batch list item from a list
*  @param[in,out] fastq_batch_list_p pointer to the list from where item will be removed
*  @return fastq_batch_list_item pointer to the removed item
*  
*  Removes the first item from the list and returns it
*/
fastq_batch_list_item_t* fastq_batch_list_remove(fastq_batch_list_t* fastq_batch_list_p);

/**
*  @brief Gets fastq batch list total length
*  @param fastq_batch_list_p pointer to the list to obtain the length
*  @return total length of the fastq batch list
*  
*  Returns fastq batch list total length
*/
int fastq_batch_list_length(fastq_batch_list_t* fastq_batch_list_p);

/**
*  @brief Gets fastq batch list length by source id
*  @param fastq_batch_list_p pointer to the list to obtain the length
*  @param source_id source id to obtain the length (pair1 or pair2)
*  @return number of batches for a given source_id in the fastq batch list
*  
*  Returns the number of batches for a given source_id in the fastq batch list
*/
int fastq_batch_list_length_by_source_id(fastq_batch_list_t* list_p, int source_id);

/**
*  @brief Gets the number of producers for the given fastq batch list
*  @param fastq_batch_list_p pointer to the list to obtain the number of producers
*  @return number of producers for the given fastq batch list
*  
*  Returns the number of producers for the given fastq batch list
*/
int fastq_batch_list_get_producers(fastq_batch_list_t* fastq_batch_list_p);

/**
*  @brief Increases by 1 the number of producers for the given fastq batch list
*  @param fastq_batch_list_p pointer to the list to increase the number of producers
*  @return number of producers for the given fastq batch list
*  
*  Increases by 1 and returns the number of producers for the given fastq batch list
*/
int fastq_batch_list_incr_producers(fastq_batch_list_t* fastq_batch_list_p);

/**
*  @brief Decreases by 1 the number of producers for the given fastq batch list
*  @param fastq_batch_list_p pointer to the list to decrease the number of producers
*  @return number of producers for the given fastq batch list
*  
*  Decreases by 1 and returns the number of producers for the given fastq batch list
*/
int fastq_batch_list_decr_producers(fastq_batch_list_t* fastq_batch_list_p);

/**
*  @brief Prints a fastq batch list in the console
*  @param fastq_batch_list_p pointer to the fastq batch list
*  @return void
*  
*  Prints the content of each fastq batch from the list in the console
*/
void fastq_batch_list_print(fastq_batch_list_t* fastq_batch_list_p);

/**
*  @brief Frees fastq batch list items
*  @param[in,out] list_p pointer to the fastq batch list
*  @return void
*  
*  Free fastq batch list items
*/
void fastq_batch_list_items_free(fastq_batch_list_t* list_p);

#endif /* FASTQ_BATCH_LIST_H */
