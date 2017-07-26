
#ifndef FASTQ_BATCH_READER_OMP_H
#define FASTQ_BATCH_READER_OMP_H

#include <stdio.h>
#include <stdio.h>
#include <pthread.h>

#include "fastq_file.h"
#include "fastq_read.h"
#include "containers/list.h"

/* **************************************
 *  		Structures		*
 * *************************************/

/**
* @brief Fastq batch reader 
* 
* Structure for implementing a fastq batch reader server
*/
typedef struct fastq_batch_reader {
  int alive;				/**< Alive indicator of the server. */
  pthread_mutex_t alive_lock;		/**< Alive variable lock. */

  int eof;				/**< End Of File indicator. */
  pthread_mutex_t eof_lock;		/**< End Of File indicator lock. */

  size_t batch_size;			/**< Fastq batch size (in bytes). */
  int batch_list_max_length;		/**< Fastq batch list maximum length. */

  pthread_t thread;			/**< pthread element. */

  int source_id;			/**< Source id of the fastq file (pair1 or pair2). */

  fastq_file_t* fastq_file_p;		/**< Pointer to the fastq file object. */
  list_t* batch_list_p;			/**< Pointer to the fastq batch list. */
} fastq_batch_reader_t;

/* **************************************
 *  		Functions		*
 * *************************************/

/**
*  @brief Creates a new fastq batch reader
*  @param filename name of the fastq file
*  @param source_id source id of the fastq file (pair1 or pair2)
*  @param batch_list_p pointer to the fastq batch list
*  @param batch_size size of the fastq batch
*  @param batch_list_max_length Fastq batch list maximum length
*  @return fastq_batch_reader_t pointer to the reader
*  
*  Creates a new fastq batch reader structure
*/
fastq_batch_reader_t* fastq_batch_reader_new(char* filename, int source_id, list_t* batch_list_p,  size_t batch_size, int batch_list_max_length);

/**
*  @brief Frees a fastq batch reader
*  @param[in,out] reader_p pointer to the fastq batch reader
*  @return void
*  
*  Frees a fastq batch reader 
*/
void fastq_batch_reader_free(fastq_batch_reader_t* reader_p);

/**
*  @brief pthread implementation of the fastq batch reader with OpenMP
*  @param param_p 
*  @return void
*  
*  Function with pthread implementation of the fastq batch reader with OpenMP
*/
void* fastq_batch_reader_thread_function(void* param_p);

#endif	/*  FASTQ_BATCH_READER_OMP_H  */

