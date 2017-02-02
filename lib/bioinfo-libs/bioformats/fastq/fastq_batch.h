#ifndef FASTQ_BATCH_H
#define FASTQ_BATCH_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "commons/string_utils.h"


/* **************************************
 *  		Structures		*
 * *************************************/

/**
* @brief Fastq batch
* 
* Structure containing fastq reads information to perform future processes over the data
*/
typedef struct fastq_batch {
    unsigned long num_reads;			/**< Number of reads in the batch. */
    unsigned long data_indices_size;		/**< Data indices size in bytes. */
    unsigned long data_size;			/**< Data size in bytes. */
    int source_id;				/**< Source id (pair1 or pair2). */
    int* header_indices;			/**< Header indices for header indexation. */
    char* header;				/**< Header with read headers. */
    int* data_indices;				/**< Data indices for quality and sequence indexation. */
    char *seq;					/**< Pointer to sequences vector. */
    char *quality;				/**< Pointer to qualities vector. */
} fastq_batch_t;

/* **************************************
 *  		Functions		*
 * *************************************/

/**
*  @brief Creates a fastq_batch_t structure
*  @return pointer to created fastq_batch_t structure
*  
*  Creates a fastq_batch_t structure
*/
fastq_batch_t* fastq_batch_new(unsigned long size);

/**
*  @brief Inits fastq_batch_t structure
*  @param[in,out] fastq_batch_p pointer to the fastq batch to initialize
*  @param size size of the fastq_batch_p data
*  @return void
*  
*  Initializes and returns fastq_batch structure
*/
//fastq_batch_t* fastq_batch_init(fastq_batch_t *fastq_batch_p, unsigned long size);

/**
*  @brief Frees fastq batch 
*  @param fastq_batch_p pointer to the fastq batch to be freed
*  @return void
*  
*  Free fastq batch structure
*/
void fastq_batch_free(fastq_batch_t* fastq_batch_p);

/**
*  @brief Prints a fastq batch to disk 
*  @param fastq_batch_p pointer to the fastq batch
*  @param fd write file descriptor
*  @return void
*  
*  Prints the content of a fastq batch to disk in the specified file
*/
void fastq_batch_print(fastq_batch_t* fastq_batch_p, FILE* fd);

/**
*  @brief Prints a read to disk
*  @param fd write file descriptor
*  @param batch_p pointer to fastq batch 
*  @param index position of the read to write within the sequence vector
*  @return void
*  
*  Prints a read to disk without preprocessing
*/
void fprintf_read(FILE* fd, fastq_batch_t* batch_p, int index);

/**
*  @brief Prints a right-side preprocessed read to disk
*  @param fd write file descriptor
*  @param batch_p pointer to fastq batch 
*  @param index position of the read to write within the sequence vector
*  @param rtrim_length number of right nucleotides to trim
*  @return void
*  
*  Prints a read to disk with right sides preprocessed
*/
void fprintf_rtrim_read(FILE* fd, fastq_batch_t* batch_p, int index, int rtrim_length);

/**
*  @brief Prints a left-side preprocessed read to disk
*  @param fd write file descriptor
*  @param batch_p pointer to fastq batch 
*  @param index position of the read to write within the sequence vector
*  @param ltrim_length number of left nucleotides to trim
*  @return void
*  
*  Prints a read to disk with left sides preprocessed
*/
void fprintf_ltrim_read(FILE* fd, fastq_batch_t* batch_p, int index, int ltrim_length);

/**
*  @brief Prints a both-side preprocessed read to disk
*  @param fd write file descriptor
*  @param batch_p pointer to fastq batch 
*  @param index position of the read to write within the sequence vector
*  @param rtrim_length number of right nucleotides to trim
*  @param ltrim_length number of left nucleotides to trim
*  @return void
*  
*  Prints a read to disk with left and right sides preprocessed
*/
void fprintf_trim_read(FILE* fd, fastq_batch_t* batch_p, int index, int rtrim_length, int ltrim_length);

#endif  /* FASTQ_BATCH_H */ 
