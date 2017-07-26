
#ifndef FASTQ_EX_BATCH_H
#define FASTQ_EX_BATCH_H

#include <stdio.h>

#include "fastq_file.h"

/* **************************************
 *  		Structures		*
 * *************************************/

/**
* @brief Fastq ex batch
* 
* Structure containing fastq reads information to perform future processes over the data
*/
typedef struct fastq_ex_batch {
    int source_id;				/**< Source id (pair1 or pair2). */ 
    unsigned long num_reads;			/**< Number of reads in the batch. */
    unsigned long data_indices_size;		/**< Data indices size in bytes (data_indices_size = seq_indices_size = quality_indices_size). */  
    unsigned long data_size;			/**< Data size in bytes (data_size = seq_size = quality_size). */
    int *data_indices;				/**< Data indices for quality and sequence indexation (data_indices = seq_indices = quality_indices). */
    int *header_indices;			/**< Header indices for header indexation. */
    char *header;				/**< Header with read headers. */
    char *seq;					/**< Pointer to sequences vector. */
    char *quality;				/**< Pointer to qualities vector. */
} fastq_ex_batch_t;

/* **************************************
 *  		Functions		*
 * *************************************/

/**
*  @brief Creates a fastq_ex_batch_t structure
*  @return pointer to created fastq_ex_batch_t structure
*  
*  Creates a fastq_ex_batch_t structure
*/
fastq_ex_batch_t* fastq_ex_batch_new(unsigned long size);

/**
*  @brief Frees fastq ex batch 
*  @param fastq_batch_p pointer to the fastq ex batch to be freed
*  @return void
*  
*  Free fastq ex batch structure
*/
void fastq_ex_batch_free(fastq_ex_batch_t* fastq_ex_batch_p);

/**
*  @brief Reads a fastq ex batch from a file 
*  @param[in,out] fastq_batch_p pointer to the fastq ex batch to fill
*  @param encode flag that indicates if nucleotides are encoded (0: not encoded, 1: encoded)
*  @param fastq_file_p pointer to the fastq file handler
*  @return void
*  
*  Reads a fastq file and fills the given fastq ex batch
*/
int fastq_ex_batch_read(fastq_ex_batch_t* batch_p, int encode, fastq_file_t* file_p);

/**
*  @brief Writes a fastq ex batch to a file 
*  @param[in,out] fastq_batch_p pointer to the fastq ex batch to write
*  @param encode flag that indicates if nucleotides are encoded (0: not encoded, 1: encoded)
*  @param fastq_file_p pointer to the fastq file handler
*  @return void
*  
*  Writes a fastq file with the content of the given fastq ex batch
*/
int fastq_ex_batch_write(fastq_ex_batch_t* batch_p, int decode, fastq_file_t* file_p);

#endif  /* FASTQ_EX_BATCH_H */ 
