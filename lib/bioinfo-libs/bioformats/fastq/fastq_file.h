#ifndef FASTQ_FILE_H
#define FASTQ_FILE_H

#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include <assert.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>

#include "commons/commons.h"
#include "commons/file_utils.h"
#include "commons/string_utils.h"

#include "containers/array_list.h"

//#include "qc_batch.h"
#include "fastq_read.h"
#include "fastq_batch.h"


#define MAX_FASTQ_FILENAME_LENGTH		1024	// Maximum filenname length
#define MAX_READ_ID_LENGTH			2048	// Maximum read ID length
#define MAX_READ_SEQUENCE_LENGTH		8192	// Maximum read sequence length

#define MAX_NUM_PRODUCERS			10

#define FQ_SEEK_BEGIN 	0
#define FQ_SEEK_CURR  	1
#define FQ_SEEK_RND   	2

#define FASTQ_FILE_PAIRED_END_MODE	1
#define FASTQ_FILE_MATE_END_MODE	2

/* **************************************
 *  		Structures		*
 * *************************************/

/**
* @brief Fastq file structure 
* 
* Structure for handling fastq files 
*/
typedef struct fastq_file {
    char *filename;				/**< Fastq file name. */
    char *mode;					/**< Opening mode ("r", "w"). */
    char *quality_encoding;		/**< Quality encoding (Illumina v1.5, Solid, ...). */

    FILE *fd;					/**< File descriptor. */

    size_t num_reads;			/**< Number of reads in the fastq file. */
    size_t num_lines;			/**< Number of lines in the fastq file. */
    //int source_id;
} fastq_file_t;

/**
* @brief Source structure 
* 
* Source structure containing id and filename
*/
typedef struct source {
    int id;					/**< Id of the source. */
    char filename[MAX_FASTQ_FILENAME_LENGTH];	/**< File name. */
} source_t;


/**
*  @brief Creates a fastq file handler and opens the file for read
*  @param filename file name of the handled fastq file
*  @return fastq_file_t pointer to a fastq handler
*  
*  Creates a fastq file handler and opens the fastq file for read
*/
fastq_file_t *fastq_fopen(char *filename);

/**
*  @brief Creates a fastq file handler and opens the file
*  @param filename file name of the handled fastq file
*  @return fastq_file_t pointer to a fastq handler
*  
*  Creates a fastq file handler and opens the fastq file in the specified mode ("r", "w")
*/
fastq_file_t *fastq_fopen_mode(char *filename, char *mode);


/*
 * SINGLE-END READ FUNCTIONS
 */
size_t fastq_fread_se(array_list_t *reads, size_t num_reads, fastq_file_t *fq_file);

size_t fastq_fread_bytes_se(array_list_t *reads, size_t bytes, fastq_file_t *fq_file);

//size_t fastq_gzread_se(array_list_t *reads, size_t num_reads, fastq_file_t *fq_file);
//
//size_t fastq_gzread_bytes_se(array_list_t *reads, size_t bytes, fastq_file_t *fq_file);


/*
 * PAIRED-END READ FUNCTIONS
 */
size_t fastq_fread_pe(array_list_t *reads, size_t num_reads, fastq_file_t *fq_file1, fastq_file_t *fq_file2, int mode);

size_t fastq_fread_bytes_pe(array_list_t *reads, size_t bytes, fastq_file_t *fq_file1, fastq_file_t *fq_file2, int mode);

size_t fastq_fread_bytes_aligner_pe(array_list_t *reads, size_t bytes, fastq_file_t *fq_file1, fastq_file_t *fq_file2);
//size_t fastq_gzread_pe(array_list_t *reads, size_t num_reads, fastq_file_t *fq_file1, fastq_file_t *fq_file2);



/**
*  @brief Reads fastq reads from a fastq file
*  @param read pointer to fastq_read_t structure
*  @param fq_file pointer to the fastq file handler
*  @return number of reads done
*  
*  Reads a single fastq read from a fastq file 
*/
int fastq_fread(fastq_read_t *read, fastq_file_t *fq_file);

/**
*  @brief Reads a given number of fastq reads from a fastq file
*  @param buffer_reads pointer to fastq_read_t structure
*  @param num_reads number of fastq reads to read
*  @param fq_file pointer to the fastq file handler
*  @return number of fastq reads read
*  
*  Reads a given number of fastq reads from a fastq file
*/
int fastq_fread_num_reads(fastq_read_t* buffer_reads, int num_reads, fastq_file_t *fq_file);

/**
*  @brief Reads fastq reads from a fastq file until a given data size is reached
*  @param buffer_reads pointer to fastq_read_t structure where reads are stored
*  @param max_size max size in bytes to read from fastq file
*  @param fq_file pointer to the fastq file handler
*  @return number of fastq reads read
*  
*  Reads fastq reads from a fastq file until a given data size is reached
*/
int fastq_fread_max_size(fastq_read_t *buffer_fq_reads, unsigned long max_size, fastq_file_t *fq_file);

/**
*  @brief Reads fastq reads from a fastq file and stores them in a batch until a given data size is reached
*  @param buffer_fq_read_batch pointer to fastq read batch where reads are stored
*  @param max_size max size in bytes to read from fastq file
*  @param fq_file pointer to the fastq file handler
*  @return number of reads stored in the fastq read batch
*  
*  Reads fastq reads from a fastq file until a given data size is reached and stores 
*  them in a fastq read batch
*/
int fastq_fread_batch_max_size(fastq_batch_t *fq_batch, unsigned long max_size, fastq_file_t *fq_file);
int fastq_fread_paired_batch_max_size(fastq_batch_t *fq_batch, unsigned long max_size, 
				      fastq_file_t *fq_file);
int fastq_fread_paired_batch_max_size2(fastq_batch_t *fq_batch, unsigned long max_size, 
				       fastq_file_t *fq_file, fastq_file_t *fq_file2);

/**
*  @brief Reads fastq reads in the given positions from a fastq file
*  @param buffer_reads pointer to fastq_read_t structure where reads are stored
*  @param index_positions vector indication positions in the fastq file to be read
*  @param fq_file pointer to the fastq file handler
*  @return number of fastq reads read
*  
*  Reads fastq reads in the given positions from a fastq file. This method is useful 
*  for read specific reads that are in known positions in the fastq file
*/
int fastq_fread_index_positions(fastq_read_t* buffer_reads, int *index_positions, fastq_file_t *fq_file);

/**
*  @brief Writes reads stored in a buffer to file 
*  @param buffer_reads pointer to fastq reads buffer
*  @param num_writes number of reads to write to file
*  @param fq_file pointer to the fastq file handler
*  @return number of written reads
*  
*  Writes reads stored in a buffer to file using the given file handler
*/
int fastq_fwrite(fastq_read_t* buffer_reads, int num_writes, fastq_file_t *fq_file);

/**
*  @brief Returns the number of reads of a fastq file
*  @param fq_file pointer to the fastq file handler
*  @return number of reads in the file
*  
*  Returns the number of reads of a fastq file. They are counted
*  during read process, not at the time of this method call
*/
unsigned int fastq_fcount(fastq_file_t *fq_file);

/**
*  @brief Removes those reads with more N than max_N_per_read parameter
*  @param buffer_reads pointer to fastq reads buffer
*  @param qc_read pointer with qc read structure with information about number of N in read
*  @param max_N_per_read maximum number of N allowed in a read
*  @return void
*  
*  Removes from the buffer_reads those reads with more N than max_N_per_read parameter.
*  N is the representation of an undefined nucleotide
*/
//void fastq_remove_Ns(fastq_read_t* buffer_reads, qc_read_t* qc_read, int max_N_per_read);

/**
*  @brief Closes a fastq file handler and the associated file
*  @param fq_file pointer to the fastq file handler
*  @return void
*  
*  Closes a fastq file handler and the associated file
*/

void fastq_fclose(fastq_file_t *fq_file);


/** @cond PRIVATE */
//size_t consume_input(int c, char **data, size_t max_len, int position_in_data);

#endif	/*  FASTQ_FILE_H  */
