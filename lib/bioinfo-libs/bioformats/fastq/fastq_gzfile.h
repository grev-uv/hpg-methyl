#ifndef FASTQ_GZFILE_H
#define FASTQ_GZFILE_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
#include <assert.h>

#include "containers/array_list.h"
#include "fastq_read.h"

#define MAX_FASTQ_FILENAME_LENGTH		64		// Maximum filenname length
#define MAX_READ_ID_LENGTH				256		// Maximum read ID length
#define MAX_READ_SEQUENCE_LENGTH		2048	// Maximum read sequence length

#define CHUNK 0x80000

/* **************************************
 *  		Structures		*
 * *************************************/

/**
* @brief Fastq file structure
*
* Structure for handling fastq files
*/
typedef struct fastq_gzfile {
    char *filename;				/**< Fastq file name. */
    char *mode;					/**< Opening mode ("r", "w"). */
    char *quality_encoding;		/**< Quality encoding (Illumina v1.5, Solid, ...). */

    FILE *fd;					/**< File descriptor. */
    z_stream strm;
    int ret;
    char *data;
    size_t data_size;

    size_t num_reads;			/**< Number of reads in the fastq file. */
    size_t num_lines;			/**< Number of lines in the fastq file. */
} fastq_gzfile_t;


/**
*  @brief Creates a fastq gzip file handler and opens the file for read
*  @param filename file name of the handled fastq gzip file
*  @return fastq_file_t pointer to a fastq gzip handler
*
*  Creates a fastq gzip file handler and opens the fastq gzip file for read
*/
fastq_gzfile_t *fastq_gzopen(char *filename);

/**
*  @brief Closes a fastq gzip file handler and the associated gzip file
*  @param fq_file pointer to the fastq gzip file handler
*  @return void
*
*  Closes a fastq gzip file handler and the associated file
*/

void fastq_gzclose(fastq_gzfile_t *fq_file);


/*
 * SINGLE-END READ FUNCTIONS
 */
size_t fastq_gzread_se(array_list_t *reads, size_t num_reads, fastq_gzfile_t *fq_file);

size_t fastq_gzread_bytes_se(array_list_t *reads, size_t bytes, fastq_gzfile_t *fq_file);


/*
 * PAIRED-END READ FUNCTIONS
 */
size_t fastq_gzread_pe(array_list_t *reads, size_t num_reads, fastq_gzfile_t *fq_file1, fastq_gzfile_t *fq_file2);

size_t fastq_gzread_bytes_pe(array_list_t *reads, size_t bytes, fastq_gzfile_t *fq_file1, fastq_gzfile_t *fq_file2);



/** @cond PRIVATE */
static size_t consume_input(int c, char **data, size_t max_len, int position_in_data);

#endif	/*  FASTQ_GZFILE_H  */
