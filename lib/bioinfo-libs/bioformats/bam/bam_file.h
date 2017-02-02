#ifndef BAM_FILE_H
#define BAM_FILE_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#if defined THRUST-GPU
  #include <thrust/host_vector.h>
  #include <thrust/device_vector.h>
  #include <thrust/sort.h>
  #include <thrust/copy.h>
#endif 



#include "commons/commons.h"
#include "commons/file_utils.h"
#include "commons/string_utils.h"
#include "commons/system_utils.h"

#include "samtools/bam.h"

#include "alignment.h"
//#include "bam_commons.h"

#define MAX_NUM_PRODUCERS  	10
#define BAM_BATCH_READ_SIZE   	50000000
#define BAM_BATCH_WRITE_SIZE   	50000000

#define SINGLE_CHROM_BATCH 	0
#define MULTIPLE_CHROM_BATCH 	1

#define TEMPORARY_HEADER_PATH	"/tmp/header.tmp.bam"

/* **************************************
 *      	Structures    		*
 * *************************************/

/**
* @brief BAM file handler
*
* BAM file handler
*/
typedef struct bam_file {
    unsigned int num_alignments;	/**< Number of alignments in the file. */
    char* filename;			/**< BAM file name. */
    char* mode;				/**< Open mode ("r", "w"). */
    bamFile bam_fd;			/**< BAM file descriptor. */
    bam_header_t* bam_header_p;		/**< Pointer to the bam_header. */
} bam_file_t;

/**
* @brief Batch of bam1_t
*
* Batch of bam1_t alignments
*/
typedef struct bam_batch {
    int type;				/**< SINGLE_CHROM_BATCH or MULTIPLE_CHROM_BATCH. */
    int num_alignments;			/**< Number of alignments in the batch. */

    int allocated_alignments;		/**< Number of allocated alignments (for reserved memory). */
    bam1_t** alignments_p;		/**< Pointers to bam1_t alignments. */
} bam_batch_t;

/* **************************************
 *      	Functions    		*
 * *************************************/

/**
*  @brief Opens a BAM file for reading and returns its handler
*  @param filename BAM file name
*  @return pointer to the BAM file handler
*  
*  Opens a BAM file and returns its handler
*/
bam_file_t* bam_fopen(char* filename);

/**
*  @brief Opens a BAM file in the specified mode and returns its handler
*  @param filename BAM file name
*  @param bam_header_p pointer to the BAM file header
*  @param mode opening mode ("r", "w")
*  @return pointer to the BAM file handler
*  
*  Opens a BAM file in the specified mode and returns its handler
*/
bam_file_t* bam_fopen_mode(char* filename, bam_header_t* bam_header_p, char* mode);

/**
*  @brief Closes a BAM file and frees its handler
*  @param bam_file[in,out] bam file handler
*  @return void
*  
*  Closes a BAM file and frees its handler
*/
void bam_fclose(bam_file_t *bam_file);

/* **********************************************
 *      	BAM read functions    		*
 * *********************************************/

/**
*  @brief Reads a bam batch of a given size
*  @param batch_p pointer to the bam batch to fill
*  @param batch_size max. size of the batch
*  @param base_quality base quality for quality values normalization
*  @param bam_file_p bam file handler
*  @return void
*  
*  Reads a bam batch of a given size
*/
int bam_fread_max_size(bam_batch_t* batch_p, size_t batch_size, int base_quality, bam_file_t* bam_file_p);

/**
*  @brief Reads a bam batch of a given size with no duplicated alignments
*  @param batch_p pointer to the bam batch to fill
*  @param batch_size max. size of the batch
*  @param base_quality base quality for quality values normalization
*  @param bam_file_p bam file handler
*  @param prev_seq last sequence of the previous batch
*  @param prev_seq_length length of the last sequence of the previous batch
*  @param prev_seq_start_coordinate start coordinate of the last sequence of the previous batch
*  @return void
*  
*  Reads a bam batch of a given size with no duplicated alignments
*  The function assumes that the file is sorted by position and then
*  only compares consecutive alignments
*/
int bam_fread_max_size_no_duplicates(bam_batch_t* batch_p, size_t batch_size, int base_quality, bam_file_t* bam_file_p, uint8_t* prev_seq, int* prev_seq_length, int* prev_seq_start_coordinate);

/**
*  @brief Reads a bam batch of a given size only with alignments from one chromosome
*  @param batch_p pointer to the bam batch to fill
*  @param batch_size max. size of the batch
*  @param chromosome chromosome of the batch (for filtering bam1_t alignments)
*  @param bam_file_p bam file handler
*  @return void
*  
*  Reads a bam batch of a given size only with alignments from one chromosome
*/
int bam_fread_max_size_by_chromosome(bam_batch_t* batch_p, size_t batch_size, int chromosome, bam_file_t* bam_file_p);

/* **********************************************
 *      	BAM write functions    		*
 * *********************************************/

/**
*  @brief Writes header content to BAM file
*  @param bam_header_p pointer to the bam_header
*  @param bam_file_p bam file handler 
*  @return void
*  
*  Writes header content to BAM file
*/
void bam_fwrite_header(bam_header_t* bam_header_p, bam_file_t* bam_file_p);

/**
*  @brief Writes header in a temporary file for further recovery
*  @param bam_header_p pointer to the bam_header
*  @return void
*  
*  Writes header in a temporary file for further recovery (TEMPORARY_HEADER_PATH)
*/
void bam_fwrite_temporary_header(bam_header_t* bam_header_p);

/**
*  @brief Reads header from a temporary file
*  @return bam_header_t* pointer to the bam_header
*  
*  Reads header form a temporary file (TEMPORARY_HEADER_PATH)
*/
bam_header_t* bam_fread_temporary_header();

/**
*  @brief Writes a bam1_t alignment to a BAM file
*  @param alignment_p pointer to bam1_t alignment
*  @param bam_file_p bam file handler 
*  @return number of bytes written to file
*  
*  Writes a bam1_t alignment to a BAM file
*/
int bam_fwrite(bam1_t* alignment_p, bam_file_t* bam_file_p);

/**
*  @brief Writes an array of bam1_t alignments to a BAM file
*  @param alignment_p pointers to bam1_t alignments
*  @param length length of the array (number of alignments)
*  @param bam_file_p bam file handler 
*  @return number of bytes written to file
*  
*  Writes an array of bam1_t alignments to a BAM file
*/
int bam_fwrite_array(bam1_t** alignment_p, int length, bam_file_t* bam_file_p);

/**
*  @brief Writes a batch of bam1_t alignments to a BAM file
*  @param batch_p pointer to the bam_batch
*  @param bam_file_p bam file handler 
*  @return number of bytes written to file
*  
*  Writes a batch of bam1_t alignments to a BAM file
*/
int bam_fwrite_batch(bam_batch_t* batch_p, bam_file_t* bam_file_p);

/**
*  @brief Writes a sorted array of bam1_t alignments to a BAM file
*  @param alignment_p pointers to bam1_t alignments
*  @param indices vector of sorted indices indicating write order
*  @param length length of the array (number of alignments)
*  @param bam_file_p bam file handler 
*  @return number of bytes written to file
*  
*  Writes a sorted array of bam1_t alignments to a BAM file
*/
#if defined THRUST-GPU

int bam_fwrite_sorted_array(bam1_t** alignment_p, thrust::host_vector<int> indices, int length, bam_file_t* bam_file_p);

#else

int bam_fwrite_sorted_array(bam1_t** alignment_p, int* indices, int length, bam_file_t* bam_file_p);

#endif

/**
*  @brief Writes an alignment to a BAM file
*  @param alignment_p pointer to the alignment
*  @param bam_file_p bam file handler 
*  @return number of bytes written to file
*  
*  Writes an alignment to a BAM file
*/
int alignment_fwrite(alignment_t* alignment_p, bam_file_t* bam_file_p);

/**
*  @brief Writes an array of alignments to a BAM file
*  @param alignment_p pointers to alignments
*  @param length length of the array (number of alignments)
*  @param bam_file_p bam file handler 
*  @return number of bytes written to file
*  
*  Writes an array of alignments to a BAM file
*/
int alignment_fwrite_array(alignment_t** alignment_p, int length, bam_file_t* bam_file_p);

/**
*  @brief Writes a batch of alignments to a BAM file
*  @param batch_p pointer to the alignment_batch
*  @param bam_file_p bam file handler 
*  @return number of bytes written to file
*  
*  Writes a batch of alignments to a BAM file
*/
int alignment_fwrite_batch(alignment_batch_t* batch_p, bam_file_t* bam_file_p);

/**
*  @brief Validates the existence of a valid header
*  @param bam_file_p bam file handler 
*  @return 1 if validated, 0 if not validated
*  
*  Validates the existence of a valid header
*/
int bam_validate_header(bam_file_t* bam_file_p);

/**
*  @brief Counts the number of alignments in the BAM file
*  @param bam_file bam file handler 
*  @return number of alignments in the BAM file
*  
*  Counts the number of alignments in the BAM file
*/
unsigned int bam_fcount(bam_file_t* bam_file);

/**
*  @brief Gets number of chromosomes/targets declared in BAM file header
*  @param filename BAM file name
*  @return number of chromosomes/targets declared in BAM file header
*  
*  Gets number of chromosomes/targets declared in BAM file header
*/
int bam_fread_num_chromosomes(char* filename);

/* **********************************************
 *      	BAM batch functions    		*
 * *********************************************/

/**
*  @brief Creates and returns a bam batch
*  @param batch_size size of the batch (for memory allocation)
*  @param type SINGLE_CHROM_BATCH or MULTIPLE_CHROM_BATCH
*  @return pointer to the created bam batch
*  
*  Creates and returns a bam batch
*/
bam_batch_t* bam_batch_new(size_t batch_size, int type);

/**
*  @brief Frees a bam batch
*  @param batch_p pointer to the bam_batch
*  @param free_alignments flag to free the inner alignments pointed by the batch
*  @return void
*  
*  Frees a bam batch
*/
void bam_batch_free(bam_batch_t* batch_p, int free_alignments);

/**
*  @brief Writes a bam batch to a file
*  @param batch_p pointer to the bam_batch
*  @param fd file descriptor
*  @return void
*  
*  Writes a bam batch to a file
*/
void bam_batch_print(bam_batch_t* batch_p, FILE* fd);

/**
*  @brief Compares two sequences
*  @param data1 sequence 1 in uint8_t format
*  @param length_seq1 length of sequence 1
*  @param start_seq1 start position of the mapped sequence 1
*  @param data2 sequence 2 in uint8_t format
*  @param length_seq2 length of sequence 2
*  @param start_seq2 start position of the mapped sequence 2 
*  @return 1 if equal, 0 if not equal
*  
*  Compares two sequences, used to prevent duplicates in the batch
*/
int bam_batch_compare_seq(uint8_t* data1, int length_seq1, int start_seq1, uint8_t* data2, int length_seq2, int start_seq2);

/* **********************************************
 *      	bam1_t functions    		*
 * *********************************************/

/**
*  @brief Frees an array of bam1_t alignments
*  @param[in,out] alignment_p pointer to the bam1_t alignments
*  @param num_alignments number of alignments in the array
*  @return void
*  
*  Frees an array of bam1_t alignments
*/
void free_bam1(bam1_t** alignment_p, int num_alignments);

/**
*  @brief Writes to disk a given bam1_t alignment
*  @param alignment_p pointer to the bam1_t alignment
*  @param fd file descriptor
*  @return void
*  
*  Writes to disk a given bam1_t alignment
*/
void print_bam1(bam1_t* alignment_p, FILE* fd);

#endif  /* BAM_FILE_H */
