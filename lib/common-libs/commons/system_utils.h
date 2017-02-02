
#ifndef SYSTEM_UTILS_H
#define SYSTEM_UTILS_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "commons.h"

#ifdef CUDA_VERSION
#include "cuda_commons.h"
#endif

#include "log.h"

#define FASTQ_QC	1
#define FASTQ_PREPRO	2
#define BAM_QC		3

#define QC_MEMORY_USAGE_FACTOR 		4
#define PREPRO_MEMORY_USAGE_FACTOR	6

#define BAM_BATCH_SIZE_TO_FREE_MEMORY_RATIO	20

#define BAM_SAM_COMPRESSION_RATIO	4
#define MEAN_COMPRESSED_ALIGNMENT_SIZE 	100

/* **************************************
 *  		Functions		*
 * *************************************/

/**
*  @brief Gets free memory from the machine (CPU memory)
*  @return unsigned long int - quantity of free memory
*  
*  Gets free memory from the machine (CPU memory)
*/
unsigned long int get_free_memory();

/**
*  @brief Gets an estimation of the memory requirement for a given process (CPU memory)
*  @param process process which is to estimate the memory
*  @param batch_size batch size of fastq reads to load from disk
*  @param max_list_length maximum length of the batch lists
*  @return unsigned long int - estimated memory
*  
*  Gets free memory from the machine (CPU memory)
*/
unsigned long int get_estimated_memory_needed(int process, int batch_size, int max_list_length);

/**
*  @brief Gets an estimation of the maximum number of alignments by chromosome for a given bam file
*  @param input_filename bam file name 
*  @return maximum number of alignments by chromosome estimated
*  
*  Gets an estimation of the maximum number of alignments by chromosome for a given bam file
*/
int get_max_estimated_alignments_by_chromosome(char* input_filename);

/**
*  @brief Estimates an optimal CPU threads number
*  @return size_t - number of threads
*  
*  Gets the optimal number of cpu threads
*/
size_t get_optimal_cpu_num_threads();

/**
*  @brief Estimates optimal GPU threads by block
*  @return optimal GPU threads by block
*  
*  Estimates optimal GPU threads by block
*  optimal_gpu_num_threads = 16 * cuda_device_warp_size = 512 (Fermi GPU)
*/
int get_optimal_gpu_num_threads();

/**
*  @brief Calculates the optimal batch size for a given process and list size
*  @param process process that will process the batch
*  @param max_list_length maximum length of the batch lists
*  @return optimal batch size
*  
*  Calculates the optimal batch size for a given process and list size
*  The optimal batch size is calculated to fill almost all free memory and 
*  taking into account GPU memory space
*/
int get_optimal_batch_size(int process, int max_list_length);

#endif	/*  SYSTEM_UTILS_H   */
