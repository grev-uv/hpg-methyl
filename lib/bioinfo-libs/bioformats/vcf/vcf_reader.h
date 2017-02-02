#ifndef VCF_READER_H
#define VCF_READER_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include <commons/file_utils.h>
#include <commons/log.h>
#include <containers/list.h>

#include "vcf_util.h"
#include "vcf_file_structure.h"

/** 
 * @cond PRIVATE
 * @brief Bytes that the ZLIB library reads in a chunk
 */
#define CHUNK 0x80000

extern int mmap_vcf;
/** @endcond */

/**
 * @file vcf_reader.h
 * @author Cristina Yenyxe Gonzalez Garcia
 * @brief Reader and parser of VCF files
 * @details This file includes functions for reading VCF files. Files are read in blocks called 
 * "batches" whose size can be specified in lines or bytes. Files can be stored plainly or compressed 
 * in GZIP format.
 * Since reading involves also parsing contents, this tasks are provided in different schemes. The
 * vcf_read_and_parse functions allow to read and parse in just one step. If these two steps must be performed 
 * separately, the vcf_light_read functions, then run_vcf_parser, must be invoked.
 */



/**
 * @brief Current status of the VCF parser
 * @details Current status of the VCF parser. It is defined by the current header entry or record it is processing. 
 * If the current line is a record, it must also keep track of the batch where it will be inserted.
 **/
typedef struct {
    vcf_header_entry_t *current_header_entry;   /**< The header entry being parsed at the moment. */
    vcf_record_t *current_record;               /**< The record being parsed at the moment. */
    vcf_batch_t *current_batch;                 /**< The batch where records will be stored until it is considered full. */
    
    size_t num_samples;                         /**< Number of sample subjects. */
    size_t num_records;                         /**< Number of records read. */
    size_t num_batches;                         /**< Number of batches read (used for error notification). */
} vcf_reader_status;


/**
 * @brief Allocates memory for a vcf_reader_status structure
 * @details Creates a new vcf_reader_status structure. The batch where the records will be inserted during parsing 
 * will have the initial length specified as argument, although it will be resized when neccessary. The 
 * user can also decide whether to store the value of each sample or ignore them for the sake of speed and 
 * memory usage.
 * 
 * @param batch_lines The initial number of lines in a batch
 * @return A new vcf_reader_status structure
 **/
vcf_reader_status *vcf_reader_status_new(size_t batch_lines, size_t current_batch_id);

/**
 * @brief Deallocates memory associated to a vcf_reader_status structure.
 * @details Deallocates memory associated to a vcf_reader_status structure.
 * 
 * @param status The structure to be freed
 **/
void vcf_reader_status_free(vcf_reader_status *status);


/**
 * @brief Read and parse blocks of the given number of lines from a VCF file
 * @details Read and parse the given number of lines from a VCF file. The data read will be stored in 
 * members of the vcf_file_t structure.
 * 
 * The user can also decide whether to parse and store the value of each sample or ignore them for 
 * the sake of speed and memory usage.
 * 
 * @param batch_lines The number of lines to read and parse
 * @param file The file the data will be read from
 * @return 0 if the file was successfully read and parsed, 1 otherwise
 **/
int vcf_read_and_parse(size_t batch_lines, vcf_file_t *file);

/**
 * @brief Read and parse blocks of the given number of bytes from a VCF file
 * @details Read and parse blocks of the given number of bytes from a VCF file. If the last byte does not correspond to a 
 * linebreak character, it keeps reading until next linebreak is found. The data read will be stored in 
 * members of the vcf_file_t structure.
 * 
 * The user can also decide whether to parse and store the value of each sample or ignore them for 
 * the sake of speed and memory usage.
 * 
 * @param batch_bytes The number of bytes to read and parse
 * @param file The file the data will be read from
 * @return 0 if the file was successfully read and parsed, 1 otherwise
 **/
int vcf_read_and_parse_bytes(size_t batch_bytes, vcf_file_t *file);

/**
 * @brief Read and parse blocks of the given number of lines from a gzipped VCF file
 * @details Read and parse blocks of the given number of lines from a gzipped VCF file. The data read will be stored in 
 * members of the vcf_file_t structure.
 * 
 * The user can also decide whether to parse and store the value of each sample or ignore them for 
 * the sake of speed and memory usage.
 * 
 * @param batch_lines The number of lines to read and parse
 * @param file The file the data will be read from
 * @return 0 if the file was successfully read and parsed, 1 otherwise
 **/
int vcf_gzip_read_and_parse(size_t batch_lines, vcf_file_t *file);

/**
 * @brief Read and parse blocks of the given number of bytes from a gzipped VCF file
 * @details Read and parse blocks of the given number of bytes from a gzipped VCF file. If the last byte does not correspond to a 
 * linebreak character, it keeps reading until next linebreak is found. The data read will be stored in 
 * members of the vcf_file_t structure.
 * 
 * The user can also decide whether to parse and store the value of each sample or ignore them for 
 * the sake of speed and memory usage.
 * 
 * @param batch_bytes The number of bytes to read and parse
 * @param file The file the data will be read from
 * @return 0 if the file was successfully read and parsed, 1 otherwise
 **/
int vcf_gzip_read_and_parse_bytes(size_t batch_bytes, vcf_file_t *file);


/**
 * @brief Read (without parsing) blocks of the given number of lines from a VCF file
 * @details Read blocks of the given number of lines from a VCF file. The data read will be stored in 
 * members of the vcf_file_t structure.
 *
 * @param batch_lines The number of lines to read
 * @param file The file the data will be read from
 * @return 0 if the file was successfully read, 1 otherwise
 **/
int vcf_light_read(size_t batch_lines, vcf_file_t *file);

/**
 * @brief Read (without parsing) blocks of the given number of bytes from a VCF file
 * @details Read blocks of the given number of bytes from a VCF file. The data read will be stored in 
 * members of the vcf_file_t structure.
 *
 * @param batch_bytes The number of bytes to read
 * @param file The file the data will be read from
 * @return 0 if the file was successfully read, 1 otherwise
 **/
int vcf_light_read_bytes(size_t batch_bytes, vcf_file_t *file);

/**
 * @brief Read (without parsing) blocks of the given number of lines from a gzipped VCF file
 * @details Read blocks of the given number of lines from a gzipped VCF file. The data read will be stored in 
 * members of the vcf_file_t structure.
 * 
 * @param batch_lines The number of lines to read
 * @param file The file the data will be read from
 * @return 0 if the file was successfully read, 1 otherwise
 **/
int vcf_gzip_light_read(size_t batch_lines, vcf_file_t *file);

/**
 * @brief Read (without parsing) blocks of the given number of bytes from a gzipped VCF file
 * @details Read blocks of the given number of bytes from a gzipped VCF file. The data read will be stored in 
 * members of the vcf_file_t structure.
 * 
 * @param batch_bytes The number of bytes to read
 * @param file The file the data will be read from
 * @return 0 if the file was successfully read, 1 otherwise
 **/
int vcf_gzip_light_read_bytes(size_t batch_bytes, vcf_file_t *file);

/**
 * @brief Read (without parsing) blocks of the given number of lines from multiple VCF files
 * @details Read blocks of the given number of lines from multiple VCF files. The data read will be stored in per-file lists 
 * that must be consumed from the outside.
 *
 * @param text_lists [out] Lists of blocks where the current text block from each file will be queued
 * @param batch_lines The number of lines to read
 * @param files The files the data will be read from
 * @param num_files The number of files to read from
 * @return 0 if the files were successfully read, 1 otherwise
 **/
int vcf_light_multiread(list_t **text_lists, size_t batch_lines, vcf_file_t **files, size_t num_files);

/**
 * @brief Parses the text from an input buffer and links the data to the VCF file given as argument
 * @details Parses the text from an input buffer and links the data to the VCF file given as argument. When finished, 
 * the VCF file will contain a list of header entries, and a list of batches, each of them with a list of 
 * records inside.
 * 
 * @param p Pointer to the beginning of the input buffer
 * @param pe Pointer to the end of the input buffer
 * @param batch_size The number of lines of a batch, or zero when its size is specified in bytes
 * @param file VCF file whose contents will be read
 * @param status Structure that stores the current status of the parser
 */
int run_vcf_parser(char *p, char *pe, size_t batch_size, vcf_file_t *file, vcf_reader_status *status);

/** @cond PRIVATE */
size_t consume_input(int c, char **data, size_t max_len, int position_in_data);
/** @endcond */

#endif
