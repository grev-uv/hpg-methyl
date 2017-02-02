#ifndef VCF_FILE_STRUCTURE_H
#define VCF_FILE_STRUCTURE_H

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include <bioformats/ped/ped_file.h>
#include <commons/log.h>
#include <containers/array_list.h>
#include <containers/khash.h>
#include <containers/list.h>

#include "vcf_util.h"

/**
 * @file vcf_file_structure.h
 * @author Cristina Yenyxe Gonzalez Garcia
 * @brief Management of VCF files structure
 * @details This file includes the structs that specify the structure of a VCF file, and functions for creating 
 * and composing a new one. A VCF file contains a header with metadata and a body with the records that correspond 
 * to each remarkable position in the sequenced genomes. For the sake of efficiency, records are grouped in blocks 
 * called "batches".
 */



/** @cond PRIVATE */
extern int mmap_vcf;
/** @endcond */

/**
 * @brief Structure that specifies a VCF file.
 * @details Structure that specifies a VCF file. The physical file is defined by its file descriptor, 
 * its filename and the mode it has been open.
 *
 * It contains a header with several entries, and a body with several records.
 */
typedef struct vcf_file {
    char* filename;     /**< Name of the file to interact with */
    char* mode;         /**< Mode the file is open (w/r/a) */
    FILE *fd;           /**< Should the file be loaded using IO functions, the file descriptor is set */
    
    char *data;         /**< Should the file be loaded using mmap, its contents are set */
    size_t data_len;    /**< Length of the contents of the file, when apply */

    char* format;       /**< Format and version (set in the first line of the VCF file) */
    int format_len;     /**< Length of the format field */
    
    array_list_t *header_entries;   /**< Entries in the header of the file (metadata) */
    array_list_t *samples_names;    /**< Names of the sequenced samples */
    
    list_t *text_batches;           /**< Blocks of text read from the file */
    list_t *record_batches;         /**< Blocks of records who constitute the body/data of the file */
} vcf_file_t;

/**
 * @brief Entry in the VCF document header.
 * @details Entry in the header of the VCF file. An entry can be of the form:
 * 
 * - '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">' or
 * - '##reference=human_b36_both.fasta' or
 * - '##Some miscellaneous information'
 */
typedef struct vcf_header_entry {
    char *name;             /**< Key of entries with pairs (key,value) */
    int name_len;           /**< Length of the name field */
    array_list_t *values;   /**< List of values of the fields describing the entry */
} vcf_header_entry_t;


/**
 * @brief Entry in the VCF document body.
 * @details Entry in the body of the VCF file. Each entry corresponds with one line in the file.
 */
typedef struct vcf_record {
    char* chromosome;           /**< Chromosome where the remarkable variation took place */
    unsigned long position;     /**< Position in chromosome where the remarkable variation took place */
    char* id;                   /**< Unique identifier of the mutation, if existing */
    char* reference;            /**< Reference allele in that position */
    char* alternate;            /**< Alternate alleles found in that position */
    float quality;              /**< Quality of the reading */
    char* filter;               /**< PASS if all filters applied were passed, another value if one of the filters failed */
    char* info;                 /**< Miscellaneous information about the reading */
    char* format;               /**< Format of the sample data stored */

    int chromosome_len;         /**< Length of the chromosome field */
    int id_len;                 /**< Length of the ID field */
    int reference_len;          /**< Length of the reference allele field */
    int alternate_len;          /**< Length of the alternate alleles field */
    int filter_len;             /**< Length of the filter field */
    int info_len;               /**< Length of the info field */
    int format_len;             /**< Length of the format field */
    
    array_list_t *samples;      /**< Values of the samples in this position */
} vcf_record_t;


/**
 * @brief Block of records of a VCF file.
 * @details Block of records of a VCF file. When reading its contents from a physical file, 
 * the input buffer with the text is also stored, so no extra allocations need to be done.
 **/
typedef struct vcf_batch {
    array_list_t *records;      /**< Records in the block */
    char *text;                 /**< Input buffer with the data for the records */
} vcf_batch_t;



/* ********************************************************
 *      (De)Allocation of header entries and records      *
 * ********************************************************/

/**
 * @brief Creates a new header entry.
 * @details Creates a new header entry.
 *
 * @return The new header entry
 **/
vcf_header_entry_t* vcf_header_entry_new();

/**
 * @brief Deallocates memory of a header entry.
 * @details Deallocates memory of a header entry.
 *
 * @param header_entry The header entry to deallocate.
 **/
void vcf_header_entry_free(vcf_header_entry_t *header_entry);

/**
 * @brief Creates a new record.
 * @details Creates a new record.
 *
 * @return The new record
 **/
vcf_record_t* vcf_record_new();

/**
 * @brief Creates a copy of an existing record.
 * @details Creates a copy of an existing record. All its fields are replicated in new memory 
 * segments, so the original record can be deallocated without risk.
 *
 * @param orig The original record
 * @return The new record
 **/
vcf_record_t *vcf_record_copy(vcf_record_t *orig);

/**
 * @brief Deallocates memory of a record.
 * @details Deallocates memory of a record. Must be used when the record is loaded from an already 
 * existing file, because its contents will be stored in the "text" field of a vcf_batch.
 *
 * @param record The record to deallocate.
 * @see vcf_batch
 **/
void vcf_record_free(vcf_record_t *record);

/**
 * @brief Deallocates memory of a record and its fields.
 * @details Deallocates memory of a record and its fields. Must be used when the record is created 
 * by the user and/or when its contents are not stored in the "text" field of a vcf_batch.
 *
 * @param record The record to deallocate.
 **/
void vcf_record_free_deep(vcf_record_t *record);



/* ********************************************************
 *           Addition of header and record entries        *
 * ********************************************************/

/**
 * @brief Adds a header entry to a VCF file.
 * @details Adds a header entry to a VCF file.
 *
 * @param header_entry The header entry to add
 * @param file The VCF file to add the entry to
 * @return If the entry was successfully added
 **/
int add_vcf_header_entry(vcf_header_entry_t *header_entry, vcf_file_t *file);

/**
 * @brief Adds a sample subject to a VCF file.
 * @details Adds a sample subject to a VCF file.
 *
 * @param name The name of the sample subject
 * @param length The length of the name
 * @param file The VCF file to add the sample to
 * @return If the sample was successfully added
 **/
int add_vcf_sample_name(char *name, int length, vcf_file_t *file);

int add_text_batch(char *batch, vcf_file_t *file);

/**
 * @brief Adds a batch of records to a VCF file.
 * @details Adds a batch of records to a VCF file.
 *
 * @param batch The batch to add
 * @param file The VCF file to add the batch to
 * @return If the batch was successfully added
 **/
int add_vcf_batch(vcf_batch_t *batch, vcf_file_t *file);

// int add_record(vcf_record_t* record, vcf_file_t *vcf_file);

char *fetch_vcf_text_batch(vcf_file_t *file);

/**
 * @brief Removes and returns the first batch in the queue of a VCF file.
 * @details Removes and returns the first batch in the queue of a VCF file.
 *
 * @param file The file to fetch the batch from
 * @return The fetched batch
 **/
vcf_batch_t *fetch_vcf_batch(vcf_file_t *file);

/**
 * @brief Removes and returns the first batch in the queue of a VCF file.
 * @details Removes and returns the first batch in the queue of a VCF file, without blocking if no 
 * batches are available. This is the recommended function to use when reading from multiple files.
 *
 * @param file The file to fetch the batch from
 * @return The fetched batch
 **/
vcf_batch_t *fetch_vcf_batch_non_blocking(vcf_file_t *file);

/**
 * @brief Returns the number of header entries in a VCF file.
 * @details Returns the number of header entries in a VCF file.
 *
 * @param file The file to query
 * @return The number of header entries in the file
 **/
size_t get_num_vcf_header_entries(vcf_file_t *file);

/**
 * @brief Returns the number of values of a header entry in a VCF file.
 * @details Returns the number of values of a header entry in a VCF file.
 *
 * @param entry The entry to query
 * @return The number of values of the header entry
 **/
size_t get_num_values_in_vcf_header_entry(vcf_header_entry_t *entry);

/**
 * @brief Returns the number of sample subjects in a VCF file.
 * @details Returns the number of sample subjects in a VCF file.
 *
 * @param file The file to query
 * @return The number of sample subjects in the file
 **/
size_t get_num_vcf_samples(vcf_file_t *file);

/**
 * @brief Returns the number of records in a VCF file.
 * @details Returns the number of records in a VCF file.
 *
 * @param file The file to query
 * @return The number of records in the file
 **/
size_t get_num_vcf_records(vcf_file_t *file);

/**
 * @brief Returns the number of batches in a VCF file.
 * @details Returns the number of batches in a VCF file.
 *
 * @param file The file to query
 * @return The number of batches in the file
 **/
size_t get_num_vcf_batches(vcf_file_t *file);



/* ********************************************************
 *                    Batch management                    *
 * ********************************************************/

/**
 * @brief Creates a new batch with the given initial capacity.
 * @details Creates a new batch with the given initial capacity. If size <= 0, a capacity of 
 * 100 elements will be set by default.
 *
 * @param size The initial capacity of the batch
 * @return The new batch of records
 **/
vcf_batch_t* vcf_batch_new(size_t size);

/**
 * @brief Deallocates memory of a batch, its records and the associated buffer.
 * @details Deallocates memory of a batch and its records. The associated buffer will be 
 * deallocated if it has been set and the file has not been loaded by using the mmap function.
 *
 * @param batch The batch to deallocate.
 * @see mmap_vcf
 **/
void vcf_batch_free(vcf_batch_t *batch);

/**
 * @brief Adds a record to a batch.
 * @details Adds a record to a batch.
 *
 * @param record The record to add
 * @param batch The batch to add the record to
 **/
void add_record_to_vcf_batch(vcf_record_t *record, vcf_batch_t *batch);


/**
 * @brief Checks if a batch is empty.
 * @details Checks if a batch is empty (has zero elements inserted).
 *
 * @param batch The batch to check
 * @return 1 if the batch is empty, 0 otherwise
 **/
int vcf_batch_is_empty(vcf_batch_t *batch);

/**
 * @brief Checks if a batch is full.
 * @details Checks if a batch is full (its size is the same as its capacity).
 *
 * @param batch The batch to check
 * @return 1 if the batch is full, 0 otherwise
 **/
int vcf_batch_is_full(vcf_batch_t *batch);


/**
 * @brief Prints basic information about the batch to the given stream.
 * @details Prints basic information about the batch to the given stream. This information 
 * consists on the size and capacity of the batch and, if at least a record has been inserted,
 * its chromosome and position.
 *
 * @param fd The stream to write the information to
 * @param batch The batch whose information will be written
 * @return The number of bytes printed
 **/
int vcf_batch_print(FILE *fd, vcf_batch_t *batch);



/* ********************************************************
 *                    Header management                   *
 * ********************************************************/

/**
 * @brief Sets the file format of a VCF file.
 * @details Sets the file format of a VCF file. The file format is specified in the first line
 * with a construct of the form <em>fileformat=VCF4.x</em>.
 *
 * @param fileformat The file format of the entry
 * @param length The length of the file format
 * @param file The entry whose file format will be set
 **/
void set_vcf_file_format(char *fileformat, int length, vcf_file_t *file);

/**
 * @brief Sets the name of a header entry, if applies.
 * @details Sets the name of a header entry, if this entry is of the form key=value.
 *
 * @param name The name of the entry
 * @param length The length of the name
 * @param entry The entry whose name will be set
 **/
void set_vcf_header_entry_name(char *name, int length, vcf_header_entry_t *entry);

/**
 * @brief Links a value with a header entry.
 * @details Links a value with a header entry.
 *
 * @param value The value
 * @param length The length of the value
 * @param entry The header entry
 **/
void add_vcf_header_entry_value(char *value, int length, vcf_header_entry_t *entry);


/* ********************************************************
 *                    Record management                   *
 * ********************************************************/

/**
 * @brief Sets the chromosome of a VCF record.
 * @details Sets the chromosome of a VCF record.
 *
 * @param chromosome The value of the chromosome
 * @param length The length of the chromosome field
 * @param record The record whose attribute will be set
 **/
void set_vcf_record_chromosome(char* chromosome, int length, vcf_record_t* record);

/**
 * @brief Sets the position of a VCF record.
 * @details Sets the position of a VCF record.
 *
 * @param position The value of the position
 * @param record The record whose attribute will be set
 **/
void set_vcf_record_position(long position, vcf_record_t* record);

/**
 * @brief Sets the ID (if existing) of a VCF record.
 * @details Sets the ID (if existing) of a VCF record.
 *
 * @param id The value of the ID
 * @param length The length of the ID field
 * @param record The record whose attribute will be set
 **/
void set_vcf_record_id(char* id, int length, vcf_record_t* record);

/**
 * @brief Sets the reference allele of a VCF record.
 * @details Sets the reference allele of a VCF record.
 *
 * @param reference The value of the reference allele
 * @param length The length of the reference allele field
 * @param record The record whose attribute will be set
 **/
void set_vcf_record_reference(char* reference, int length, vcf_record_t* record);

/**
 * @brief Sets the alternate allele of a VCF record.
 * @details Sets the alternate allele of a VCF record.
 *
 * @param alternate The value of the alternate allele
 * @param length The length of the alternate allele field
 * @param record The record whose attribute will be set
 **/
void set_vcf_record_alternate(char* alternate, int length, vcf_record_t* record);

/**
 * @brief Sets the quality of a VCF record.
 * @details Sets the quality of a VCF record.
 *
 * @param quality The value of the quality
 * @param record The record whose attribute will be set
 **/
void set_vcf_record_quality(float quality, vcf_record_t* record);

/**
 * @brief Specifies if a VCF record passed some filters.
 * @details Specifies if a VCF record passed some filters.
 *
 * @param filter The value of the filter description
 * @param length The length of the filter field
 * @param record The record whose attribute will be set
 **/
void set_vcf_record_filter(char* filter, int length, vcf_record_t* record);

/**
 * @brief Sets the miscellaneous information of a VCF record.
 * @details Sets the miscellaneous information of a VCF record.
 *
 * @param info The value of the miscellaneous information
 * @param length The length of the info field
 * @param record The record whose attribute will be set
 **/
void set_vcf_record_info(char* info, int length, vcf_record_t* record);

/**
 * @brief Sets the format of the samples of a VCF record.
 * @details Sets the format of the samples of a VCF record.
 *
 * @param format The value of the format of the samples
 * @param length The length of the format field
 * @param record The record whose attribute will be set
 **/
void set_vcf_record_format(char* format, int length, vcf_record_t* record);

/**
 * @brief Adds the value of a sample in a VCF record.
 * @details Adds the value of a sample in a VCF record.
 *
 * @param sample The value of the samples
 * @param length The length of the sample value
 * @param record The record to add the sample to
 **/
void add_vcf_record_sample(char* sample, int length, vcf_record_t* record);



KHASH_MAP_INIT_STR(ids, int);

individual_t **sort_individuals(vcf_file_t *vcf, ped_file_t *ped);

khash_t(ids)* associate_samples_and_positions(vcf_file_t* file);

#endif
