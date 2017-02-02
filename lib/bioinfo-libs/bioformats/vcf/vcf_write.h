#ifndef VCF_WRITE_H
#define VCF_WRITE_H

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <containers/array_list.h>

#include "vcf_file_structure.h"

/**
 * @file vcf_write.h
 * @author Cristina Yenyxe Gonzalez Garcia
 * @brief Writing a whole VCF file or any of its parts.
 * @details This file includes functions for writing a whole VCF file or any of its parts, including:
 * the file format, the header entries, the line that separates header from body (aka delimiter)
 * and the body itself, which is divided in batches that group several lines/records.
 */


/* **********************************************
 *                Whole file writing            *
 * **********************************************/

/**
 * @brief Serialize a VCF file given its structure in main memory.
 * @details Writes a whole VCF file given a structure stored in main memory. This process consists on 
 * 2 steps: writing the header (which contains the file format, the header entries and the 
 * line separating header and body) and writing the records (grouped in batches).
 * 
 * @param file The structure to be serialized to file
 * @param fd The descriptor of the file to serialize to
 * @return 0 if the data is successfully written, 1 otherwise
 */
int write_vcf_file(vcf_file_t *file, FILE *fd);


/* **********************************************
 *                  Header writing              *
 * **********************************************/

/**
 * @brief Writes a VCF header to the given file descriptor
 * @details Writes a VCF header to the given file descriptor. A header contains a declaration of the 
 * file format, several entries (INFO, FORMAT, etc.) and the line that separates it from the file 
 * body.
 * 
 * @param file The file whose header will be serialized
 * @param fd The descriptor of the file to serialize to
 * @return 0 if the data is successfully written, 1 otherwise
 **/
int write_vcf_header(vcf_file_t *file, FILE *fd);

/**
 * @brief Writes the information about file format to the given file descriptor
 * @details Writes the information about the format of a VCF file.
 * 
 * @param file The file whose file format will be serialized
 * @param fd The descriptor of the file to serialize to
 * @return 0 if the data is successfully written, 1 otherwise
 **/
int write_vcf_fileformat(vcf_file_t *file, FILE *fd);

/**
 * @brief Writes an entry of the header of a VCF file
 * @details Writes an entry of the header of a VCF file. An entry can be of the form:
 * 
 * - '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">' or
 * - '##reference=human_b36_both.fasta' or
 * - '##Some miscellaneous information'
 *
 * @param entry The entry to be serialized
 * @param fd The descriptor of the file to serialize to
 * @return 0 if the data is successfully written, 1 otherwise
 **/
int write_vcf_header_entry(vcf_header_entry_t *entry, FILE *fd);

/**
 * @brief Writes the line that delimits the end of the header and the beginning of the body of a VCF file
 * @details Writes the line that delimits the end of the header and the beginning of the body of a VCF file. 
 * This line contains the name of the different fields, as well as the name of the samples from the 
 * sequencing process.
 * 
 * @param file The file whose delimiter will be serialized
 * @param fd The descriptor of the file to serialize to
 * @return 0 if the data is successfully written, 1 otherwise
 **/
int write_vcf_delimiter(vcf_file_t *file, FILE *fd);

/**
 * @brief Given a list of sample names, writes the line that separates the header and the body of a VCF file
 * @details Given a list of sample names, writes the line that delimits the end of the header and the beginning 
 * of the body of a VCF file. This line contains the name of the different fields, as well as the name of the 
 * samples from the sequencing process.
 * 
 * @param sample_names The names of the samples
 * @param num_sample The number of samples
 * @param fd The descriptor of the file to serialize to
 * @return 0 if the data is successfully written, 1 otherwise
 **/
int write_vcf_delimiter_from_samples(char **sample_names, int num_samples, FILE *fd);


/* **********************************************
 *                  Body writing                *
 * **********************************************/

/**
 * @brief Writes a batch of VCF records to the given file descriptor
 * @details Writes a batch of VCF records to the given file descriptor. The fact that the records are grouped 
 * in a batch doesn't add any special contents to the resulting output.
 * 
 * @param batch The batch whose records will be serialized
 * @param fd The descriptor of the file to serialize to
 * @return 0 if the data is successfully written, 1 otherwise
 **/
int write_vcf_batch(vcf_batch_t *batch, FILE *fd);


/**
 * @brief Writes an entry/record of a VCF file to the given file descriptor
 * @details Writes a record to the given file descriptor, including all its basic fields and the data 
 * about that position in the genome for all the sample subjects.
 * 
 * @param record The record that will be serialized
 * @param fd The descriptor of the file to serialize to
 * @return 0 if the data is successfully written, 1 otherwise
 **/
int write_vcf_record(vcf_record_t* record, FILE *fd);


#endif 
