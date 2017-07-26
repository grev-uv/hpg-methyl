#ifndef UTIL_H
#define UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include <commons/log.h>


/**
 * @file vcf_util.h
 * @author Cristina Yenyxe Gonzalez Garcia
 * @brief Getting diverse information from a VCF file
 * @details This file includes functions for getting information from a VCF file that can make 
 * easier information retrieval, such as a certain value in the FORMAT field, the value of the 
 * alleles of a sample, and so on.
 */


/**
 * @brief Status when reading a genotype alleles
 * @details Status when reading a genotype alleles (all available, something missing, haploid...)
 */
enum alleles_code { ALLELES_OK, FIRST_ALLELE_MISSING, SECOND_ALLELE_MISSING, ALL_ALLELES_MISSING, HAPLOID };

/**
 * @brief Flag which defines whether VCF files will be memory-mapped instead of using the I/O API.
 * @details Flag which defines whether VCF files will be memory-mapped instead of using the I/O API.
 * 
 * @see http://www.kernel.org/doc/man-pages/online/pages/man2/mmap.2.html
 */
int mmap_vcf;


/**
 * @brief Returns the number of comma-separated regions represented in a string
 * @details Given a string containing comma-separated regions, returns the number of them. A region could 
 * be represented as chromosome:position:ref_allele:alt_allele.
 *
 * @param regions_string The string that contains the regions
 * @return size_t The number of regions
 **/
size_t count_regions(char *regions_string);


/**
 * @brief Given the acronym of a field of the INFO column, retrieves its value
 * @details The INFO column in a VCF file is composed of pairs (acronym, value). The structure is 
 * the following: acronym1=value1,acronym2=value2,...,acronymN=valueN. Given one of the 
 * acronyms, this function returns its associated value.
 *
 * @param field The acronym of the field
 * @param info The contents of the INFO column
 * @return The value associated to the acronym
 **/
char *get_field_value_in_info(const char *field, char *info);

char *set_field_value_in_info(char *key, char *value, int append, char *info_in, int info_len);


/**
 * @brief Given the acronym of a field in the FORMAT column, retrives its position
 * @details The FORMAT column in a VCF file is of the form acronym1:acronym2:...:acronymN. Given one 
 * of the acronyms, this function returns its position.
 *
 * @param field The acronym of the field
 * @param format The contents of the FORMAT column
 * @return The position associated to the acronym
 **/
int get_field_position_in_format(const char *field, char *format);


/**
 * @brief Retrieves the value of the i-th field of a sample
 * @details Retrieves the value of the i-th colon-separated field in a sample.
 *
 * @param sample The sample to extract the value of the field from
 * @param position The position of the field whose value we are interested in
 * @return The text of the value (must be cast by the user)
 **/
char *get_field_value_in_sample(char *sample, int position);

void set_field_value_in_sample(char **sample, int position, char* value);


/**
 * @brief Given the value of a sample and the position of its genotype, returns the value of its alleles
 * @details A sample in a VCF file is described in the FORMAT column, so its values are of the form 
 * value1:value2:value3. Given the position where the genotype field is, this function returns
 * the value of its alleles. In a haploid position, only the first allele is retrieved.
 *
 * @param sample The contents of a sample column
 * @param genotype_position Position of the genotype field
 * @param allele1 [out] Value of the first allele
 * @param allele2 [out] Value of the second allele
 * @return ALLELES_OK if both alleles are present, FIRST_ALLELE_MISSING if the first one is missing, 
 * SECOND_ALLELE_MISSING if the second one is missing, ALL_ALLELES_MISSING if both are missing, 
 * HAPLOID if the position is haploid
 * @see get_field_position_in_format()
 **/
enum alleles_code get_alleles(char *sample, int genotype_position, int *allele1, int *allele2);


#ifdef __cplusplus
}
#endif

#endif
