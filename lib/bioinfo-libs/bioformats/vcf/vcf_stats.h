#ifndef VCF_STATS_H
#define VCF_STATS_H

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <bioformats/family/checks_family.h>
#include <bioformats/family/family.h>
#include <commons/log.h>
#include <containers/array_list.h>
#include <containers/khash.h>
#include <containers/list.h>

#include "vcf_file_structure.h"
#include "vcf_util.h"

/**
 * @file vcf_stats.h
 * @author Cristina Yenyxe Gonzalez Garcia
 * @brief Getting statistics from VCF files
 * @details This file includes functions for getting statistics from a VCF file. These statistics can 
 * be obtained one per variant or after analyzing the whole file.
 */


/**
 * @brief Statistics global to a VCF file
 * @details Statistics global to a VCF file. These are, among others, the number of variants, SNPs 
 * and indels it contains, and the mean quality of the measures.
 **/
typedef struct file_stats {
    int variants_count;         /**< Number of variants of the file. */
    int samples_count;          /**< Number of samples of the file. */
    
    int snps_count;             /**< Number of SNPs of the file. */
    int indels_count;           /**< Number of insertions/deletions of the file. */
    
    int transitions_count;      /**< Number of transitions (A <-> G or C <-> T) of the file. */
    int transversions_count;    /**< Number of transversions (A <-> C, A <-> T, C <-> G, G <-> T) of the file. */
    
    int biallelics_count;       /**< Number of variants with only one alternate allele of the file. */
    int multiallelics_count;    /**< Number of variants with more than one alternate allele of the file. */
    
    int pass_count;             /**< Number of variants that passed a filter. */
    float accum_quality;        /**< Sum of all values of the QUALITY column. */
    float mean_quality;         /**< Mean of all values of the QUALITY column. */
} file_stats_t;


/**
 * @brief Statistics of a variant of a VCF file
 * 
 * @details Statistics of a variant of a VCF file. These are, among others, the count and frequency of 
 * each allele and genotype that can be obtained by combining them.
 **/
typedef struct variant_stats {
    char *chromosome;           /**< Chromosome of the variant. */
    unsigned long position;     /**< Position of the variant. */
    
    char *ref_allele;           /**< Reference allele of the variant. */
    char **alternates;          /**< List of alternate alleles of the variant. */
    char *maf_allele;           /**< Allele with MAF. */
    char *mgf_genotype;         /**< Genotype with MGF. */
    
    int num_alleles;            /**< Number of alleles of the variant (1 reference + N alternates). */
    int *alleles_count;         /**< Times each allele has been counted. */
    int *genotypes_count;       /**< Times each possible genotype has been counted. */
    float *alleles_freq;        /**< Frequency of each allele in relation to the total. */
    float *genotypes_freq;      /**< Frequency of each genotype in relation to the total. */
    float maf;                  /**< Minimum allele frequency. */
    float mgf;                  /**< Minimum genotype frequency. */
    
    int missing_alleles;        /**< Number of alleles whose information is missing. */
    int missing_genotypes;      /**< Number of genotypes with at least one allele missing. */
    int mendelian_errors;       /**< Number of mendelian errors found. */
    int is_indel;               /**< Whether this variant is an indel or not. */
    
    float cases_percent_dominant;       /**< Percentage of cases that follow a dominant inheritance pattern */
    float controls_percent_dominant;    /**< Percentage of controls that follow a dominant inheritance pattern */
    float cases_percent_recessive;      /**< Percentage of cases that follow a recessive inheritance pattern */
    float controls_percent_recessive;   /**< Percentage of controls that follow a recessive inheritance pattern */
} variant_stats_t;

/**
 * @brief Statistics of a sample of a VCF file
 * 
 * @details Statistics of a samples of a VCF file. These are, among others, the number of missing 
 * genotypes and mendelian errors.
 **/
typedef struct sample_stats {
    char *name;                 /**< Name of the sample. */
    
    size_t mendelian_errors;    /**< Number of mendelian errors. */
    size_t missing_genotypes;   /**< Number of genotypes with at least one allele missing. */
} sample_stats_t;


/* ********************************
 * Initialization and destruction *
 * ********************************/

/**
 * @brief Allocates memory for a file_stats_t structure
 * @details Allocates memory for a file_stats_t structure
 * 
 * @return A new file_stats_t structure
 **/
file_stats_t *file_stats_new();

/**
 * @brief Deallocates memory associated to a file_stats_t structure
 * @details Deallocates memory associated to a file_stats_t structure
 * 
 * @param stats The structure to be freed
 */
void file_stats_free(file_stats_t *stats);


/**
 * @brief Initializes a variant_stats_t structure mandatory fields
 * @details Initializes a variant_stats_t structure mandatory fields, which are its chromosome, position and 
 * reference allele.
 * 
 * @param chromosome Chromosome of the variant in the genome
 * @param position Position of the variant in the chromosome
 * @param ref_allele Reference allele
 * @return A new variant_stats_t structure
 */
variant_stats_t *variant_stats_new(char *chromosome, unsigned long position, char *ref_allele);

/**
 * @brief Deallocates memory associated to a variant_stats_t structure
 * @details Deallocates memory associated to a variant_stats_t structure
 * 
 * @param stats The structure to be freed
 */
void variant_stats_free(variant_stats_t *stats);


/**
 * @brief Initializes a sample_stats_t structure mandatory fields
 * @details Initializes a sample_stats_t structure mandatory field, the sample name.
 * 
 * @param name Name of the sample
 * @return A new sample_stats_t structure
 */
sample_stats_t *sample_stats_new(char *name);

/**
 * @brief Deallocates memory associated to a sample_stats_t structure
 * @details Deallocates memory associated to a sample_stats_t structure
 * 
 * @param stats The structure to be freed
 */
void sample_stats_free(sample_stats_t *stats);



/* ******************************
 *           Execution          *
 * ******************************/

/**
 * @brief Given a list of variants, gets their statistics and also the ones that apply to the VCF file
 * @details Given a list of variants, gets their statistics and also the ones that apply to the VCF file. The statistics 
 * per variant are queued in the output_list argument, and the statistics of the whole file are stored in the 
 * file_stats structure.
 *
 * @param variants The list of variants whose statistics will be got
 * @param num_variants The number of variants
 * @param individuals The list of samples that will be used to retrieve some statistics
 * @param sample_ids Relationship between the name of a sample and its position in the VCF file
 * @param output_list [out] The list where the statistics per variant will be stored
 * @param file_stats [in,out] The statistics of the VCF file
 * @return Whether the statistics were successfully retrieved
 **/
int get_variants_stats(vcf_record_t **variants, int num_variants, individual_t **individuals, khash_t(ids) *sample_ids, 
                       list_t *output_list, file_stats_t *file_stats);

/**
 * @brief Given a list of variants, gets the statistics related to their samples and also the ones that apply to the VCF file
 * @details Given a list of variants, gets the statistics that apply to samples and to the whole VCF file. The statistics 
 * per sample are stored in a sample_stats_t structure, and the one about the file are in the file_stats structure.
 * The individuals and sample_ids arguments must be retrieved by using the 'sort_individuals' and 'associate_samples_and_positions'
 * functions in vcf_file_structure.h, respectively.
 * 
 * @param variants The list of variants whose samples statistics will be got
 * @param num_variants The number of variants
 * @param individuals The list of samples whose statistics will be got
 * @param sample_ids Relationship between the name of a sample and its position in the VCF file
 * @param sample_stats [in,out] The statistics of the samples
 * @param file_stats [in,out] The statistics of the VCF file
 * @return Whether the statistics were successfully retrieved
 **/
int get_sample_stats(vcf_record_t **variants, int num_variants, individual_t **individuals, khash_t(ids) *sample_ids, 
                     sample_stats_t **sample_stats, file_stats_t *file_stats);

/**
 * @brief Given the statistics of a file (supposedly not fully-processed), updates the value of its statistics
 * @details Given the statistics of a file (supposedly not fully-processed), updates the value of its statistics. The variants, 
 * samples, SNPs, transition, transversion, indels, biallelic, multiallelic and PASS count will be accumulated to the 
 * previous values. The accumulated quality will also be summed, and as a result its mean will be recalculated.
 *
 * @param variants_count Number of new variants found
 * @param samples_count Number of samples
 * @param snps_count Number of new SNPs found
 * @param transitions_count Number of new transitions found
 * @param transversions_count Number of new transversions found
 * @param indels_count Number of new insertions/deletions found
 * @param biallelics_count Number of new biallelic variants found
 * @param multiallelics_count Number of new multiallelic variants found
 * @param pass_count Number of new variants that passed a filter found
 * @param accum_quality Sum of all new values of the QUALITY column found
 * @param stats [in,out] The statistics of the VCF file
 **/
void update_file_stats(int variants_count, int samples_count, int snps_count, int transitions_count, int transversions_count, 
                       int indels_count, int biallelics_count, int multiallelics_count, int pass_count, float accum_quality, 
                       file_stats_t *stats);


#endif
