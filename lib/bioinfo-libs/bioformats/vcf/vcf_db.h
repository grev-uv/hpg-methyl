#ifndef VCF_DB_H
#define VCF_DB_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <commons/log.h>
#include <commons/sqlite/sqlite3.h>
#include <containers/array_list.h>

//------------------------------------------------------------------------

#define VCF_CHUNKSIZE 10000

//------------------------------------------------------------------------

typedef struct vcf_query_fields {
    char *chromosome;
    unsigned long position;
    
    char *allele_ref;                   /**< Reference allele of the variant. */
    char *allele_maf;                   /**< Allele with MAF of the variant. */
    char *genotype_maf;                 /**< Genotype with MAF of the variant. */
    
    float allele_maf_freq;              /**< Frequency of the MAF allele. */
    float genotype_maf_freq;            /**< Frequency of the MAF genotype. */
    
    int missing_alleles;                /**< Number of alleles whose information is missing. */
    int missing_genotypes;              /**< Number of genotypes with at least one allele missing. */
    int mendelian_errors;               /**< Number of mendelian errors found. */
    int is_indel;                       /**< Whether this variant is an indel or not. */
    
    float cases_percent_dominant;       /**< Percentage of cases that follow a dominant inheritance pattern */
    float controls_percent_dominant;    /**< Percentage of controls that follow a dominant inheritance pattern */
    float cases_percent_recessive;      /**< Percentage of cases that follow a recessive inheritance pattern */
    float controls_percent_recessive;   /**< Percentage of controls that follow a recessive inheritance pattern */
} vcf_query_fields_t;

//------------------------------------------------------------------------

vcf_query_fields_t *vcf_query_fields_new(char *chr, unsigned long position, char *allele_ref, char *allele_maf, float allele_maf_freq, 
                                         char *genotype_maf, float genotype_maf_freq, int missing_alleles, int missing_genotypes,
                                         int mendelian_errors, int is_indel, float cases_percent_dominant, float controls_percent_dominant,
                                         float cases_percent_recessive, float controls_percent_recessive);

void vcf_query_fields_free(vcf_query_fields_t *p);

void print_vcf_query_fields(vcf_query_fields_t *p);

//------------------------------------------------------------------------
// 
//------------------------------------------------------------------------

int create_vcf_query_fields(sqlite3 *db);
int create_vcf_index(sqlite3 *db);

//------------------------------------------------------------------------

int insert_vcf_query_fields(void *custom_fields, sqlite3 *db);

int insert_vcf_query_fields_list(array_list_t *list, sqlite3 *db);

int prepare_statement_vcf_query_fields(sqlite3 *db, sqlite3_stmt **stmt);

int insert_statement_vcf_query_fields(void *custom_fields, 
				      sqlite3_stmt *stmt, sqlite3 *db);


//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif // end of VCF_DB_H
