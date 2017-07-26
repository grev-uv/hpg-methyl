#ifndef VCF_DB_H
#define VCF_DB_H

#ifdef __cplusplus
extern "C" {
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <commons/log.h>
#include <sqlite/sqlite3.h>
#include <containers/array_list.h>


#define VCF_CHUNKSIZE 2000

typedef struct variant_stats_db_fields {
    char *chromosome;
    unsigned long position;
    
    char *allele_ref;                   /**< Reference allele of the variant. */
    char *allele_alt;                   /**< Alternative allele(s) of the variant. */
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
} variant_stats_db_fields_t;

typedef struct sample_stats_db_fields {
    char *name;
    int missing_genotypes;
    int mendelian_errors;
} sample_stats_db_fields_t;


variant_stats_db_fields_t *variant_stats_db_fields_new(
                char *chr, unsigned long position, char *allele_ref, char *allele_alt, 
                char *allele_maf, float allele_maf_freq, char *genotype_maf, float genotype_maf_freq, 
                int missing_alleles, int missing_genotypes, int mendelian_errors, int is_indel, 
                float cases_percent_dominant, float controls_percent_dominant,
                float cases_percent_recessive, float controls_percent_recessive);

void variant_stats_db_fields_free(variant_stats_db_fields_t *p);

void print_variant_stats_db_fields(variant_stats_db_fields_t *p);


sample_stats_db_fields_t *sample_stats_db_fields_new(char *name, int missing_genotypes, int mendelian_errors);

void sample_stats_db_fields_free(sample_stats_db_fields_t* p);

void print_sample_stats_db_fields(sample_stats_db_fields_t* p);



int pre_variant_stats_db(sqlite3 *db);

int post_variant_stats_db(sqlite3 *db);



int insert_variant_stats_db_fields(void *custom_fields, sqlite3 *db);

int insert_variant_stats_db_fields_list(array_list_t *list, sqlite3 *db);

int prepare_statement_variant_stats_db_fields(sqlite3 *db, sqlite3_stmt **stmt);

int insert_statement_variant_stats_db_fields(void *custom_fields, sqlite3_stmt *stmt, sqlite3 *db);

int prepare_statement_sample_stats_db_fields(sqlite3 *db, sqlite3_stmt **stmt);

int insert_sample_stats_db_fields_list(array_list_t *list, sqlite3 *db);


#ifdef __cplusplus
}
#endif


#endif // end of VCF_DB_H
