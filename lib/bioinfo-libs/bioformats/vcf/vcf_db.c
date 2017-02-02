#include "vcf_db.h"

//------------------------------------------------------------------------

vcf_query_fields_t *vcf_query_fields_new(char *chr, unsigned long position, char *allele_ref, char *allele_maf, float allele_maf_freq, 
                                         char *genotype_maf, float genotype_maf_freq, int missing_alleles, int missing_genotypes,
                                         int mendelian_errors, int is_indel, float cases_percent_dominant, float controls_percent_dominant,
                                         float cases_percent_recessive, float controls_percent_recessive) {
    vcf_query_fields_t *p = malloc(sizeof(vcf_query_fields_t));

    p->chromosome = strdup(chr);
    p->position = position;

    p->allele_ref = strdup(allele_ref);
    p->allele_maf = strdup(allele_maf);
    p->allele_maf_freq = allele_maf_freq;
    p->genotype_maf = strdup(genotype_maf);
    p->genotype_maf_freq = genotype_maf_freq;

    p->missing_alleles = missing_alleles;
    p->missing_genotypes = missing_genotypes;
    p->mendelian_errors = mendelian_errors;
    p->is_indel = is_indel;

    p->cases_percent_dominant = cases_percent_dominant;
    p->controls_percent_dominant = controls_percent_dominant;
    p->cases_percent_recessive = cases_percent_recessive;
    p->controls_percent_recessive = controls_percent_recessive;

    return p;
}

void vcf_query_fields_free(vcf_query_fields_t *p) {
    assert(p);
    free(p->chromosome);
    free(p->allele_ref);
    free(p->allele_maf);
    free(p->genotype_maf);
}

void print_vcf_query_fields(vcf_query_fields_t *fields) {
    printf("chr = %s, pos = %ld, AL_ref = %s, AL_maf = %s (%f), GT_maf = %s (%f), miss %d|%d, Mendel = %d, indel = %d, dominant %f|%f, recessive %f|%f\n", 
            fields->chromosome, fields->position, fields->allele_ref, fields->allele_maf, fields->allele_maf_freq, 
            fields->genotype_maf, fields->genotype_maf_freq, fields->missing_alleles, fields->missing_genotypes, 
            fields->mendelian_errors, fields->is_indel, fields->cases_percent_dominant, fields->controls_percent_dominant, 
            fields->cases_percent_recessive, fields->controls_percent_recessive);
}

//------------------------------------------------------------------------
//
//------------------------------------------------------------------------

int create_vcf_query_fields(sqlite3 *db) {
    // create record_stats table for vcf files
    int rc;
    char *error_msg;
    char sql[1000];
    sprintf(sql, "CREATE TABLE record_query_fields (\
                               chromosome TEXT, position INT64, allele_ref TEXT, \
                               allele_maf TEXT, genotype_maf TEXT, \
                               allele_maf_freq DOUBLE, genotype_maf_freq DOUBLE, \
                               miss_allele INT, miss_gt INT, \
                               mendel_err INT, is_indel INT, \
                               cases_percent_dominant DOUBLE, \
                               controls_percent_dominant DOUBLE, \
                               cases_percent_recessive DOUBLE, \
                               controls_percent_recessive DOUBLE)");

    if (rc = sqlite3_exec(db, sql, NULL, NULL, &error_msg)) {
        LOG_FATAL_F("Stats database failed: %s\n", error_msg);
    }
    return 0;
}

//------------------------------------------------------------------------

int create_vcf_index(sqlite3 *db) {
    int rc;
    char sql[128];
    char *error_msg;

    // create chunks index
    sprintf(sql, "CREATE INDEX record_query_fields_chromosome_start_end_idx ON chunk (chromosome, start, end)");
    if (rc = sqlite3_exec(db, sql, NULL, NULL, &error_msg)) {
        LOG_ERROR_F("Stats database failed creating VCF index: %s\n", error_msg);
    }
    return rc;
}

//------------------------------------------------------------------------

int insert_vcf_query_fields(void *custom_fields, sqlite3 *db) {
    vcf_query_fields_t *fields = (vcf_query_fields_t *) custom_fields;

    //  print_vcf_query_fields(custom_fields);

    // insert into record_query_fields
    int rc;
    char sql[500 + strlen(fields->chromosome) + strlen(fields->allele_ref) * 4]; // AL_ref + AL_maf + GT_maf
    sprintf(sql, "INSERT INTO record_query_fields VALUES('%s', %ld, '%s', '%s', '%s', %f, %f, %i, %i, %i, %i, %f, %f, %f, %f)", 
            fields->chromosome, fields->position, fields->allele_ref, fields->allele_maf, fields->genotype_maf, 
            fields->allele_maf_freq, fields->genotype_maf_freq, fields->missing_alleles, fields->missing_genotypes, 
            fields->mendelian_errors, fields->is_indel, fields->cases_percent_dominant, fields->controls_percent_dominant, 
            fields->cases_percent_recessive, fields->controls_percent_recessive);

    char *error_msg;
    if (rc = sqlite3_exec(db, sql, NULL, NULL, &error_msg)) {
        LOG_DEBUG_F("Stats database failed: %s\n", error_msg);
    }

    return rc;
}

//------------------------------------------------------------------------

int insert_vcf_query_fields_list(array_list_t *list, sqlite3 *db) {

    int rc;
    sqlite3_stmt *stmt;
    vcf_query_fields_t *fields;
    char *error_message;

    prepare_statement_vcf_query_fields(db, &stmt);

    if (rc = sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &error_message)) {
        LOG_DEBUG_F("Stats databases failed: %s (%d)\n", rc, error_message);
    }

    int num_items = array_list_size(list);
    for (int i = 0; i < num_items; i++) {
        fields = array_list_get(i, list);

        sqlite3_bind_text(stmt, 1, fields->chromosome, strlen(fields->chromosome), SQLITE_STATIC);
        sqlite3_bind_int64(stmt, 2, fields->position);
        sqlite3_bind_text(stmt, 3, fields->allele_ref, strlen(fields->allele_ref), SQLITE_STATIC);
        sqlite3_bind_text(stmt, 4, fields->allele_maf, strlen(fields->allele_maf), SQLITE_STATIC);
        sqlite3_bind_text(stmt, 5, fields->genotype_maf, strlen(fields->genotype_maf), SQLITE_STATIC);
        sqlite3_bind_double(stmt, 6, fields->allele_maf_freq);
        sqlite3_bind_double(stmt, 7, fields->genotype_maf_freq);
        sqlite3_bind_int(stmt, 8, fields->missing_alleles);
        sqlite3_bind_int(stmt, 9, fields->missing_genotypes);
        sqlite3_bind_int(stmt, 10, fields->mendelian_errors);
        sqlite3_bind_int(stmt, 11, fields->is_indel);
        sqlite3_bind_double(stmt, 12, fields->cases_percent_dominant);
        sqlite3_bind_double(stmt, 13, fields->controls_percent_dominant);
        sqlite3_bind_double(stmt, 14, fields->cases_percent_recessive);
        sqlite3_bind_double(stmt, 15, fields->controls_percent_recessive);

        if (rc = sqlite3_step(stmt) != SQLITE_DONE) {
            LOG_DEBUG_F("Stats databases failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
        }

        sqlite3_reset(stmt);
    }

    if (rc = sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &error_message)) {
        LOG_DEBUG_F("Stats databases failed: %s (%d)\n", rc, error_message);
    }

    sqlite3_finalize(stmt);
}

//------------------------------------------------------------------------

int prepare_statement_vcf_query_fields(sqlite3 *db, sqlite3_stmt **stmt) {
  char sql[] = "INSERT INTO record_query_fields VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14, ?15)";
  return sqlite3_prepare_v2(db, sql, strlen(sql) + 300, stmt, NULL);
}

//------------------------------------------------------------------------

int insert_statement_vcf_query_fields(void *custom_fields, sqlite3_stmt *stmt, sqlite3 *db) {  
    int rc;
    vcf_query_fields_t *fields = (vcf_query_fields_t *) custom_fields;

    sqlite3_bind_text(stmt, 1, fields->chromosome, strlen(fields->chromosome), SQLITE_STATIC);
    sqlite3_bind_int64(stmt, 2, fields->position);
    sqlite3_bind_text(stmt, 3, fields->allele_ref, strlen(fields->allele_ref), SQLITE_STATIC);
    sqlite3_bind_text(stmt, 4, fields->allele_maf, strlen(fields->allele_maf), SQLITE_STATIC);
    sqlite3_bind_text(stmt, 5, fields->genotype_maf, strlen(fields->genotype_maf), SQLITE_STATIC);
    sqlite3_bind_double(stmt, 6, fields->allele_maf_freq);
    sqlite3_bind_double(stmt, 7, fields->genotype_maf_freq);
    sqlite3_bind_int(stmt, 8, fields->missing_alleles);
    sqlite3_bind_int(stmt, 9, fields->missing_genotypes);
    sqlite3_bind_int(stmt, 10, fields->mendelian_errors);
    sqlite3_bind_int(stmt, 11, fields->is_indel);
    sqlite3_bind_double(stmt, 12, fields->cases_percent_dominant);
    sqlite3_bind_double(stmt, 13, fields->controls_percent_dominant);
    sqlite3_bind_double(stmt, 14, fields->cases_percent_recessive);
    sqlite3_bind_double(stmt, 15, fields->controls_percent_recessive);

    if (rc = sqlite3_step(stmt) != SQLITE_DONE) {
        LOG_DEBUG_F("Stats databases failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
    }

    sqlite3_reset(stmt);

    return rc;
}

