#ifndef VCF_STATS_REPORT_H
#define VCF_STATS_REPORT_H

#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include <bioformats/db/db_utils.h>
#include <commons/file_utils.h>
#include <commons/log.h>
#include <commons/sqlite/sqlite3.h>
#include <containers/array_list.h>

#include "vcf_db.h"
#include "vcf_stats.h"


char *get_vcf_stats_filename_prefix(char *vcf_filename, char *out_filename, char *outdir);


char *get_vcf_file_stats_output_filename(char *prefix);

void report_vcf_summary_stats(FILE *stats_fd, void *db, file_stats_t *stats);


char *get_variant_stats_output_filename(char *prefix);

void report_vcf_variant_stats(FILE *stats_fd, void *db, khash_t(stats_chunks) *hash, int num_variants, variant_stats_t **stats);

void report_vcf_variant_stats_header(FILE *stats_fd);


char *get_sample_stats_output_filename(char *prefix);

void report_vcf_sample_stats(FILE *stats_fd, void *db, size_t num_samples, sample_stats_t **stats);

void report_vcf_sample_stats_header(FILE *stats_fd);

#endif

