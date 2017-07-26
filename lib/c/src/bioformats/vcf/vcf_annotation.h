#ifndef VCF_ANNOTATION_H
#define	VCF_ANNOTATION_H

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>

#include <containers/array_list.h>
#include <containers/list.h>

#include "vcf_stats.h"

void precalculate_aux_values_for_annotation(int calculate_stats, int calculate_dp, int calculate_mq, vcf_record_t *record,
                                            variant_stats_t **variant_stats, file_stats_t *file_stats, list_t *stats_list, 
                                            int *dp, int *mq0, double *mq);

char *get_annotation_allele_count(variant_stats_t *variant_stats, int *output_len);

char *get_annotation_allele_freq(variant_stats_t *variant_stats, int *output_len);

char *get_annotation_allele_number(variant_stats_t *variant_stats, int *output_len);


char *get_annotation_read_depth(int dp, int *output_len);

char *get_annotation_quality_by_depth(vcf_record_t *record, int dp, int *output_len);


char *get_annotation_mapping_quality(double mq, int *output_len);

char *get_annotation_mapping_quality_zero(int mq0, int *output_len);


char *get_annotation_non_missing_samples(vcf_record_t *record, char *empty_sample, int *output_len);


#ifdef __cplusplus
}
#endif

#endif
