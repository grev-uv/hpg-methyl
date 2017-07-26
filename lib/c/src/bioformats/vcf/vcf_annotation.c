#include "vcf_annotation.h"
 
void precalculate_aux_values_for_annotation(int calculate_stats, int calculate_dp, int calculate_mq, vcf_record_t *record,
                                            variant_stats_t **variant_stats, file_stats_t *file_stats, list_t *stats_list, 
                                            int *dp, int *mq0, double *mq) {
    // Variant statistics
    if (calculate_stats) {
        get_variants_stats(&record, 1, NULL, NULL,0, stats_list, file_stats);
        list_item_t *item = list_remove_item(stats_list);
        *variant_stats = item->data_p;
        list_item_free(item);
    }
    
    // Read depth
    if (calculate_dp) {
        int dp_pos = -1;
        *dp = 0;
        for (int j = 0; j < record->samples->size; j++) {
            char *aux = strndup(record->format, record->format_len);
            dp_pos = get_field_position_in_format("DP", aux);
            free(aux);
            if (dp_pos >= 0) {
                aux = strdup((char*) array_list_get(j, record->samples));
                *dp += atoi(get_field_value_in_sample(aux, dp_pos));
                free(aux);
            }
        }
    }
    
    // Mapping quality
    if (calculate_mq) {
        for (int j = 0; j < record->samples->size; j++) {
            char *aux = strndup(record->format, record->format_len);
            int mq_pos = get_field_position_in_format("GQ", aux);
            free(aux);
            if (mq_pos < 0) {
                continue;
            }

            aux = strdup((char*) array_list_get(j, record->samples));
            int cur_gq = atoi(get_field_value_in_sample(aux, mq_pos));
            free(aux);
//                 printf("sample = %s\tmq_pos = %d\tvalue = %d\n", array_list_get(j, record->samples), mq_pos,
//                        atoi(get_field_value_in_sample(strdup((char*) array_list_get(j, record->samples)), mq_pos)));
            if (cur_gq == 0) {
                (*mq0)++;
            } else {
                *mq += cur_gq * cur_gq;
            }
        }
        *mq = sqrt(*mq / record->samples->size);
    }
}

char *get_annotation_allele_count(variant_stats_t *variant_stats, int *output_len) {
    char *result = calloc(3 + variant_stats->num_alleles * 10, sizeof(char));
    int len = 3;
    strncat(result, "AC=", 3);
    
    for (int j = 1; j < variant_stats->num_alleles - 1; j++) {
        sprintf(result + len, "%d,", variant_stats->alleles_count[j]);
        len = strlen(result);
    }
    sprintf(result + len, "%d", variant_stats->alleles_count[variant_stats->num_alleles - 1]);
    len = strlen(result);
    
    *output_len = len;
    return result;
}

char *get_annotation_allele_freq(variant_stats_t *variant_stats, int *output_len) {
    // For each ALT, AF = alt_freq / (1 - ref_freq)
    char *result = calloc(3 + variant_stats->num_alleles * 7, sizeof(char));
    int len = 3;
    strncat(result, "AF=", 3);
    
    for (int j = 1; j < variant_stats->num_alleles - 1; j++) {
        sprintf(result + len, "%.3f,", variant_stats->alleles_freq[j]);
        len = strlen(result);
    }
    sprintf(result+len, "%.3f", variant_stats->alleles_freq[variant_stats->num_alleles - 1]);
    len = strlen(result);
    
    *output_len = len;
    return result;
}

char *get_annotation_allele_number(variant_stats_t *variant_stats, int *output_len) {
    char *result = calloc(8, sizeof(char));
    sprintf(result, "AN=%d", variant_stats->num_alleles);
    *output_len = strlen(result);
    return result;
}


char *get_annotation_read_depth(int dp, int *output_len) {
    char *result = calloc(16, sizeof(char));
    sprintf(result, "DP=%d", dp);
    *output_len = strlen(result);
    return result;
}

char *get_annotation_quality_by_depth(vcf_record_t *record, int dp, int *output_len) {
    char *result = calloc(16, sizeof(char));
    int len = 0;
    
    if (record->quality < 0) {
        strncat(result, "QD=.;", 1);
        len = 5;
    } else if (record->quality > 0) {
        sprintf(result, "QD=%.3f;", record->quality / dp);
        len = strlen(result);
    } else {
        strncat(result, "QD=0.000;", 1);
        len = 5;
    }
    
    *output_len = len;
    return result;
}


char *get_annotation_mapping_quality(double mq, int *output_len) {
    char *result = calloc(16, sizeof(char));
    sprintf(result, "MQ=%.3f", mq);
    *output_len = strlen(result);
    return result;
}

char *get_annotation_mapping_quality_zero(int mq0, int *output_len) {
    char *result = calloc(16, sizeof(char));
    sprintf(result, "MQ0=%d", mq0);
    *output_len = strlen(result);
    return result;
}


char *get_annotation_non_missing_samples(vcf_record_t *record, char *empty_sample, int *output_len) {
    char *result = calloc(16, sizeof(char));
    int ns = 0;
    for (int j = 0; j < record->samples->size; j++) {
        if (strcmp(array_list_get(j, record->samples), empty_sample)) {
            ns++;
        }
    }
    sprintf(result, "NS=%d", ns);
    *output_len = strlen(result);
    return result;
}
