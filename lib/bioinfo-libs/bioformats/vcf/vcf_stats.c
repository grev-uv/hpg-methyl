#include "vcf_stats.h"


static int is_mendelian_error(individual_t *father, individual_t *mother, individual_t *child, 
                              int child_allele1, int child_allele2, int gt_position, 
                              vcf_record_t *record, khash_t(ids) *sample_ids);

/* ******************************
 *      Whole file statistics   *
 * ******************************/
 
file_stats_t* file_stats_new() {
    return calloc (1, sizeof(file_stats_t));
}

void file_stats_free(file_stats_t* stats) {
    assert(stats);
    free(stats);
}

void update_file_stats(int variants_count, int samples_count, int snps_count, int transitions_count, int transversions_count, 
                       int indels_count, int biallelics_count, int multiallelics_count, int pass_count, float accum_quality, 
                       file_stats_t *stats) {
    assert(stats);
    stats->samples_count = samples_count;
    stats->variants_count += variants_count;
    stats->snps_count += snps_count;
    stats->transitions_count += transitions_count;
    stats->transversions_count += transversions_count;
    stats->indels_count += indels_count;
    stats->biallelics_count += biallelics_count;
    stats->multiallelics_count += multiallelics_count;
    stats->pass_count += pass_count;
    stats->accum_quality += accum_quality;
    stats->mean_quality = accum_quality / variants_count;
}


/* ******************************
 *     Per variant statistics   *
 * ******************************/
 
variant_stats_t* variant_stats_new(char *chromosome, unsigned long position, char *ref_allele) {
    assert(chromosome);
    assert(ref_allele);
    
    variant_stats_t *stats = (variant_stats_t*) malloc (sizeof(variant_stats_t));
    stats->chromosome = chromosome;
    stats->position = position;
    stats->ref_allele = ref_allele;
    stats->alternates = NULL;
    stats->maf_allele = NULL;
    
    stats->num_alleles = 1;
    stats->alleles_count = NULL;
    stats->genotypes_count = NULL;
    stats->alleles_freq = NULL;
    stats->genotypes_freq = NULL;
    stats->maf = 0.0f;
    stats->mgf = 0.0f;
    
    stats->missing_alleles = 0;
    stats->missing_genotypes = 0;
    stats->mendelian_errors = 0;
    stats->is_indel = 0;
    
    stats->cases_percent_dominant = 0.0f;
    stats->controls_percent_dominant = 0.0f;
    stats->cases_percent_recessive = 0.0f;
    stats->controls_percent_recessive = 0.0f;
    
    return stats;
}

void variant_stats_free(variant_stats_t* stats) {
    assert(stats);
    
    if (stats->chromosome) { free(stats->chromosome); }
    if (stats->ref_allele) { free(stats->ref_allele); }
    if (stats->alternates) {
        for (int i = 0; i < stats->num_alleles-1; i++) {
            free(stats->alternates[i]);
        }
        free(stats->alternates);
    }
    if (stats->alleles_count) { free(stats->alleles_count); }
    if (stats->genotypes_count) { free(stats->genotypes_count); }
    if (stats->alleles_freq) { free(stats->alleles_freq); }
    if (stats->genotypes_freq) { free(stats->genotypes_freq); }
    free(stats);
}

int get_variants_stats(vcf_record_t **variants, int num_variants, individual_t **individuals, khash_t(ids) *sample_ids, 
                       list_t *output_list, file_stats_t *file_stats) {
    assert(variants);
    assert(output_list);
    assert(file_stats);
    
    char *copy_buf = NULL, *sample;
    
    int num_alternates, gt_position, curr_position;
    int allele1, allele2, alleles_code;
    
    // Temporary variables for file stats updating
    int variants_count = 0, samples_count = 0, snps_count = 0, indels_count = 0, pass_count = 0;
    int transitions_count = 0, transversions_count = 0, biallelics_count = 0, multiallelics_count = 0;
    float accum_quality = 0;
    // Temporary variables for variant stats updating
    int total_alleles_count = 0, total_genotypes_count = 0;
    int cases_dominant = 0;     // Number of cases that follow a dominant inheritance pattern
    int controls_dominant = 0;  // Number of controls that follow a dominant inheritance pattern
    int cases_recessive = 0;    // Number of cases that follow a recessive inheritance pattern
    int controls_recessive = 0; // Number of controls that follow a recessive inheritance pattern
    float maf = 1.0f, mgf = 1.0f;
    float cur_gt_freq;
    
    // Variant stats management
    vcf_record_t *record;
    variant_stats_t *stats;
    for (int i = 0; i < num_variants; i++) {
        record = variants[i];
        stats = variant_stats_new(strndup(record->chromosome, record->chromosome_len), 
                                  record->position, 
                                  strndup(record->reference, record->reference_len));
        
        // Reset counters
        total_alleles_count = total_genotypes_count = 0;
        cases_dominant = controls_dominant = cases_recessive = controls_recessive = 0;
        maf = mgf = 1.0f;
        
        // Create list of alternates
        copy_buf = strndup(record->alternate, record->alternate_len);
        stats->alternates = split(copy_buf, ",", &num_alternates);
        if (copy_buf) {
            free(copy_buf);
            copy_buf = NULL;
        }
        
        stats->num_alleles = num_alternates + 1;
        LOG_DEBUG_F("num alternates = %d\tnum_alleles = %d\n", num_alternates, stats->num_alleles);
        
        // Create lists of allele and genotypes counters and frequencies
        stats->alleles_count = (int*) calloc (stats->num_alleles, sizeof(int));
        stats->genotypes_count = (int*) calloc (stats->num_alleles * stats->num_alleles, sizeof(int));
        stats->alleles_freq = (float*) calloc (stats->num_alleles, sizeof(float));
        stats->genotypes_freq = (float*) calloc (stats->num_alleles * stats->num_alleles, sizeof(float));
        
        // Get position where GT is in sample
        copy_buf = strndup(record->format, record->format_len);
        gt_position = get_field_position_in_format("GT", copy_buf);
        if (copy_buf) {
            free(copy_buf);
            copy_buf = NULL;
        }
        
        LOG_DEBUG_F("Genotype position = %d\n", gt_position);
        if (gt_position < 0) { continue; }   // This variant has no GT field
        
        /* Traverse samples and find:
         * - Present and missing alleles
         * - Mendelian errors
         * - (Dis)agreements with inheritance models
         */
        for (int j = 0; j < record->samples->size; j++) {
            sample = (char*) array_list_get(j, record->samples);
            
            // Get to GT position
            copy_buf = strdup(sample);
            alleles_code = get_alleles(copy_buf, gt_position, &allele1, &allele2);
            if (copy_buf) {
                free(copy_buf);
                copy_buf = NULL;
            }
            LOG_DEBUG_F("sample = %s, alleles = %d/%d\n", sample, allele1, allele2);
            
            // Check missing alleles and genotypes
            if (!alleles_code) {
                // Both alleles set
                curr_position = allele1 * (stats->num_alleles) + allele2;
                assert(allele1 <= stats->num_alleles);
                assert(allele2 <= stats->num_alleles);
                assert(curr_position <= stats->num_alleles * stats->num_alleles);
                
                stats->alleles_count[allele1] += 1;
                stats->alleles_count[allele2] += 1;
                stats->genotypes_count[curr_position] += 1;
                total_alleles_count += 2;
                total_genotypes_count++;
            } else {
                // Missing genotype (one or both alleles missing)
                stats->missing_genotypes++;
                if (allele1 < 0) { 
                    stats->missing_alleles++; 
                } else {
                    stats->alleles_count[allele1]++;
                    total_alleles_count++;
                }
                    
                if (allele2 < 0) { 
                    stats->missing_alleles++;
                } else {
                    stats->alleles_count[allele2]++;
                    total_alleles_count++;
                }
            }
            
            // Include statistics that depend on pedigree information
            if (individuals) {
                // Check mendelian errors (pedigree data must be given)
                if (sample_ids && !alleles_code) {
                    if (is_mendelian_error(individuals[j]->father, individuals[j]->mother, individuals[j], 
                                        allele1, allele2, gt_position, record, sample_ids) > 0) {
                        (stats->mendelian_errors)++;
                    }
                }

                // Check inheritance models
                if (individuals[j]->condition == UNAFFECTED) {
                    if (!allele1 && !allele2) { // 0|0
                        controls_dominant++;
                        controls_recessive++;
                    } else if ((!allele1 && allele2) || (allele1 && !allele2)) { // 0|1 or 1|0
                        controls_recessive++;
                    }
                } else if (individuals[j]->condition == AFFECTED) {
                    if (allele1 && allele2 && allele1 == allele2) { // 1|1, 2|2, and so on
                        cases_recessive++;
                        cases_dominant++;
                    } else if (allele1 || allele2) { // 0|1, 1|0, 1|2, 2|1, 1|3, and so on
                        cases_dominant++;
                    }
                }
            }
        }
        
        assert(cases_dominant >= cases_recessive);
        assert(controls_recessive >= controls_dominant);
        
        // Once all samples have been traverse, calculate % that follow inheritance model
        stats->controls_percent_dominant = (float) controls_dominant * 100 / (record->samples->size - stats->missing_genotypes);
        stats->cases_percent_dominant = (float) cases_dominant * 100 / (record->samples->size - stats->missing_genotypes);
        stats->controls_percent_recessive = (float) controls_recessive * 100 / (record->samples->size - stats->missing_genotypes);
        stats->cases_percent_recessive = (float) cases_recessive * 100 / (record->samples->size - stats->missing_genotypes);
        
        // Get allele and genotype frequencies, as well as MAF and MGF
        for (int j = 0; j < stats->num_alleles; j++) {
            stats->alleles_freq[j] = (float) stats->alleles_count[j] / total_alleles_count;
            if (stats->alleles_freq[j] < maf) {
                maf = stats->alleles_freq[j];
                stats->maf_allele = (j == 0) ? stats->ref_allele : stats->alternates[j-1];
            }
        }
        stats->maf = maf;
        
        for (int j = 0; j < stats->num_alleles * stats->num_alleles; j++) {
            printf("");
            stats->genotypes_freq[j] = (float) stats->genotypes_count[j] / total_genotypes_count;
        }
        
        for (int j = 0; j < stats->num_alleles; j++) {
            for (int k = j; k < stats->num_alleles; k++) {
                int idx1 = j * stats->num_alleles + k;
                if (j == k) {
                    cur_gt_freq = stats->genotypes_freq[idx1];
                } else {
                    int idx2 = k * stats->num_alleles + j;
                    cur_gt_freq = stats->genotypes_freq[idx1] + stats->genotypes_freq[idx2];
                }

                if (cur_gt_freq < mgf) {
                    if (copy_buf) { 
                        free(copy_buf); 
                        copy_buf = NULL;
                    }
                    char *first_allele = (j == 0) ? stats->ref_allele : stats->alternates[j-1];
                    char *second_allele = (k == 0) ? stats->ref_allele : stats->alternates[k-1];
                    copy_buf = malloc((strlen(first_allele) + strlen(second_allele) + 2) * sizeof(char));
                    sprintf(copy_buf, "%s|%s", first_allele, second_allele);
                    mgf = cur_gt_freq;
                }
            }
        }
        
        stats->mgf = mgf;
        stats->mgf_genotype = copy_buf;
        
        
        // Update variables finally used to update file_stats_t structure
        variants_count++;
        if (i == 0) { samples_count = record->samples->size; }  // Just once per batch
        if (strncmp(record->id, ".", record->id_len)) { snps_count++; }
        if (!strncmp(record->filter, "PASS", record->filter_len)) { pass_count++; }
        if (record->quality >= 0) { accum_quality += record->quality; } // -1 = N/A
        if (stats->num_alleles > 2) {
            multiallelics_count++; 
        } else if (stats->num_alleles > 1) {
            biallelics_count++;
        }
        
        /* 
            * 3 possibilities for being an INDEL:
            * - The value of the ALT field is <DEL> or <INS>
            * - The REF allele is not . but the ALT is
            * - The REF allele is . but the ALT is not
            * - The REF field length is different than the ALT field length
            */
        if ((strncmp(".", stats->ref_allele, 1) && !strncmp(".", record->alternate, 1)) ||
            (strncmp(".", record->alternate, 1) && !strncmp(".", stats->ref_allele, 1)) ||
            !strncmp("<INS>", record->alternate, record->alternate_len) ||
            !strncmp("<DEL>", record->alternate, record->alternate_len) ||
             record->reference_len != record->alternate_len) {
            stats->is_indel = 1;
            indels_count++;
        } else {
            stats->is_indel = 0;
        }
            
        int ref_len = strlen(stats->ref_allele);
        int alt_len;
        for (int j = 0; j < num_alternates; j++) {
            alt_len = strlen(stats->alternates[j]);
            
            // Transitions and transversions
            if (ref_len == 1 && alt_len == 1) {
                switch (stats->ref_allele[0]) {
                    case 'C':
                        if (stats->alternates[j][0] == 'T') {
                            transitions_count++;
                        } else {
                            transversions_count++;
                        }
                        break;
                    case 'T':
                        if (stats->alternates[j][0] == 'C') {
                            transitions_count++;
                        } else {
                            transversions_count++;
                        }
                        break;
                    case 'A':
                        if (stats->alternates[j][0] == 'G') {
                            transitions_count++;
                        } else {
                            transversions_count++;
                        }
                        break;
                    case 'G':
                        if (stats->alternates[j][0] == 'A') {
                            transitions_count++;
                        } else {
                            transversions_count++;
                        }
                        break;
                }
            }
        }
        
        // Insert results in output list
        list_item_t *variant_result = list_item_new(i, 0, stats);
        list_insert_item(variant_result, output_list);
    }
    
    // Update file_stats_t structure
#pragma omp critical 
    {
        update_file_stats(variants_count, samples_count, snps_count, transitions_count, transversions_count, 
                          indels_count, biallelics_count, multiallelics_count, pass_count, accum_quality, file_stats);
    }
    
    return 0;
}


/* ******************************
 *     Per sample statistics    *
 * ******************************/

sample_stats_t* sample_stats_new(char* name) {
    assert(name);
    sample_stats_t *stats = malloc (sizeof(sample_stats_t));
    stats->name = strdup(name);
    stats->mendelian_errors = 0;
    stats->missing_genotypes = 0;
    return stats;
}

void sample_stats_free(sample_stats_t* stats) {
    assert(stats);
    free(stats->name);
    free(stats);
}

int get_sample_stats(vcf_record_t **variants, int num_variants, individual_t **individuals, khash_t(ids) *sample_ids, 
                     sample_stats_t **sample_stats, file_stats_t *file_stats) {
    assert(variants);
    assert(sample_stats);
    assert(file_stats);
    
    int gt_position;
    int allele1, allele2, alleles_code;
    char *copy_buf, *sample;
    
    // Sample stats management
    vcf_record_t *record;
    for (int i = 0; i < num_variants; i++) {
        record = variants[i];
        for(int j = 0; j < record->samples->size; j++) {
            sample = (char*) array_list_get(j, record->samples);
            char *format_dup = strndup(record->format, record->format_len);
            gt_position = get_field_position_in_format("GT", format_dup);
            free(format_dup);
            
            // Get to GT position
            copy_buf = strdup(sample);
            alleles_code = get_alleles(copy_buf, gt_position, &allele1, &allele2);
            if (copy_buf) {
                free(copy_buf);
            }
            LOG_DEBUG_F("sample = %s, alleles = %d/%d\n", sample, allele1, allele2);
            
            // Find the missing alleles
            if (alleles_code > 0) {
                // Missing genotype (one or both alleles missing)
                #pragma omp atomic
                (sample_stats[j]->missing_genotypes)++;
                continue;
            }
            
            assert(individuals[j]);
            
            // Check mendelian errors
            if (individuals && sample_ids && !alleles_code &&
                    is_mendelian_error(individuals[j]->father, individuals[j]->mother, individuals[j], 
                                   allele1, allele2, gt_position, record, sample_ids) > 0) {
                #pragma omp atomic
                (sample_stats[j]->mendelian_errors)++;
            }
        }
        
    }
    
    return 0;
}

static int is_mendelian_error(individual_t *father, individual_t *mother, individual_t *child, 
                              int child_allele1, int child_allele2, int gt_position, 
                              vcf_record_t *record, khash_t(ids) *sample_ids) {
    int is_error = 0;
    int father_allele1, father_allele2;
    int mother_allele1, mother_allele2;
    
    if (!father || !mother) {
        return -1;
    }
            
    int father_pos = -1, mother_pos = -1;
    khiter_t iter = kh_get(ids, sample_ids, father->id);
    if (iter != kh_end(sample_ids)) {
        father_pos = kh_value(sample_ids, iter);
    }
    iter = kh_get(ids, sample_ids, mother->id);
    if (iter != kh_end(sample_ids)) {
        mother_pos = kh_value(sample_ids, iter);
    }

    if (father_pos < 0 || mother_pos < 0) {
        return -1;
    }

    char **sample_data = (char**) record->samples->items;
    char *father_sample = strdup(sample_data[father_pos]);
    char *mother_sample = strdup(sample_data[mother_pos]);

    //LOG_DEBUG_F("Samples: Father = %s\tMother = %s\tChild = %s\n", sample_data[father_pos], sample_data[mother_pos], sample);

    // If any parent's alleles can't be read or is missing, can't decide
    if (get_alleles(father_sample, gt_position, &father_allele1, &father_allele2) ||
        get_alleles(mother_sample, gt_position, &mother_allele1, &mother_allele2)) {
        free(father_sample);
        free(mother_sample);
        return -1;
    }

    // Increment mendelian errors counter when impossible combination is found
    char *aux_chromosome = strndup(record->chromosome, record->chromosome_len);
    if (check_mendel(aux_chromosome, father_allele1, father_allele2, mother_allele1, mother_allele2, 
                     child_allele1, child_allele2, child->sex)) {
        is_error = 1;
    }

    free(father_sample);
    free(mother_sample);
    free(aux_chromosome);
    
    return is_error;
}