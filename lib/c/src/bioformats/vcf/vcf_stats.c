#include "vcf_stats.h"

static int is_mendelian_error(individual_t *father, individual_t *mother, individual_t *child, 
                              int child_allele1, int child_allele2, int gt_position, 
                              vcf_record_t *record, khash_t(ids) *sample_ids);
void hardy_weinberg_test(hardy_weinberg_stats_t *hw);

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
 
variant_stats_t* variant_stats_new(char *chromosome, unsigned long position, char *ref_allele, char *alt_alleles, int num_phenotypes) {
    assert(chromosome);
    assert(ref_allele);
    
    variant_stats_t *stats = (variant_stats_t*) malloc (sizeof(variant_stats_t));
    stats->chromosome = chromosome;
    stats->position = position;
    stats->ref_allele = ref_allele;
    stats->alt_alleles = alt_alleles;
    stats->alternates = NULL;
    stats->maf_allele = NULL;
    stats->mgf_genotype = NULL;
    
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
    
    stats->num_phenotypes = num_phenotypes;
    memset(&(stats->hw_all),0, sizeof(hardy_weinberg_stats_t));
    stats->pheno_stats = (phenotype_stats_t*)calloc(num_phenotypes , sizeof(phenotype_stats_t));

    return stats;
}

void variant_stats_free(variant_stats_t* stats) {
    assert(stats);
    
    if (stats->chromosome) { free(stats->chromosome); }
    if (stats->ref_allele) { free(stats->ref_allele); }
    if (stats->alt_alleles) { free(stats->alt_alleles); }
    if (stats->mgf_genotype) { free(stats->mgf_genotype); }
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
    if (stats->pheno_stats) {
        for(int i = 0; i < stats->num_phenotypes; i++) {
            if (stats->pheno_stats[i].alleles_count)   { free(stats->pheno_stats[i].alleles_count); }
            if (stats->pheno_stats[i].genotypes_count) { free(stats->pheno_stats[i].genotypes_count); }
            if (stats->pheno_stats[i].alleles_freq)    { free(stats->pheno_stats[i].alleles_freq); }
            if (stats->pheno_stats[i].genotypes_freq)  { free(stats->pheno_stats[i].genotypes_freq); }
        }
        free(stats->pheno_stats); 
    }
    free(stats);
}

int get_variants_stats_old(vcf_record_t **variants, int num_variants, individual_t **individuals, khash_t(ids) *sample_ids,
                        list_t *output_list, file_stats_t *file_stats) {
    get_variants_stats(variants, num_variants, individuals, sample_ids, 0, output_list, file_stats);
    return 0;
}

int get_variants_stats(vcf_record_t **variants, int num_variants, individual_t **individuals, khash_t(ids) *sample_ids,
                        int num_variables, list_t *output_list, file_stats_t *file_stats) {
    assert(variants);
    assert(output_list);
    assert(file_stats);
    
    char *copy_buf = NULL, *sample;
    
    int num_alternates, gt_position, curr_position = 0;
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
    float maf = INT_MAX, mgf = INT_MAX;
    float cur_gt_freq;
    
    // Struct for temporary variables for each phenotype stats
    struct phenotype_stats_var_count {
        int samples_num;
        
        int total_alleles_count;   
        int total_genotypes_count;  
    
        float cases_dominant;
        float controls_dominant;
        float cases_recessive;
        float controls_recessive;
    }   *pheno_count = NULL, *aux_pheno_count;
    
    if (num_variables > 0) {
        pheno_count = malloc(sizeof(struct phenotype_stats_var_count)*num_variables);
    }
    
    // Variant stats management
    vcf_record_t *record;
    variant_stats_t *stats;
    phenotype_stats_t *aux_pheno_stats;
    for (int i = 0; i < num_variants; i++) {
        record = variants[i];
        stats = variant_stats_new(strndup(record->chromosome, record->chromosome_len), 
                                  record->position, 
                                  strndup(record->reference, record->reference_len),
                                  strndup(record->alternate, record->alternate_len),
                                  num_variables);
        // Reset counters
        total_alleles_count = total_genotypes_count = 0;
        cases_dominant = controls_dominant = cases_recessive = controls_recessive = 0;
        maf = mgf = INT_MAX;
        memset(pheno_count, 0, sizeof(struct phenotype_stats_var_count)*(num_variables));
        
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
        
        for (int j = 0; j < num_variables; j++) {
            aux_pheno_stats = &stats->pheno_stats[j];
            
            aux_pheno_stats->num_alleles = stats->num_alleles;
            aux_pheno_stats->alleles_count = (int*) calloc (stats->num_alleles, sizeof(int));
            aux_pheno_stats->genotypes_count = (int*) calloc (stats->num_alleles * stats->num_alleles, sizeof(int));
            aux_pheno_stats->alleles_freq = (float*) calloc (stats->num_alleles, sizeof(float));
            aux_pheno_stats->genotypes_freq = (float*) calloc (stats->num_alleles * stats->num_alleles, sizeof(float));
        }
        
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
            if (alleles_code == ALLELES_OK) {
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
                
                // Counting genotypes for Hardy-Weinberg (all phenotypes)
                if (!allele1 && !allele2) { // 0|0
                    stats->hw_all.n_AA++;
                } else if ((!allele1 && allele2==1) || (allele1==1 && !allele2)) { // 0|1, 1|0
                    stats->hw_all.n_Aa++;
                } else if(allele1==1 && allele2==1){ // 1|1
                    stats->hw_all.n_aa++;
                }
                
            } else if (alleles_code == HAPLOID) {
                // Haploid (chromosomes X/Y)
                stats->alleles_count[allele1]++;
                total_alleles_count++;
                
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
                // Check mendelian errors (pedigree data must be provided)
                if (sample_ids && (alleles_code == ALLELES_OK || alleles_code == HAPLOID)) {
                    if (is_mendelian_error(individuals[j]->father, individuals[j]->mother, individuals[j], 
                                           allele1, allele2, gt_position, record, sample_ids) > 0) {
                        (stats->mendelian_errors)++;
                    }
                }

                if (alleles_code == ALLELES_OK) {
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
                
                // Stats only for variables in variable_group
                if(num_variables > 0 && individuals[j]->variable >= 0) {
                    aux_pheno_stats = &(stats->pheno_stats[individuals[j]->variable]);
                    aux_pheno_count = &(pheno_count[individuals[j]->variable]);
                    aux_pheno_count->samples_num++;
                    if(alleles_code == ALLELES_OK) {
                        aux_pheno_count->total_alleles_count += 2;
                        aux_pheno_stats->alleles_count[allele1] += 1;
                        aux_pheno_stats->alleles_count[allele2] += 1;
                        aux_pheno_stats->genotypes_count[curr_position] += 1;
                        aux_pheno_count->total_genotypes_count++;
                        
                        // Check inheritance models
                        if (individuals[j]->condition == UNAFFECTED) {
                            if (!allele1 && !allele2) { // 0|0
                                aux_pheno_count->controls_dominant++;
                                aux_pheno_count->controls_recessive++;
                            } else if ((!allele1 && allele2) || (allele1 && !allele2)) { // 0|1 or 1|0
                                aux_pheno_count->controls_recessive++;
                            }
                        } else if (individuals[j]->condition == AFFECTED) {
                            if (allele1 && allele2 && allele1 == allele2) { // 1|1, 2|2, and so on
                                aux_pheno_count->cases_recessive++;
                                aux_pheno_count->cases_dominant++;
                            } else if (allele1 || allele2) { // 0|1, 1|0, 1|2, 2|1, 1|3, and so on
                                aux_pheno_count->cases_dominant++;
                            }
                        }
                        
                        // Counting genotypes for Hardy-Weinberg
                        if (!allele1 && !allele2) { // 0|0
                            aux_pheno_stats->hw.n_AA++;
                        } else if ((!allele1 && allele2==1) || (allele1==1 && !allele2)) { // 0|1, 1|0
                            aux_pheno_stats->hw.n_Aa++;
                        } else if(allele1==1 && allele2==1){ // 1|1
                            aux_pheno_stats->hw.n_aa++;
                        }
                        
                    } else if (alleles_code == HAPLOID) {
                        // Haploid (chromosomes X/Y)
                        aux_pheno_stats->alleles_count[allele1]++;
                        aux_pheno_count->total_alleles_count++;
                        aux_pheno_count->total_genotypes_count++;
                        
                    } else {
                        // Missing genotype (one or both alleles missing)
                        aux_pheno_stats->missing_genotypes++;
                        if (allele1 < 0) { 
                            aux_pheno_stats->missing_alleles++; 
                        } else {
                            aux_pheno_count->total_alleles_count++;
                            aux_pheno_stats->alleles_count[allele1]++;
                        }
                            
                        if (allele2 < 0) { 
                            aux_pheno_stats->missing_alleles++;
                        } else {
                            aux_pheno_count->total_alleles_count++;
                            aux_pheno_stats->alleles_count[allele2]++;
                        }
                    }
                }
            }
        }   // Finish all samples loop
        
        assert(cases_dominant >= cases_recessive);
        assert(controls_recessive >= controls_dominant);
        assert(controls_dominant <= record->samples->size - stats->missing_genotypes);
        assert(cases_dominant <= record->samples->size - stats->missing_genotypes);
        assert(controls_recessive <= record->samples->size - stats->missing_genotypes);
        assert(cases_recessive <= record->samples->size - stats->missing_genotypes);
        
        // Once all samples have been traverse, calculate % that follow inheritance model
        stats->controls_percent_dominant = (float) controls_dominant * 100 / (record->samples->size - stats->missing_genotypes);
        stats->cases_percent_dominant = (float) cases_dominant * 100 / (record->samples->size - stats->missing_genotypes);
        stats->controls_percent_recessive = (float) controls_recessive * 100 / (record->samples->size - stats->missing_genotypes);
        stats->cases_percent_recessive = (float) cases_recessive * 100 / (record->samples->size - stats->missing_genotypes);
        
        // Get allele and genotype frequencies, as well as MAF and MGF
        for (int j = 0; j < stats->num_alleles; j++) {
            stats->alleles_freq[j] = (total_alleles_count > 0) ? (float) stats->alleles_count[j] / total_alleles_count : 0;
            if (stats->alleles_freq[j] < maf) {
                maf = stats->alleles_freq[j];
                stats->maf_allele = (j == 0) ? stats->ref_allele : stats->alternates[j-1];
            }
        }
        stats->maf = maf;
        
        for (int j = 0; j < stats->num_alleles * stats->num_alleles; j++) {
            stats->genotypes_freq[j] = (total_genotypes_count > 0) ? (float) stats->genotypes_count[j] / total_genotypes_count : 0;
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
        
        if(individuals) {
            for(int pheno_iter = 0; pheno_iter < num_variables; pheno_iter++) {
                aux_pheno_stats = &(stats->pheno_stats[pheno_iter]);
                aux_pheno_count = &(pheno_count[pheno_iter]);
                
                aux_pheno_stats->controls_percent_dominant = (float) aux_pheno_count->controls_dominant * 100 / (aux_pheno_count->samples_num - aux_pheno_stats->missing_genotypes);
                aux_pheno_stats->cases_percent_dominant = (float) aux_pheno_count->cases_dominant * 100 / (aux_pheno_count->samples_num - aux_pheno_stats->missing_genotypes);
                aux_pheno_stats->controls_percent_recessive = (float) aux_pheno_count->controls_recessive * 100 / (aux_pheno_count->samples_num - aux_pheno_stats->missing_genotypes);
                aux_pheno_stats->cases_percent_recessive = (float) aux_pheno_count->cases_recessive * 100 / (aux_pheno_count->samples_num - aux_pheno_stats->missing_genotypes);
                
                maf = mgf = INT_MAX;
                for (int j = 0; j < stats->num_alleles; j++) {
                    aux_pheno_stats->alleles_freq[j] = (aux_pheno_count->total_alleles_count > 0) ? (float) aux_pheno_stats->alleles_count[j] / aux_pheno_count->total_alleles_count : 0;
                    if (aux_pheno_stats->alleles_freq[j] < maf) {
                        maf = aux_pheno_stats->alleles_freq[j];
                    }
                }
                aux_pheno_stats->maf = maf;
            
                if(aux_pheno_count->total_genotypes_count > 0)
                    for (int j = 0; j < stats->num_alleles * stats->num_alleles; j++) {
                        aux_pheno_stats->genotypes_freq[j] = (float) aux_pheno_stats->genotypes_count[j] / aux_pheno_count->total_genotypes_count;
                    }
                
                for (int j = 0; j < aux_pheno_stats->num_alleles; j++) {
                    for (int k = j; k < aux_pheno_stats->num_alleles; k++) {
                        int idx1 = j * aux_pheno_stats->num_alleles + k;
                        if (j == k) {
                            cur_gt_freq = aux_pheno_stats->genotypes_freq[idx1];
                        } else {
                            int idx2 = k * stats->num_alleles + j;
                            cur_gt_freq = aux_pheno_stats->genotypes_freq[idx1] + aux_pheno_stats->genotypes_freq[idx2];
                        }
                        if (cur_gt_freq < mgf) 
                            mgf = cur_gt_freq;
                        
                    }
                }
                aux_pheno_stats->mgf = mgf;
            }
        }
        
        
        // Testing for Hardy-Weinberg Equilibrium (HWE)
        //printf("Start hwe for the %d variant\n",i);
        hardy_weinberg_test(&stats->hw_all);
        for (int j = 0; j < num_variables; j++) {
            hardy_weinberg_test(&(stats->pheno_stats[j].hw));
        }//printf("\n\n");
        
        
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
         * 4 possibilities for being an INDEL:
         * - The value of the ALT field is <DEL> or <INS>
         * - The REF allele is not . but the ALT is
         * - The REF allele is . but the ALT is not
         * - The REF field length is different than the ALT field length
         * When checking REF vs ALT, all alternatives must be traversed
         */
        stats->is_indel = 0;
        if ((strncmp(".", stats->ref_allele, 1) && !strncmp(".", record->alternate, 1)) ||
            (strncmp(".", record->alternate, 1) && !strncmp(".", stats->ref_allele, 1)) ||
            !strncmp("<INS>", record->alternate, record->alternate_len) ||
            !strncmp("<DEL>", record->alternate, record->alternate_len)) {
            stats->is_indel = 1;
            indels_count++;
        } else {
            for (int j = 0; j < stats->num_alleles - 1; j++) {
                size_t alt_len = strlen(stats->alternates[j]);
                if (record->reference_len != alt_len) {
                    stats->is_indel = 1;
                    indels_count++;
                    break;
                }
            }
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
    
    free(pheno_count);
    
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
    stats->homozygotes_number = 0;
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
            if (alleles_code != ALLELES_OK) {
                // Missing genotype (one or both alleles missing)
                #pragma omp atomic
                (sample_stats[j]->missing_genotypes)++;
                continue;
            }
            
            // Check mendelian errors
            if (individuals && sample_ids && alleles_code == ALLELES_OK &&
                    is_mendelian_error(individuals[j]->father, individuals[j]->mother, individuals[j], 
                                   allele1, allele2, gt_position, record, sample_ids) > 0) {
                #pragma omp atomic
                (sample_stats[j]->mendelian_errors)++;
            }
            
            //Count homozygotes
            if(allele1 == allele2)
            {
                #pragma omp atomic
                (sample_stats[j]->homozygotes_number)++;
            }
        }
        
    }
    
    return 0;
}


/* ******************************
 *      Auxiliary functions     *
 * ******************************/

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

    //LOG_DEBUG_F("Samples: Father = %s\tMother = %s\tChild = %d/%d\n", sample_data[father_pos], sample_data[mother_pos], child_allele1, child_allele2);

    // If any parent's alleles can't be read or is missing, can't decide
    if (get_alleles(father_sample, gt_position, &father_allele1, &father_allele2) != ALLELES_OK ||
        get_alleles(mother_sample, gt_position, &mother_allele1, &mother_allele2) != ALLELES_OK) {
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

void hardy_weinberg_test(hardy_weinberg_stats_t *hw) {
    int n = hw->n = hw->n_AA + hw->n_Aa + hw->n_aa;
    int n_AA = hw->n_AA;
    int n_Aa = hw->n_Aa;
    int n_aa = hw->n_aa;
    if (n) {
        float p = hw->p = (2.0*n_AA+n_Aa)/(2*n);
        float q = hw->q = 1-hw->p;
        
        hw->e_AA = (p*p*n);
        hw->e_Aa = (2*p*q*n);
        hw->e_aa = (q*q*n);
        
        /* * *///printf("Observed Values: %d/%d/%d\n", n_AA, n_Aa, n_aa,n);

        //printf("\ni=%d, n = %d, p= %f, q= %f\n",i, n, p, q);
             //   printf("Expected Values: %.2f,\t   %.2f,\t   %.2f,\te_n = %.4f\n", hw->e_AA, hw->e_Aa, hw->e_aa,hw->e_AA+hw->e_Aa+hw->e_aa);
        LOG_DEBUG_F("O(HET) %f,\tE(HET) %f\n", ((float)n_Aa)/n, (float)hw->e_Aa/n);
        
        if( hw->e_AA == n_AA) hw->e_AA = n_AA = 1;
        if( hw->e_Aa == n_Aa) hw->e_Aa = n_Aa = 1;
        if( hw->e_aa == n_aa) hw->e_aa = n_aa = 1;
        
        hw->chi2 = (n_AA - hw->e_AA)*(n_AA - hw->e_AA)/hw->e_AA 
                 + (n_Aa - hw->e_Aa)*(n_Aa - hw->e_Aa)/hw->e_Aa
                 + (n_aa - hw->e_aa)*(n_aa - hw->e_aa)/hw->e_aa; 
        hw->p_value = 1-gsl_cdf_chisq_P(hw->chi2,1);
        
        //printf("Expected Values: %d, %d, %d, e_n = %d\n", e_AA, e_Aa, e_aa,e_AA+e_Aa+e_aa);
        /*printf("CHI %f\t", hw->chi2);
        printf("p-val %f\n", hw->p_value);*/
    }
    
}
