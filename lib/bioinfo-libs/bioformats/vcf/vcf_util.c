#include "vcf_util.h"

size_t count_regions(char *regions_string) {
    assert(regions_string);
    
    size_t num_regions = 0;
    char *aux = regions_string;
    while (*aux) {
        if (*aux++ == ',') ++num_regions;
    }
    return ++num_regions;
}

char *get_field_value_in_info(const char *field, char *info) {
    assert(field);
    assert(info);
    
    char *save_strtok, *token;
    char *value;
    
    // Search for field in info (must begin with '<field>=')
    token = strtok_r(info, ";", &save_strtok);
    while (token != NULL && 
	   (!starts_with(token, field) || // Field is completely different than the requested one
	   (strlen(token) > strlen(field) && token[strlen(field)] != '='))) { // Same beginning, but longer than the requested one
        token = strtok_r(NULL, ";", &save_strtok);
    }
    
    if (token == NULL) {
        return NULL;  // Field not found
    } else {
        value = strtok_r(token, ";", &save_strtok);
    }
    
    // Search for the field value
    value = strtok_r(token, "=", &save_strtok);
    if (!strcmp(value, field)) {
//         free(value);
        value = strtok_r(NULL, "=", &save_strtok);
    }
    
    return value;
}

int get_field_position_in_format(const char *field, char *format) {
    assert(field);
    assert(format);
    
    int field_pos = 0;
    char *save_strtok, *token;
    token = strtok_r(format, ":", &save_strtok);
    while (token != NULL && strcmp(token, field)) {
        token = strtok_r(NULL, ":", &save_strtok);
        field_pos++;
    }
    
    return (token == NULL) ? -1 : field_pos;
}

char *get_field_value_in_sample(char *sample, int position) {
    assert(sample);
    assert(position >= 0);
    
    int field_pos = 0;
    char *save_strtok, *token = NULL;
    token = strtok_r(sample, ":", &save_strtok);
    while (token != NULL && field_pos < position) {
        token = strtok_r(NULL, ":", &save_strtok);
        field_pos++;
    }
    
    return token;
}

int get_alleles(char* sample, int genotype_position, int* allele1, int* allele2) {
    assert(sample);
    assert(allele1);
    assert(allele2);
    
    char *aux_buffer, *allele, *genotype;
    int ret_code = 0, cur_pos = -1;
    
    // If the sample is not in a form that can contain a genotype, mark both alleles as missing
    if (strnlen(sample, 3) < 3) {
        *allele1 = -1;
        *allele2 = -1;
        ret_code = 3;
        return ret_code;
    }
    
    // Process a correct sample, at least in length terms
    
    while (cur_pos < genotype_position) {
        genotype = strtok_r(sample, ":", &aux_buffer);
        sample = NULL;
        cur_pos++;
    }
    
//     LOG_DEBUG_F("genotype = %s\n", genotype);
    
    allele = strtok_r(genotype, "/|", &aux_buffer);
    if (strcmp(".", allele) == 0) {
        *allele1 = -1;
        ret_code += 1;
    } else {
        *allele1 = atoi(allele);
    }
    
//     LOG_DEBUG_F("allele 1 = %s\n", allele);
    
    allele = strtok_r(NULL, "/|", &aux_buffer);
    if (strcmp(".", allele) == 0) {
        *allele2 = -1;
        ret_code += 2;
    } else {
        *allele2 = atoi(allele);
    }
    
//     LOG_DEBUG_F("allele2 = %s\n", allele);
    
    return ret_code;
}
