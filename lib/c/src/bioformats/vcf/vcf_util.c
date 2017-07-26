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

char *set_field_value_in_info(char *key, char *value, int append, char *info_in, int info_len) {
    assert(info_in);
    assert(key);
    assert(value);
    
    char **splits_info, **splits_key;
    char *aux_string, *key_string, *val_string;
    char *copy_buf;
    int num_splits, num_key_val;
    int find = 0;

    char *info = strndup(info_in, info_len);
    if(strcmp(info, ".") == 0) {
        copy_buf = (char*) calloc(strlen(key) + strlen(value) + 1 + 1, sizeof(char));
        strcpy(copy_buf, key);
        strcat(copy_buf, "=");
        strcat(copy_buf, value);
        
        free(info);
        return copy_buf;
    } else {
        splits_info = split(info, ";", &num_splits);
        copy_buf = (char*) calloc(info_len + strlen(key) + strlen(value) + 1 + 1 + 1, sizeof(char));
        
        for (int i = 0; i < num_splits; i++) {
            aux_string = strdup(splits_info[i]);
            splits_key = split(aux_string, "=", &num_key_val);
            key_string = strdup(splits_key[0]);
            val_string = strdup(splits_key[1]);
            
            if (strcmp(key, key_string) == 0) {
                strcat(copy_buf, key_string);
                strcat(copy_buf, "=");
                if (append) { strcat(copy_buf, val_string); }
                
                if (!strstr(val_string, value)) {
                    if(append) { strcat(copy_buf, ","); }
                    strcat(copy_buf, value);    
                }
                find = 1;
            } else {
                strcat(copy_buf, key_string);
                strcat(copy_buf, "=");
                strcat(copy_buf, val_string);
            }
            if (i < (num_splits - 1)) {
                strcat(copy_buf, ";");
            }

            free(splits_key[0]);
            free(splits_key[1]);
            free(splits_key);
            free(key_string);
            free(val_string);
            free(aux_string);
        }
        if (find == 0) {
            strcat(copy_buf, ";");
            strcat(copy_buf, key);
            strcat(copy_buf, "=");
            strcat(copy_buf, value);
        }
        
        for (int i = 0; i < num_splits; i++) {
            free(splits_info[i]);
        }
        free(splits_info);
        
        free(info);
        return copy_buf;
    }
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

void set_field_value_in_sample(char **sample, int position, char* value) {
    assert(sample);
    assert(value);
    assert(position >= 0);


    int num_splits;

    char **splits = split(*sample, ":", &num_splits);
    *sample = realloc(*sample, (strlen(*sample) + strlen(value) + 1) * sizeof(char));
    strcpy(*sample, "");
    for (int i = 0; i < num_splits; i++) {
        if (i == position) {
            strcat(*sample, value);
        } else {
            strcat(*sample, splits[i]);
        }

        if (i < (num_splits - 1)) {
            strcat(*sample, ":");
        }
    }
    for (int i = 0; i < num_splits; i++) {
        free(splits[i]);
    }
    free(splits);
}


enum alleles_code get_alleles(char* sample, int genotype_position, int* allele1, int* allele2) {
    assert(sample);
    assert(allele1);
    assert(allele2);
    
    char *aux_buffer, *allele, *genotype = NULL;
    enum alleles_code ret_code = ALLELES_OK;
    int cur_pos = -1;
    
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
    
    if (genotype == NULL) {     // calling strtok with NULL would not mean what we want
        *allele1 = -1;
        *allele2 = -1;
        ret_code = ALL_ALLELES_MISSING;
        return ret_code;    //if genotype is NULL, the previous while() was skipped, and we don't have the string with the alleles
    }
    
//     LOG_DEBUG_F("genotype = %s\n", genotype);
    
    allele = strtok_r(genotype, "/|", &aux_buffer);
    if (strcmp(".", allele) == 0) {
        *allele1 = -1;
        ret_code = FIRST_ALLELE_MISSING;
    } else {
        *allele1 = atoi(allele);
    }
    
//     LOG_DEBUG_F("allele 1 = %s\n", allele);
    
    if (strlen(aux_buffer) == 0) { // Haploid
        *allele2 = -1;
        ret_code = HAPLOID;
        return ret_code;
    }
    
    allele = strtok_r(NULL, "/|", &aux_buffer);
    if (strcmp(".", allele) == 0) {
        *allele2 = -1;
        ret_code = (ret_code == FIRST_ALLELE_MISSING) ? ALL_ALLELES_MISSING : SECOND_ALLELE_MISSING;
    } else {
        *allele2 = atoi(allele);
    }
    
//     LOG_DEBUG_F("allele2 = %s\n", allele);
    
    return ret_code;
}

