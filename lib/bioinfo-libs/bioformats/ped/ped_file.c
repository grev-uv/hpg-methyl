#include "ped_file.h"


//====================================================================================
//  ped_file.c
//  ped file management functions
//====================================================================================

//-----------------------------------------------------
// ped_open
//-----------------------------------------------------


ped_file_t *ped_open(char *filename) {
    size_t len;
    char *data = mmap_file(&len, filename);

    ped_file_t *ped_file = (ped_file_t *) malloc(sizeof(ped_file_t));
    ped_file->filename = filename;
    ped_file->data = data;
    ped_file->data_len = len;
    
    ped_file->families = cp_hashtable_create_by_option(COLLECTION_MODE_DEEP,
                                                       10,
                                                       cp_hash_istring,         // Hash function
                                                       (cp_compare_fn) strcasecmp,     // Key comparison function
                                                       NULL,                    // Key copy function
                                                       NULL,                    // Key destructor function
                                                       NULL,                    // Value copy function
                                                       (cp_destructor_fn) family_free // Value destructor function
                                                      );
    return ped_file;
}


//-----------------------------------------------------
// ped_close and memory freeing
//-----------------------------------------------------

void ped_close(ped_file_t *ped_file, int free_families) {
    // Free families if asked to
    if (free_families) {
        cp_hashtable_destroy(ped_file->families);
    }
    
    munmap((void*) ped_file->data, ped_file->data_len);
    free(ped_file);
}

void ped_record_free(ped_record_t* ped_record) {
    free(ped_record->family_id);
    free(ped_record->individual_id);
    if (ped_record->father_id) { free(ped_record->father_id); }
    if (ped_record->mother_id) { free(ped_record->mother_id); }
    free(ped_record);
}

//-----------------------------------------------------
// I/O operations (read and write) in various ways
//-----------------------------------------------------

int ped_read(ped_file_t *ped_file) {
    int ret_code = 0;
    list_t *ped_batches = (list_t*) malloc (sizeof(list_t));
    list_init("batches", 1, 100, ped_batches);
    
#pragma omp parallel sections
{
#pragma omp section
    {
        ret_code = ped_read_batches(ped_batches, 1000, ped_file);
        if (ret_code) {
            LOG_FATAL_F("Error %d while reading the file %s\n", ret_code, ped_file->filename);
        }
        
        list_decr_writers(ped_batches);
    }
    
#pragma omp section
    {
        list_item_t *item = NULL;
        list_item_t *record_item = NULL;
        while ((item = list_remove_item(ped_batches)) != NULL) {
            ped_batch_t *batch = (ped_batch_t*) item->data_p;
            record_item = batch->first_p;
            while (record_item) {
                ret_code &= add_ped_record(record_item->data_p, ped_file);
                record_item = record_item->next_p;
            }
            ped_batch_free(batch);
            list_item_free(item);
        }
    }
}
    list_free_deep(ped_batches, NULL);
    
    return ret_code;
}

int ped_read_batches(list_t *batches_list, size_t batch_size, ped_file_t *ped_file) {
    return ped_ragel_read(batches_list, batch_size, ped_file);
}

int ped_write(ped_file_t *ped_file, char *filename) {
    FILE *fd = fopen(filename, "w");
    if (fd < 0) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return 1;
    }
    
    if (ped_write_to_file(ped_file, fd)) {
        fprintf(stderr, "Error writing file: %s\n", filename);
        fclose(fd);
        return 2;
    }
    
    fclose(fd);
    return 0;
}

//-----------------------------------------------------
// load data into the ped_file_t
//-----------------------------------------------------

int add_family(family_t* family, ped_file_t* ped_file) {
    return cp_hashtable_put(ped_file->families, family->id, family) == NULL;
}

int get_num_families(ped_file_t* ped_file) {
    assert(ped_file);
    return cp_hashtable_count(ped_file->families);
}

int add_ped_record(ped_record_t* record, ped_file_t *ped_file) {
    assert(record);
    assert(ped_file);
    
    int result = 0;
    individual_t *father = NULL, *mother = NULL, *individual = NULL;
    enum Condition condition = MISSING_CONDITION;
    
    // Get family or, should it not exist yet, create it
    family_t *family = cp_hashtable_get(ped_file->families, record->family_id);
    if (family == NULL) {
        family = family_new(strdup(record->family_id));
        if (add_family(family, ped_file)) {
            return ALREADY_EXISTING_FAMILY;
        }
    }
    
    LOG_DEBUG_F("family id = %s\tindiv id = %s\tfather id = %s\tmother id = %s\n", record->family_id, record->individual_id, record->father_id, record->mother_id);
    
    // If it is an ancestor with no sex defined, add to the list of unknown members
    if (!record->father_id && !record->mother_id && record->sex == UNKNOWN_SEX) {
        condition = get_condition_from_phenotype(record->phenotype);
        individual = individual_new(strdup(record->individual_id), record->phenotype, record->sex, condition, NULL, NULL, family);
        return family_add_unknown(individual, family);
    }
    
    // Get parents from family or, should they not exist yet, create them
    if (!family->father) {
        // Non-existing father, set his ID from the record (if available)
        if (record->father_id) {
            LOG_DEBUG_F("Set family %s father", family->id);
            father = individual_new(strdup(record->father_id), -9, MALE, MISSING_CONDITION, NULL, NULL, family);
            family_set_parent(father, family);
        }
    
    } else if (record->father_id && strcasecmp(family->father->id, record->father_id)) {
        // Father already exists and IDs do not match, error!
        return FATHER_APPEARS_MORE_THAN_ONCE;
    
    } else if (!strcasecmp(family->father->id, record->individual_id)) {
        // The father was created while reading one of his children
        // Should his status be 'missing', fill the phenotypical information
        father = family->father;
        LOG_DEBUG_F("Father already found, condition = %d\n", father->condition);
        
        // If the father struct members are missing, fill them
        if (father->condition == MISSING_CONDITION) {
            father->phenotype = record->phenotype;
            father->condition = get_condition_from_phenotype(father->phenotype);
            LOG_DEBUG_F("Father modified, condition = %d\n", father->condition);
        }
        return 0;   // Nothing more to do, he already belongs to the family
    }
    
    if (!family->mother) {
        // Non-existing mother, set his ID from the record (if available)
        if (record->mother_id) {
            LOG_DEBUG_F("Set family %s mother", family->id);
            mother = individual_new(strdup(record->mother_id), -9, FEMALE, MISSING_CONDITION, NULL, NULL, family);
            family_set_parent(mother, family);
        }
    
    } else if (record->mother_id && strcasecmp(family->mother->id, record->mother_id)) {
        // Mother already exists and IDs do not match, error!
        return MOTHER_APPEARS_MORE_THAN_ONCE;
    
    } else if (!strcasecmp(family->mother->id, record->individual_id)) {
        // The mother was created while reading one of his children
        // Should his status be 'missing', fill the phenotypical information
        mother = family->mother;
        LOG_DEBUG_F("Mother already found, condition = %d\n", mother->condition);
        
        // If the mother struct members are missing, fill them
        if (mother->condition == MISSING_CONDITION) {
            mother->phenotype = record->phenotype;
            mother->condition = get_condition_from_phenotype(mother->phenotype);
            LOG_DEBUG_F("Mother modified, condition = %d\n", mother->condition);
        }
        return 0;   // Nothing more to do, he already belongs to the family
    }
    
    // Create individual with the information extracted from the PED record
    condition = get_condition_from_phenotype(record->phenotype);
    individual = individual_new(strdup(record->individual_id), record->phenotype, record->sex, condition, father, mother, family);
    if (father || mother) {
        LOG_DEBUG_F("** add family %s child (id %s)\n", family->id, individual->id);
        family_add_child(individual, family);
    } else {
        LOG_DEBUG_F("** set family %s parent of sex %d (id %s)\n", family->id, individual->sex, individual->id);
        result = family_set_parent(individual, family);
        if (result == 1) {
            result = FATHER_APPEARS_MORE_THAN_ONCE;
        } else if (result == 2) {
            result = MOTHER_APPEARS_MORE_THAN_ONCE;
        }
    }

    return result;
}

