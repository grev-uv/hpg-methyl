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
    ped_file->people = kh_init(family_members);
    
    ped_file->variables = kh_init(str);
    ped_file->num_variables = 0;
    ped_file->accept_new_values = 1;
    
    ped_file->affected = NULL;
    ped_file->unaffected = NULL;
    set_unaffected_phenotype("1", ped_file);
    set_affected_phenotype("2", ped_file);
    
    ped_file->variable_field = NULL;
    ped_file->num_field = 0;
    
    return ped_file;
}

//-----------------------------------------------------
// ped_close and memory freeing
//-----------------------------------------------------

void ped_close(ped_file_t *ped_file, int free_families, int free_phenotype) {
    // Free families if asked to
    if (free_families) {
        cp_hashtable_destroy(ped_file->families);
        kh_destroy(family_members, ped_file->people);
    }
    // Free phenotype hash if asked to
    if (free_phenotype) {
        kh_destroy(str,ped_file->variables);
    } 
    
    munmap((void*) ped_file->data, ped_file->data_len);
    if(ped_file->affected) free(ped_file->affected);
    if(ped_file->unaffected) free(ped_file->unaffected);
    if(ped_file->variable_field) free(ped_file->variable_field);
    free(ped_file);
}

void ped_record_free(ped_record_t* ped_record) {
    free(ped_record->family_id);
    free(ped_record->individual_id);
    if (ped_record->father_id) { free(ped_record->father_id); }
    if (ped_record->mother_id) { free(ped_record->mother_id); }
    if (ped_record->phenotype) { free(ped_record->phenotype); }
    if (ped_record->custom_field) { free(ped_record->custom_field); }
    free(ped_record);
}

//-----------------------------------------------------
// I/O operations (read and write) in various ways
//-----------------------------------------------------

int ped_read(ped_file_t *ped_file) {
    int ret_code = 0;
    list_t *ped_batches = (list_t*) malloc (sizeof(list_t));
    list_init("batches", 1, 100, ped_batches);
    
    // Read the whole file in batches
    // Do not link children with their parents, only sets their IDs
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
                ret_code = add_ped_record(record_item->data_p, ped_file);
                if (ret_code > 0) {
                    LOG_FATAL_F("%s - %s\n", ((ped_record_t*) record_item->data_p)->family_id, get_ped_semantic_error_msg(ret_code));
                }
                record_item = record_item->next_p;
            }
            ped_batch_free(batch);
            list_item_free(item);
        }
    }
}
    
    // Link children with their parents
    family_t **families = (family_t**) cp_hashtable_get_values(ped_file->families);
    int num_families = get_num_families(ped_file);
    
    for (int f = 0; f < num_families; f++) {
        family_t *family = families[f];
        for (int k = kh_begin(family->members); k < kh_end(family->members); k++) {
            if (!kh_exist(family->members, k)) {
                continue;
            }

            individual_t *individual = kh_value(family->members, k);
            if (strcmp(individual->father_id, "0")) {
                khiter_t iter = kh_get(family_members, family->members, individual->father_id);
                // Check the father is also included in the PED file
                if (iter != kh_end(family->members)) {
                    individual_t *father = kh_value(family->members, iter);
                    if (father->sex != MALE) {  // Check the father is, in fact, a man
                        LOG_FATAL_F("Father %s of individual %s is not a man\n", father->id, individual->id);
                    }
                    
                    individual->father = father;
                    individual_add_child(individual, father);
                } else {
                    LOG_WARN_F("Child %s references father %s who does not exist in the PED file\n",
                               individual->id, individual->father_id);
                }
            }
            
            if (strcmp(individual->mother_id, "0")) {
                khiter_t iter = kh_get(family_members, family->members, individual->mother_id);
                // Check the mother is also included in the PED file
                if (iter != kh_end(family->members)) {
                    individual_t *mother = kh_value(family->members, iter);
                    if (mother->sex != FEMALE) {  // Check the mother is, in fact, a woman
                        LOG_FATAL_F("Mother %s of individual %s is not a woman\n", mother->id, individual->id);
                    }
                    
                    individual->mother = mother;
                    individual_add_child(individual, mother);
                } else {
                    LOG_WARN_F("Child %s references mother %s who does not exist in the PED file\n",
                               individual->id, individual->mother_id);
                }
            }
        }
    }
    
    // Only for testing purposes
    // ped_write(ped_file, "/tmp/variant/pedoutput.ped");
    
    free(families);
    
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
    assert(family);
    assert(ped_file);
    return cp_hashtable_put(ped_file->families, family->id, family) == NULL;
}

int get_num_families(ped_file_t* ped_file) {
    assert(ped_file);
    return cp_hashtable_count(ped_file->families);
}

int ped_add_individual(individual_t *individual, ped_file_t *ped_file) {
    assert(individual);
    assert(ped_file);
    int ret = 0;
    khiter_t iter = kh_put(family_members, ped_file->people, strdup(individual->id), &ret);
    if (ret) {
        kh_value(ped_file->people, iter) = individual;
        return 0;
    } else {
        // Do not accept repeated entries
        LOG_FATAL_F("Individual %s appears more than once in the PED file\n", individual->id);
    }
    return 1;
}

int add_ped_record(ped_record_t* record, ped_file_t *ped_file) {
    assert(record);
    assert(ped_file);
    
    // Get family or, should it not exist yet, create it
    family_t *family = cp_hashtable_get(ped_file->families, record->family_id);
    if (family == NULL) {
        family = family_new(strdup(record->family_id));
        if (add_family(family, ped_file)) {
            return ALREADY_EXISTING_FAMILY;
        }
    }
    
    LOG_DEBUG_F("family id = %s\tindiv id = %s\tfather id = %s\tmother id = %s\n", 
                record->family_id, record->individual_id, record->father_id, record->mother_id);

    if (!strcmp(record->individual_id, record->father_id)) {
        LOG_FATAL_F("Individual %s has the same ID as his/her father\n", record->individual_id);
    } else if (!strcmp(record->individual_id, record->mother_id)) {
        LOG_FATAL_F("Individual %s has the same ID as his/her mother\n", record->individual_id);
    }
    
    // Create individual
    enum Condition condition = get_condition_from_phenotype(record->phenotype, ped_file);
    individual_t *individual = individual_new_ids_only(strdup(record->individual_id), record->var_index, record->sex, condition, 
                                                       strdup(record->father_id), strdup(record->mother_id), family);
    
    // Add individual to pedigree file
    ped_add_individual(individual, ped_file);
    
    // If it is an ancestor with no sex defined, add to the list of unknown members in the family
    if (!strcmp(record->father_id, "0") && !strcmp(record->mother_id, "0") && record->sex == UNKNOWN_SEX) {
        return family_add_unknown(individual, family);
    }
    
    // Otherwise, add as a normal member
    return family_add_member(individual, family);
}


family_t **ped_flatten_families(ped_file_t *ped_file, int *num_families) {
    khash_t(family) *hash_families = kh_init(family);
    
    family_t **families = (family_t**) cp_hashtable_get_values(ped_file->families);
    for (int f = 0; f < get_num_families(ped_file); f++) {
        family_t *family = families[f];
        for (int k = kh_begin(family->members); k < kh_end(family->members); k++) {
            if (!kh_exist(family->members, k)) { continue; }

            individual_t *individual = kh_value(family->members, k);
            if (individual->father && individual->mother) {
                // Compose the key in the hashtable as "fatherid_motherid"
                char *family_key = calloc(strlen(individual->father_id) + strlen(individual->mother_id) + 2, sizeof(char));
                sprintf(family_key, "%s_%s", individual->father_id, individual->mother_id);
                
                // If the family has already been flattened, add the child as new member
                // Otherwise, create the family and add both parents and child
                khiter_t fam_iter = kh_get(family, hash_families, family_key);
                family_t* flat_family = NULL;
                
                if (fam_iter != kh_end(hash_families)) {
                    flat_family = kh_value(hash_families, fam_iter);
                    family_add_member(individual, flat_family);
                } else {
                    // Configure family
                    flat_family = family_new(individual->family->id);
                    family_add_member(individual->father, flat_family);
                    family_add_member(individual->mother, flat_family);
                    family_add_member(individual, flat_family);
                    
                    // Add new family
                    int ret = 0;
                    khiter_t iter = kh_put(family, hash_families, family_key, &ret);
                    if (ret) {
                        kh_value(hash_families, iter) = flat_family;
                    } else {
                        free(flat_family); // Do not free members!
                        LOG_ERROR_F("Could not process family with pair of founders '%s'\n", family_key);
                    }
                }
            }
        }
    }
    
    // Retrieve the families in an array
    family_t **flat_families = malloc(kh_n_buckets(hash_families) * sizeof(family_t*));
    int curr_family = 0;
    for (int k = kh_begin(hash_families); k < kh_end(hash_families); k++) {
        if (!kh_exist(hash_families, k)) { continue; }
        
        flat_families[curr_family] = kh_value(hash_families, k);
        
/*
        family_t *family = flat_families[curr_family];
        assert(family->id);
        printf("Family %s\t\t", family->id);
        for (int k = kh_begin(family->members); k < kh_end(family->members); k++) {
            if (!kh_exist(family->members, k)) { continue; }
            individual_t *individual = kh_value(family->members, k);
            printf("%s\t", individual->id);
        }
        printf("\n");
*/
        
        curr_family++;
    }
    *num_families = curr_family;
    
    kh_destroy(family, hash_families);
    
    return flat_families;
}


//-----------------------------------------------------
// phenotypes management
//-----------------------------------------------------

khash_t(str)* get_phenotypes(ped_file_t *ped_file){
    assert(ped_file);
    return ped_file->variables;
}

int get_num_variables(ped_file_t* ped_file) {
    assert(ped_file);
    return ped_file->num_variables;
}

void set_unaffected_phenotype(const char* id, ped_file_t *ped_file){
    assert(ped_file);
    if(ped_file->unaffected)
        free(ped_file->unaffected);
    ped_file->unaffected = strdup(id);
}

void set_affected_phenotype(const char* id, ped_file_t *ped_file){
    assert(ped_file);
    if(ped_file->affected)
        free(ped_file->affected);
    ped_file->affected = strdup(id);
}

void set_variable_field(const char* id, int num_field, ped_file_t *ped_file){
    assert(ped_file);
    if(ped_file->variable_field){
        free(ped_file->variable_field);
    }
    ped_file->variable_field = strdup(id);
    ped_file->num_field = num_field;
}


int set_phenotype_group(char** ids, int n , ped_file_t *ped_file){
    int ret, fail = 0;
    int k;
    
    for(int i = 0; i < n; i++){
        //printf("Element %d from group %d : %s\n", i, ped_file->num_phenotypes, ids[i]);
        k = kh_put(str, ped_file->variables, ids[i], &ret);
        if(!ret) {        //Phenotype already inserted. 
            fail = 1;
        }
        kh_value(ped_file->variables, k) = ped_file->num_variables; //Overwritten value 
    }
    ped_file->num_variables++;
    if(fail){
        return -1;
    } else {
        return ped_file->num_variables;
    }
}
