#include "ped_write.h"

int ped_write_to_file(ped_file_t *ped_file, FILE *fd) {
    assert(ped_file);
    assert(fd);
    
    family_t **families = (family_t**) cp_hashtable_get_values(ped_file->families);
    int num_families = get_num_families(ped_file);
    
    // Write members of each family
    family_t *family;
    LOG_DEBUG_F("Number of families read: %d\n", num_families);
    for (int i = 0; i < num_families; i++) {
        family = families[i];
        
        for (int k = kh_begin(family->members); k < kh_end(family->members); k++) {
            if (!kh_exist(family->members, k)) { continue; }
            individual_t *individual = kh_value(family->members, k);
            write_ped_individual(individual, fd);
        }
    }
    
    return 0;
}


void write_ped_individual(individual_t* individual, FILE* fd) {
    assert(individual);
    assert(fd);
    
    char *cond_str = NULL;
    switch(individual->condition) {
        case UNAFFECTED:
            cond_str = "1";
            break;
        case AFFECTED:
            cond_str = "2";
            break;
        default:
            cond_str = "0";
            break;
    }
    
    char *sex_str = NULL;
    switch(individual->sex) {
        case MALE:
            sex_str = "1";
            break;
        case FEMALE:
            sex_str = "2";
            break;
        default:
            sex_str = "0";
            break;
    }
    
    
    fprintf(fd, "%s\t%s\t%s\t%s\t%s\t%s\n", individual->family->id, individual->id, 
            individual->father_id, individual->mother_id, sex_str, cond_str);

}



void write_ped_batch(ped_batch_t *ped_batch, FILE *fd) {
    assert(ped_batch);
    assert(fd);
    
    for (list_item_t *i = ped_batch->first_p; i != NULL; i = i->next_p) {
        write_ped_record(i->data_p, fd);
    }
}

void write_ped_record(ped_record_t* ped_record, FILE *fd) {
    assert(ped_record);
    assert(fd);
    
    fprintf(fd, "%s\t%s\t%s\t%s\t%d\t%s", ped_record->family_id, ped_record->individual_id, 
            ped_record->father_id, ped_record->mother_id, ped_record->sex,
            ped_record->phenotype); //Doesn't print the custom field value
}
