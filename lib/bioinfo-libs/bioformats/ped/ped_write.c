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
        // Write mother and father
        write_ped_individual(family->father, fd);
        write_ped_individual(family->mother, fd);
        // Write children
        LOG_DEBUG_F("Family %s has %ld children\n", family->id, family->children->size);
            
        linked_list_iterator_t *iterator = linked_list_iterator_new(family->children);
        individual_t *child = NULL;
        while (child = linked_list_iterator_curr(iterator)) {
            write_ped_individual(child, fd);
            linked_list_iterator_next(iterator);
        }
        linked_list_iterator_free(iterator);
    }
    
    return 0;
}


void write_ped_individual(individual_t* individual, FILE* fd) {
    assert(individual);
    assert(fd);
    
    fprintf(fd, "%s\t%s\t", individual->family->id, individual->id);
    fprintf(fd, "%s\t%s\t%d\t", (individual->father == NULL) ? "0" : individual->father->id,
                                (individual->mother == NULL) ? "0" : individual->mother->id,
                                individual->sex);
    if (individual->phenotype - ((int) individual->phenotype) > 1e3) {
        fprintf(fd, "%f\t\n", individual->phenotype);
    } else {
        fprintf(fd, "%d\t\n", (int) individual->phenotype);
    }
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
    
    fprintf(fd, "%s\t%s\t%s\t%s\t%d\t", ped_record->family_id, ped_record->individual_id, 
            ped_record->father_id, ped_record->mother_id, ped_record->sex);
    if (ped_record->phenotype - ((int) ped_record->phenotype) > 1e3) {
        fprintf(fd, "%f\t\n", ped_record->phenotype);
    } else {
        fprintf(fd, "%d\t\n", (int) ped_record->phenotype);
    }
}
