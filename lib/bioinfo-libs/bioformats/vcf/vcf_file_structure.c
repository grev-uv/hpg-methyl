#include "vcf_file_structure.h"


/* ********************************************************
 *      (De)Allocation of header entries and records      *
 * ********************************************************/

vcf_header_entry_t* vcf_header_entry_new() {
    vcf_header_entry_t *entry = (vcf_header_entry_t*) malloc (sizeof(vcf_header_entry_t));
    entry->name = NULL;
    entry->name_len = 0;
    entry->values = array_list_new(4, 1.5, COLLECTION_MODE_ASYNCHRONIZED);
    return entry;
}

void vcf_header_entry_free(vcf_header_entry_t *header_entry) {
    assert(header_entry);
    free(header_entry->name);
    array_list_free(header_entry->values, free);
    free(header_entry);
}

vcf_record_t* vcf_record_new() {
    vcf_record_t *record = (vcf_record_t*) calloc (1, sizeof(vcf_record_t));
    record->samples = array_list_new(16, 1.5, COLLECTION_MODE_ASYNCHRONIZED);
    return record;
}

vcf_record_t *vcf_record_copy(vcf_record_t *orig) {
    vcf_record_t *record = (vcf_record_t*) calloc (1, sizeof(vcf_record_t));
    record->chromosome = strndup(orig->chromosome, orig->chromosome_len);
    record->chromosome_len = orig->chromosome_len;
    record->position = orig->position;
    record->id = strndup(orig->id, orig->id_len);
    record->id_len = orig->id_len;
    record->reference = strndup(orig->reference, orig->reference_len);
    record->reference_len = orig->reference_len;
    record->alternate = strndup(orig->alternate, orig->alternate_len);
    record->alternate_len = orig->alternate_len;
    record->quality = orig->quality;
    record->filter = strndup(orig->filter, orig->filter_len);
    record->filter_len = orig->filter_len;
    record->info = strndup(orig->info, orig->info_len);
    record->info_len = orig->info_len;
    record->format = strndup(orig->format, orig->format_len);
    record->format_len = orig->format_len;
    record->samples = array_list_new(orig->samples->size + 1, 1.5, COLLECTION_MODE_ASYNCHRONIZED);
    for (int i = 0; i < orig->samples->size; i++) {
        array_list_insert(strdup(array_list_get(i, orig->samples)), record->samples);
    }
    return record;
}

void vcf_record_free(vcf_record_t *record) {
    assert(record);
    array_list_free(record->samples, free);
    free(record);
}

void vcf_record_free_deep(vcf_record_t *record) {
    assert(record);
    free(record->chromosome);
    free(record->id);
    free(record->reference);
    free(record->alternate);
    free(record->filter);
    free(record->info);
    free(record->format);
    array_list_free(record->samples, free);
    free(record);
}


/* ********************************************************
 *           Addition of header and record entries        *
 * ********************************************************/

int add_vcf_header_entry(vcf_header_entry_t *header_entry, vcf_file_t *file) {
    assert(header_entry);
    assert(file);
    return array_list_insert(header_entry, file->header_entries);
}

int add_vcf_sample_name(char *name, int length, vcf_file_t *file) {
    assert(name);
    assert(file);
    return array_list_insert(strndup(name, length), file->samples_names);
}

int add_text_batch(char *batch, vcf_file_t *file) {
    assert(batch);
    assert(file);
    list_item_t *item = list_item_new(rand(), 1, batch); 
    return list_insert_item(item, file->text_batches);
}

int add_vcf_batch(vcf_batch_t *batch, vcf_file_t *file) {
    assert(batch);
    assert(file);
    list_item_t *item = list_item_new(rand(), 1, batch); 
    return list_insert_item(item, file->record_batches);
}

char *fetch_vcf_text_batch(vcf_file_t *file) {
    assert(file);
    list_item_t *item = list_remove_item(file->text_batches);
    if (item) {
        char *batch = item->data_p;
        list_item_free(item);
        return batch;
    }
    return NULL;
}

vcf_batch_t *fetch_vcf_batch(vcf_file_t *file) {
    assert(file);
    list_item_t *item = list_remove_item(file->record_batches);
    if (item) {
        vcf_batch_t *batch = item->data_p;
        list_item_free(item);
        return batch;
    }
    return NULL;
}

vcf_batch_t *fetch_vcf_batch_non_blocking(vcf_file_t *file) {
    assert(file);
    list_item_t *item = list_remove_item_async(file->record_batches);
    if (item) {
        vcf_batch_t *batch = item->data_p;
        list_item_free(item);
        return batch;
    }
    return NULL;
}


size_t get_num_vcf_header_entries(vcf_file_t *file) {
    assert(file);
    return file->header_entries->size;
}

size_t get_num_values_in_vcf_header_entry(vcf_header_entry_t *entry) {
    assert(entry);
    return entry->values->size;
}

size_t get_num_vcf_samples(vcf_file_t *file) {
    assert(file);
    return file->samples_names->size;
}

size_t get_num_vcf_records(vcf_file_t *file) {
    assert(file);
    size_t num_records = 0;
    for (list_item_t *i = file->record_batches->first_p; i != NULL; i = i->next_p) {
        num_records += file->record_batches->length;
    }
    return num_records;
}

size_t get_num_vcf_batches(vcf_file_t *file) {
    assert(file);
    return file->record_batches->length;
}


/* ********************************************************
 *                    Batch management                    *
 * ********************************************************/

vcf_batch_t* vcf_batch_new(size_t size) {
    vcf_batch_t *vcf_batch = malloc(sizeof(vcf_batch_t));
    vcf_batch->text = NULL;
    
    if (size < 1) {
        size = 100;
    }
    vcf_batch->records = array_list_new(size, 1.2, COLLECTION_MODE_ASYNCHRONIZED);
    
    return vcf_batch;
}

void vcf_batch_free(vcf_batch_t* batch) {
    assert(batch);
    
    if (batch->text && !mmap_vcf) {
//         printf("text to free = '%.*s'\n", 50, batch->text);
        free(batch->text);
    }
    array_list_free(batch->records, vcf_record_free);
    free(batch);
}

void add_record_to_vcf_batch(vcf_record_t *record, vcf_batch_t *batch) {
    assert(record);
    assert(batch);
    array_list_insert(record, batch->records);
}

inline int vcf_batch_is_empty(vcf_batch_t *batch) {
    assert(batch);
    return batch->records->size == 0;
}

inline int vcf_batch_is_full(vcf_batch_t *batch) {
    assert(batch);
//     printf("batch size = %zu\tcapacity = %zu\n", vcf_batch->size, vcf_batch->capacity);
    return batch->records->size == batch->records->capacity;
}

int vcf_batch_print(FILE *fd, vcf_batch_t *batch) {
    assert(fd);
    assert(batch);
    
    vcf_record_t *first_record = array_list_get(0, batch->records);
    if (first_record) {
        return fprintf(fd, "Batch with %zu/%zu records: begin in chr%s:%ld\n", 
                      batch->records->size, batch->records->capacity, 
                      first_record->chromosome, first_record->position);
    } else {
        return fprintf(fd, "Batch with %zu/%zu records (empty)\n", 
                      batch->records->size, batch->records->capacity); 
    }
}


/* ********************************************************
 *                    Header management                   *
 * ********************************************************/

void set_vcf_file_format(char *fileformat, int length, vcf_file_t *file) {
    assert(fileformat);
    assert(file);
    file->format = strndup(fileformat, length);
    file->format_len = length;
//     LOG_DEBUG_F("set format = %.*s\n", file->format_len, file->format);
}

void set_vcf_header_entry_name(char *name, int length, vcf_header_entry_t *entry) {
    assert(name);
    assert(entry);
    entry->name = strndup(name, length);
    entry->name_len = length;
//     LOG_DEBUG_F("set entry name: %.*s\n", entry->name_len, entry->name);
}

void add_vcf_header_entry_value(char *value, int length, vcf_header_entry_t *entry) {
    assert(value);
    assert(entry);
    int result = array_list_insert(strndup(value, length), entry->values);
}


/* ********************************************************
 *                    Record management                   *
 * ********************************************************/

void set_vcf_record_chromosome(char* chromosome, int length, vcf_record_t* record) {
    assert(chromosome);
    assert(record);
    
    if (starts_with_n(chromosome, "chrom", 5)) {
        record->chromosome = chromosome + 5;
        record->chromosome_len = length - 5;
    } else if (starts_with_n(chromosome, "chr", 3)) {
        record->chromosome = chromosome + 3;
        record->chromosome_len = length - 3;
    } else {
        record->chromosome = chromosome;
        record->chromosome_len = length;
    }
//     LOG_DEBUG_F("set chromosome: '%.*s'\n", record->chromosome_len, record->chromosome);
}

void set_vcf_record_position(long position, vcf_record_t* record) {
    assert(record);
    record->position = position;
//     LOG_DEBUG_F("set position: %ld\n", record->position);
}

void set_vcf_record_id(char* id, int length, vcf_record_t* record) {
    assert(id);
    assert(record);
    record->id = id;
    record->id_len = length;
//     LOG_DEBUG_F("set id: %s\n", record->id);
}

void set_vcf_record_reference(char* reference, int length, vcf_record_t* record) {
    assert(reference);
    assert(record);
    record->reference = reference;
    record->reference_len = length;
//     LOG_DEBUG_F("set reference: %s\n", record->reference);
}

void set_vcf_record_alternate(char* alternate, int length, vcf_record_t* record) {
    assert(alternate);
    assert(record);
    record->alternate = alternate;
    record->alternate_len = length;
//     LOG_DEBUG_F("set alternate: %s\n", record->alternate);
}

void set_vcf_record_quality(float quality, vcf_record_t* record) {
    assert(record);
    record->quality = quality;
//     LOG_DEBUG_F("set quality: %f\n", record->quality);
}

void set_vcf_record_filter(char* filter, int length, vcf_record_t* record) {
    assert(filter);
    assert(record);
    record->filter = filter;
    record->filter_len = length;
//     LOG_DEBUG_F("set filter: %s\n", record->filter);
}

void set_vcf_record_info(char* info, int length, vcf_record_t* record) {
    assert(info);
    assert(record);
    record->info = info;
    record->info_len = length;
//     LOG_DEBUG_F("set info: %s\n", record->info);
}

void set_vcf_record_format(char* format, int length, vcf_record_t* record) {
    assert(format);
    assert(record);
    record->format = format;
    record->format_len = length;
//     LOG_DEBUG_F("set format: %s\n", record->format);
}

void add_vcf_record_sample(char* sample, int length, vcf_record_t* record) {
    assert(sample);
    assert(record);
//     int result = array_list_insert(sample, record->samples);
    int result = array_list_insert(strndup(sample, length), record->samples);
//     if (result) {
//         LOG_DEBUG_F("sample %s inserted\n", sample);
//     } else {
//         LOG_DEBUG_F("sample %s not inserted\n", sample);
//     }
}


/* *******************
 *      Sorting      *
 * *******************/

individual_t **sort_individuals(vcf_file_t *vcf, ped_file_t *ped) {
    family_t *family;
    family_t **families = (family_t**) cp_hashtable_get_values(ped->families);
    int num_families = get_num_families(ped);

    individual_t **individuals = calloc (get_num_vcf_samples(vcf), sizeof(individual_t*));
    khash_t(ids) *positions = associate_samples_and_positions(vcf);
    int pos = 0;

    for (int f = 0; f < num_families; f++) {
        family = families[f];
        individual_t *father = family->father;
        individual_t *mother = family->mother;

        if (father != NULL) {
            pos = 0;
            LOG_DEBUG_F("father ID = %s\n", father->id);
            khiter_t iter = kh_get(ids, positions, father->id);
            if (iter != kh_end(positions)) {
                pos = kh_value(positions, iter);
                individuals[pos] = father;
            }
        }

        if (mother != NULL) {
            pos = 0;
            LOG_DEBUG_F("mother ID = %s\n", mother->id);
            khiter_t iter = kh_get(ids, positions, mother->id);
            if (iter != kh_end(positions)) {
                pos = kh_value(positions, iter);
                individuals[pos] = mother;
            }
        }

        linked_list_iterator_t *iterator = linked_list_iterator_new(family->children);
        individual_t *child = NULL;
        while (child = linked_list_iterator_curr(iterator)) {
            pos = 0;
            LOG_DEBUG_F("child ID = %s\n", child->id);
            khiter_t iter = kh_get(ids, positions, child->id);
            if (iter != kh_end(positions)) {
                pos = kh_value(positions, iter);
                individuals[pos] = child;
            }
            linked_list_iterator_next(iterator);
        }
        linked_list_iterator_free(iterator);
        
        iterator = linked_list_iterator_new(family->unknown);
        individual_t *unknown = NULL;
        while (unknown = linked_list_iterator_curr(iterator)) {
            pos = 0;
            LOG_DEBUG_F("unknown ID = %s\n", unknown->id);
            khiter_t iter = kh_get(ids, positions, unknown->id);
            if (iter != kh_end(positions)) {
                pos = kh_value(positions, iter);
                individuals[pos] = unknown;
            }
            linked_list_iterator_next(iterator);
        }
        linked_list_iterator_free(iterator);
        
        assert(father || mother || linked_list_size(family->unknown) > 0);
    }

    kh_destroy(ids, positions);
    free(families);

    return individuals;
}

khash_t(ids)* associate_samples_and_positions(vcf_file_t* file) {
    LOG_DEBUG_F("** %zu sample names read\n", file->samples_names->size);
    array_list_t *sample_names = file->samples_names;
    khash_t(ids) *sample_ids = kh_init(ids);
    
    for (int i = 0; i < sample_names->size; i++) {
        char *name = sample_names->items[i];
        int ret;
        khiter_t iter = kh_get(ids, sample_ids, name);
        if (iter != kh_end(sample_ids)) {
            LOG_FATAL_F("Sample %s appears more than once. File can not be analyzed.\n", name);
        } else {
            iter = kh_put(ids, sample_ids, name, &ret);
            if (ret) {
                kh_value(sample_ids, iter) = i;
            }
        }
    }
    
//     char **keys = (char**) cp_hashtable_get_keys(sample_names);
//     int num_keys = cp_hashtable_count(sample_names);
//     for (int i = 0; i < num_keys; i++) {
//         printf("%s\t%d\n", keys[i], *((int*) cp_hashtable_get(sample_ids, keys[i])));
//     }
    
    return sample_ids;
}
