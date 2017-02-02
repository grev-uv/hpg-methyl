#include "bed_batch.h"
#include "bed_file.h"

bed_batch_t* bed_batch_new(size_t size) {
    bed_batch_t *bed_batch = malloc(sizeof(bed_batch_t));
    bed_batch->text = NULL;
    
    if (size < 1) {
        size = 100;
    }
    bed_batch->records = array_list_new(size, 1.4, COLLECTION_MODE_ASYNCHRONIZED);
    
    return bed_batch;
}

void bed_batch_free(bed_batch_t* batch) {
    assert(batch);
    
    if (batch->text) { free(batch->text); }
    array_list_free(batch->records, bed_record_free);
    free(batch);
}

void add_record_to_bed_batch(bed_record_t *record, bed_batch_t *batch) {
    assert(record);
    assert(batch);
    array_list_insert(record, batch->records);
}

inline int bed_batch_is_empty(bed_batch_t *batch) {
    assert(batch);
    return batch->records->size == 0;
}

inline int bed_batch_is_full(bed_batch_t *batch) {
    assert(batch);
    return batch->records->size == batch->records->capacity;
}

void bed_batch_print(FILE *fd, bed_batch_t *batch) {
    bed_record_t *first_record = (bed_record_t*) batch->records->items[0];
    fprintf(fd, "Batch with %zu/%zu records - %s in %ld\n", 
            batch->records->size, batch->records->capacity, 
            first_record->sequence, first_record->start);
}
