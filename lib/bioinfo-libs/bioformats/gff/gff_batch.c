#include "gff_batch.h"
#include "gff_file.h"

gff_batch_t* gff_batch_new(size_t size) {
    gff_batch_t *gff_batch = malloc(sizeof(gff_batch_t));
    gff_batch->text = NULL;
    
    if (size < 1) {
        size = 100;
    }
    gff_batch->records = array_list_new(size, 1.4, COLLECTION_MODE_ASYNCHRONIZED);
    
    return gff_batch;
}

void gff_batch_free(gff_batch_t* batch) {
    assert(batch);
    
    if (batch->text) { free(batch->text); }
    array_list_free(batch->records, gff_record_free);
    free(batch);
}

void add_record_to_gff_batch(gff_record_t *record, gff_batch_t *batch) {
    assert(record);
    assert(batch);
    array_list_insert(record, batch->records);
}

inline int gff_batch_is_empty(gff_batch_t *batch) {
    assert(batch);
    return batch->records->size == 0;
}

inline int gff_batch_is_full(gff_batch_t *batch) {
    assert(batch);
    return batch->records->size == batch->records->capacity;
}

void gff_batch_print(FILE *fd, gff_batch_t *batch) {
    gff_record_t *first_record = (gff_record_t*) batch->records->items[0];
    fprintf(fd, "Batch with %zu/%zu records - %s in %ld\n", 
            batch->records->size, batch->records->capacity, 
            first_record->sequence, first_record->start);
}
