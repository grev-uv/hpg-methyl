#include "ped_batch.h"
#include "ped_file.h"

int freed = 0;
ped_batch_t* ped_batch_new(size_t size) {
    ped_batch_t *ped_batch_p = (ped_batch_t*) malloc (sizeof(ped_batch_t));
    list_init("batch", 1, size, ped_batch_p);
    return ped_batch_p;
}

void ped_batch_free(ped_batch_t* ped_batch_p) {
    list_item_t *item;
    while ((item = list_remove_item_async(ped_batch_p))) {
        ped_record_free(item->data_p);
        list_item_free(item);
    }
    free(ped_batch_p);
}

void add_record_to_ped_batch(ped_record_t *record, ped_batch_t *ped_batch) {
    list_item_t *item = list_item_new(ped_batch->length, 1, record);
    list_insert_item(item, ped_batch);
}

inline int ped_batch_is_empty(ped_batch_t *ped_batch) {
    return ped_batch->length == 0;
}

inline int ped_batch_is_full(ped_batch_t *ped_batch) {
    return ped_batch->length == ped_batch->max_length;
}

void ped_batch_print(FILE *fd, ped_batch_t *ped_batch) {
    fprintf(fd, "Batch with %zu/%zu records\n", ped_batch->length, ped_batch->max_length);
}
