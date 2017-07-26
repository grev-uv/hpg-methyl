#include "bed_write.h"

int bed_write_to_file(bed_file_t *bed_file, FILE *fd) {
    assert(bed_file);
    assert(fd);

    // Write header entries
    linked_list_iterator_t *headers_iter = linked_list_iterator_new(bed_file->header_entries);
    bed_header_entry_t *entry = NULL;
    while ((entry = linked_list_iterator_curr(headers_iter))) {
        write_bed_header_entry(entry, fd);
        linked_list_iterator_next(headers_iter);
    }
    linked_list_iterator_free(headers_iter);
    
    // Write records
    linked_list_iterator_t *records_iter = linked_list_iterator_new(bed_file->records);
    bed_record_t *record = NULL;
    while ((entry = linked_list_iterator_curr(records_iter))) {
        write_bed_record(record, fd);
        linked_list_iterator_next(records_iter);
    }
    linked_list_iterator_free(records_iter);

    return 0;
}

void write_bed_header_entry(bed_header_entry_t *entry, FILE *fd) {
    assert(entry);
    assert(fd);

    fprintf(fd, "##%s\n", entry->text);
}

void write_bed_batch(bed_batch_t *bed_batch, FILE *fd) {
    assert(bed_batch);
    assert(fd);
    
    for (int i = 0; i < bed_batch->records->size; i++) {
        write_bed_record(bed_batch->records->items[i], fd);
    }
}

void write_bed_record(bed_record_t* bed_record, FILE *fd) {
    assert(bed_record);
    assert(fd);
    
    fprintf(fd, "%.*s\t%ld\t%ld", bed_record->sequence_len, bed_record->sequence, bed_record->start, bed_record->end);
    if (bed_record->name) {
        fprintf(fd, "\t%.*s", bed_record->name_len, bed_record->name);
    } else {
        fprintf(fd, "\t.");
    }
    if (bed_record->score >= 0) {
        fprintf(fd, "\t%.3f", bed_record->score);
    } else {
        fprintf(fd, "\t.");
    }
    if (bed_record->strand) {
        fprintf(fd, "\t%c", bed_record->strand);
    } else {
        fprintf(fd, "\t.");
    }
    if (bed_record->thickstart >= 0) {
        fprintf(fd, "\t%ld", bed_record->thickstart);
    } else {
        fprintf(fd, "\t.");
    }
    if (bed_record->thickend >= 0) {
        fprintf(fd, "\t%ld", bed_record->thickend);
    } else {
        fprintf(fd, "\t.");
    }
    if (bed_record->item_rgb) {
        fprintf(fd, "\t%.*s", bed_record->item_rgb_len, bed_record->item_rgb);
    } else {
        fprintf(fd, "\t.");
    }
    if (bed_record->block_count >= 0) {
        fprintf(fd, "\t%d", bed_record->block_count);
    } else {
        fprintf(fd, "\t.");
    }
    if (bed_record->block_sizes) {
        fprintf(fd, "\t%.*s", bed_record->block_sizes_len, bed_record->block_sizes);
    } else {
        fprintf(fd, "\t.");
    }
    if (bed_record->block_starts) {
        fprintf(fd, "\t%.*s", bed_record->block_starts_len, bed_record->block_starts);
    } else {
        fprintf(fd, "\t.");
    }
}
