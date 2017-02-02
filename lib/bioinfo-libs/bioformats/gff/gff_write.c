#include "gff_write.h"

int gff_write_to_file(gff_file_t *gff_file, FILE *fd) {
    assert(gff_file);
    assert(fd);

    // Write header entries
    linked_list_iterator_t *headers_iter = linked_list_iterator_new(gff_file->header_entries);
    gff_header_entry_t *entry = NULL;
    while ((entry = linked_list_iterator_curr(headers_iter))) {
        write_gff_header_entry(entry, fd);
        linked_list_iterator_next(headers_iter);
    }
    linked_list_iterator_free(headers_iter);
    
    // Write records
    linked_list_iterator_t *records_iter = linked_list_iterator_new(gff_file->records);
    gff_record_t *record = NULL;
    while ((entry = linked_list_iterator_curr(records_iter))) {
        write_gff_record(record, fd);
        linked_list_iterator_next(records_iter);
    }
    linked_list_iterator_free(records_iter);

    return 0;
}

void write_gff_header_entry(gff_header_entry_t *entry, FILE *fd) {
    assert(entry);
    assert(fd);

    fprintf(fd, "##%s\n", entry->text);
}

void write_gff_batch(gff_batch_t *gff_batch, FILE *fd) {
    assert(gff_batch);
    assert(fd);
    
    for (int i = 0; i < gff_batch->records->size; i++) {
        write_gff_record(gff_batch->records->items[i], fd);
    }
}

void write_gff_record(gff_record_t* gff_record, FILE *fd) {
    assert(gff_record);
    assert(fd);
    
    fprintf(fd, "%s\t%s\t%s\t%ld\t%ld\t", gff_record->sequence, gff_record->source, gff_record->feature, gff_record->start, gff_record->end);
    if (gff_record->score < 0) {
        fprintf(fd, ".\t");
    } else {
        fprintf(fd, "%d\t", gff_record->score);
    }
    fprintf(fd, "%c\t", gff_record->strand);
    if (gff_record->frame < 0) {
        fprintf(fd, ".\t");
    } else {
        fprintf(fd, "%d\t", gff_record->frame);
    }
    fprintf(fd, "%s\t\n", gff_record->attribute);
}
