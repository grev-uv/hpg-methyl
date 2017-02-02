#include "bed_file.h"


//====================================================================================
//  bed_file.c
//  bed file management functions
//====================================================================================

//-----------------------------------------------------
// bed_open
//-----------------------------------------------------


bed_file_t *bed_open(char *filename) {
    size_t len;
    char *data = mmap_file(&len, filename);

    bed_file_t *bed_file = (bed_file_t *) malloc(sizeof(bed_file_t));
    bed_file->filename = filename;
    bed_file->data = data;
    bed_file->data_len = len;
    bed_file->header_entries = linked_list_new(COLLECTION_MODE_SYNCHRONIZED);
    bed_file->records = linked_list_new(COLLECTION_MODE_SYNCHRONIZED);
    return bed_file;
}


//-----------------------------------------------------
// bed_close and memory freeing
//-----------------------------------------------------

void bed_close(bed_file_t *bed_file, int free_records) {
    // Free header entries
    linked_list_free(bed_file->header_entries, NULL);   // TODO doesn't work! :(
    
    // Free records list if asked to
    if (free_records) {
        linked_list_free(bed_file->records, bed_record_free);
    }
    
    munmap((void*) bed_file->data, bed_file->data_len);
    free(bed_file);
}

void bed_header_entry_free(bed_header_entry_t *bed_header_entry) {
    assert(bed_header_entry);
    free(bed_header_entry->text);
    free(bed_header_entry);
}

void bed_record_free(bed_record_t *bed_record) {
    assert(bed_record);
    free(bed_record);
}

void bed_record_free_deep(bed_record_t *bed_record) {
    assert(bed_record);
    if (bed_record->sequence)   { free(bed_record->sequence); }
    if (bed_record->name)       { free(bed_record->name); }
    if (bed_record->item_rgb)   { free(bed_record->item_rgb); }
    if (bed_record->block_sizes) { free(bed_record->block_sizes); }
    if (bed_record->block_starts) { free(bed_record->block_starts); }
    free(bed_record);
}


//-----------------------------------------------------
// I/O operations (read and write) in various ways
//-----------------------------------------------------

int bed_read(bed_file_t *bed_file) {
    return bed_ragel_read(NULL, 1, bed_file);
}

int bed_read_batches(list_t *batches_list, size_t batch_size, bed_file_t *bed_file) {
    return bed_ragel_read(batches_list, batch_size, bed_file);
}

int bed_write(bed_file_t *bed_file, char *filename) {
    FILE *fd = fopen(filename, "w");
    if (fd < 0) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return 1;
    }
    
    if (bed_write_to_file(bed_file, fd)) {
        fprintf(stderr, "Error writing file: %s\n", filename);
        fclose(fd);
        return 2;
    }
    
    fclose(fd);
    return 0;
}

//-----------------------------------------------------
// load data into the bed_file_t
//-----------------------------------------------------

int add_bed_header_entry(bed_header_entry_t *header_entry, bed_file_t *bed_file) {
    if (linked_list_insert(header_entry, bed_file->header_entries)) {
        LOG_DEBUG_F("header entry %zu\n", bed_file->header_entries->size);
        return 1;
    } else {
        LOG_WARN_F("header entry %zu not inserted\n", bed_file->header_entries->size);
        return 0;
    }
}

int add_bed_record(bed_record_t* record, bed_file_t *bed_file) {
    if (linked_list_insert(record, bed_file->records)) {
        LOG_DEBUG_F("record %zu\n", bed_file->records->size);
        return 1;
    } else {
        LOG_DEBUG_F("record %zu not inserted\n", bed_file->records->size);
        return 0;
    }
}

