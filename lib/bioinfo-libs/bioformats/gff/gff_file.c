#include "gff_file.h"
#include "containers/linked_list.h"


//====================================================================================
//  gff_file.c
//  gff file management functions
//====================================================================================

//-----------------------------------------------------
// gff_open
//-----------------------------------------------------


gff_file_t *gff_open(char *filename) {
    size_t len;
    char *data = mmap_file(&len, filename);

    gff_file_t *gff_file = (gff_file_t *) malloc(sizeof(gff_file_t));
    gff_file->filename = filename;
    gff_file->data = data;
    gff_file->data_len = len;
    gff_file->header_entries = linked_list_new(COLLECTION_MODE_SYNCHRONIZED);
    gff_file->records = linked_list_new(COLLECTION_MODE_SYNCHRONIZED);
    return gff_file;
}


//-----------------------------------------------------
// gff_close and memory freeing
//-----------------------------------------------------

void gff_close(gff_file_t *gff_file, int free_records) {
    // Free header entries
    linked_list_free(gff_file->header_entries, NULL);   // TODO doesn't work! :(
    
    // Free records list if asked to
    if (free_records) {
        linked_list_free(gff_file->records, gff_record_free);
    }
    
    munmap((void*) gff_file->data, gff_file->data_len);
    free(gff_file);
}

void gff_header_entry_free(gff_header_entry_t *gff_header_entry) {
    assert(gff_header_entry);
    free(gff_header_entry->text);
    free(gff_header_entry);
}

void gff_record_free(gff_record_t *gff_record) {
    assert(gff_record);
    free(gff_record);
}

void gff_record_free_deep(gff_record_t *gff_record) {
    assert(gff_record);
    if (gff_record->sequence) { free(gff_record->sequence); }
    if (gff_record->source) { free(gff_record->source); }
    if (gff_record->feature) { free(gff_record->feature); }
    if (gff_record->attribute) { free(gff_record->attribute); }
    free(gff_record);
}


//-----------------------------------------------------
// I/O operations (read and write) in various ways
//-----------------------------------------------------

int gff_read(gff_file_t *gff_file) {
    return gff_ragel_read(NULL, 1, gff_file);
}

int gff_read_batches(list_t *batches_list, size_t batch_size, gff_file_t *gff_file) {
    return gff_ragel_read(batches_list, batch_size, gff_file);
}

int gff_write(gff_file_t *gff_file, char *filename) {
    FILE *fd = fopen(filename, "w");
    if (fd < 0) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return 1;
    }
    
    if (gff_write_to_file(gff_file, fd)) {
        fprintf(stderr, "Error writing file: %s\n", filename);
        fclose(fd);
        return 2;
    }
    
    fclose(fd);
    return 0;
}

//-----------------------------------------------------
// load data into the gff_file_t
//-----------------------------------------------------

int add_gff_header_entry(gff_header_entry_t *header_entry, gff_file_t *gff_file) {
    if (linked_list_insert(header_entry, gff_file->header_entries)) {
        LOG_DEBUG_F("header entry %zu\n", gff_file->header_entries->size);
        return 1;
    } else {
        LOG_WARN_F("header entry %zu not inserted\n", gff_file->header_entries->size);
        return 0;
    }
}

int add_gff_record(gff_record_t* record, gff_file_t *gff_file) {
    if (linked_list_insert(record, gff_file->records)) {
        LOG_DEBUG_F("record %zu\n", gff_file->records->size);
        return 1;
    } else {
        LOG_DEBUG_F("record %zu not inserted\n", gff_file->records->size);
        return 0;
    }
}

