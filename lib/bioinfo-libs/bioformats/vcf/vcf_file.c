#include "vcf_file.h"


//====================================================================================
//  vcf_file.c
//  vcf file management functions
//====================================================================================

//-----------------------------------------------------
// vcf_open
//-----------------------------------------------------


vcf_file_t *vcf_open(char *filename, size_t max_simultaneous_batches) {
    if (!exists(filename)) {
        return NULL;
    }
    
    vcf_file_t *vcf_file = (vcf_file_t *) malloc(sizeof(vcf_file_t));
    vcf_file->filename = filename;
    
    // Initialize file descriptor or mmap'd buffers
    if (mmap_vcf) {
        size_t len;
        char *data = mmap_file(&len, filename);
        vcf_file->fd = fopen(filename, "r");
        vcf_file->data = data;
        vcf_file->data_len = len;
    } else {
        vcf_file->fd = fopen(filename, "r");
        vcf_file->data = NULL;
        vcf_file->data_len = 0;
    }

    // Initialize header
    vcf_file->header_entries = array_list_new(10, 1.5, COLLECTION_MODE_SYNCHRONIZED);

    // Initialize samples names list
    vcf_file->samples_names = array_list_new(10, 1.5, COLLECTION_MODE_SYNCHRONIZED);

    // Initialize text batches
    vcf_file->text_batches = (list_t*) malloc(sizeof(list_t));
    if (max_simultaneous_batches <= 0) {
        list_init("text batches", 1, INT_MAX, vcf_file->text_batches);
    } else {
        list_init("text batches", 1, max_simultaneous_batches, vcf_file->text_batches);
    }
    
    // Initialize record batches
    vcf_file->record_batches = (list_t*) malloc(sizeof(list_t));
    if (max_simultaneous_batches <= 0) {
        list_init("vcf_file batches", 1, INT_MAX, vcf_file->record_batches);
    } else {
        list_init("vcf_file batches", 1, max_simultaneous_batches, vcf_file->record_batches);
    }

    return vcf_file;
}

vcf_file_t *vcf_file_new(char *filename, size_t max_simultaneous_batches) {
    vcf_file_t *vcf_file = (vcf_file_t *) malloc(sizeof(vcf_file_t));
    vcf_file->filename = filename;
    
    vcf_file->data = NULL;
    vcf_file->data_len = 0;

    // Initialize header
    vcf_file->header_entries = array_list_new(10, 1.5, COLLECTION_MODE_SYNCHRONIZED);

    // Initialize samples names list
    vcf_file->samples_names = array_list_new(10, 1.5, COLLECTION_MODE_SYNCHRONIZED);

    // Initialize text batches
    vcf_file->text_batches = (list_t*) malloc(sizeof(list_t));
    if (max_simultaneous_batches <= 0) {
        list_init("text batches", 1, INT_MAX, vcf_file->text_batches);
    } else {
        list_init("text batches", 1, max_simultaneous_batches, vcf_file->text_batches);
    }
    
    // Initialize record batches
    vcf_file->record_batches = (list_t*) malloc(sizeof(list_t));
    if (max_simultaneous_batches <= 0) {
        list_init("vcf_file batches", 1, INT_MAX, vcf_file->record_batches);
    } else {
        list_init("vcf_file batches", 1, max_simultaneous_batches, vcf_file->record_batches);
    }

    return vcf_file;
}


//-----------------------------------------------------
// vcf_close and memory freeing
//-----------------------------------------------------

void vcf_close(vcf_file_t *vcf_file) {
    free(vcf_file->format);
    array_list_free(vcf_file->samples_names, free);
    array_list_free(vcf_file->header_entries, vcf_header_entry_free);
    list_free_deep(vcf_file->text_batches, free);
    list_free_deep(vcf_file->record_batches, vcf_batch_free);

    if (mmap_vcf) {
        munmap((void*) vcf_file->data, vcf_file->data_len);
    }
    
    fclose(vcf_file->fd);
    free(vcf_file);
}


//-----------------------------------------------------
// I/O operations (read and write) in various ways
//-----------------------------------------------------

int vcf_read(vcf_file_t *file, bool also_parse, size_t batch_size, bool size_in_lines) {
    int ret_code = 0;
    if (also_parse) {
        if (size_in_lines) {
            ret_code = vcf_parse_batches(batch_size, file);
        } else {
            ret_code = vcf_parse_batches_in_bytes(batch_size, file);
        }
    } else {
        if (size_in_lines) {
            ret_code = vcf_read_batches(batch_size, file);
        } else {
            ret_code = vcf_read_batches_in_bytes(batch_size, file);
        }
    }
    return ret_code;
}

int vcf_parse_batches(size_t batch_lines, vcf_file_t *vcf_file) {
    if (ends_with(vcf_file->filename, ".vcf")) {
        return vcf_read_and_parse(batch_lines, vcf_file);
    } else if (ends_with(vcf_file->filename, ".gz")) {
        return vcf_gzip_read_and_parse(batch_lines, vcf_file);
    }
    LOG_FATAL_F("The format of file %s can't be processed\n", vcf_file->filename);
}

int vcf_parse_batches_in_bytes(size_t batch_bytes, vcf_file_t *vcf_file) {
    if (ends_with(vcf_file->filename, ".vcf")) {
        return vcf_read_and_parse_bytes(batch_bytes, vcf_file);
    } else if (ends_with(vcf_file->filename, ".gz")) {
        return vcf_gzip_read_and_parse_bytes(batch_bytes, vcf_file);
    }
    LOG_FATAL_F("The format of file %s can't be processed\n", vcf_file->filename);
}

int vcf_read_batches(size_t batch_lines, vcf_file_t *vcf_file) {
    if (ends_with(vcf_file->filename, ".vcf")) {
        return vcf_light_read(batch_lines, vcf_file);
    } else if (ends_with(vcf_file->filename, ".gz")) {
        return vcf_gzip_light_read(batch_lines, vcf_file);
    }
    LOG_FATAL_F("The format of file %s can't be processed\n", vcf_file->filename);
}

int vcf_read_batches_in_bytes(size_t batch_bytes, vcf_file_t *vcf_file) {
    if (ends_with(vcf_file->filename, ".vcf")) {
        return vcf_light_read_bytes(batch_bytes, vcf_file);
    } else if (ends_with(vcf_file->filename, ".gz")) {
        return vcf_gzip_light_read_bytes(batch_bytes, vcf_file);
    }
    LOG_FATAL_F("The format of file %s can't be processed\n", vcf_file->filename);
}

int vcf_multiread_batches(list_t **text_lists, size_t batch_lines, vcf_file_t **vcf_files, int num_files) {
    return vcf_light_multiread(text_lists, batch_lines, vcf_files, num_files);
}

void notify_end_reading(vcf_file_t *vcf_file) {
    list_decr_writers(vcf_file->text_batches);
}

void notify_end_parsing(vcf_file_t *vcf_file) {
    list_decr_writers(vcf_file->record_batches);
}




int vcf_write(vcf_file_t *vcf_file, char *filename) {
    FILE *fd = fopen(filename, "w");
    if (fd < 0) {
        LOG_FATAL_F("Error opening file: %s\n", filename);
    }

    if (write_vcf_file(vcf_file, fd)) {
        fclose(fd);
        LOG_FATAL_F("Error writing file: %s\n", filename);
    }

    fclose(fd);
    return 0;
}
