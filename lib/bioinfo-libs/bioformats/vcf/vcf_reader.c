#include "vcf_reader.h"

vcf_reader_status *vcf_reader_status_new(size_t batch_lines, size_t current_batch_id) {
    vcf_reader_status *status = (vcf_reader_status *) malloc (sizeof(vcf_reader_status));
    status->current_record = NULL;
    status->current_header_entry = vcf_header_entry_new();
    if (batch_lines == 0) {
        status->current_batch = vcf_batch_new(500);
    } else if (batch_lines > 0) {
        status->current_batch = vcf_batch_new(batch_lines);
    }

    status->num_samples = 0;
    status->num_records = 0;
    status->num_batches = current_batch_id;
    
    return status;
}

void vcf_reader_status_free(vcf_reader_status *status) {
    assert(status);
    
    if (status->current_header_entry->name == NULL && status->current_header_entry->values->size == 0) {
        vcf_header_entry_free(status->current_header_entry);
    }
    if (vcf_batch_is_empty(status->current_batch)) {
        vcf_batch_free(status->current_batch);
    }
    if (status->current_record->chromosome == NULL && status->current_record->samples->size == 0) {
        vcf_record_free(status->current_record);
    }
    free(status);
}


/* **********************************************
 *              Reading and parsing             *
 * **********************************************/

int vcf_read_and_parse(size_t batch_lines, vcf_file_t *file) {
    assert(file);
    assert(batch_lines > 0);
    
    int cs = 0;
    char *p, *pe;

    vcf_reader_status *status = vcf_reader_status_new(batch_lines, 0);
    
    if (mmap_vcf) {
        LOG_DEBUG("Using mmap for file loading\n");
        p = file->data;
        pe = p + file->data_len;
        cs = run_vcf_parser(p, pe, batch_lines, file, status);
    } else {
        LOG_DEBUG("Using file-IO functions for file loading\n");
        size_t max_len = 256;
        int eof_found = 0;
        char *aux;

        // Read text of a batch and call ragel parser in a loop
        while (!eof_found) {
            char *data = (char*) calloc (max_len, sizeof(char));
            int c = 0;
            int lines = 0;

            for (int i = 0; !eof_found && lines < batch_lines; i++) {
                c = fgetc(file->fd);
                
                if (c != EOF) {
                    max_len = consume_input(c, &data, max_len, i);
                    if (c == '\n') {
                        lines++;
                    }
                    (file->data_len)++;
                } else {
                    eof_found = 1;
                }

            }

            data[file->data_len] = '\0';

            p = data;
            pe = p + file->data_len;
            cs |= run_vcf_parser(p, pe, batch_lines, file, status);
            file->data_len = 0;

            // Prepare status for next batch
            status->current_batch = vcf_batch_new(batch_lines);
        }
    }

    // Insert the last batch
    // TODO this should not be neccessary because it is now inserted in execute_vcf_ragel_machine
//     if (!vcf_batch_is_empty(status->current_batch))
//     {
//         list_item_t *item = list_item_new(get_num_vcf_batches(file), 1, status->current_batch);
//         list_insert_item(item, batches_list);
//         printf("Batch added - %zu records (last)\n", status->current_batch->records->size);
//     }

    if ( cs ) {
        LOG_INFO("Last state was not the expected");
    } 

    LOG_INFO_F("Records read = %zu\n", status->num_records);
    LOG_INFO_F("Samples per record = %zu\n", get_num_vcf_samples(file));

    // Free status->current_xxx pointers if not needed in another module
    vcf_reader_status_free(status);

    return cs ;
}

int vcf_read_and_parse_bytes(size_t batch_bytes, vcf_file_t *file) {
    assert(file);
    assert(batch_bytes > 0);
    
    int cs = 0;
    char *p, *pe;

    vcf_reader_status *status = vcf_reader_status_new(0, 0);
    
    if (mmap_vcf) {
        LOG_DEBUG("Using mmap for file loading\n");
        p = file->data;
        pe = p + file->data_len;
        cs = run_vcf_parser(p, pe, 0, file, status);
    } else {
        LOG_DEBUG("Using file-IO functions for file loading\n");
        size_t max_len = 256;
        int eof_found = 0;
        char *aux;

        // Read text of a batch and call ragel parser in a loop
        while (!eof_found) {
            char *data = (char*) calloc (max_len, sizeof(char));
            int c = 0;
            int lines = 0;

            for (int i = 0; !eof_found; i++) {
                c = fgetc(file->fd);
                
                if (c != EOF) {
                    max_len = consume_input(c, &data, max_len, i);
                    if (c == '\n') {
                        lines++;
                        if (i >= batch_bytes) {
                            break;
                        }
                    }
                    (file->data_len)++;
                } else {
                    eof_found = 1;
                }

            }

            data[file->data_len] = '\0';

            p = data;
            pe = p + file->data_len;
            cs |= run_vcf_parser(p, pe, 0, file, status);
            file->data_len = 0;

            // Prepare status for next batch
            status->current_batch = vcf_batch_new(500);
        }
    }

    // Insert the last batch
    // TODO this should not be neccessary because it is now inserted in execute_vcf_ragel_machine
//     if (!vcf_batch_is_empty(status->current_batch))
//     {
//         list_item_t *item = list_item_new(get_num_vcf_batches(file), 1, status->current_batch);
//         list_insert_item(item, batches_list);
//         printf("Batch added - %zu records (last)\n", status->current_batch->records->size);
//     }

    if ( cs ) {
        LOG_INFO("Last state was not the expected");
    } 

    LOG_INFO_F("Records read = %zu\n", status->num_records);
    LOG_INFO_F("Samples per record = %zu\n", get_num_vcf_samples(file));

    // Free status->current_xxx pointers if not needed in another module
    vcf_reader_status_free(status);

    return cs ;
}

int vcf_gzip_read_and_parse(size_t batch_lines, vcf_file_t *file) {
    assert(file);
    assert(batch_lines > 0);
    
    int cs = 0;
    char *p, *pe;

    vcf_reader_status *status = vcf_reader_status_new(batch_lines, 0);
    
    LOG_DEBUG("Using file-IO functions for file loading\n");

    size_t max_len = 256;
    int eof_found = 0;
    int c = 0, i = 0, lines = 0;
    char *aux;
    char *data = (char*) calloc (max_len, sizeof(char));

    // ZLIB variables
    int ret;
    unsigned have = 0, consumed = 0;
    z_stream strm;
    unsigned char in[CHUNK];
    unsigned char out[CHUNK];

    // ZLIB stream initialization
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit2 (&strm, 15 + 32);    // Using inflateInit2 for GZIP support
    if (ret != Z_OK) {
        LOG_ERROR("gzipped file could not be decompressed");
        return 1;
    }


    do {
        strm.avail_in = fread(in, 1, CHUNK, file->fd);
        if (ferror(file->fd)) {
            (void)inflateEnd(&strm);
            return Z_ERRNO;
        }
        if (strm.avail_in == 0)
            break;
        strm.next_in = in;

        /* run inflate() on input until output buffer not full */
        do {
            strm.avail_out = CHUNK;
            strm.next_out = out;
            ret = inflate(&strm, Z_NO_FLUSH);
            assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
            switch (ret) {
            case Z_NEED_DICT:
                ret = Z_DATA_ERROR;     /* and fall through */
            case Z_DATA_ERROR:
            case Z_MEM_ERROR:
                (void)inflateEnd(&strm);
                return ret;
            }
            have = CHUNK - strm.avail_out;
            
            for (consumed = 0; consumed < have && !eof_found; consumed++) {
                c = out[consumed];

                if (c != EOF) {
                    max_len = consume_input(c, &data, max_len, i);
                    if (c == '\n') {
                        lines++;
                    }
                    i++;
                    (file->data_len)++;
                } else {
                    eof_found = 1;
                }

                // Process batch
                if (lines == batch_lines) {
                    data[file->data_len] = '\0';
                    p = data;
                    pe = p + file->data_len;
                    cs |= run_vcf_parser(p, pe, batch_lines, file, status);
                    file->data_len = 0;

                    // Setup for next batch
                    status->current_batch = vcf_batch_new(batch_lines);
                    i = 0;
                    lines = 0;
                    data = (char*) calloc (max_len, sizeof(char));
                }
            }

        } while (strm.avail_out == 0);

        /* done when inflate() says it's done */
    } while (ret != Z_STREAM_END);

    // Consume last batch
    if (lines > 0 && lines < batch_lines) {
        data[file->data_len] = '\0';
        p = data;
        pe = p + file->data_len;
        cs |= run_vcf_parser(p, pe, batch_lines, file, status);
        file->data_len = 0;
    } else {
        // Empty batch, the data buffer must be free'd
        free(data);
    }

    if ( cs ) {
        LOG_INFO("Last state was not the expected");
    } 

    LOG_INFO_F("Records read = %zu\n", status->num_records);
    LOG_INFO_F("Samples per record = %zu\n", get_num_vcf_samples(file));

    // Free status->current_xxx pointers if not needed in another module
    vcf_reader_status_free(status);

    /* clean up and return */
    (void)inflateEnd(&strm);

    return cs ;
}

int vcf_gzip_read_and_parse_bytes(size_t batch_bytes, vcf_file_t *file) {
    assert(file);
    assert(batch_bytes > 0);
    
    int cs = 0;
    char *p, *pe;

    vcf_reader_status *status = vcf_reader_status_new(0, 0);
    
    LOG_DEBUG("Using file-IO functions for file loading\n");

    size_t max_len = 256;
    int eof_found = 0;
    int c = 0, i = 0, lines = 0;
    char *aux;
    char *data = (char*) calloc (max_len, sizeof(char));

    // ZLIB variables
    int ret;
    unsigned have = 0, consumed = 0;
    z_stream strm;
    unsigned char in[CHUNK];
    unsigned char out[CHUNK];

    // ZLIB stream initialization
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit2 (&strm, 15 + 32);    // Using inflateInit2 for GZIP support
    if (ret != Z_OK) {
        LOG_ERROR("gzipped file could not be decompressed");
        return 1;
    }


    do {
        strm.avail_in = fread(in, 1, CHUNK, file->fd);
        if (ferror(file->fd)) {
            (void)inflateEnd(&strm);
            return Z_ERRNO;
        }
        if (strm.avail_in == 0)
            break;
        strm.next_in = in;

        /* run inflate() on input until output buffer not full */
        do {
            strm.avail_out = CHUNK;
            strm.next_out = out;
            ret = inflate(&strm, Z_NO_FLUSH);
            assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
            switch (ret) {
            case Z_NEED_DICT:
                ret = Z_DATA_ERROR;     /* and fall through */
            case Z_DATA_ERROR:
            case Z_MEM_ERROR:
                (void)inflateEnd(&strm);
                return ret;
            }
            have = CHUNK - strm.avail_out;
            
            for (consumed = 0; consumed < have && !eof_found; consumed++) {
                c = out[consumed];

                if (c != EOF) {
                    max_len = consume_input(c, &data, max_len, i);
                    i++;
                    (file->data_len)++;
                    if (c == '\n') {
                        lines++;

                        // Process batch
                        if (i >= batch_bytes) {
                            data[i+1] = '\0';
                            p = data;
                            pe = p + file->data_len;
                            cs |= run_vcf_parser(p, pe, 0, file, status);
                            file->data_len = 0;

                            // Setup for next batch
                            status->current_batch = vcf_batch_new(500);
                            i = 0;
                            lines = 0;
                            data = (char*) calloc (max_len, sizeof(char));
                        }
                    }

                } else {
                    eof_found = 1;
                }
            }

        } while (strm.avail_out == 0);

        /* done when inflate() says it's done */
    } while (ret != Z_STREAM_END);

    // Consume last batch
    if (i > 0) {
        data[i+1] = '\0';
        p = data;
        pe = p + file->data_len;
        cs |= run_vcf_parser(p, pe, 0, file, status);
    }

    if ( cs ) {
        LOG_INFO("Last state was not the expected");
    } 

    LOG_INFO_F("Records read = %zu\n", status->num_records);
    LOG_INFO_F("Samples per record = %zu\n", get_num_vcf_samples(file));

    // Free status->current_xxx pointers if not needed in another module
    vcf_reader_status_free(status);

    /* clean up and return */
    (void)inflateEnd(&strm);

    return cs ;
}


/* **********************************************
 *                  Only reading                *
 * **********************************************/

int vcf_light_read(size_t batch_lines, vcf_file_t *file) {
    assert(file);
    assert(batch_lines > 0);
    
    LOG_DEBUG("Using file-IO functions for file loading\n");

    size_t max_len = 256;
    int eof_found = 0;
    char *aux;

    // Read text of a batch and call ragel parser in a loop
    while (!eof_found) {
        char *data = (char*) calloc (max_len, sizeof(char));
        int c = 0;
        int lines = 0;

        for (int i = 0; !eof_found && lines < batch_lines; i++) {
            c = fgetc(file->fd);

            if (c != EOF) {
                max_len = consume_input(c, &data, max_len, i);
                if (c == '\n') {
                    lines++;
                }
            } else {
                eof_found = 1;
            }
        }

        list_item_t *item = list_item_new(get_num_vcf_batches(file), 1, data);
        list_insert_item(item, file->text_batches);
//         printf("Text batch inserted = '%.*s'\n", 200, data);
    }

    return 0;
}

int vcf_light_read_bytes(size_t batch_bytes, vcf_file_t *file) {
    assert(file);
    assert(batch_bytes > 0);
    
    LOG_DEBUG("Using file-IO functions for file loading\n");

    size_t max_len = 256;
    int eof_found = 0;
    int last_idx = 0;
    int i = 0;

    // Read text of a batch and call ragel parser in a loop
    while (!eof_found) {
        char *data = (char*) calloc (max_len, sizeof(char));
        int c = 0;
        int lines = 0;

        for (i = 0; !eof_found; i++) {
            c = fgetc(file->fd);

            if (c != EOF) {
                max_len = consume_input(c, &data, max_len, i);
                if (c == '\n') {
                    lines++;
                    if (i >= batch_bytes) {
                        break;
                    }
                }
            } else {
                eof_found = 1;
            }
        }

        data[i+1] = '\0';

        // Enqueue current batch
        list_item_t *item = list_item_new(get_num_vcf_batches(file), 1, data);
        list_insert_item(item, file->text_batches);
//             printf("Text batch inserted = '%s'\n", data);
    }

    return 0;
}

int vcf_gzip_light_read(size_t batch_lines, vcf_file_t *file) {
    assert(file);
    assert(batch_lines > 0);
    
    LOG_DEBUG("Using file-IO functions for file loading\n");

    size_t max_len = 256;
    int eof_found = 0;
    int c = 0, i = 0, lines = 0;
    char *aux;
    char *data = (char*) calloc (max_len, sizeof(char));

    // ZLIB variables
    int ret;
    unsigned have;
    z_stream strm;
    unsigned char in[CHUNK];
    unsigned char out[CHUNK];

    // ZLIB stream initialization
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit2 (&strm, 15 + 32);    // Using inflateInit2 for GZIP support
    if (ret != Z_OK) {
        LOG_ERROR("gzipped file could not be decompressed");
        return 1;
    }


    do {
        strm.avail_in = fread(in, 1, CHUNK, file->fd);
        if (ferror(file->fd)) {
            (void)inflateEnd(&strm);
            return Z_ERRNO;
        }
        if (strm.avail_in == 0)
            break;
        strm.next_in = in;

        /* run inflate() on input until output buffer not full */
        do {
            strm.avail_out = CHUNK;
            strm.next_out = out;
            ret = inflate(&strm, Z_NO_FLUSH);
            assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
            switch (ret) {
            case Z_NEED_DICT:
                ret = Z_DATA_ERROR;     /* and fall through */
            case Z_DATA_ERROR:
            case Z_MEM_ERROR:
                (void)inflateEnd(&strm);
                return ret;
            }
            have = CHUNK - strm.avail_out;
            
            for (int j = 0; j < have && !eof_found; j++) {
                c = out[j];

                if (c != EOF) {
                    max_len = consume_input(c, &data, max_len, i);
                    if (c == '\n') {
                        lines++;
                    }
                    i++;
                    (file->data_len)++;
                } else {
                    eof_found = 1;
                }

                // Process batch
                if (lines == batch_lines) {
                    list_item_t *item = list_item_new(get_num_vcf_batches(file), 1, data);
                    list_insert_item(item, file->text_batches);

                    // Setup for next batch
                    i = 0;
                    lines = 0;
                    data = (char*) calloc (max_len, sizeof(char));
                }
            }

        } while (strm.avail_out == 0);

        /* done when inflate() says it's done */
    } while (ret != Z_STREAM_END);

    // Consume last batch
    if (lines > 0 && lines < batch_lines) {
        list_item_t *item = list_item_new(get_num_vcf_batches(file), 1, data);
        list_insert_item(item, file->text_batches);
    }

    /* clean up and return */
    (void)inflateEnd(&strm);
    return ret == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
}

int vcf_gzip_light_read_bytes(size_t batch_bytes, vcf_file_t *file) {
    assert(file);
    assert(batch_bytes > 0);
    
    LOG_DEBUG("Using file-IO functions for file loading\n");

    size_t max_len = 256;
    int eof_found = 0;
    int c = 0, i = 0, lines = 0;
    char *aux;
    char *data = (char*) calloc (max_len, sizeof(char));

    // ZLIB variables
    int ret;
    unsigned have;
    z_stream strm;
    unsigned char in[CHUNK];
    unsigned char out[CHUNK];

    // ZLIB stream initialization
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit2 (&strm, 15 + 32);    // Using inflateInit2 for GZIP support
    if (ret != Z_OK) {
        LOG_ERROR("gzipped file could not be decompressed");
        return 1;
    }


    do {
        strm.avail_in = fread(in, 1, CHUNK, file->fd);
        if (ferror(file->fd)) {
            (void)inflateEnd(&strm);
            return Z_ERRNO;
        }
        if (strm.avail_in == 0)
            break;
        strm.next_in = in;

        /* run inflate() on input until output buffer not full */
        do {
            strm.avail_out = CHUNK;
            strm.next_out = out;
            ret = inflate(&strm, Z_NO_FLUSH);
            assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
            switch (ret) {
            case Z_NEED_DICT:
                ret = Z_DATA_ERROR;     /* and fall through */
            case Z_DATA_ERROR:
            case Z_MEM_ERROR:
                (void)inflateEnd(&strm);
                return ret;
            }
            have = CHUNK - strm.avail_out;
            
            for (int j = 0; j < have && !eof_found; j++) {
                c = out[j];

                if (c != EOF) {
                    max_len = consume_input(c, &data, max_len, i);
                    i++;
                    (file->data_len)++;
                    if (c == '\n') {
                        lines++;

                        // Process batch
                        if (i >= batch_bytes) {
                            data[i+1] = '\0';
                            list_item_t *item = list_item_new(get_num_vcf_batches(file), 1, data);
                            list_insert_item(item, file->text_batches);

                            // Setup for next batch
                            i = 0;
                            lines = 0;
                            data = (char*) calloc (max_len, sizeof(char));
                        }
                    }

                } else {
                    eof_found = 1;
                }
            }

        } while (strm.avail_out == 0);

        /* done when inflate() says it's done */
    } while (ret != Z_STREAM_END);

    // Consume last batch
//     if (lines > 0 && lines < batch_size) {
    if (i > 0) {
        data[i+1] = '\0';
        list_item_t *item = list_item_new(get_num_vcf_batches(file), 1, data);
        list_insert_item(item, file->text_batches);
    }

    /* clean up and return */
    (void)inflateEnd(&strm);
    return ret == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
}


/* **********************************************
 *      Only reading from multiple files        *
 * **********************************************/

int vcf_light_multiread(list_t **text_lists, size_t batch_lines, vcf_file_t **files, size_t num_files) {
    LOG_DEBUG("Using file-IO functions for file loading\n");

    // Initialize file-private variables
    size_t max_len[num_files];
    
    for (int i = 0; i < num_files; i++) {
        max_len[i] = 256;
        files[i]->data_len = 0;
    }
    
    int num_eof_found = 0;
    int eof_found[num_files];
    memset(eof_found, 0, num_files * sizeof(int));

    // Read text of a batch from each file and call ragel parser in a loop
    while (num_eof_found < num_files) {
        for (int f = 0; f < num_files; f++) {
            if (eof_found[f]) {
                continue;
            }

            char *data = (char*) calloc (max_len[f], sizeof(char));
            int c = 0;
            int lines = 0;

            int i;
            for (i = 0; !eof_found[f] && lines < batch_lines; i++) {
                c = fgetc(files[f]->fd);

                if (c != EOF) {
                    max_len[f] = consume_input(c, &data, max_len[f], i);
                    if (c == '\n') {
                        lines++;
                    }
                } else {
                    printf("EOF found in file %d\n", f);
                    eof_found[f] = 1;
                    num_eof_found++;
                    list_decr_writers(text_lists[f]);
                }
            }
            
            data[i] = '\0';
            
            // Enqueue current batch
            list_item_t *item = list_item_new(get_num_vcf_batches(files[f]), 1, data);
            list_insert_item(item, text_lists[f]);
//             printf("Text batch inserted = '%.*s'\n", 200, data + strlen(data) - 200);
        }
    }

    return 0;
}



/* **********************************************
 *              Auxiliary functions             *
 * **********************************************/

size_t consume_input(int c, char **data, size_t max_len, int position_in_data) {
    assert(data);
    
    (*data)[position_in_data] = c;
    // Text too long to be stored in 'data', realloc
    if (position_in_data == max_len - 1) {
        char *aux = realloc(*data, max_len + 10240);
        if (aux) {
            *data = aux;
            return max_len + 10240;
        } else {
            LOG_FATAL("Could not allocate enough memory for reading input VCF file\n");
        }
    }
    return max_len;
}
