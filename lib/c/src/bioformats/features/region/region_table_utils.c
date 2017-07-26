#include "region_table_utils.h"

region_table_t *parse_regions(char *input_regions, int as_positions, const char *url, const char *species, const char *version) {
    region_table_t *regions_table = new_region_table_from_ws(url, species, version);

    char *saveptr = NULL, *token;
    size_t token_len;

    int num_regions;
    char **regions_data = split(input_regions, ",", &num_regions);
    
    region_t *regions[num_regions];
    
    for (int i = 0; i < num_regions; i++) {
        // Set chromosome
        token = strtok_r(regions_data[i], ":", &saveptr);
        token_len = strlen(token);
        char *chromosome = strndup(token, token_len);
        
        // Set start position
        token = strtok_r(NULL, "-", &saveptr);
        size_t start_position, end_position;
        start_position = (token != NULL) ? atol(token) : 1;
        
        // Set end position
        token = strtok_r(NULL, "-", &saveptr);
        if (token != NULL) {
            end_position = atol(token);
        } else {
            if (as_positions) {
                end_position = start_position;
            } else {
                end_position = UINT_MAX;
            }
        }
        
        regions[i] = region_new(chromosome, start_position, end_position, NULL, NULL);
        
        LOG_DEBUG_F("region '%s:%u-%u'\n", regions[i]->chromosome, regions[i]->start_position, regions[i]->end_position);
    }
    
    insert_regions(regions, num_regions, regions_table);
    finish_region_table_loading(regions_table);
    
    for (int i = 0; i < num_regions; i++) {
        free(regions_data[i]);
        free(regions[i]);
    }
    free(regions_data);
    
    return regions_table;
}

region_table_t *parse_regions_from_gff_file(char *filename, const char *url, const char *species, const char *version) {
    gff_file_t *file = gff_open(filename);
    if (file == NULL) {
        return NULL;
    } 
    
    region_table_t *regions_table = new_region_table_from_ws(url, species, version);
    
    int ret_code = 0;
    size_t max_batches = 20, batch_size = 2000;
    list_t *read_list = (list_t*) malloc (sizeof(list_t));
    list_init("batches", 1, max_batches, read_list);
    
    #pragma omp parallel sections
    {
        // The producer reads the GFF file
        #pragma omp section
        {
            LOG_DEBUG_F("Thread %d reads the GFF file\n", omp_get_thread_num());
            ret_code = gff_read_batches(read_list, batch_size, file);
            list_decr_writers(read_list);
            
            if (ret_code) {
                LOG_FATAL_F("Error while reading GFF file %s (%d)\n", filename, ret_code);
            }
        }
        
        // The consumer inserts regions in the structure 
        #pragma omp section
        {
            list_item_t *item = NULL;
            gff_batch_t *batch;
            gff_record_t *record;
            
            region_t *regions_batch[REGIONS_CHUNKSIZE];
            int avail_regions = 0;
            
            while (( item = list_remove_item(read_list) )) {
                batch = item->data_p;
                // For each record in the batch, generate a new region
                for (int i = 0; i < batch->records->size; i++) {
                    record = batch->records->items[i];
                    
                    region_t *region = region_new(strndup(record->sequence, record->sequence_len), 
                                                  record->start, record->end,
                                                  record->strand ? strndup(&record->strand, 1) : NULL, 
                                                  record->feature ? strndup(record->feature, record->feature_len) : NULL);
                    
                    LOG_DEBUG_F("region '%s:%u-%u'\n", region->chromosome, region->start_position, region->end_position);
                    
                    regions_batch[avail_regions++] = region;
                    
                    // Save when the recommended size is reached
                    if (avail_regions == REGIONS_CHUNKSIZE) {
                        insert_regions(regions_batch, avail_regions, regions_table);
                        for (int i = 0; i < avail_regions; i++) {
                            free(regions_batch[i]);
                        }
                        avail_regions = 0;
                    }
                }
               
                gff_batch_free(batch);
                list_item_free(item);
            }
            
            // Save the remaining regions that did not fill a batch
            if (avail_regions > 0) {
                insert_regions(regions_batch, avail_regions, regions_table);
                for (int i = 0; i < avail_regions; i++) {
                    free(regions_batch[i]);
                }
                avail_regions = 0;
            }
        }
    }
    
    finish_region_table_loading(regions_table);
    
    list_free_deep(read_list, NULL);
    
    gff_close(file, 1);
    
    return regions_table;
}

region_table_t *parse_regions_from_bed_file(char *filename, const char *url, const char *species, const char *version) {
    bed_file_t *file = bed_open(filename);
    if (file == NULL) {
        return NULL;
    } 
    
    region_table_t *regions_table = new_region_table_from_ws(url, species, version);
    
    int ret_code = 0;
    size_t max_batches = 20, batch_size = 2000;
    list_t *read_list = (list_t*) malloc (sizeof(list_t));
    list_init("batches", 1, max_batches, read_list);
    
    #pragma omp parallel sections
    {
        // The producer reads the bed file
        #pragma omp section
        {
            LOG_DEBUG_F("Thread %d reads the BED file\n", omp_get_thread_num());
            ret_code = bed_read_batches(read_list, batch_size, file);
            list_decr_writers(read_list);
            
            if (ret_code) {
                LOG_FATAL_F("Error while reading BED file %s (%d)\n", filename, ret_code);
            }
        }
        
        // The consumer inserts regions in the structure 
        #pragma omp section
        {    
            list_item_t *item = NULL;
            bed_batch_t *batch;
            bed_record_t *record;
            
            region_t *regions_batch[REGIONS_CHUNKSIZE];
            int avail_regions = 0;
            
            while (( item = list_remove_item(read_list) )) {
                batch = item->data_p;
                // For each record in the batch, generate a new region
                for (int i = 0; i < batch->records->size; i++) {
                    record = batch->records->items[i];
                    
                    region_t *region = region_new(strndup(record->sequence, record->sequence_len), 
                                                  record->start, record->end, strndup(&record->strand, 1),
                                                  NULL);
                    
                    LOG_DEBUG_F("region '%s:%u-%u'\n", region->chromosome, region->start_position, region->end_position);
                    
                    regions_batch[avail_regions++] = region;
                    
                    // Save when the recommended size is reached
                    if (avail_regions == REGIONS_CHUNKSIZE) {
                        insert_regions(regions_batch, avail_regions, regions_table);
                        for (int i = 0; i < avail_regions; i++) {
                            free(regions_batch[i]);
                        }
                        avail_regions = 0;
                    }
                }
               
                bed_batch_free(batch);
                list_item_free(item);
            }
            
            // Save the remaining regions that did not fill a batch
            if (avail_regions > 0) {
                insert_regions(regions_batch, avail_regions, regions_table);
                for (int i = 0; i < avail_regions; i++) {
                    free(regions_batch[i]);
                }
                avail_regions = 0;
            }
        }
    }
    
    list_free_deep(read_list, NULL);
    
    bed_close(file, 1);
    
    return regions_table;
}

int region_table_parse_from_string(char *input_regions, region_table_t *regions_table) {
    int as_positions = 1;
    char *str_1 = input_regions;
    char *str_2 = (char*) malloc (64 * sizeof(char));
    char *saveptr1 = NULL, *saveptr2 = NULL;
    char *token, *subtoken;
    size_t token_len, subtoken_len;

    int i = 0;
    while ((token = strtok_r(str_1, ",", &saveptr1)) != NULL) {

        region_t *region = (region_t*) malloc (sizeof(region_t));

        token_len = strlen(token);
        
        LOG_DEBUG_F("token = %s, len = %zu\n", token, token_len);
        
        strncpy(str_2, token, 63);
        str_2[token_len] = '\0';
        
        // Set chromosome
        subtoken = strtok_r(str_2, ":", &saveptr2);
        subtoken_len = strlen(subtoken);
        region->chromosome = (char*) malloc ((subtoken_len+1) * sizeof(char));
        strncpy(region->chromosome, subtoken, subtoken_len);
        region->chromosome[subtoken_len] = '\0';
        
        // Set start position
        subtoken = strtok_r(NULL, "-", &saveptr2);
        region->start_position = (subtoken != NULL) ? atol(subtoken) : 1;
        
        // Set end position
        subtoken = strtok_r(NULL, "-", &saveptr2);
        if (subtoken != NULL) {
            region->end_position = atol(subtoken);
        } else {
            if (as_positions) {
                region->end_position = region->start_position;
            } else {
                region->end_position = UINT_MAX;
            }
        }
        
        LOG_DEBUG_F("region '%s:%u-%u'\n", region->chromosome, region->start_position, region->end_position);

	region->strand = NULL;
	region->type = NULL;
        insert_region(region, regions_table);
        
        str_1 = NULL;
        
        i++;
    }

    free(str_1); 
    free(str_2);

    return 1;
}

int region_table_parse_from_gff_file(char *filename, region_table_t *regions_table) {
    gff_file_t *file = gff_open(filename);
    if (file == NULL) {
        return 0;
    } 
    
    int ret_code = 0;
    size_t max_batches = 20, batch_size = 2000;
    list_t *read_list = (list_t*) malloc (sizeof(list_t));
    list_init("batches", 1, max_batches, read_list);
    
    #pragma omp parallel sections
    {
        // The producer reads the GFF file
        #pragma omp section
        {
            LOG_DEBUG_F("Thread %d reads the GFF file\n", omp_get_thread_num());
            ret_code = gff_read_batches(read_list, batch_size, file);
            list_decr_writers(read_list);
            
            if (ret_code) {
                LOG_FATAL_F("Error while reading GFF file %s (%d)\n", filename, ret_code);
            }
        }
        
        // The consumer inserts regions in the structure 
        #pragma omp section
        {
            list_item_t *item = NULL;
            gff_batch_t *batch;
            gff_record_t *record;
            
            region_t *regions_batch[REGIONS_CHUNKSIZE];
            int avail_regions = 0;
            
            while (( item = list_remove_item(read_list) )) {
                batch = item->data_p;
                // For each record in the batch, generate a new region
                for (int i = 0; i < batch->records->size; i++) {
                    record = batch->records->items[i];
                    
                    region_t *region = region_new(strndup(record->sequence, record->sequence_len), 
                                                  record->start, record->end,
                                                  record->strand ? strndup(&record->strand, 1) : NULL, 
                                                  record->feature ? strndup(record->feature, record->feature_len) : NULL);
                    
                    LOG_DEBUG_F("region '%s:%u-%u'\n", region->chromosome, region->start_position, region->end_position);
                    
                    regions_batch[avail_regions++] = region;
                    
                    // Save when the recommended size is reached
                    if (avail_regions == REGIONS_CHUNKSIZE) {
                        insert_regions(regions_batch, avail_regions, regions_table);
                        for (int i = 0; i < avail_regions; i++) {
                            free(regions_batch[i]);
                        }
                        avail_regions = 0;
                    }
                }
               
                gff_batch_free(batch);
                list_item_free(item);
            }
            
            // Save the remaining regions that did not fill a batch
            if (avail_regions > 0) {
                insert_regions(regions_batch, avail_regions, regions_table);
                for (int i = 0; i < avail_regions; i++) {
                    free(regions_batch[i]);
                }
                avail_regions = 0;
            }
        }
    }
    
    finish_region_table_loading(regions_table);
    
    list_free_deep(read_list, NULL);
    
    gff_close(file, 1);
    
    return 1;
}

