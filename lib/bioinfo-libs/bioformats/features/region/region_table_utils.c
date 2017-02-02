#include "region_table_utils.h"


region_table_t *parse_regions(char *input_regions, int as_positions, const char *url, const char *species, const char *version) {
    region_table_t *regions_table = create_table(url, species, version);

    char *str_1 = input_regions;
    char *str_2 = (char*) malloc (64 * sizeof(char));
    char *saveptr1, *saveptr2;
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
        
        insert_region(region, regions_table);
        
        str_1 = NULL;
        
        i++;
    }

    free(str_1); 
    free(str_2);

    return regions_table;
}

region_table_t *parse_regions_from_gff_file(char *filename, const char *url, const char *species, const char *version) {
    gff_file_t *file = gff_open(filename);
    if (file == NULL) {
        return NULL;
    } 
    
    region_table_t *regions_table = create_table(url, species, version);
    
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
            while ( (item = list_remove_item(read_list)) != NULL ) {
                batch = item->data_p;
                // For each record in the batch, generate a new region
                for (int i = 0; i < batch->records->size; i++) {
                    record = batch->records->items[i];
                    
                    region_t *region = (region_t*) malloc (sizeof(region_t));
                    region->chromosome = strndup(record->sequence, record->sequence_len);
                    region->start_position = record->start;
                    region->end_position = record->end;
                    LOG_DEBUG_F("region '%s:%u-%u'\n", region->chromosome, region->start_position, region->end_position);
                    
                    insert_region(region, regions_table);
                }
               
                gff_batch_free(batch);
                list_item_free(item);
            }
        }
    }
    
    list_free_deep(read_list, NULL);
    
    gff_close(file, 1);
    
    return regions_table;
}

region_table_t *parse_regions_from_bed_file(char *filename, const char *url, const char *species, const char *version) {
    bed_file_t *file = bed_open(filename);
    if (file == NULL) {
        return NULL;
    } 
    
    region_table_t *regions_table = create_table(url, species, version);
    
    int ret_code = 0;
    size_t max_batches = 20, batch_size = 2000;
    list_t *read_list = (list_t*) malloc (sizeof(list_t));
    list_init("batches", 1, max_batches, read_list);
    
    #pragma omp parallel sections
    {
        // The producer reads the bed file
        #pragma omp section
        {
            LOG_DEBUG_F("Thread %d reads the bed file\n", omp_get_thread_num());
            ret_code = bed_read_batches(read_list, batch_size, file);
            list_decr_writers(read_list);
            
            if (ret_code) {
                LOG_FATAL_F("Error while reading bed file %s (%d)\n", filename, ret_code);
            }
        }
        
        // The consumer inserts regions in the structure 
        #pragma omp section
        {    
            list_item_t *item = NULL;
            bed_batch_t *batch;
            bed_record_t *record;
            while ( (item = list_remove_item(read_list)) != NULL ) {
                batch = item->data_p;
                // For each record in the batch, generate a new region
                for (int i = 0; i < batch->records->size; i++) {
                    record = batch->records->items[i];
                    
                    region_t *region = (region_t*) malloc (sizeof(region_t));
                    region->chromosome = strndup(record->sequence, record->sequence_len);
                    region->start_position = record->start;
                    region->end_position = record->end;
                    LOG_DEBUG_F("region '%s:%u-%u'\n", region->chromosome, region->start_position, region->end_position);
                    
                    insert_region(region, regions_table);
                }
               
                bed_batch_free(batch);
                list_item_free(item);
            }
        }
    }
    
    list_free_deep(read_list, NULL);
    
    bed_close(file, 1);
    
    return regions_table;
}
