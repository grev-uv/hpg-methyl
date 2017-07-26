#include "region_table.h"

static void prepare_region_table_statements(region_table_t *table) {
    sqlite3 *db = table->storage;

    // Insert regions
    char *sql_insert = "INSERT INTO regions VALUES (?1, ?2, ?3, ?4, ?5)";
    if (sqlite3_prepare_v2(table->storage, sql_insert, strlen(sql_insert), &(table->insert_region_stmt), NULL) != SQLITE_OK) {
        LOG_FATAL_F("Could not create regions database: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
    }

    // Find exact regions
    char *sql_find_exact = "SELECT COUNT(*) FROM regions WHERE chromosome = ?1 AND start = ?2 AND end = ?3";
    if (sqlite3_prepare_v2(table->storage, sql_find_exact, strlen(sql_find_exact), &(table->find_exact_region_stmt), NULL) != SQLITE_OK) {
        LOG_FATAL_F("Could not create regions database: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
    }

    char *sql_find_exact_type = "SELECT COUNT(*) FROM regions WHERE chromosome = ?1 AND start = ?2 AND end = ?3 AND type = ?4";
    if (sqlite3_prepare_v2(table->storage, sql_find_exact_type, strlen(sql_find_exact_type), &(table->find_exact_region_type_stmt), NULL) != SQLITE_OK) {
        LOG_FATAL_F("Could not create regions database: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
    }

    // Find regions
    char *sql_find = "SELECT COUNT(*) FROM regions WHERE chromosome = ?1 AND start <= ?3 AND end >= ?2";
    if (sqlite3_prepare_v2(table->storage, sql_find, strlen(sql_find), &(table->find_region_stmt), NULL) != SQLITE_OK) {
        LOG_FATAL_F("Could not create regions database: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
    }

    char *sql_find_type = "SELECT COUNT(*) FROM regions WHERE chromosome = ?1 AND start <= ?3 AND end >= ?2 AND type = ?4";
    if (sqlite3_prepare_v2(table->storage, sql_find_type, strlen(sql_find_type), &(table->find_region_type_stmt), NULL) != SQLITE_OK) {
        LOG_FATAL_F("Could not create regions database: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
    }

    // Remove regions
    char *sql_remove_exact = "DELETE FROM regions WHERE chromosome = ?1 AND start = ?2 AND end = ?3";
    if (sqlite3_prepare_v2(table->storage, sql_remove_exact, strlen(sql_remove_exact), &(table->remove_exact_region_stmt), NULL) != SQLITE_OK) {
        LOG_FATAL_F("Could not create regions database: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
    }

    char *sql_remove = "DELETE FROM regions WHERE chromosome = ?1 AND start <= ?3 AND end >= ?2";
    if (sqlite3_prepare_v2(table->storage, sql_remove, strlen(sql_remove), &(table->remove_region_stmt), NULL) != SQLITE_OK) {
        LOG_FATAL_F("Could not create regions database: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
    }

    // Query chromosomes
    char *sql_get_chromosome = "SELECT * FROM regions WHERE chromosome = ?1";
    if (sqlite3_prepare_v2(table->storage, sql_get_chromosome, strlen(sql_get_chromosome), &(table->get_chromosome_stmt), NULL) != SQLITE_OK) {
        LOG_FATAL_F("Could not create regions database: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
    }

    char *sql_count_in_chromosome = "SELECT COUNT(*) FROM regions WHERE chromosome = ?1";
    if (sqlite3_prepare_v2(table->storage, sql_count_in_chromosome, strlen(sql_count_in_chromosome), &(table->count_in_chromosome_stmt), NULL) != SQLITE_OK) {
        LOG_FATAL_F("Could not create regions database: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
    }

}

region_table_t *new_region_table_from_ws(const char *url, const char *species, const char *version) {
    int num_chromosomes;
    
    // Initialize structure
    region_table_t *table = (region_table_t*) malloc (sizeof(region_table_t));
    table->ordering = get_chromosome_order(url, species, version, &num_chromosomes);
    table->max_chromosomes = num_chromosomes;
    table->chunks = kh_init(stats_chunks);
    table->is_ready = 0;

    srand( time(NULL) );
    int suffix = rand();
    char db_name[32];
    sprintf(db_name, "/tmp/regions_%d.db", suffix);
    create_regions_db(db_name, REGIONS_CHUNKSIZE, &(table->storage));
    prepare_region_table_statements(table);
    
    return table;
}

region_table_t *new_region_table(int num_chromosomes, char **chromosomes) {
    // Initialize structure
    region_table_t *table = (region_table_t*) malloc (sizeof(region_table_t));
    table->ordering = chromosomes;
    table->max_chromosomes = num_chromosomes;
    table->chunks = kh_init(stats_chunks);
    table->is_ready = 0;

    srand( time(NULL) );
    int suffix = rand();
    char db_name[32];
    sprintf(db_name, "/tmp/regions_%d.db", suffix);
    create_regions_db(db_name, REGIONS_CHUNKSIZE, &(table->storage));
    prepare_region_table_statements(table);
    
    return table;
}

void free_region_table(region_table_t* regions) {
    // Destroy prepared statements
    sqlite3_finalize(regions->insert_region_stmt);
    sqlite3_finalize(regions->find_exact_region_stmt);
    sqlite3_finalize(regions->find_exact_region_type_stmt);
    sqlite3_finalize(regions->find_region_stmt);
    sqlite3_finalize(regions->find_region_type_stmt);
    sqlite3_finalize(regions->remove_exact_region_stmt);
    sqlite3_finalize(regions->remove_region_stmt);
    sqlite3_finalize(regions->get_chromosome_stmt);
    sqlite3_finalize(regions->count_in_chromosome_stmt);

    // Close database
    sqlite3_close_v2(regions->storage);
    
    // Free ordering array
    char **ordering = regions->ordering;
    for (int i = 0; i < regions->max_chromosomes; i++) {
        free(ordering[i]);
    }
    free(ordering);

    // Destroy the chunks hashtable
    kh_destroy(stats_chunks, regions->chunks);

    free(regions);
}


void finish_region_table_loading(region_table_t *table) {
    // Save chunks from hashtable
    insert_chunk_hash(REGIONS_CHUNKSIZE, table->chunks, table->storage);
    // Index contents
    create_regions_index(table->storage);
    // Mark as ready!
    table->is_ready = 1;
}


/* ******************************
 *  		Regions		*
 * ******************************/

int insert_region(region_t *region, region_table_t *table) {
    return insert_regions(&region, 1, table);
}

int insert_regions(region_t **regions, int num_regions, region_table_t *table) {
    int rc;
    sqlite3_stmt *stmt = table->insert_region_stmt;
    sqlite3* db = table->storage;

    char *sql_begin = "BEGIN TRANSACTION";
    rc = exec_sql(sql_begin, db);
    if (rc != SQLITE_OK) {
        LOG_ERROR_F("Could not insert regions: %s (%d)\n", 
                    sqlite3_errmsg(db), sqlite3_errcode(db));
        return rc;
    }
    
    for (int i = 0; i < num_regions; i++) {
        region_t *region = regions[i];

        sqlite3_bind_text(stmt, 1, region->chromosome, strlen(region->chromosome), SQLITE_STATIC);
        sqlite3_bind_int64(stmt, 2, region->start_position);
        sqlite3_bind_int64(stmt, 3, region->end_position);

        if (region->strand) { 
            sqlite3_bind_text(stmt, 4, region->strand, strlen(region->strand), SQLITE_STATIC); 
        } else {
            sqlite3_bind_null(stmt, 4);
        }
        if (region->type) {
            sqlite3_bind_text(stmt, 5, region->type, strlen(region->type), SQLITE_STATIC);
        } else {
            sqlite3_bind_null(stmt, 5);
        }

        if ((rc = sqlite3_step(stmt) == SQLITE_DONE)) {
            // Update value in chunks hashtable
            update_chunks_hash(region->chromosome, UINT_MAX, REGIONS_CHUNKSIZE, 
                               region->start_position, region->end_position, table->chunks);
        } else {
            LOG_ERROR_F("Could not insert region %s:%ld-%ld: %s (%d)\n", 
                        region->chromosome, region->start_position, region->end_position, 
                        sqlite3_errmsg(db), sqlite3_errcode(db));
        }

        sqlite3_reset(stmt);
    }
    
    char *sql_end = "END TRANSACTION";
    rc = exec_sql(sql_end, db);
    if (rc != SQLITE_OK) {
        LOG_ERROR_F("Could not insert regions: %s (%d)\n", 
                    sqlite3_errmsg(db), sqlite3_errcode(db));
        return rc;
    }

    return rc;
}


int find_exact_region(region_t *region, region_table_t *table) {
    assert(region);
    assert(table);
    
    if (!table->is_ready) { // Don't allow queries over an un-indexed DB
        finish_region_table_loading(table);
    }
    
    int count = 0;
    
    sqlite3_stmt *query_stmt = table->find_exact_region_stmt;
    sqlite3_bind_text(query_stmt, 1, region->chromosome, strlen(region->chromosome), SQLITE_STATIC);
    sqlite3_bind_int64(query_stmt, 2, region->start_position);
    sqlite3_bind_int64(query_stmt, 3, region->end_position);
    
    if (sqlite3_step(query_stmt) == SQLITE_ROW) {
        count = sqlite3_column_int(query_stmt, 0);
    } else {
        // Retry once
        sqlite3_reset(query_stmt);
        sqlite3_bind_text(query_stmt, 1, region->chromosome, strlen(region->chromosome), SQLITE_STATIC);
        sqlite3_bind_int64(query_stmt, 2, region->start_position);
        sqlite3_bind_int64(query_stmt, 3, region->end_position);
        
        if (sqlite3_step(query_stmt) == SQLITE_ROW) {
            count = sqlite3_column_int(query_stmt, 0);
        } else {
            sqlite3* db = table->storage;
            LOG_ERROR_F("Regions table failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
        }
    }
    
    sqlite3_reset(query_stmt);
    
    return count > 0;
}

int find_exact_region_by_type(region_t *region, region_table_t *table) {
    assert(region);
    assert(region->type);
    assert(table);
    
    if (!table->is_ready) { // Don't allow queries over an un-indexed DB
        finish_region_table_loading(table);
    }
    
    int count = 0;
    
    sqlite3_stmt *query_stmt = table->find_exact_region_type_stmt;
    sqlite3_bind_text(query_stmt, 1, region->chromosome, strlen(region->chromosome), SQLITE_STATIC);
    sqlite3_bind_int64(query_stmt, 2, region->start_position);
    sqlite3_bind_int64(query_stmt, 3, region->end_position);
    sqlite3_bind_text(query_stmt, 4, region->type, strlen(region->type), SQLITE_STATIC);
    
    if (sqlite3_step(query_stmt) == SQLITE_ROW) {
        count = sqlite3_column_int(query_stmt, 0);
    } else {
        // Retry once
        sqlite3_reset(query_stmt);
        sqlite3_bind_text(query_stmt, 1, region->chromosome, strlen(region->chromosome), SQLITE_STATIC);
        sqlite3_bind_int64(query_stmt, 2, region->start_position);
        sqlite3_bind_int64(query_stmt, 3, region->end_position);
    sqlite3_bind_text(query_stmt, 4, region->type, strlen(region->type), SQLITE_STATIC);
        
        if (sqlite3_step(query_stmt) == SQLITE_ROW) {
            count = sqlite3_column_int(query_stmt, 0);
        } else {
            sqlite3* db = table->storage;
            LOG_ERROR_F("Regions table failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
        }
    }
    
    sqlite3_reset(query_stmt);
    
    return count > 0;
}

int find_region(region_t *region, region_table_t *table) {
    assert(region);
    assert(table);
    
    if (!table->is_ready) { // Don't allow queries over an un-indexed DB
        finish_region_table_loading(table);
    }
    
    int count = 0;
    
    sqlite3_stmt *query_stmt = table->find_region_stmt;
    sqlite3_bind_text(query_stmt, 1, region->chromosome, strlen(region->chromosome), SQLITE_STATIC);
    sqlite3_bind_int64(query_stmt, 2, region->start_position);
    sqlite3_bind_int64(query_stmt, 3, region->end_position);
    
    if (sqlite3_step(query_stmt) == SQLITE_ROW) {
        count = sqlite3_column_int(query_stmt, 0);
    } else {
        // Retry once
        sqlite3_reset(query_stmt);
        sqlite3_bind_text(query_stmt, 1, region->chromosome, strlen(region->chromosome), SQLITE_STATIC);
        sqlite3_bind_int64(query_stmt, 2, region->start_position);
        sqlite3_bind_int64(query_stmt, 3, region->end_position);
        
        if (sqlite3_step(query_stmt) == SQLITE_ROW) {
            count = sqlite3_column_int(query_stmt, 0);
        } else {
            sqlite3* db = table->storage;
            LOG_ERROR_F("Regions table failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
        }
    }

    sqlite3_reset(query_stmt);
    
    return count > 0;
}

int find_region_by_type(region_t *region, region_table_t *table) {
    assert(region);
    assert(region->type);
    assert(table);
    
    if (!table->is_ready) { // Don't allow queries over an un-indexed DB
        finish_region_table_loading(table);
    }
    
    int count = 0;
    sqlite3* db = table->storage;
    
    sqlite3_stmt *query_stmt = table->find_region_type_stmt;
    sqlite3_bind_text(query_stmt, 1, region->chromosome, strlen(region->chromosome), SQLITE_STATIC);
    sqlite3_bind_int64(query_stmt, 2, region->start_position);
    sqlite3_bind_int64(query_stmt, 3, region->end_position);
    sqlite3_bind_text(query_stmt, 4, region->type, strlen(region->type), SQLITE_STATIC);
    
    if (sqlite3_step(query_stmt) == SQLITE_ROW) {
        count = sqlite3_column_int(query_stmt, 0);
    } else {
        // Retry once
        sqlite3_reset(query_stmt);
        sqlite3_bind_text(query_stmt, 1, region->chromosome, strlen(region->chromosome), SQLITE_STATIC);
        sqlite3_bind_int64(query_stmt, 2, region->start_position);
        sqlite3_bind_int64(query_stmt, 3, region->end_position);
        sqlite3_bind_text(query_stmt, 4, region->type, strlen(region->type), SQLITE_STATIC);
        
        if (sqlite3_step(query_stmt) == SQLITE_ROW) {
            count = sqlite3_column_int(query_stmt, 0);
        } else {
            LOG_ERROR_F("Regions table failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
        }
    }
    
    sqlite3_reset(query_stmt);
    
    return count > 0;
}

int remove_exact_region(region_t *region, region_table_t *table) {
    int rc;
    sqlite3_stmt *stmt1 = table->remove_exact_region_stmt, *stmt2;
    sqlite3* db = table->storage;
    
    char *sql_begin = "BEGIN TRANSACTION";
    rc = exec_sql(sql_begin, db);
    if (rc != SQLITE_OK) {
        LOG_ERROR_F("Could not remove region: %s (%d)\n", 
                    sqlite3_errmsg(db), sqlite3_errcode(db));
        return rc;
    }
    
    sqlite3_bind_text(stmt1, 1, region->chromosome, strlen(region->chromosome), SQLITE_STATIC);
    sqlite3_bind_int64(stmt1, 2, region->start_position);
    sqlite3_bind_int64(stmt1, 3, region->end_position);

    if ((rc = sqlite3_step(stmt1)) == SQLITE_DONE) {
        sqlite3_reset(stmt1);
        
        int chunk_start = region->start_position / REGIONS_CHUNKSIZE;
        int chunk_end = region->end_position / REGIONS_CHUNKSIZE;
        
        // Decrement features_count in chunk table
        char sql_chunks[] = "UPDATE chunk SET features_count = features_count - 1 WHERE chromosome = ?1 AND \
                             (chunk_id >= ?2 AND chunk_id <= ?3)";
        rc = sqlite3_prepare_v2(db, sql_chunks, strlen(sql_chunks), &stmt2, NULL);
        
        if (rc != SQLITE_OK) {
            LOG_ERROR_F("Could not remove region %s:%ld-%ld: %s (%d)\n", 
                        region->chromosome, region->start_position, region->end_position, 
                        sqlite3_errmsg(db), sqlite3_errcode(db));
            return rc;
        }
        
        sqlite3_bind_text(stmt2, 1, region->chromosome, strlen(region->chromosome), SQLITE_STATIC);
        sqlite3_bind_int64(stmt2, 2, chunk_start);
        sqlite3_bind_int64(stmt2, 3, chunk_end);
    
        if ((rc = sqlite3_step(stmt2)) != SQLITE_DONE) {
            LOG_ERROR_F("Could not remove region %s:%ld-%ld: %s (%d)\n", 
                        region->chromosome, region->start_position, region->end_position, 
                        sqlite3_errmsg(db), sqlite3_errcode(db));
            return rc;
        }
        
        sqlite3_finalize(stmt2);
        
    } else {
        LOG_ERROR_F("Could not remove region %s:%ld-%ld: %s (%d)\n", 
                    region->chromosome, region->start_position, region->end_position, 
                    sqlite3_errmsg(db), sqlite3_errcode(db));
        return rc;
    }
    
    char *sql_end = "END TRANSACTION";
    rc = exec_sql(sql_end, db);
    if (rc != SQLITE_OK) {
        LOG_ERROR_F("Could not remove regions: %s (%d)\n", 
                    sqlite3_errmsg(db), sqlite3_errcode(db));
        return rc;
    }

    return rc;
}


int remove_region(region_t *region, region_table_t *table) {
    int rc;
    sqlite3_stmt *stmt1 = table->remove_region_stmt, *stmt2;
    sqlite3* db = table->storage;
    
    char *sql_begin = "BEGIN TRANSACTION";
    rc = exec_sql(sql_begin, db);
    if (rc != SQLITE_OK) {
        LOG_ERROR_F("Could not remove region: %s (%d)\n", 
                    sqlite3_errmsg(db), sqlite3_errcode(db));
        return rc;
    }

    sqlite3_bind_text(stmt1, 1, region->chromosome, strlen(region->chromosome), SQLITE_STATIC);
    sqlite3_bind_int64(stmt1, 2, region->start_position);
    sqlite3_bind_int64(stmt1, 3, region->end_position);

    if ((rc = sqlite3_step(stmt1)) == SQLITE_DONE) {
        sqlite3_reset(stmt1);
        
        int affected_rows = sqlite3_changes(db);
        int chunk_start = region->start_position / REGIONS_CHUNKSIZE;
        int chunk_end = region->end_position / REGIONS_CHUNKSIZE;
        
        // Decrement features_count in chunk table depending on the number of affected rows
        // If zero is reached, don't decrement any more
        char sql_chunks[] = "UPDATE chunk SET features_count = MAX(0, features_count - ?1) WHERE chromosome = ?2 AND \
                             (chunk_id >= ?3 AND chunk_id <= ?4)";
        rc = sqlite3_prepare_v2(db, sql_chunks, strlen(sql_chunks), &stmt2, NULL);
        
        if (rc != SQLITE_OK) {
            LOG_ERROR_F("Could not remove region %s:%ld-%ld: %s (%d)\n", 
                        region->chromosome, region->start_position, region->end_position, 
                        sqlite3_errmsg(db), sqlite3_errcode(db));
            return rc;
        }
        
        sqlite3_bind_int(stmt2, 1, affected_rows);
        sqlite3_bind_text(stmt2, 2, region->chromosome, strlen(region->chromosome), SQLITE_STATIC);
        sqlite3_bind_int64(stmt2, 3, chunk_start);
        sqlite3_bind_int64(stmt2, 4, chunk_end);
    
        if ((rc = sqlite3_step(stmt2)) != SQLITE_DONE) {
            LOG_ERROR_F("Could not remove region %s:%ld-%ld: %s (%d)\n", 
                        region->chromosome, region->start_position, region->end_position, 
                        sqlite3_errmsg(db), sqlite3_errcode(db));
            return rc;
        }
        
        sqlite3_finalize(stmt2);
        
    } else {
        LOG_ERROR_F("Could not remove region %s:%ld-%ld: %s (%d)\n", 
                    region->chromosome, region->start_position, region->end_position, 
                    sqlite3_errmsg(db), sqlite3_errcode(db));
        return rc;
    }
    
    char *sql_end = "END TRANSACTION";
    rc = exec_sql(sql_end, db);
    if (rc != SQLITE_OK) {
        LOG_ERROR_F("Could not remove regions: %s (%d)\n", 
                    sqlite3_errmsg(db), sqlite3_errcode(db));
        return rc;
    }

    return rc;
}


/* ******************************
 *  Chromosomes (region trees)	*
 * ******************************/

array_list_t *get_chromosome(const char *key, region_table_t *table) {
    if (!table->is_ready) { // Don't allow queries over an un-indexed DB
        finish_region_table_loading(table);
    }

    array_list_t *regions_in_chromosome = array_list_new(100, 1.5, COLLECTION_MODE_ASYNCHRONIZED);
    
    sqlite3_stmt *query_stmt = table->get_chromosome_stmt;
    sqlite3_bind_text(query_stmt, 1, key, strlen(key), SQLITE_STATIC);
    
    if (sqlite3_step(query_stmt) == SQLITE_ROW) {
        do {
	  char *chromosome = strdup((char *)sqlite3_column_text(query_stmt, 0));
            size_t start = sqlite3_column_int64(query_stmt, 1);
            size_t end = sqlite3_column_int64(query_stmt, 2);
            const char *strand = (char *)sqlite3_column_text(query_stmt, 3);
            const char *type = (char *)sqlite3_column_text(query_stmt, 4);
            
            region_t *region = region_new(strdup(chromosome), start, end,
                                          strand ? strdup(strand) : NULL, 
                                          type ? strdup(type) : NULL);
            array_list_insert(region, regions_in_chromosome);
        } while (sqlite3_step(query_stmt) == SQLITE_ROW);
    } else {
        // Retry once
        sqlite3_reset(query_stmt);
        sqlite3_bind_text(query_stmt, 1, key, strlen(key), SQLITE_STATIC);
        if (sqlite3_step(query_stmt) == SQLITE_ROW) {
            do {
                const char *chromosome = (char *)sqlite3_column_text(query_stmt, 0);
                size_t start = sqlite3_column_int64(query_stmt, 1);
                size_t end = sqlite3_column_int64(query_stmt, 2);
                const char *strand = (char *)sqlite3_column_text(query_stmt, 3);
                const char *type = (char *)sqlite3_column_text(query_stmt, 4);
                
                region_t *region = region_new(strdup(chromosome), start, end,
                                              strand ? strdup(strand) : NULL, 
                                              type ? strdup(type) : NULL);
                array_list_insert(region, regions_in_chromosome);
            } while (sqlite3_step(query_stmt) == SQLITE_ROW);
        } else {
            sqlite3* db = table->storage;
            LOG_ERROR_F("Regions table failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
        }
    }
    
    sqlite3_reset(query_stmt);
    
    return regions_in_chromosome;
}

int count_regions_in_chromosome(const char *key, region_table_t *table) {
    if (!table->is_ready) { // Don't allow queries over an un-indexed DB
        finish_region_table_loading(table);
    }
    
    int count = -1;
    
    sqlite3_stmt *query_stmt = table->count_in_chromosome_stmt;
    sqlite3_bind_text(query_stmt, 1, key, strlen(key), SQLITE_STATIC);
    
    if (sqlite3_step(query_stmt) == SQLITE_ROW) {
        count = sqlite3_column_int(query_stmt, 0);
    } else {
        // Retry once
        sqlite3_reset(query_stmt);
        sqlite3_bind_text(query_stmt, 1, key, strlen(key), SQLITE_STATIC);
        if (sqlite3_step(query_stmt) == SQLITE_ROW) {
            count = sqlite3_column_int(query_stmt, 0);
        } else {
            sqlite3* db = table->storage;
            LOG_ERROR_F("Regions table failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
        }
    }
    
    sqlite3_reset(query_stmt);
    
    return count;
}

