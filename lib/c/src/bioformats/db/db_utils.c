/*
 * db_utils.c
 *
 *  Created on: Apr 22, 2013
 *      Author: jtarraga
 */

#include "db_utils.h"

//------------------------------------------------------------------------

inline int exec_sql(char *sql, sqlite3* db) {
  int rc;
  char *error_msg;
  if ((rc = sqlite3_exec(db, sql, NULL, NULL, &error_msg))) {
    LOG_ERROR_F("Stats database failed (%s): %s\n", sql, error_msg);
  }
  return rc;
}

//------------------------------------------------------------------------
//             S T A T I S T I C S        D A T A B A S E
//------------------------------------------------------------------------

//-------------------------- CREATION FUNCTIONS --------------------------

int create_stats_db(const char *db_name, int chunksize, 
		    int (*create_custom_fields)(sqlite3 *), sqlite3** db) {

  // create sqlite db
  if (sqlite3_open(db_name, db)) {
    LOG_FATAL_F("Could not open stats database (%s): %s\n", 
		db_name, sqlite3_errmsg(*db));
  }

  int rc;
  char sql[128];
  
  sprintf(sql, "BEGIN TRANSACTION");
  rc = exec_sql(sql, *db);

  // create global stats table and index, and insert the chunksize
  sprintf(sql, "CREATE TABLE IF NOT EXISTS global_stats (name TEXT PRIMARY KEY, title TEXT, value TEXT)");
  rc = exec_sql(sql, *db);


  sprintf(sql, "%i", chunksize);
  rc = insert_global_stats("CHUNK_SIZE", "Chunk size", sql, *db);
  rc = insert_global_stats("CHR_PREFIX", "Chromosome prefix", "", *db);

  // create chunks table
  sprintf(sql, "CREATE TABLE IF NOT EXISTS chunk (chromosome TEXT, chunk_id INT, start INT, end INT, features_count INT)");
  rc = exec_sql(sql, *db);

  // create record_query_fields table for bam, vcf.. files
  if (create_custom_fields) {
    rc = create_custom_fields(*db);
  }
  
  sprintf(sql, "END TRANSACTION");
  rc = exec_sql(sql, *db);

  return rc;
}

int create_stats_index(int (*create_custom_index)(sqlite3 *), sqlite3* db) {
  int rc;
  char sql[128];

  // create chunks index
  sprintf(sql, "CREATE INDEX chunk_chromosome_chunk_id_idx ON chunk (chromosome, chunk_id)");
  rc = exec_sql(sql, db);

  // create custom index (for the record_query_fields table)
  if (create_custom_index) {
    rc = create_custom_index(db);
  }

  return rc;
}

int close_stats_db(sqlite3* db, khash_t(stats_chunks) *hash) {
    // Close database connection
    sqlite3_close_v2(db);
    // Destroy the chunks hashtable
    kh_destroy(stats_chunks, hash);

    return 0;

}


//-------------------------- GLOBAL STATS TABLE --------------------------

int insert_global_stats(const char *name, const char *title, 
			const char *value, sqlite3 *db) {
  int rc;

  // IMPORTANT !!!
  //
  // The sqlite3_exec() interface is a convenience wrapper around 
  // sqlite3_prepare_v2(), sqlite3_step(), and sqlite3_finalize(), 
  // that allows an application to run multiple statements of SQL 
  // without having to use a lot of C code.


  char sql[strlen(name) + strlen(title) + strlen(value) + 100];
  sprintf(sql, "INSERT INTO global_stats VALUES('%s', '%s', '%s')", 
	  name, title, value);
  rc = exec_sql(sql, db);

  return rc;

}

//------------------------------------------------------------------------

int prepare_statement_global_stats(sqlite3 *db, sqlite3_stmt **stmt) {
  char sql[] = "INSERT INTO global_stats VALUES(?1, ?2, ?3)"; 
  return sqlite3_prepare_v2(db, sql, strlen(sql), stmt, NULL);
}

//------------------------------------------------------------------------

int finalize_statement_global_stats(sqlite3_stmt *stmt) {
  return sqlite3_finalize(stmt);
}

//------------------------------------------------------------------------

int insert_statement_global_stats(const char *name, const char *title, 
				  const char *value, sqlite3_stmt *stmt, 
				  sqlite3 *db) {  
  int rc;

  sqlite3_bind_text(stmt, 1, name, strlen(name), SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt, 2, title, strlen(title), SQLITE_TRANSIENT);
  sqlite3_bind_text(stmt, 3, value, strlen(value), SQLITE_TRANSIENT);

  if ((rc = sqlite3_step(stmt) != SQLITE_DONE)) {
    LOG_ERROR_F("Stats databases failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
  } else {
    rc = SQLITE_OK;
  }
  sqlite3_reset(stmt);

  return rc;
}



//------------------------------------------------------------------------
//                  R E G I O N S    D A T A B A S E
//------------------------------------------------------------------------

int create_regions_db(const char *db_name, int chunksize, sqlite3** db) {
  // create sqlite db
  if (sqlite3_open(db_name, db)) {
    LOG_FATAL_F("Could not open regions database (%s): %s\n", 
		db_name, sqlite3_errmsg(*db));
  }

  int rc;
  char sql[128];
  
  sprintf(sql, "BEGIN TRANSACTION");
  rc = exec_sql(sql, *db);

  // create regions table
  sprintf(sql, "CREATE TABLE IF NOT EXISTS regions (chromosome TEXT, start INT, end INT, strand CHARACTER(1), type TEXT)");
  rc = exec_sql(sql, *db);

  // create chunks table
  sprintf(sql, "CREATE TABLE IF NOT EXISTS chunk (chromosome TEXT, chunk_id INT, start INT, end INT, features_count INT)");
  rc = exec_sql(sql, *db);

  sprintf(sql, "END TRANSACTION");
  rc = exec_sql(sql, *db);

  return rc;
}

//------------------------------------------------------------------------

int create_regions_index(sqlite3* db) {
  int rc;
  char sql[128];

  // create chunks index
  sprintf(sql, "CREATE INDEX chunk_chromosome_chunk_id_idx ON chunk (chromosome, chunk_id)");
  rc = exec_sql(sql, db);

  // create regions index
  sprintf(sql, "CREATE INDEX regions_chromosome_start_end_idx ON chunk (chromosome, start, end)");
  rc = exec_sql(sql, db);
  return rc;
}

//------------------------------------------------------------------------

static int prepare_statement_regions_query_fields(sqlite3 *db, sqlite3_stmt **stmt) {
  char sql[] = "INSERT INTO regions VALUES (?1, ?2, ?3, ?4, ?5)";
  return sqlite3_prepare_v2(db, sql, strlen(sql), stmt, NULL);
}

//------------------------------------------------------------------------

int insert_regions_query_fields_list(array_list_t *list, sqlite3 *db) {

    int rc;
    sqlite3_stmt *stmt;
    region_t *fields;
    char *error_message;

    prepare_statement_regions_query_fields(db, &stmt);

    if ((rc = sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &error_message))) {
        LOG_ERROR_F("Regions database failed: %s (%d)\n", rc, error_message);
    }

    int num_items = array_list_size(list);
    for (int i = 0; i < num_items; i++) {
        fields = array_list_get(i, list);

        sqlite3_bind_text(stmt, 1, fields->chromosome, strlen(fields->chromosome), SQLITE_STATIC);
        sqlite3_bind_int64(stmt, 2, fields->start_position);
        sqlite3_bind_int64(stmt, 3, fields->end_position);
        sqlite3_bind_text(stmt, 4, fields->strand, strlen(fields->strand), SQLITE_STATIC);
        sqlite3_bind_text(stmt, 5, fields->type, strlen(fields->type), SQLITE_STATIC);

        if ((rc = sqlite3_step(stmt) != SQLITE_DONE)) {
            LOG_ERROR_F("Regions database failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
        }

        sqlite3_reset(stmt);
    }

    if ((rc = sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &error_message))) {
        LOG_ERROR_F("Regions database failed: %s (%d)\n", rc, error_message);
    }

    sqlite3_finalize(stmt);

    return 0;

}



//------------------------------------------------------------------------
//                     C H U N K       T A B L E
//------------------------------------------------------------------------

int insert_chunk(const char *chr, int chunk_id, int start, int end, 
		 int features_count, sqlite3 *db) {

  int rc;
  sqlite3_stmt* stmt;
  char sql[] = "INSERT INTO chunk VALUES(?1, ?2, ?3, ?4, ?5)";

  sqlite3_prepare_v2(db, sql, strlen(sql), &stmt, NULL);
  sqlite3_bind_text(stmt, 1, chr, strlen(chr), SQLITE_STATIC);
  sqlite3_bind_int(stmt, 2, chunk_id);
  sqlite3_bind_int(stmt, 3, start);
  sqlite3_bind_int(stmt, 4, end);
  sqlite3_bind_int(stmt, 5, features_count);

  if ((rc = sqlite3_step(stmt) != SQLITE_DONE)) {
    LOG_ERROR_F("Database failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
  } else {
    rc = SQLITE_OK;
  }
  sqlite3_finalize(stmt);

  return rc;
}

//------------------------------------------------------------------------

int inc_chunk(const char *chr, int chunk_id, int chunk_start, int chunk_end, sqlite3 *db) {
  int rc;
  sqlite3_stmt *stmt;
  char sql[strlen(chr) + 200];

  strcpy(sql, "SELECT features_count FROM chunk WHERE chromosome = ?1 AND chunk_id = ?2");

  sqlite3_prepare_v2(db, sql, strlen(sql), &stmt, NULL);
  sqlite3_bind_text(stmt, 1, chr, strlen(chr), SQLITE_STATIC);
  sqlite3_bind_int(stmt, 2, chunk_id);
  
  rc = sqlite3_step(stmt);

  if (rc == SQLITE_ROW) {
    sprintf(sql, "UPDATE chunk SET features_count = features_count + 1 WHERE chromosome = '%s' AND chunk_id = %i", 
	    chr, chunk_id);
    rc = exec_sql(sql, db);
  } else if (rc == SQLITE_DONE) {
    sprintf(sql, "INSERT INTO chunk VALUES('%s', %i, %i, %i, 1)", chr, chunk_id,
	    chunk_start, chunk_end);
    rc = exec_sql(sql, db);
  }

  sqlite3_finalize(stmt);

  return 0;

}

//------------------------------------------------------------------------

int update_chunks_hash(const char *chr, int chr_length, int chunksize, 
		       int start, int end, khash_t(stats_chunks) *hash) {
  
  // update features counter for "touched" chunks
  int ret;
  khiter_t k;
  char *key;

  int chunk_start, chunk_end;
  int chunk_id_start = start / chunksize;
  int chunk_id_end = end / chunksize;
  
  for (int chunk_id = chunk_id_start; chunk_id <= chunk_id_end; chunk_id++) {
    chunk_start = chunk_id * chunksize;
    chunk_end = chunk_start + chunksize - 1;
    if (chunk_end > chr_length) {
      chunk_end = chr_length - 1;
    }

    key = calloc(strlen(chr) + 32,  sizeof(char));
    sprintf(key, "%i::%s", chunk_id, chr);

    k = kh_put(stats_chunks, hash, key, &ret);
    if (ret == 0) {
      kh_value(hash, k) = (kh_value(hash, k) + 1);
      free(key);
    } else if (ret == 1) {
      kh_value(hash, k) = 1;
    }

    //    LOG_DEBUG_F("key = %s, value = %i\n", key, kh_value(hash, k));
  }
  return 0;
}

//------------------------------------------------------------------------

int insert_chunk_hash(int chunksize, khash_t(stats_chunks) *hash, sqlite3 *db) {
  int rc;
  sqlite3_stmt *stmt;
  char *errorMessage;

  char chr[1024];
  char *key;
  int chunk_id, chunk_start, features_count, counter = 0;

  // prepared statement
  char sql[] = "INSERT INTO chunk VALUES (?1, ?2, ?3, ?4, ?5)";
  sqlite3_prepare_v2(db, sql, strlen(sql), &stmt, NULL);

  if ((rc = sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage))) {
    LOG_ERROR_F("Database failed: %s (%d)\n", errorMessage, rc);
  }

  for (khiter_t k = kh_begin(hash); k != kh_end(hash); ++k) {

    if (kh_exist(hash, k)) {
      counter++;
           
      key = (char *)kh_key(hash, k);
      features_count = kh_value(hash, k);
      //      LOG_DEBUG_F("key %s -> value = %i\n", key, features_count);
      sscanf(key, "%i::%s", &chunk_id, chr);
      //      LOG_DEBUG_F("\tchunk_id = %i, chr = %s\n", chunk_id, chr);

      chunk_start = chunk_id * chunksize;

      sqlite3_bind_text(stmt, 1, chr, strlen(chr), SQLITE_TRANSIENT);
      sqlite3_bind_int(stmt, 2, chunk_id);
      sqlite3_bind_int(stmt, 3, chunk_start);
      sqlite3_bind_int(stmt, 4, chunk_start + chunksize - 1);
      sqlite3_bind_int(stmt, 5, features_count);

      if ((rc = sqlite3_step(stmt) != SQLITE_DONE)) {
	LOG_ERROR_F("Database failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
      }

      sqlite3_reset(stmt);

      // we must free the key, it was allocated previously 
      free(key);
      kh_del(stats_chunks, hash, k);


      if (counter % 100000 == 0) {
	// commit the current transaction
	if ((rc = sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage))) {
	  LOG_ERROR_F("Database failed: %s (%d)\n", errorMessage, rc);
	}

	// start a new transaction
	if ((rc = sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage))) {
	  LOG_ERROR_F("Database failed: %s (%d)\n", errorMessage, rc);
	}
	counter = 0;
      }
    }
  }

  if (counter) {
    if ((rc = sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage))) {
      LOG_ERROR_F("Database failed: %s (%d)\n", errorMessage, rc);
    }
  }

  sqlite3_finalize(stmt);
  
  return rc;
}
