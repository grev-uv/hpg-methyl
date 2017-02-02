#ifndef DB_UTILS_H
#define DB_UTILS_H

/*
 * db_utils.h
 *
 *  Created on: Apr 22, 2013
 *      Author: jtarraga
 */

//------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>

#include <commons/log.h>
#include <commons/sqlite/sqlite3.h>
#include <containers/khash.h>

//------------------------------------------------------------------------

KHASH_MAP_INIT_STR(stats_chunks, int)

//------------------------------------------------------------------------

int create_stats_db(const char *db_name, int chunksize, 
		    int (*record_stats_creation)(sqlite3 *), sqlite3** db);

int create_stats_index(int (*create_custom_index)(sqlite3 *), sqlite3* db);

//------------------------------------------------------------------------

int insert_global_stats(const char *id, const char *title, 
			const char *value, sqlite3 *db);

int prepare_statement_global_stats(sqlite3 *db, sqlite3_stmt **stmt);

int finalize_statement_global_stats(sqlite3_stmt *stmt);

int insert_statement_global_stats(const char *id, const char *title, 
				  const char *value, sqlite3_stmt *stmt, 
				  sqlite3 *db);

//------------------------------------------------------------------------

int insert_chunk(const char *chromosome, int chunk_id, int start, int end, 
		 int features_count, sqlite3 *db);

int inc_chunk(const char *chr, int chunk_id, int chunk_start, int chunk_end, 
	      sqlite3 *db);

int update_chunks_hash(const char *chr, int chr_length, int chunksize, 
		       int start, int end, khash_t(stats_chunks) *hash);

int insert_chunk_hash(int chunksize, khash_t(stats_chunks) *hash, sqlite3 *db);

//------------------------------------------------------------------------

int insert_record_query_fields(const char *chr, int chr_length, int chunksize, 
			       int start, int end, void *fields,
			       int (*insert_custom_fields)(void *, sqlite3 *), 
			       sqlite3* db);

int insert_statement_record_query_fields(const char *chr, int chr_length, int chunksize, 
					 int start, int end, void *fields,
					 int (*insert_statement_custom_fields)(void *, sqlite3_stmt *, sqlite3 *), 
					 sqlite3_stmt *custom_stmt,
					 sqlite3* db);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif // end of DB_UTILS_H
