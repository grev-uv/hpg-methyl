/*
 * bam_db.c
 *
 *  Created on: Apr 23, 2013
 *      Author: jtarraga
 */

#include "bam_db.h"

//------------------------------------------------------------------------

bam_query_fields_t *bam_query_fields_new(char *id, char *chr, int chr_length,
					 int strand, int start, int end, 
					 int flag, int mapping_quality,
					 int num_errors, int num_indels,
					 int indels_length, int template_length) {

  bam_query_fields_t *p = calloc(1, sizeof(bam_query_fields_t));

  p->id = strdup(id);
  p->chr = strdup(chr);
  p->chr_length = chr_length;
  p->strand = strand;
  p->start = start;
  p->end = end;
  p->flag = flag;
  p->mapping_quality = mapping_quality;
  p->num_errors = num_errors;
  p->num_indels = num_indels;
  p->indels_length = indels_length;
  p->template_length = template_length;

  return p;
}

void bam_query_fields_free(bam_query_fields_t *p) {
  if (p) {
    if (p->chr) free(p->chr);
    if (p->id) free(p->id);
    free(p);
  }
}

void print_bam_query_fields(bam_query_fields_t *p) {
  printf("id = %s, chr = %s, chr_length = %i, strand = %i, start = %i, end = %i, flag = %i, mquality = %u, num_erros = %i, num_indels = %i, indels_length = %i, template_length = %i\n", 
          p->id, p->chr, p->chr_length, p->strand, p->start, p->end, p->flag, p->mapping_quality, p->num_errors, p->num_indels, p->indels_length, p->template_length);
}

//------------------------------------------------------------------------
//
//------------------------------------------------------------------------

int create_bam_query_fields(sqlite3 *db) {

  // create record_stats table for bam files
  int rc;
  char *error_msg;
  char sql[1000];
  sprintf(sql, "CREATE TABLE record_query_fields (chromosome TEXT, strand CHAR, start INT, end INT, flag CHAR, mapping_quality CHAR, num_errors INT, num_indels INT, indels_length INT, template_length INT, id TEXT)");

  if (rc = sqlite3_exec(db, sql, NULL, NULL, &error_msg)) {
    LOG_FATAL_F("Stats database failed: %s\n", error_msg);
  }
  return 0;
}

//------------------------------------------------------------------------

int create_bam_index(sqlite3 *db) {
  int rc;
  char sql[128];
  char *error_msg;

  // create chunks index
  sprintf(sql, "CREATE INDEX record_query_fields_chromosome_start_end_idx ON record_query_fields (chromosome, start, end)");
  if (rc = sqlite3_exec(db, sql, NULL, NULL, &error_msg)) {
    LOG_DEBUG_F("Stats database failed creating BAM index: %s\n", error_msg);
  }
  return rc;
}

//------------------------------------------------------------------------

int insert_bam_query_fields(void *custom_fields, sqlite3 *db) {
  bam_query_fields_t *fields = (bam_query_fields_t *) custom_fields;

  //  print_bam_query_fields(custom_fields);

  // insert into record_query_fields
  int rc;
  char sql[500 + strlen(fields->id) + strlen(fields->chr)];
  sprintf(sql, "INSERT INTO record_query_fields VALUES('%s', %i, %i, %i, %i, %i, %i, %i, %i, %i, '%s')", 
	  fields->chr, fields->strand, fields->start, fields->end, fields->flag, fields->mapping_quality, fields->num_errors, fields->num_indels, fields->indels_length, fields->template_length, fields->id);

  char *error_msg;
  if (rc = sqlite3_exec(db, sql, NULL, NULL, &error_msg)) {
    LOG_DEBUG_F("Stats database failed: %s\n", error_msg);
  }

  return rc;
}

//------------------------------------------------------------------------

int insert_bam_query_fields_list(array_list_t *list, sqlite3 *db) {

  int rc;
  sqlite3_stmt *stmt;
  bam_query_fields_t *fields;
  char *errorMessage;

  prepare_statement_bam_query_fields(db, &stmt);

  if (rc = sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage)) {
    LOG_DEBUG_F("Stats databases failed: %s (%d)\n", rc, errorMessage);
  }

  int num_items = array_list_size(list);
  for (int i = 0; i < num_items; i++) {
    fields = array_list_get(i, list);

    sqlite3_bind_text(stmt, 1, fields->chr, strlen(fields->chr), SQLITE_STATIC);
    sqlite3_bind_int(stmt, 2, fields->strand);
    sqlite3_bind_int(stmt, 3, fields->start);
    sqlite3_bind_int(stmt, 4, fields->end);
    sqlite3_bind_int(stmt, 5, fields->flag);
    sqlite3_bind_int(stmt, 6, fields->mapping_quality);
    sqlite3_bind_int(stmt, 7, fields->num_errors);
    sqlite3_bind_int(stmt, 8, fields->num_indels);
    sqlite3_bind_int(stmt, 9, fields->indels_length);
    sqlite3_bind_int(stmt, 10, fields->template_length);
    sqlite3_bind_text(stmt, 11, fields->id, strlen(fields->id), SQLITE_STATIC);

    if (rc = sqlite3_step(stmt) != SQLITE_DONE) {
      LOG_DEBUG_F("Stats databases failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
    }

    sqlite3_reset(stmt);
  }

  if (rc = sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage)) {
    LOG_DEBUG_F("Stats databases failed: %s (%d)\n", rc, errorMessage);
  }

  sqlite3_finalize(stmt);
}

//------------------------------------------------------------------------

int prepare_statement_bam_query_fields(sqlite3 *db, sqlite3_stmt **stmt) {
  //char sql[] = "INSERT INTO record_query_fields (chromosome) VALUES (?1)";
  char sql[] = "INSERT INTO record_query_fields VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11)";
  return sqlite3_prepare_v2(db, sql, strlen(sql) + 200, stmt, NULL);
}

//------------------------------------------------------------------------

int insert_statement_bam_query_fields(void *custom_fields, 
				      sqlite3_stmt *stmt, sqlite3 *db) {  
  int rc;
  bam_query_fields_t *fields = (bam_query_fields_t *) custom_fields;

  sqlite3_bind_text(stmt, 1, fields->chr, strlen(fields->chr), SQLITE_STATIC);
  sqlite3_bind_int(stmt, 2, fields->strand);
  sqlite3_bind_int(stmt, 3, fields->start);
  sqlite3_bind_int(stmt, 4, fields->end);
  sqlite3_bind_int(stmt, 5, fields->flag);
  sqlite3_bind_int(stmt, 6, fields->mapping_quality);
  sqlite3_bind_int(stmt, 7, fields->num_errors);
  sqlite3_bind_int(stmt, 8, fields->num_indels);
  sqlite3_bind_int(stmt, 9, fields->indels_length);
  sqlite3_bind_int(stmt, 10, fields->template_length);
  sqlite3_bind_text(stmt, 11, fields->id, strlen(fields->id), SQLITE_STATIC);

  if (rc = sqlite3_step(stmt) != SQLITE_DONE) {
    LOG_DEBUG_F("Stats databases failed: %s (%d)\n", sqlite3_errmsg(db), sqlite3_errcode(db));
  }

  sqlite3_reset(stmt);

  return rc;
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------

