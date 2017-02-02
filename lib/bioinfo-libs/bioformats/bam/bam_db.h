#ifndef BAM_DB_H
#define BAM_DB_H

/*
 * bam_db.h
 *
 *  Created on: Apr 23, 2013
 *      Author: jtarraga
 */

//------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>

#include "commons/log.h"
#include "commons/sqlite/sqlite3.h"
#include "containers/array_list.h"

//------------------------------------------------------------------------

#define BAM_CHUNKSIZE 10000

//------------------------------------------------------------------------

typedef struct bam_query_fields {
  int chr_length;
  int strand;
  int start;
  int end;
  int flag;
  int mapping_quality;
  int num_errors;
  int num_indels;
  int indels_length;
  int template_length;
  char *chr;
  char *id;
} bam_query_fields_t;

//------------------------------------------------------------------------

bam_query_fields_t *bam_query_fields_new(char *id, char *chr, int chr_length,
					 int strand, int start, int end, 
					 int flag, int mapping_quality,
					 int num_errors, int num_indels,
					 int indels_length, int template_length);

void bam_query_fields_free(bam_query_fields_t *p);

void print_bam_query_fields(bam_query_fields_t *p);

//------------------------------------------------------------------------
// 
//------------------------------------------------------------------------

int create_bam_query_fields(sqlite3 *db);
int create_bam_index(sqlite3 *db);

//------------------------------------------------------------------------

int insert_bam_query_fields(void *custom_fields, sqlite3 *db);

int insert_bam_query_fields_list(array_list_t *list, sqlite3 *db);

int prepare_statement_bam_query_fields(sqlite3 *db, sqlite3_stmt **stmt);
int insert_statement_bam_query_fields(void *custom_fields, 
				      sqlite3_stmt *stmt, sqlite3 *db);


//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif // end of BAM_DB_H
