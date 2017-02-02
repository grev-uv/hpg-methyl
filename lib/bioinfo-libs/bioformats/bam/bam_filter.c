/*
 * bam_filter.c
 *
 *  Created on: Feb 25, 2013
 *      Author: jtarraga
 */

#include "bam_filter.h"

//------------------------------------------------------------------------

bam_filter_input_t *new_bam_filter_input(char *in_filename, char *out_dirname,
					 int by_mapped, int by_unmapped, 
					 int by_proper_pairs, int by_unique,
					 int by_num_errors, int min_num_errors, int max_num_errors,
					 int by_quality, int min_quality, int max_quality,
					 int by_length, int min_length, int max_length,
					 region_table_t *region_table,
					 int num_threads,int batch_size) {
  bam_filter_input_t *input = (bam_filter_input_t *) calloc(1, sizeof(bam_filter_input_t));
  
  input->in_filename = strdup(in_filename);
  input->out_dirname = strdup(out_dirname);
  input->by_mapped = by_mapped;
  input->by_unmapped = by_unmapped;
  input->by_proper_pairs = by_proper_pairs;
  input->by_unique = by_unique;
  input->by_num_errors = by_num_errors;
  input->min_num_errors = min_num_errors;
  input->max_num_errors = max_num_errors;
  input->by_quality = by_quality;
  input->min_quality = min_quality;
  input->max_quality = max_quality;
  input->by_length = by_length;
  input->min_length = min_length;
  input->max_length = max_length;

  input->region_table = region_table;
  input->num_threads = num_threads;
  input->batch_size = batch_size;

  return input;
}

void free_bam_filter_input(bam_filter_input_t *input) {
  if (input) {

    if (input->in_filename) free(input->in_filename);
    if (input->out_dirname) free(input->out_dirname);

    free(input);
  }
}

//====================================================================
// W O R K F L O W     F O R      S T A T I S T I C S
//====================================================================

#define CONSUMER_STAGE   -1

//--------------------------------------------------------------------
// structure between the different workflow stages
//--------------------------------------------------------------------

typedef struct bam_filter_wf_batch {
  bam_file_t *out_file;
  bam_filter_input_t *in_filter;
  array_list_t *bam1_list;
  char **sequence_labels;
  sqlite3 *db;
} bam_filter_wf_batch_t;

bam_filter_wf_batch_t *new_bam_filter_wf_batch(bam_file_t *out_file,
					       bam_filter_input_t *in_filter,
					       char **sequence_labels,
					       array_list_t *bam1_list,
					       sqlite3 *db) {
  
  bam_filter_wf_batch_t *b = (bam_filter_wf_batch_t *) calloc(1, sizeof(bam_filter_wf_batch_t));
  
  b->out_file = out_file;
  b->in_filter = in_filter;
  b->bam1_list = bam1_list;
  b->sequence_labels = sequence_labels;
  b->db = db;
  
  return b;
}

void free_bam_filter_wf_batch(bam_filter_wf_batch_t *b) {
  if (b) free(b);
}

//--------------------------------------------------------------------
// workflow input
//--------------------------------------------------------------------

typedef struct bam_filter_wf_input {
  bam_file_t *in_file;
  bam_file_t *out_file;
  char **sequence_labels;
  bam_filter_input_t *in_filter;
  sqlite3 *db;
} bam_filter_wf_input_t;

bam_filter_wf_input_t *new_bam_filter_wf_input(bam_file_t *in_file,
					       bam_file_t *out_file,
					       char **sequence_labels,
					       bam_filter_input_t *in_filter,
					       sqlite3 *db) {
  
  bam_filter_wf_input_t *wfi = (bam_filter_wf_input_t *) calloc(1, sizeof(bam_filter_wf_input_t));
  
  wfi->in_file = in_file;
  wfi->out_file = out_file;
  wfi->sequence_labels = sequence_labels;
  wfi->in_filter = in_filter;
  wfi->db = db;

  return wfi;
}

void free_bam_filter_wf_input(bam_filter_wf_input_t *wfi) {
  if (wfi) free(wfi);
}

//--------------------------------------------------------------------
// workflow producer
//--------------------------------------------------------------------

void *bam_filter_producer(void *input) {
  
  bam_filter_wf_input_t *wf_input = (bam_filter_wf_input_t *) input;
  bam_filter_wf_batch_t *new_batch = NULL;
  int max_num_bam1s = wf_input->in_filter->batch_size;

  bam1_t *bam1;
  array_list_t *bam1_list = array_list_new(max_num_bam1s, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);  
  for (int i = 0; i < max_num_bam1s; i++) {
    bam1 = bam_init1();    
    if (bam_read1(wf_input->in_file->bam_fd, bam1) > 0) {
      array_list_insert(bam1, bam1_list);
    } else {
      bam_destroy1(bam1);
      break;
    }
  }
  
  size_t num_items = array_list_size(bam1_list);
  
  if (num_items == 0) {
    array_list_free(bam1_list, NULL);
  } else {
    new_batch = new_bam_filter_wf_batch(wf_input->out_file,
					wf_input->in_filter,
					wf_input->sequence_labels,
					bam1_list,
					wf_input->db);				      
  }
    
  return new_batch;
}

//--------------------------------------------------------------------
// workflow consumer
//--------------------------------------------------------------------

int bam_filter_consumer(void *data) {
  bam_filter_wf_batch_t *batch = (bam_filter_wf_batch_t *) data;
  
  bam_file_t *out_file = batch->out_file;

  bam1_t *bam1;
  array_list_t *bam1_list = batch->bam1_list;
  size_t num_items = array_list_size(bam1_list);

  for (int i = 0; i < num_items; i++) {
    bam1 = array_list_get(i, bam1_list);
    //    printf("write bam1\n");
    bam_write1(out_file->bam_fd, bam1);
    bam_destroy1(bam1);
  }

  array_list_free(bam1_list, NULL);

  // free memory
  free_bam_filter_wf_batch(batch);
}

//--------------------------------------------------------------------
// workflow worker
//--------------------------------------------------------------------
void filtered_bam(bam1_t *bam1, sqlite3 *db);

int bam_filter_worker(void *data) {
  bam_filter_wf_batch_t *batch = (bam_filter_wf_batch_t *) data;

  array_list_t *bam1_list = batch->bam1_list;
  bam_filter_input_t *in_filter = batch->in_filter;
  region_table_t * region_table = in_filter->region_table;

  region_t region;
  char **sequence_labels = batch->sequence_labels;
  
  bam1_t *bam1;
  int bam_seq_len, num_cigar_ops, quality, length;
  uint8_t* bam_seq;
  uint32_t bam_flag, cigar_int, *cigar;
  int isize, strand, num_errors;
  size_t num_items = array_list_size(bam1_list);

  int match;

  array_list_t *out_bam1_list = array_list_new(num_items, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);  

  for (int i = 0; i < num_items; i++) {
    
    bam1 = array_list_get(i, bam1_list);
    bam_flag = (uint32_t) bam1->core.flag;

    // by unmapped
    if (in_filter->by_unmapped) {
      if (bam_flag & BAM_FUNMAP) {
	array_list_insert(bam1, out_bam1_list);
      } else {
	filtered_bam(bam1, batch->db);
	//	bam_destroy1(bam1);
      }
      continue;
    }

    if (region_table) {
      region.chromosome = sequence_labels[bam1->core.tid];
      region.start_position = bam1->core.pos;
      region.end_position = region.start_position + bam1->core.l_qseq;
    }

    if ((region_table == NULL) || find_region(&region, region_table)) {
      /*
      // by mapped
      if (in_filter->by_mapped && (bam_flag & BAM_FUNMAP)) {
	bam_destroy1(bam1);
	continue;
      }
      */
      // by proper pair
      if (in_filter->by_proper_pairs && (bam_flag & BAM_FPAIRED)) {
	if ( !(bam_flag & BAM_FPROPER_PAIR)) {
	  filtered_bam(bam1, batch->db);
	  //	bam_destroy1(bam1);
	  continue;
	}
      }

      // by unique
      if (in_filter->by_unique && (bam_flag & BAM_FSECONDARY)) {
	filtered_bam(bam1, batch->db);
	//	bam_destroy1(bam1);
	continue;
      }

      // by num. errors
      if (in_filter->by_num_errors) {
	num_errors = bam_aux2i(bam_aux_get(bam1, "NM"));
	if ( ((in_filter->max_num_errors != -1) && (num_errors > in_filter->max_num_errors)) ||
	     ((in_filter->min_num_errors != -1) && (num_errors < in_filter->min_num_errors))   ) {
	  filtered_bam(bam1, batch->db);
	  //	  bam_destroy1(bam1);
	  continue;
	}
      }

      // by quality
      if (in_filter->by_quality) {
	quality = bam1->core.qual;
	if ( ((in_filter->max_quality != -1) && (quality > in_filter->max_quality)) ||
	     ((in_filter->min_quality != -1) && (quality < in_filter->min_quality))   ) {
	  filtered_bam(bam1, batch->db);
	  //	  bam_destroy1(bam1);
	  continue;
	}
      }

      // by length
      if (in_filter->by_length) {
	length = bam1->core.l_qseq;
	if ( ((in_filter->max_length != -1) && (length > in_filter->max_length)) ||
	     ((in_filter->min_length != -1) && (length < in_filter->min_length))   ) {
	  filtered_bam(bam1, batch->db);
	  //	  bam_destroy1(bam1);
	  continue;
	}
      }

      // finally, this bam1 passed all the filters,
      // insert it in the output list
      array_list_insert(bam1, out_bam1_list);
    } else {
      filtered_bam(bam1, batch->db);
      //      bam_destroy1(bam1);
    }
  }

  // free the input list and update the output list
  array_list_free(bam1_list, NULL);
  batch->bam1_list = out_bam1_list;
  
  return CONSUMER_STAGE;
}

//--------------------------------------------------------------------

inline void filtered_bam(bam1_t *bam1, sqlite3 *db) {

  if (bam1->core.flag & BAM_FPAIRED) {
    char *qname = bam1_qname(bam1);
    char sql[strlen(qname) + 100];
    sprintf(sql, "INSERT INTO ids VALUES('%s')", qname);
    printf("insert  in db: %s\n", sql);
    sqlite3_exec(db, sql, 0, 0, 0);
  }
  bam_destroy1(bam1);
}

//--------------------------------------------------------------------
// workflow description
//--------------------------------------------------------------------

void filter_bam(bam_filter_input_t *input) {
  bam_file_t *in_file = bam_fopen(input->in_filename);
  
  char path[strlen(input->in_filename) + strlen(input->out_dirname) + 100];
  sprintf(path, "%s/%s.filtered", input->out_dirname, input->in_filename);

  size_t ref_length = 0;
  int num_targets = in_file->bam_header_p->n_targets;

  char **sequence_labels = (char **) calloc(num_targets, sizeof(char *));
  bam_header_t *bam_header = (bam_header_t *) calloc(1, sizeof(bam_header_t));
  bam_header->n_targets = num_targets;
  bam_header->target_name = (char **) calloc(num_targets, sizeof(char *));
  bam_header->target_len = (uint32_t*) calloc(num_targets, sizeof(uint32_t));
  for (int i = 0; i < num_targets; i++) {
    sequence_labels[i] = strdup(in_file->bam_header_p->target_name[i]);

    bam_header->target_name[i] = strdup(in_file->bam_header_p->target_name[i]);
    bam_header->target_len[i] = in_file->bam_header_p->target_len[i];
  }
  bam_header->text = strdup("@PG\tID:hpg-bam filter\tVN:1.0\n");
  bam_header->l_text = strlen(bam_header->text);

  bam_file_t *out_file = bam_fopen_mode(path, bam_header, "w");
  bam_fwrite_header(bam_header, out_file);

  // create a db for unmapped reads
  sqlite3 *db;
  create_unmapped_read_db(&db);

  //------------------------------------------------------------------
  // workflow management
  //
  bam_filter_wf_input_t *wf_input = new_bam_filter_wf_input(in_file,
							    out_file,
							    sequence_labels,
							    input,
							    db);
  
  // create and initialize workflow
  workflow_t *wf = workflow_new();
  
  workflow_stage_function_t stage_functions[] = {bam_filter_worker};
  char *stage_labels[] = {"BAM filter worker"};
  workflow_set_stages(1, &stage_functions, stage_labels, wf);
  
  // optional producer and consumer functions
  workflow_set_producer(bam_filter_producer, "BAM filter producer", wf);
  workflow_set_consumer(bam_filter_consumer, "BAM filter consumer", wf);
  
  workflow_run_with(input->num_threads, wf_input, wf);
  
  // free memory
  workflow_free(wf);
  free_bam_filter_wf_input(wf_input);
  //
  // end of workflow management
  //------------------------------------------------------------------

  // free memory
  for (int i = 0; i < num_targets; i++) {
    free(sequence_labels[i]);
  }
  free(sequence_labels);

  // close files
  bam_fclose(in_file);
  out_file->bam_header_p = NULL;
  bam_fclose(out_file);


  // check if the output file must be post-processed
  int count = count_rows("unmapped", db);
  if (count) {
    int rc;
    char sql[1024];
    sqlite3_stmt *stmt;

    char tmp[strlen(path) + 100];
    sprintf(tmp, "%s.tmp", path);
    
    // move the output file to the temporary file
    char cmd[2 * strlen(tmp) + 100];
    sprintf(cmd, "mv %s %s", path, tmp);
    system(cmd);
    
    in_file = bam_fopen(tmp);
    out_file = bam_fopen_mode(path, bam_header, "w");
    bam_fwrite_header(bam_header, out_file);
    
    // read the whole temporary file and
    // copy to the output file only those bam1_t
    // that do not appear in the sqlite3 database
    char *id;
    bam1_t *bam1;
    bam1 = bam_init1();    
    //sprintf(sql, "SELECT id FROM ids WHERE id = '%s'", bam1_qname(bam1));
    sprintf(sql, "SELECT id FROM ids WHERE id = ?");
    rc = sqlite3_prepare_v2(db, sql, -1, &stmt, 0);
    while (bam_read1(in_file->bam_fd, bam1) > 0) {
      id = bam1_qname(bam1);
      printf("id = %s\n", id);
      rc = sqlite3_bind_text(stmt, 1, id, strlen(id), 0);
      printf("result form bind_text rc = %i\n", rc);
      if (rc == SQLITE_OK) {
	printf("step\n");
	rc = sqlite3_step(stmt);
	if (rc == SQLITE_DONE) {
	  printf("not found in the db: id = %s\n", id);
	  // not found in the db,
	  // then write it to the output file
	  bam_write1(out_file->bam_fd, bam1);
	}
      }
      sqlite3_reset(stmt);
      
      // to be sure that all bam1_t pointers are free
      bam_destroy1(bam1);
      bam1 = bam_init1();    
    }
    bam_destroy1(bam1);
    
    // remove temporary file
    //      sprintf(cmd, "rm -rf %s", tmp);
    //      system(cmd);
    
    // close files
    bam_fclose(in_file);
    out_file->bam_header_p = NULL;
    bam_fclose(out_file);
  }

  bam_header_destroy(bam_header);
  
  // close sqlite db
  sqlite3_close(db);
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
