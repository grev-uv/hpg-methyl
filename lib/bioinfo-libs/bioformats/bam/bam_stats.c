/*
 * bam_stats.c
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#include "bam_stats.h"

//------------------------------------------------------------------------

bam_stats_input_t *bam_stats_input_new(char *in_filename, region_table_t *region_table,
				       int num_threads,int batch_size, void *db, void *hash) {
  bam_stats_input_t *input = (bam_stats_input_t *) calloc(1, sizeof(bam_stats_input_t));
  
  input->in_filename = strdup(in_filename);
  input->region_table = region_table;
  input->num_threads = num_threads;
  input->batch_size = batch_size;
  input->db = db;
  input->hash = hash;
  
  return input;
}

void bam_stats_input_free(bam_stats_input_t *input) {
  if (input) {
    if (input->in_filename) free(input->in_filename);

    free(input);
  }
}

//------------------------------------------------------------------------

bam_stats_output_t *bam_stats_output_new() {
  bam_stats_output_t *output = (bam_stats_output_t *) calloc(1, sizeof(bam_stats_output_t));
  
  // global statistics
  output->ref_length = 0;
  output->num_sequences = 0;
  
  output->single_end = 1;
  
  output->num_reads = 0;
  output->num_unique_alignments = 0;
  output->num_mapped_reads = 0;
  output->num_unmapped_reads = 0;

  output->min_alignment_length = 999999999;
  output->max_alignment_length = 0;

  // stats per strand
  output->num_unique_alignments_strand[0] = 0;
  output->num_unique_alignments_strand[1] = 0;
  output->num_mapped_reads_strand[0] = 0;
  output->num_mapped_reads_strand[1] = 0;

  // errors stats
  output->num_indels = 0;
  output->indels_acc = 0;
  for (int i = 0; i <= NUM_ERRORS_STATS ; i++) {
    output->num_errors[i] = 0;
  }

  // nucleotide content
  output->num_As = 0;
  output->num_Cs = 0;
  output->num_Gs = 0;
  output->num_Ts = 0;
  output->num_Ns = 0;
  for (int i = 0; i < 100 ; i++) {
    output->GC_content[i] = 0;
  }

  // insert
  output->min_insert_size = 99999999;
  output->max_insert_size = 0;
  output->insert_size_acc = 0;

  // quality
  output->min_quality = 999999;
  output->max_quality = 0;
  output->quality_acc = 0;
  for (int i = 0; i < QUALITY_STATS ; i++) {
    output->quality[i] = 0;
  }

  // coverage (depth): global and per chromosome
  output->unmapped_nts = 0;
  output->sequence_labels = NULL;
  output->sequence_depths_per_nt = NULL;
  output->sequence_lengths = NULL;
  output->depth_per_sequence = NULL;
  output->depth = 0.0f;

  //  output->gc_hash = kh_init(32);

  return output;
}

void bam_stats_output_free(bam_stats_output_t *output) {
  if (output) {
    
    if (output->sequence_labels) {
      for(int i = 0; i < output->num_sequences; i++) {
	if (output->sequence_labels[i]) free(output->sequence_labels[i]);
      }
      free(output->sequence_labels);
    } 
    
    if (output->sequence_depths_per_nt) {
      for(int i = 0; i < output->num_sequences; i++) {
	if (output->sequence_depths_per_nt[i]) free(output->sequence_depths_per_nt[i]);
      }
      free(output->sequence_depths_per_nt);
    }
    
    if (output->sequence_lengths) {
      free(output->sequence_lengths);
    }

    if (output->depth_per_sequence) {
      free(output->depth_per_sequence);
    }

    free(output);
  }
}

//====================================================================
// W O R K F L O W     F O R      S T A T I S T I C S
//====================================================================

#define CONSUMER_STAGE   -1

//--------------------------------------------------------------------
// structure between the different workflow stages
//--------------------------------------------------------------------

typedef struct bam_stats_wf_batch {
  bam_stats_input_t *in_stats;
  array_list_t *bam1_list;
  array_list_t *stats_list;
  bam_stats_output_t *tmp_stats;
  bam_stats_output_t *out_stats;
} bam_stats_wf_batch_t;

bam_stats_wf_batch_t *bam_stats_wf_batch_new(bam_stats_input_t *in_stats,
					     array_list_t *bam1_list,
					     bam_stats_output_t *tmp_stats,
					     bam_stats_output_t *out_stats) {
  
  bam_stats_wf_batch_t *b = (bam_stats_wf_batch_t *) calloc(1, sizeof(bam_stats_wf_batch_t));
  
  b->in_stats = in_stats;
  b->bam1_list = bam1_list;
  b->stats_list = NULL;
  b->tmp_stats = tmp_stats;
  b->out_stats = out_stats;
  
  return b;
}

void bam_stats_wf_batch_free(bam_stats_wf_batch_t *b) {
  if (b) free(b);
}

//--------------------------------------------------------------------
// workflow input
//--------------------------------------------------------------------

typedef struct bam_stats_wf_input {
  bam_file_t *in_file;
  bam_stats_input_t *in_stats;
  bam_stats_output_t *out_stats;
} bam_stats_wf_input_t;

bam_stats_wf_input_t *bam_stats_wf_input_new(bam_file_t *in_file,
					     bam_stats_input_t *in_stats,
					     bam_stats_output_t *out_stats) {
  
  bam_stats_wf_input_t *wfi = (bam_stats_wf_input_t *) calloc(1, sizeof(bam_stats_wf_input_t));

  wfi->in_file = in_file;
  wfi->in_stats = in_stats;
  wfi->out_stats = out_stats;

  return wfi;
}

void bam_stats_wf_input_free(bam_stats_wf_input_t *wfi) {
  if (wfi) free(wfi);
}

//--------------------------------------------------------------------
// workflow producer
//--------------------------------------------------------------------

void *bam_stats_producer(void *input) {
  
  bam_stats_wf_input_t *wf_input = (bam_stats_wf_input_t *) input;
  bam_stats_wf_batch_t *new_batch = NULL;
  int max_num_bam1s = wf_input->in_stats->batch_size;

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
    new_batch = bam_stats_wf_batch_new(wf_input->in_stats,
				       bam1_list,  
				       bam_stats_output_new(),
				       wf_input->out_stats);				      
  }
    
  return new_batch;
}

//--------------------------------------------------------------------
// workflow consumer
//--------------------------------------------------------------------
int bam_progress = 0;

int bam_stats_consumer(void *data) {
  bam_stats_wf_batch_t *batch = (bam_stats_wf_batch_t *) data;
  
  bam_stats_output_t *tmp = batch->tmp_stats;
  bam_stats_output_t *out = batch->out_stats;

  bam_progress += tmp->num_reads;

  // merge tmp stats into "final" stats
  out->num_reads += tmp->num_reads;
  out->num_unique_alignments += tmp->num_unique_alignments;
  out->num_mapped_reads += tmp->num_mapped_reads;
  out->num_unmapped_reads += tmp->num_unmapped_reads;
  out->num_mapped_reads_1 += tmp->num_mapped_reads_1;
  out->num_unmapped_reads_1 += tmp->num_unmapped_reads_1;
  out->num_mapped_reads_2 += tmp->num_mapped_reads_2;
  out->num_unmapped_reads_2 += tmp->num_unmapped_reads_2;

  out->num_As += tmp->num_As;
  out->num_Cs += tmp->num_Cs;
  out->num_Gs += tmp->num_Gs;
  out->num_Ts += tmp->num_Ts;
  out->num_Ns += tmp->num_Ns;
  for (int i = 0; i < 100 ; i++) {
    out->GC_content[i] += tmp->GC_content[i];
  }

  out->single_end = tmp->single_end;

  out->num_unique_alignments_strand[0] += tmp->num_unique_alignments_strand[0];
  out->num_unique_alignments_strand[1] += tmp->num_unique_alignments_strand[1];
  out->num_mapped_reads_strand[0] += tmp->num_mapped_reads_strand[0];
  out->num_mapped_reads_strand[1] += tmp->num_mapped_reads_strand[1];

  out->num_indels += tmp->num_indels;
  out->indels_acc += tmp->indels_acc;
  for (int i = 0; i <= NUM_ERRORS_STATS ; i++) {
    out->num_errors[i] += tmp->num_errors[i];
  }

  if (tmp->max_alignment_length > out->max_alignment_length) out->max_alignment_length = tmp->max_alignment_length;
  if (tmp->min_alignment_length < out->min_alignment_length) out->min_alignment_length = tmp->min_alignment_length;

  if (tmp->max_insert_size > out->max_insert_size) out->max_insert_size = tmp->max_insert_size;
  if (tmp->min_insert_size < out->min_insert_size) out->min_insert_size = tmp->min_insert_size;
  out->insert_size_acc += tmp->insert_size_acc;

  if (tmp->max_quality > out->max_quality) out->max_quality = tmp->max_quality;
  if (tmp->min_quality < out->min_quality) out->min_quality = tmp->min_quality;
  out->quality_acc += tmp->quality_acc;
  for (int i = 0; i < QUALITY_STATS ; i++) {
    out->quality[i] += tmp->quality[i];
  }

  region_t region;
  char **sequence_labels = out->sequence_labels;
  region_table_t * region_table = batch->in_stats->region_table;

  bam1_t *bam1;
  uint32_t bam_flag;
  size_t seq_id, start, end;
  array_list_t *bam1_list = batch->bam1_list;
  size_t num_items = array_list_size(bam1_list);
  for (int i = 0; i < num_items; i++) {
    bam1 = array_list_get(i, bam1_list);

    bam_flag = (uint32_t) bam1->core.flag;
    if (!(bam_flag & BAM_FUNMAP)) {
      seq_id = bam1->core.tid;
      start = bam1->core.pos;
      end = start + bam1->core.l_qseq;
      if (region_table) {
	region.chromosome = sequence_labels[seq_id];
	region.start_position = start;
	region.end_position = end;
      }

      if ((region_table == NULL) || find_region(&region, region_table)) {
	for (int p = start; p < end; p++) {
	  out->sequence_depths_per_nt[seq_id][p]++;
	}
      }
    }

    bam_destroy1(bam1);
  }
  array_list_free(bam1_list, NULL);

  // stats
  if (batch->in_stats->db && batch->stats_list) {
    khash_t(stats_chunks) *hash = batch->in_stats->hash;
    sqlite3 *db = batch->in_stats->db;
    array_list_t *stats_list = batch->stats_list;

    // insert stats-list into the db record-query-fields table
    insert_bam_query_fields_list(stats_list, db);

    bam_query_fields_t *fields;
    size_t num_items = array_list_size(stats_list);
    for (int i = 0; i < num_items; i++) {
      fields = array_list_get(i, stats_list);

      // update chunks hash and free bam query fields
      update_chunks_hash(fields->chr, fields->chr_length, BAM_CHUNKSIZE,
			 fields->start, fields->end, hash);
      bam_query_fields_free(fields);
    }

    array_list_free(batch->stats_list, NULL);
  }

  if (bam_progress % 500000 == 0) {
    LOG_INFO_F("%i reads processed !\n", bam_progress);
  }

  // free memory
  bam_stats_output_free(tmp);
  bam_stats_wf_batch_free(batch);
}

//--------------------------------------------------------------------
// workflow worker
//--------------------------------------------------------------------

int bam_stats_worker(void *data) {
  bam_stats_wf_batch_t *batch = (bam_stats_wf_batch_t *) data;

  array_list_t *bam1_list = batch->bam1_list;
  bam_stats_output_t *tmp = batch->tmp_stats;
  region_table_t * region_table = batch->in_stats->region_table;

  region_t region;
  char **sequence_labels = batch->out_stats->sequence_labels;
  size_t *sequence_lengths = batch->out_stats->sequence_lengths;
  
  bam1_t *bam1;
  int bam_seq_len, num_cigar_ops, gc_content;
  uint8_t* bam_seq;
  uint32_t bam_flag, cigar_int, *cigar;
  int isize, strand, num_errors, num_indels, indels_length;
  size_t chr_length, num_items = array_list_size(bam1_list);
  char *chr, *qname;


  sqlite3 *db;
  bam_query_fields_t *fields;
  array_list_t *stats_list = NULL;
  if (batch->in_stats->db) {
    //    db = (sqlite3 *) batch->in_stats->db;
    stats_list = array_list_new(num_items, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);  
  }

  for (int i = 0; i < num_items; i++) {
    
    bam1 = array_list_get(i, bam1_list);
    bam_flag = (uint32_t) bam1->core.flag;
    isize = 0;
    num_errors = 0;
    num_indels = 0;
    indels_length = 0;
    chr = sequence_labels[bam1->core.tid];
    chr_length = sequence_lengths[bam1->core.tid];
    qname = bam1_qname(bam1);
    strand = 0;
    /*
    printf("qname = %s\n", qname);
    printf("\tcore.tid = %i\n", bam1->core.tid);
    printf("\tchr = %s\n", chr);
    */

    if (region_table) {
      region.chromosome = chr;
      region.start_position = bam1->core.pos;
      region.end_position = region.start_position + bam1->core.l_qseq;
    }

    if ((region_table == NULL) || find_region(&region, region_table)) {
	tmp->num_reads++;

	if (bam_flag & BAM_FUNMAP) {
	  tmp->num_unmapped_reads++;
	} else {
	  
	  strand = (int) ((bam_flag & BAM_FREVERSE) > 0);
	  
	  num_errors = bam_aux2i(bam_aux_get(bam1, "NM"));
	  if (num_errors < NUM_ERRORS_STATS) {
	    tmp->num_errors[num_errors]++;
	  } else {
	    tmp->num_errors[NUM_ERRORS_STATS]++;
	  }

	  cigar = bam1_cigar(bam1);
	  num_cigar_ops = (int) bam1->core.n_cigar; 
	  for (int j = 0; j < num_cigar_ops; j++) {
	    cigar_int = cigar[j];
	    switch (cigar_int & BAM_CIGAR_MASK) {
	    case BAM_CINS:  //I: insertion to the reference
	    case BAM_CDEL:  //D: deletion from the reference
	      num_indels++;
	      indels_length += (cigar_int >> BAM_CIGAR_SHIFT);
	      break;
	    }
	  }
	  tmp->num_indels += num_indels;
	  tmp->indels_acc += indels_length;

	  if (!(bam_flag & BAM_FSECONDARY)) {
	    tmp->num_unique_alignments++;
	    tmp->num_unique_alignments_strand[strand]++;
	  }
	  
	  tmp->num_mapped_reads++;
	  tmp->num_mapped_reads_strand[strand]++;
	}
	
	if (bam_flag & BAM_FPAIRED) {
	  if (bam_flag & BAM_FUNMAP) {
	    if (bam_flag & BAM_FREAD1) {
	      tmp->num_unmapped_reads_1++;
	    } else {
	      tmp->num_unmapped_reads_2++;
	    }
	  } else {
	    if (bam_flag & BAM_FREAD1) {
	      tmp->num_mapped_reads_1++;
	    } else {
	      tmp->num_mapped_reads_2++;
	    }
	  }
	  
	  if (!(bam_flag & BAM_FUNMAP) && !(bam_flag & BAM_FMUNMAP) && (bam_flag & BAM_FPROPER_PAIR)) { 
	    //      if (bam_flag & BAM_FPROPER_PAIR) { 
	    isize = abs(bam1->core.isize);
	    if (isize > tmp->max_insert_size) tmp->max_insert_size = isize;
	    if (isize < tmp->min_insert_size) tmp->min_insert_size = isize;
	    tmp->insert_size_acc += isize;
	  }
	}
	
	if (!(bam_flag & BAM_FUNMAP)) {
	  if (bam1->core.qual > tmp->max_quality) tmp->max_quality = bam1->core.qual;
	  if (bam1->core.qual < tmp->min_quality) tmp->min_quality = bam1->core.qual;
	  tmp->quality_acc += bam1->core.qual;
	  tmp->quality[bam1->core.qual]++;
	}
	
	bam_seq = bam1_seq(bam1);
	bam_seq_len = bam1->core.l_qseq;
	
	if (bam_seq_len > tmp->max_alignment_length) tmp->max_alignment_length = bam_seq_len;
	if (bam_seq_len < tmp->min_alignment_length) tmp->min_alignment_length = bam_seq_len;
	
	gc_content = 0;
	for (int i = 0; i < bam_seq_len; i++) {
	  switch (bam1_seqi(bam_seq, i)) {
	  case 1:
	    tmp->num_As++;
	    break;
	  case 2:
	    tmp->num_Cs++;
	    gc_content++;
	    break;
	  case 4:
	    tmp->num_Gs++;
	    gc_content++;
	    break;
	  case 8:
	    tmp->num_Ts++;
	    break;
	  case 15:
	    tmp->num_Ns++;
	    break;
	  }
	}
	gc_content = round(100.0 * gc_content / bam_seq_len);
	tmp->GC_content[gc_content]++;

	if (stats_list) {
	  fields = bam_query_fields_new(qname, chr, chr_length, strand, bam1->core.pos + 1, 
					bam1->core.pos + bam1->core.l_qseq + 1, bam_flag,
					bam1->core.qual, num_errors, num_indels, indels_length,
					isize);

	  array_list_insert(fields, stats_list);
	}
    }
	
    tmp->single_end = (bam_flag & BAM_FPAIRED ? 0 : 1);
  } // end for num_items

  if (stats_list) {
    batch->stats_list = stats_list;
  }

  return CONSUMER_STAGE;
}

//--------------------------------------------------------------------
// workflow description
//--------------------------------------------------------------------

void bam_stats(bam_stats_input_t *input, bam_stats_output_t *output) {
  bam_file_t *bam_file = bam_fopen(input->in_filename);
  
  size_t ref_length = 0;
  int num_targets = bam_file->bam_header_p->n_targets;
  
  output->sequence_labels = (char **) calloc(num_targets, sizeof(char *));
  output->sequence_depths_per_nt = (int **) calloc(num_targets, sizeof(int *));
  output->sequence_lengths = (size_t *) calloc(num_targets, sizeof(size_t));
  output->depth_per_sequence = (double *) calloc(num_targets, sizeof(size_t));

  for (int i = 0; i < num_targets; i++) {
    ref_length += bam_file->bam_header_p->target_len[i];

    output->sequence_labels[i] = strdup(bam_file->bam_header_p->target_name[i]);
    output->sequence_depths_per_nt[i] = (int *) calloc(bam_file->bam_header_p->target_len[i], sizeof(int));
    output->sequence_lengths[i] = bam_file->bam_header_p->target_len[i];
  }

  output->ref_length = ref_length;
  output->num_sequences = num_targets;

  //------------------------------------------------------------------
  // workflow management
  //
  bam_stats_wf_input_t *wf_input = bam_stats_wf_input_new(bam_file,
							  input, output);
  
  // create and initialize workflow
  workflow_t *wf = workflow_new();
  
  workflow_stage_function_t stage_functions[] = {bam_stats_worker};
  char *stage_labels[] = {"BAM stats worker"};
  workflow_set_stages(1, &stage_functions, stage_labels, wf);
  
  // optional producer and consumer functions
  workflow_set_producer(bam_stats_producer, "BAM stats producer", wf);
  workflow_set_consumer(bam_stats_consumer, "BAM stats consumer", wf);
  
  workflow_run_with(input->num_threads, wf_input, wf);
  
  // free memory
  workflow_free(wf);
  bam_stats_wf_input_free(wf_input);
  //
  // end of workflow management
  //------------------------------------------------------------------

  output->num_nucleotides = output->num_As + output->num_Cs + 
                            output->num_Gs + output->num_Ts + output->num_Ns;
  
  size_t nt_depth, unmapped_nts = 0;
  size_t seq_len, acc_per_sequence = 0, acc = 0;
  for (int i = 0; i < num_targets; i++) {
    seq_len = output->sequence_lengths[i];
    acc_per_sequence = 0;
    for (int j = 0; j < seq_len; j++) {
      nt_depth = ((int) output->sequence_depths_per_nt[i][j]);
      if (nt_depth) {
	acc_per_sequence += nt_depth;
      } else {
	unmapped_nts++;
      }
    }
    acc += acc_per_sequence;
    output->depth_per_sequence[i] = 1.0f * acc_per_sequence / output->sequence_lengths[i];
  }

  output->depth = 1.0f * acc / ref_length;
  output->unmapped_nts = unmapped_nts;

  bam_fclose(bam_file);
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
