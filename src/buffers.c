
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "buffers.h"

//------------------------------------------------------------------------------------

batch_t *batch_new(bwt_server_input_t *bwt_input,
                   region_seeker_input_t *region_input,
                   cal_seeker_input_t *cal_input,
                   pair_server_input_t *pair_input,
                   sw_server_input_t *sw_input,
                   batch_writer_input_t *writer_input,
		               int mapping_mode,
                   mapping_batch_t *mapping_batch,
				           bool write_mcontext) {
  batch_t *b = (batch_t *) calloc(1, sizeof(batch_t));
  b->bwt_input = bwt_input;
  b->region_input = region_input;
  b->cal_input = cal_input;
  b->pair_input = pair_input;
  b->sw_input = sw_input;
  b->writer_input = writer_input;
  b->mapping_batch = mapping_batch;
  b->mapping_mode = mapping_mode;
  b->write_mcontext = write_mcontext;

  return b;
}

//------------------------------------------------------------------------------------

void batch_free(batch_t *b) {
  if (b) {
    free(b);
  }
}

//------------------------------------------------------------------------------------

report_optarg_t *report_optarg_new(int all, int n_best, int n_hits, int only_paired, int best) {
  report_optarg_t *p = (report_optarg_t*) calloc(1, sizeof(report_optarg_t));

  p->all = all;
  p->n_best = n_best;
  p->n_hits = n_hits;
  p->only_paired = only_paired;
  p->best = best;

  return p;
}

//------------------------------------------------------------------------------------

void report_optarg_free(report_optarg_t *p) {
  if (p != NULL) {
    free(p);
  }
}

//------------------------------------------------------------------------------------

pair_mng_t *pair_mng_new(int pair_mode, size_t min_distance, size_t max_distance, 
                         int report_only_paired) {
  pair_mng_t *p = (pair_mng_t*) calloc(1, sizeof(pair_mng_t));

  p->pair_mode = pair_mode;
  p->min_distance = min_distance;
  p->max_distance = max_distance;
  p->report_only_paired = report_only_paired;

  return p;
}

//------------------------------------------------------------------------------------

void pair_mng_free(pair_mng_t *p) {
  if (p != NULL) {
    free(p);
  }
}

//------------------------------------------------------------------------------------

mapping_batch_t *mapping_batch_new(array_list_t *fq_batch, pair_mng_t *pair_mng) {
  mapping_batch_t *p = (mapping_batch_t *) calloc(1, sizeof(mapping_batch_t));
  size_t num_reads = array_list_size(fq_batch);

  p->action = BWT_ACTION;
  p->num_targets = 0;
  p->num_extra_targets = 0;
  p->num_allocated_targets = num_reads;
  p->extra_stage_do = 0;

  if (!pair_mng) { 
    p->pair_mng = pair_mng_new(SINGLE_END_MODE, 0, 0, 0); 
  } else {
    p->pair_mng = pair_mng_new(pair_mng->pair_mode, pair_mng->min_distance, 
			       pair_mng->max_distance, pair_mng->report_only_paired); 
  }

  p->num_gaps = 0;
  p->num_sws = 0;
  p->num_ext_sws = 0;

  p->num_to_do = 0;
  p->fq_batch = fq_batch;
  p->targets = (size_t *) calloc(num_reads, sizeof(size_t));
  p->extra_targets = (size_t *) calloc(num_reads, sizeof(size_t));
  p->extra_stage_id = (unsigned char *) calloc(num_reads, sizeof(unsigned char));
  p->mapping_lists = (array_list_t **) calloc(num_reads, sizeof(array_list_t*));
  p->bwt_mappings = (unsigned char *)calloc(num_reads, sizeof(unsigned char));

  for (size_t i = 0; i < num_reads; i++) {
    p->mapping_lists[i] = array_list_new(500, 
					 1.25f, 
					 COLLECTION_MODE_ASYNCHRONIZED); 
  }
  
  p->margin = BS_HIST_MARGIN;
  p->histogram_sw = (float *)calloc(num_reads, sizeof(float));

  // Added by PP for bisulfite
  p->num_targets2 = 0;
  p->num_to_do2 = 0;
  p->targets2 = (size_t *) calloc(num_reads, sizeof(size_t));

  return p;
}

//------------------------------------------------------------------------------------

mapping_batch_t *mapping_batch_new_2(size_t num_reads, array_list_t *fq_batch, pair_mng_t *pair_mng) {
  mapping_batch_t *p = (mapping_batch_t *) calloc(1, sizeof(mapping_batch_t));

  p->action = BWT_ACTION;
  p->num_targets = 0;
  p->num_extra_targets = 0;
  p->num_allocated_targets = num_reads;
  p->extra_stage_do = 0;

  if (!pair_mng) { 
    p->pair_mng = pair_mng_new(SINGLE_END_MODE, 0, 0, 0); 
  } else {
    p->pair_mng = pair_mng_new(pair_mng->pair_mode, pair_mng->min_distance, 
			       pair_mng->max_distance, pair_mng->report_only_paired); 
  }

  p->num_gaps = 0;
  p->num_sws = 0;
  p->num_ext_sws = 0;
  p->num_to_do = 0;
  p->fq_batch = fq_batch;
  p->targets = (size_t *) calloc(num_reads, sizeof(size_t));
  p->extra_targets = (size_t *) calloc(num_reads, sizeof(size_t));
  p->extra_stage_id = (unsigned char *) calloc(num_reads, sizeof(unsigned char));
  p->mapping_lists = (array_list_t **) calloc(num_reads, sizeof(array_list_t*));
  p->bwt_mappings = (unsigned char *)calloc(num_reads, sizeof(unsigned char));

  return p;
}

//------------------------------------------------------------------------------------

void mapping_batch_free(mapping_batch_t *p) {
  if (p == NULL) return;
  
  if (p->fq_batch) { array_list_free(p->fq_batch, (void *) fastq_read_free); }
  if (p->targets) { free(p->targets); }
  if (p->mapping_lists) { free(p->mapping_lists); }
  if (p->pair_mng) { free(p->pair_mng); }
  if (p->extra_stage_id) { free(p->extra_stage_id); }
  if (p->extra_targets) { free(p->extra_targets); }

  if (p->old_mapping_lists) { free(p->old_mapping_lists); }
  if (p->bwt_mappings) free(p->bwt_mappings);

  // Added by PP
  if (p->CT_fq_batch) { array_list_free(p->CT_fq_batch, (void *) fastq_read_free); }
  if (p->CT_rev_fq_batch) { array_list_free(p->CT_rev_fq_batch, (void *) fastq_read_free); }
  if (p->GA_fq_batch) { array_list_free(p->GA_fq_batch, (void *) fastq_read_free); }
  if (p->GA_rev_fq_batch) { array_list_free(p->GA_rev_fq_batch, (void *) fastq_read_free); }
  if (p->mapping_lists2) { free(p->mapping_lists2); }
  if (p->targets2) { free(p->targets2); }
  if (p->bs_status) {free(p->bs_status); }
  if (p->histogram_sw) {free(p->histogram_sw); }
  
  free(p);
}

//------------------------------------------------------------------------------------

void alignment_aux_init(alignment_t* alignment, alignment_aux_t *alignment_aux) {
  alignment_aux->seq_strand = alignment->seq_strand;
  alignment_aux->chromosome = alignment->chromosome; 
  alignment_aux->position = alignment->position;
  alignment_aux->num_cigar_operations = alignment->num_cigar_operations;
  alignment_aux->map_quality = alignment->map_quality;
  alignment_aux->optional_fields_length = alignment->optional_fields_length;
  alignment_aux->mapping_len = strlen(alignment->sequence);
  alignment_aux->cigar_len = strlen(alignment->cigar);
}

//------------------------------------------------------------------------------------

void file_write_alignments(fastq_read_t *fq_read, array_list_t *items, FILE *fd) {
  size_t num_items = array_list_size(items);
  
  if (num_items <= 0) { 
    return; 
  }

  int tot_len_cigar = 0, tot_len_of = 0;

  //[HEAD][SEQUENCE][QUALITY]
  size_t max_len_cigar = num_items*1024*2;
  char *buffer_cigar = (char *)malloc(sizeof(char)*max_len_cigar);

  size_t max_len_of = num_items*1024*2;
  uint8_t *buffer_of = (uint8_t *)malloc(sizeof(uint8_t)*max_len_of);
  alignment_aux_t alignment_aux[num_items];
  alignment_aux_t *alignment_a;

  memset(alignment_aux, 0, sizeof(alignment_aux_t)*num_items);  

  for (int i = 0; i < num_items; i++) {
    alignment_a = &alignment_aux[i];
    alignment_t *alignment = array_list_get(i, items);
    alignment_aux_init(alignment, alignment_a);     

    int cigar_len = strlen(alignment->cigar);    
    memcpy(&buffer_cigar[tot_len_cigar], alignment->cigar, cigar_len);
    tot_len_cigar += cigar_len;
    
    if (tot_len_cigar >= max_len_cigar) { 
      max_len_cigar = max_len_cigar * 2;
      buffer_cigar = realloc(buffer_cigar, max_len_cigar); 
    }

    int of_len = alignment->optional_fields_length;
    memcpy(&buffer_of[tot_len_of], alignment->optional_fields, of_len);
    tot_len_of += of_len;
    
    if (tot_len_of >= max_len_of) {
      max_len_of = max_len_of * 2;
      buffer_of = realloc(buffer_of, max_len_cigar); 
    }
  }

  fwrite(alignment_aux, sizeof(alignment_aux_t), num_items, fd);
  fwrite(buffer_cigar, sizeof(char), tot_len_cigar, fd);  
  fwrite(buffer_of, sizeof(uint8_t), tot_len_of, fd);  

  free(buffer_cigar);
}

//------------------------------------------------------------------------------------

void file_write_meta_alignments(fastq_read_t *fq_read, array_list_t *items, FILE *fd) {
  size_t seq_size  = fq_read->length;
  size_t num_items = array_list_size(items);

  if (!num_items) { 
    return; 
  }

  size_t max_len = num_items * 1024;
  char *cigar_buffer = (char *)calloc(max_len, sizeof(char));
  size_t tot_len = 0;

  simple_alignment_t simple_alignment[num_items];
  simple_alignment_t *simple_a;

  memset(simple_alignment, 0, sizeof(simple_alignment_t)*num_items);
  
  //[type][size head][size seq][num items][HEAD][SEQUENCE][QUALITY][CAL 0][CAL n]
  for (int i = 0; i < num_items; i++) {
    meta_alignment_t *meta_alignment = array_list_get(i, items);
    cal_t *first_cal = array_list_get(0, meta_alignment->cals_list);
    cal_t *last_cal = array_list_get(meta_alignment->cals_list->size - 1, meta_alignment->cals_list);
    seed_region_t *first_seed = linked_list_get_first(first_cal->sr_list);
    seed_region_t *last_seed  = linked_list_get_last(last_cal->sr_list);
    cigar_code_t *cigar_code = meta_alignment->cigar_code;    
    char *cigar_str = new_cigar_code_string(cigar_code);
    int cigar_len = strlen(cigar_str);

    simple_a = &simple_alignment[i];

    if (meta_alignment->cigar_left != NULL) {
      simple_a->gap_start = 0;
    } else {
      simple_a->gap_start = first_seed->read_start;
    }

    if (meta_alignment->cigar_right != NULL) {
      simple_a->gap_end = 0;
    } else {
      simple_a->gap_end = seq_size - last_seed->read_end - 1;
    }

    simple_a->map_strand = first_cal->strand;
    simple_a->map_chromosome = first_cal->chromosome_id;
    simple_a->map_start = first_cal->start;
    simple_a->map_distance = cigar_code->distance;
    simple_a->cigar_len = cigar_len;
    
    memcpy(&cigar_buffer[tot_len], cigar_str, cigar_len);
    tot_len += cigar_len;

    if (tot_len >= max_len) { 
      max_len = max_len * 2;
      cigar_buffer = realloc(cigar_buffer, max_len); 
    }
  }

  // Write binary file 
  //[size head][size seq][num items][HEAD][SEQUENCE][QUALITY][CAL 0][CAL n]
  fwrite(simple_alignment, sizeof(simple_alignment_t), num_items, fd);
  fwrite(cigar_buffer, sizeof(char), tot_len, fd);  

  free(cigar_buffer);
}

//------------------------------------------------------------------------------------

void file_write_cals(fastq_read_t *fq_read, array_list_t *items, FILE *fd) {
  size_t num_items = array_list_size(items);
  if (!num_items) { return; }
  
  bwt_anchor_t bwt_anchor[num_items];
  memset(bwt_anchor, 0, sizeof(bwt_anchor_t)*num_items);

  for (int i = 0; i < num_items; i++) {
    cal_t *cal = array_list_get(i, items);
    
    bwt_anchor[i].strand     = cal->strand;
    bwt_anchor[i].chromosome = cal->chromosome_id - 1;
    bwt_anchor[i].start      = cal->start;
    bwt_anchor[i].end        = cal->end;

    seed_region_t *seed = linked_list_get_first(cal->sr_list);

    if (seed->read_start == 0) {
      bwt_anchor[i].type = FORWARD_ANCHOR;
    } else {
      bwt_anchor[i].type = BACKWARD_ANCHOR;
    }
  }

  fwrite(bwt_anchor, sizeof(bwt_anchor_t), num_items, fd);
}

void file_write_fastq_read(fastq_read_t *fq_read, size_t num_items, FILE *fd) {
  size_t head_size = strlen(fq_read->id);
  size_t seq_size  = fq_read->length;

  // Write binary file 
  //[type][size head][size seq][num items][HEAD][SEQUENCE][QUALITY][CAL 0][CAL n]
  size_t items_sizes[3] = {head_size, seq_size, num_items};

  //[size head][size seq][num items]
  fwrite(items_sizes, sizeof(size_t), 3, fd);
  
  //[HEAD][SEQUENCE][QUALITY]
  size_t total_size = head_size + 2*seq_size;
  char *buffer = (char *)malloc(sizeof(char)*total_size);
  
  memcpy(buffer, fq_read->id, head_size);
  memcpy(&buffer[head_size], fq_read->sequence, seq_size);
  memcpy(&buffer[head_size + seq_size], fq_read->quality, seq_size);

  fwrite(buffer, sizeof(char), total_size, fd);

  free(buffer);
}

void file_write_items(fastq_read_t *fq_read, array_list_t *items, 
		      unsigned char data_type, FILE *fd1, FILE *fd2,
		      int mode) {
  FILE *fd;

  if (mode == 0) {
    fd = fd1;
  } else {
    fd = fd2;
  }

  fwrite(&data_type, sizeof(unsigned char), 1, fd);
  file_write_fastq_read(fq_read, array_list_size(items), fd);

  if (data_type == CAL_TYPE) {
    file_write_cals(fq_read, items, fd);
  } else if (data_type == META_ALIGNMENT_TYPE) {
    file_write_meta_alignments(fq_read, items, fd);    
  } else {
    file_write_alignments(fq_read, items, fd);    
  }
}

//------------------------------------------------------------------------------------

bs_context_t *bs_context_new(size_t num_reads) {
  bs_context_t *p = (bs_context_t*) calloc(1, sizeof(bs_context_t));

  p->context_CpG = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  p->context_CHG = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  p->context_CHH = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  p->context_MUT = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  return p;
}

//------------------------------------------------------------------------------------

void bs_context_free(bs_context_t *p) {
  if (p) {
    free(p);
  }
}
