
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

bs_context_t *bs_context_new(size_t num_reads, size_t num_chromosomes) {
  bs_context_t *p = (bs_context_t*) calloc(1, sizeof(bs_context_t));

  p->context_CpG = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  p->context_CHG = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  p->context_CHH = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  p->context_MUT = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  p->methyl_reads = calloc(num_chromosomes, sizeof(uint32_t));

  return p;
}

//------------------------------------------------------------------------------------

void bs_context_free(bs_context_t *p) {
  if (p) {
    free(p->methyl_reads);
    free(p);
  }
}
