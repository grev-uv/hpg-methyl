#include "workflow_functions.h"

int limit_of_reads = 0;
int n_insert = 0;
int tot_reads2 = 0;

//====================================================================================

wf_input_t *wf_input_new(fastq_batch_reader_input_t *fq_reader_input,
                         batch_t *batch) {

  wf_input_t *wfi = (wf_input_t *) calloc(1, sizeof(wf_input_t));
  wfi->fq_reader_input = fq_reader_input;
  wfi->batch = batch;

  return wfi;
}

//--------------------------------------------------------------------------------------

void wf_input_free(wf_input_t *wfi) {
  if (wfi) free(wfi);
}

//--------------------------------------------------------------------
// workflow producer
//--------------------------------------------------------------------

void *fastq_reader(void *input) {
  struct timeval start, end;
  double time;

  if (time_on) { 
    start_timer(start); 
  }

  wf_input_t *wf_input = (wf_input_t *) input;
  batch_t *new_batch = NULL;
  batch_t *batch = wf_input->batch;

  fastq_batch_reader_input_t *fq_reader_input = wf_input->fq_reader_input;
  array_list_t *reads = array_list_new(10000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  if (fq_reader_input->flags == SINGLE_END_MODE) {
    fastq_fread_bytes_se(reads, fq_reader_input->batch_size, fq_reader_input->fq_file1);
  } else {
    fastq_fread_bytes_aligner_pe(reads, fq_reader_input->batch_size, fq_reader_input->fq_file1, fq_reader_input->fq_file2);
  }

  size_t num_reads = array_list_size(reads);

  if (num_reads == 0) {
    array_list_free(reads, (void *)fastq_read_free);
  } else {
    mapping_batch_t *mapping_batch = mapping_batch_new(reads, batch->pair_input->pair_mng);

    new_batch = batch_new(batch->bwt_input, batch->region_input, batch->cal_input, 
    batch->pair_input, batch->sw_input, batch->writer_input, 
    batch->mapping_mode, mapping_batch, batch->write_mcontext);
  }

  if (time_on) { 
    stop_timer(start, end, time); 
    timing_add(time, FASTQ_READER, timing); 
  }

  return new_batch;
}

//--------------------------------------------------------------------
// stage functions
//--------------------------------------------------------------------

int bwt_stage_bs(work_item_t *item, int thread_id) {
  work_item_t *wf = (work_item_t*)item;

  batch_t *batch = wf->data;
  workflow_t *workflow = (workflow_t*)wf->context;

  if (workflow->thread_workspace[thread_id][BS_BWT_STAGE] == NULL) {
    workflow->thread_workspace[thread_id][BS_BWT_STAGE] = calloc(1, sizeof(bwt_stage_bs_workspace_t));
  }

  return apply_bwt_bs(batch->bwt_input, batch, workflow->thread_workspace[thread_id][BS_BWT_STAGE]);
}

//--------------------------------------------------------------------

int cal_stage_bs(work_item_t *item, int thread_id) {
  work_item_t *wf = (work_item_t*)item;

  batch_t *batch = wf->data;
  workflow_t *workflow = (workflow_t*)wf->context;

  if (workflow->thread_workspace[thread_id][BS_CAL_STAGE] == NULL) {
    workflow->thread_workspace[thread_id][BS_CAL_STAGE] = calloc(1, sizeof(caling_bs_stage_workspace_t));
  }

  return apply_caling_bs(batch->cal_input, batch, workflow->thread_workspace[thread_id][BS_CAL_STAGE]);
}

//---------------------------------------------------------------------

int pre_pair_stage(work_item_t *item, int thread_id) {
  work_item_t *wf = (work_item_t*)item;

  batch_t *batch = wf->data;
  workflow_t *workflow = (workflow_t*)wf->context;

  if (workflow->thread_workspace[thread_id][BS_PRE_PAIR_STAGE] == NULL) {
    workflow->thread_workspace[thread_id][BS_PRE_PAIR_STAGE] = calloc(1, sizeof(apply_pair_bs_stage_workspace_t));
  }

  return apply_pair(batch->pair_input, batch, workflow->thread_workspace[thread_id][BS_PRE_PAIR_STAGE]);
}

//--------------------------------------------------------------------

int sw_stage_bs(work_item_t *item, int thread_id) {
  work_item_t *wf = (work_item_t*)item;

  batch_t *batch = wf->data;
  workflow_t *workflow = (workflow_t*)wf->context;

  if (workflow->thread_workspace[thread_id][BS_SW_STAGE] == NULL) {
    workflow->thread_workspace[thread_id][BS_SW_STAGE] = calloc(1, sizeof(apply_sw_bs_stage_workspace_t));
  }

  return apply_sw_bs(batch->sw_input, batch, workflow->thread_workspace[thread_id][BS_SW_STAGE]);
}

//--------------------------------------------------------------------

int post_pair_stage_bs(work_item_t *item, int thread_id) {
  work_item_t *wf = (work_item_t*)item;

  batch_t *batch = wf->data;
  workflow_t *workflow = (workflow_t*)wf->context;

  if (workflow->thread_workspace[thread_id][BS_POST_PAIR_STAGE] == NULL) {
    workflow->thread_workspace[thread_id][BS_POST_PAIR_STAGE] = calloc(1, sizeof(prepare_alignments_bs_workspace_t));
  }

  return prepare_alignments_bs(batch->pair_input, batch, workflow->thread_workspace[thread_id][BS_POST_PAIR_STAGE]);
}

//--------------------------------------------------------------------

int bs_status_stage(work_item_t *item, int thread_id) {
  work_item_t *wf = (work_item_t*)item;

  batch_t *batch = wf->data;
  workflow_t *workflow = (workflow_t*)wf->context;

  if (workflow->thread_workspace[thread_id][BS_STATUS_STAGE] == NULL) {
    workflow->thread_workspace[thread_id][BS_STATUS_STAGE] = calloc(1, sizeof(methylation_stage_workspace_t));
  }

  return methylation_status_report(batch->sw_input, batch, workflow->thread_workspace[thread_id][BS_STATUS_STAGE]);
}
