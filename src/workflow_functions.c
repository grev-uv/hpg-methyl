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
    fastq_fread_bytes_aligner_pe(reads, fq_reader_input->batch_size, 
    fq_reader_input->fq_file1, fq_reader_input->fq_file2);
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

int bwt_stage_bs(void *data) {
  batch_t *batch = (batch_t *) data;
  return apply_bwt_bs(batch->bwt_input, batch);
}

//--------------------------------------------------------------------

int cal_stage_bs(void *data) {
  batch_t *batch = (batch_t *) data;
  return apply_caling_bs(batch->cal_input, batch);
}

//---------------------------------------------------------------------

int pre_pair_stage(void *data) {
  batch_t *batch = (batch_t *) data;
  return apply_pair(batch->pair_input, batch);
}

//--------------------------------------------------------------------

int sw_stage_bs(void *data) {
  batch_t *batch = (batch_t *) data;
  return apply_sw_bs(batch->sw_input, batch);
}

//--------------------------------------------------------------------

int post_pair_stage_bs(void *data) {
  batch_t *batch = (batch_t *) data;
  return prepare_alignments_bs(batch->pair_input, batch);
}

//--------------------------------------------------------------------

int bs_status_stage(void *data) {
  batch_t *batch = (batch_t *) data;
  return methylation_status_report(batch->sw_input, batch);
}
