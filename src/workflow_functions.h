#ifndef WORKFLOW_FUNCTIONS_H
#define WORKFLOW_FUNCTIONS_H

#include <stdio.h>

#include "commons/log.h"
#include "commons/file_utils.h"

#include "commons/workflow_scheduler.h"

#include "bioformats/fastq/fastq_batch_reader.h"
#include "bioformats/bam/bam_file.h"

#include "buffers.h"
#include "bwt_server.h"
#include "region_seeker.h"
#include "cal_seeker.h"
#include "pair_server.h"
#include "sw_server.h"
#include "batch_writer.h"

extern pthread_mutex_t mutex_sp;

//--------------------------------------------------------------------
// Workflow producers
//--------------------------------------------------------------------

void *fastq_reader(void *input);

//--------------------------------------------------------------------
// Stage bs functions
//--------------------------------------------------------------------

// Burrows-Wheeler Transform
int bwt_stage_bs(void *data);

// Candidate Alignment Location
int cal_stage_bs(void *data);

// Pre-pair
int pre_pair_stage(void *data);

// Smith-Waterman Algorithm
int sw_stage_bs(void *data);

// Post pair
int post_pair_stage_bs(void *data);

// Methylation status
int bs_status_stage(void *data);

//---------------------------------------------------------------------
// Workflow input
//--------------------------------------------------------------------- 

typedef struct wf_input {
  fastq_batch_reader_input_t *fq_reader_input;
  batch_t *batch;
} wf_input_t;

wf_input_t *wf_input_new(fastq_batch_reader_input_t *fq_reader_input, batch_t *batch);
void wf_input_free(wf_input_t *wfi);


#endif // WORKFLOW_FUNCTIONS_H
