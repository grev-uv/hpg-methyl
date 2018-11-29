#ifndef WORKFLOW_SCHEDULER_H
#define WORKFLOW_SCHEDULER_H

#ifdef __cplusplus
extern "C" {
#endif

#include <sys/syscall.h>
#include <sys/types.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

#include "commons/commons.h"
#include "containers/array_list.h"

//----------------------------------------------------------------------------------------

#define WORKFLOW_STATUS_RUNNING  1
#define WORKFLOW_STATUS_FINISHED 2

#define NUM_MAX_CPUS 128

//----------------------------------------------------------------------------------------
// work_item
//----------------------------------------------------------------------------------------

typedef struct work_item {
  int stage_id;
  void *data;
  void *context;
} work_item_t;

work_item_t *work_item_new(int stage_id, void *data, void *context);
void work_item_free(work_item_t *item);

typedef int (*workflow_stage_function_t) (work_item_t *item, int thread_id);
typedef void (*workflow_stage_workspace_cleanup_function_t) (void *workspace);

typedef void* (*workflow_producer_function_t) (void *data);
typedef int (*workflow_consumer_function_t) (void *data);

//----------------------------------------------------------------------------------------
// Workflow
//----------------------------------------------------------------------------------------

typedef struct workflow workflow_t;

struct workflow {
  int num_threads;
  int max_num_work_items;
  int num_stages;
  int completed_producer;
  int num_pending_items;     
  int running_producer;
  int running_consumer;
  
  int complete_extra_stage;
  
  pthread_cond_t producer_cond;
  pthread_cond_t consumer_cond;
  pthread_cond_t workers_cond;
  
  pthread_mutex_t producer_mutex;
  pthread_mutex_t consumer_mutex;     
  
  pthread_mutex_t main_mutex;
  
  double workflow_time;
  double producer_time;
  double consumer_time;
  double *stage_times;

  pthread_mutex_t *stage_times_mutex;     

  array_list_t **pending_items;
  array_list_t *completed_items;
  
  workflow_stage_function_t *stage_functions;
  char** stage_labels;
  
  workflow_producer_function_t *producer_function;
  char* producer_label;
  
  workflow_consumer_function_t *consumer_function;
  char* consumer_label;

  void*** thread_workspace;
  workflow_stage_workspace_cleanup_function_t *stage_cleanup_functions;
};

//----------------------------------------------------------------------------------------

typedef struct workflow_context {
  void *input;
  workflow_t *wf;
} workflow_context_t;

//----------------------------------------------------------------------------------------

typedef struct thread_context {
  workflow_context_t *wf;
  int thread_id;
} thread_context_t;

//----------------------------------------------------------------------------------------

workflow_t *workflow_new();
void workflow_free(workflow_t *wf);

void workflow_set_stages(int num_stages, workflow_stage_function_t *functions, char **labels, workflow_t *wf, 
                          workflow_stage_workspace_cleanup_function_t *cleanup_functions);
void workflow_set_producer(workflow_producer_function_t *function, char *label, workflow_t *wf);
void workflow_set_consumer(workflow_consumer_function_t *function, char *label, workflow_t *wf);

int workflow_get_num_items(workflow_t *wf);
int workflow_get_num_items_at(int stage_id, workflow_t *wf);
int workflow_get_num_completed_items(workflow_t *wf);
int workflow_is_producer_finished(workflow_t *wf);

void workflow_insert_item(void *data, workflow_t *wf);
void workflow_insert_item_at(int stage_id, void *data, workflow_t *wf);
void workflow_insert_stage_item_at(void *data, int new_stage, workflow_t *wf);
void *workflow_remove_item(workflow_t *wf);

int workflow_get_status(workflow_t *wf);
void workflow_producer_finished(workflow_t *wf);
void workflow_change_status_function(workflow_t *wf);

void workflow_run(void *input, workflow_t *wf);
void workflow_run_with(int num_threads, void *input, workflow_t *wf);

thread_context_t* thread_context_new(workflow_context_t *wf, int tid);
void thread_context_free(thread_context_t* th);

//----------------------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

#endif
