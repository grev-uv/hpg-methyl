#ifndef PAIR_SERVER_H
#define PAIR_SERVER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "commons/commons.h"
#include "commons/system_utils.h"
#include "containers/list.h"

#include "buffers.h"
#include "timing.h"
#include "sw_server.h"

//====================================================================================
//  structures and prototypes
//====================================================================================

struct pair_server_input {
  pair_mng_t *pair_mng;
  report_optarg_t *report_optarg;

  list_t* pair_list;
  list_t* sw_list;
  list_t* write_list;
};

//--------------------------------------------------------------------------------------

typedef struct pair {
  int index1;
  int index2;
  float score;
} pair_t;

//--------------------------------------------------------------------------------------

typedef struct apply_pair_bs_stage_workspace {
  size_t *mapped1;
  size_t *mapped2;

  size_t allocated_mapped1;
  size_t allocated_mapped2;
} apply_pair_bs_stage_workspace_t;

//--------------------------------------------------------------------------------------

typedef struct prepare_alignments_bs_workspace {

} prepare_alignments_bs_workspace_t;

//--------------------------------------------------------------------------------------
// In case the build is being compiled with all optimizations disabled (ie: debug build)
// inline functions must be forward declared to enable the compiler to find the appropriate
// simbols to link into them on all the translation units.
// - Date: 14 / 11 / 2016
// - Who: Cesar
#ifdef __GNUC__
#ifdef __NO_INLINE__
array_list_t *create_new_list(size_t*, size_t, array_list_t*);
void update_mispaired_pairs(size_t, size_t, array_list_t*, array_list_t*);
size_t select_n_hits(array_list_t*, size_t);
size_t select_best_hits(array_list_t*, size_t);
void select_best (array_list_t*);
void filter_alignments(char, size_t, size_t, int, array_list_t*);
#endif
#endif

//------------------------------------------------------------------------------------

void pair_server_input_init(pair_mng_t *pair_mng, report_optarg_t *report_optarg,
			    list_t* pair_list, list_t *sw_list,
			    list_t *write_list, pair_server_input_t* input);

//====================================================================================

void pair_server(pair_server_input_t* input);

//------------------------------------------------------------------------------------

int apply_pair(pair_server_input_t* input, batch_t *batch, apply_pair_bs_stage_workspace_t *workspace);
void clean_prepare_alignments_bs_stage_workspace(void *workspace);

int prepare_alignments_bs(pair_server_input_t* input, batch_t *batch, prepare_alignments_bs_workspace_t *workspace);
void clean_apply_pair_stage_workspace(void* workspace);

pair_t *pair_new(int index1, int index2, float score);
void pair_free(pair_t *p);

#endif // PAIR_SERVER_H
