#include "pair_server.h"

//------------------------------------------------------------------------------------
// this functions writes reads to disk, reads came from
// a list
//------------------------------------------------------------------------------------

inline void filter_alignments(char report_all, 
			      size_t report_n_best, 
			      size_t report_n_hits,
			      int report_best,
			      array_list_t *mapping_list);

//------------------------------------------------------------------------------------

inline array_list_t *create_new_list(size_t *valid_items, size_t num_valids, array_list_t *list) {
  size_t num_items = array_list_size(list);
  int flag = array_list_get_flag(list);

  array_list_t *new_list = array_list_new(num_valids, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_set_flag(flag, new_list);

  for (int k = 0; k < num_items; k++) {
    if (valid_items[k] == 1) {
      array_list_insert(array_list_get(k, list), new_list);
      array_list_set(k, NULL, list);
    }
  }

  if (flag == ALIGNMENTS_FOUND) {
    array_list_free(list, (void *) alignment_free);
  } else {
    array_list_free(list, (void *) cal_free);
  }

  return new_list;
}

//------------------------------------------------------------------------------------

void update_mispaired_pair(int pair_num, size_t num_items, array_list_t *list) {
  alignment_t *alig;

  for (size_t i = 0; i < num_items; i++) {
    alig = (alignment_t *) array_list_get(i, list);

    // set pair fields
    alig->mate_position = 0;
    alig->mate_chromosome = 0;
    alig->template_length = 0;
     
    alig->is_paired_end = 1;
    alig->is_paired_end_mapped = 0;
    alig->is_mate_mapped = 0;
    alig->mate_strand = 0;
    alig->pair_num = pair_num;
  }
}

//------------------------------------------------------------------------------------

inline void update_mispaired_pairs(size_t num_items1, size_t num_items2,
				   array_list_t *list1, array_list_t *list2) {
  alignment_t *alig;
  alignment_t *first1 = array_list_get(0, list1);
  alignment_t *first2 = array_list_get(0, list2);

  for (size_t i = 0; i < num_items1; i++) {
    alig = (alignment_t *) array_list_get(i, list1);

    // set pair1 fields
    alig->mate_position = first2->position;
    alig->mate_chromosome = first2->chromosome;  
    alig->template_length = 0;
     
    alig->is_paired_end = 1;
    alig->is_paired_end_mapped = 0;
    alig->is_mate_mapped = 1;
    alig->mate_strand = first2->seq_strand;
    alig->pair_num = 1;
  }

  for (size_t i = 0; i < num_items2; i++) {
    alig = (alignment_t *) array_list_get(i, list2);

    // set pair2 fields
    alig->mate_position = first1->position;
    alig->mate_chromosome = first1->chromosome;
    alig->template_length = 0;
     
    alig->is_paired_end = 1;
    alig->is_paired_end_mapped = 0;
    alig->is_mate_mapped = 1;
    alig->mate_strand = first1->seq_strand;
    alig->pair_num = 2;
  }
}

//------------------------------------------------------------------------------------

pair_t *pair_new(int index1, int index2, float score) {
  pair_t *p = (pair_t *) calloc(1, sizeof(pair_t));
  p->index1 = index1;
  p->index2 = index2;
  p->score = score;
  return p;
}

//------------------------------------------------------------------------------------

void pair_free(pair_t *p) {
  if (p) 
    free(p);
}

//------------------------------------------------------------------------------------

inline size_t select_n_hits(array_list_t *mapping_list, size_t report_n_hits) {
  alignment_t *aux_alignment;
  size_t num_mappings = array_list_size(mapping_list);
  int i;

  for (i = num_mappings - 1; i >= report_n_hits; i--) {
    aux_alignment = array_list_remove_at(i, mapping_list);
    alignment_free(aux_alignment);
  }

  return array_list_size(mapping_list);
}

//-----------------------------------------------------------------------------
extern unsigned int alignmentcmp(alignment_t *alignment_1, alignment_t *alignment_2);

inline size_t select_best_hits(array_list_t *mapping_list, size_t report_n_best) {
  int j, i;
  int primary_delete = 0;
  size_t best_pos, array_size;
  size_t num_mappings = array_list_size(mapping_list);
  alignment_t *best_alignment, *aux_alignment;
  array_list_t *mapping_list_filter;

  mapping_list_filter = array_list_new(num_mappings + 1, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  
  for (j = 0; j < report_n_best; j++) {     
    best_pos = 0;
    best_alignment = array_list_get( 0, mapping_list);
    array_size = array_list_size(mapping_list);

    for (i = 1; i < array_size; i++) {
      aux_alignment = array_list_get(i, mapping_list);

      if (alignmentcmp(best_alignment, aux_alignment) == 2) {
        best_alignment = aux_alignment;
        best_pos = i;
      }
    }

    array_list_insert(array_list_remove_at(best_pos, mapping_list), mapping_list_filter);
  }
  
  // free all mapings discarded
  array_size = array_list_size(mapping_list);
  
  for (j = array_size - 1; j >= 0; j--) {
    aux_alignment = array_list_remove_at(j, mapping_list);

    if (!is_secondary_alignment(aux_alignment)) { 
      primary_delete = 1; 
    } 

    alignment_free(aux_alignment);
  }

  for (j = report_n_best - 1; j >= 0; j--) {
    aux_alignment = array_list_remove_at(j, mapping_list_filter);
    array_list_insert(aux_alignment, mapping_list);
  }
  
  array_list_free(mapping_list_filter, NULL);

  if (primary_delete) {
    aux_alignment = array_list_get(0, mapping_list);

    if (aux_alignment) {
      set_secondary_alignment(0, aux_alignment);
    }
  }

  return array_list_size(mapping_list);
}

//-----------------------------------------------------------------------------

inline void select_best (array_list_t *mapping_list) {
  size_t num_items = array_list_size(mapping_list);
  int max_quality = 0;
  alignment_t *alig;

  for (int i = 0; i < num_items; i++) {
    alig = array_list_get(i, mapping_list);

    if (alig->map_quality > max_quality) {
      max_quality = alig->map_quality;
    }
  }

  for (int i = num_items - 1; i >= 0; i--) {
    alig = array_list_get(i, mapping_list);
    
    if (alig->map_quality != max_quality) {
      alig = array_list_remove_at(i, mapping_list);
      alignment_free(alig);
    } 
  }
}

//-----------------------------------------------------------------------------

inline void filter_alignments(char report_all, size_t report_n_best, size_t report_n_hits,
			      int report_best, array_list_t *mapping_list) {
  size_t num_mappings = array_list_size(mapping_list);

  if (!report_all && num_mappings) {
    if (report_best) {
      // max-score
      select_best(mapping_list);
    } else if (report_n_best > 0) {
      // n-best
      if (num_mappings > report_n_best) { 
	      select_best_hits(mapping_list, report_n_best);
      }
    } else if (report_n_hits > 0) {
      // n-hits
      if (num_mappings > report_n_hits) { 	
	      select_n_hits(mapping_list, report_n_hits);        
      }
    } else {
    }
  }
}

//-----------------------------------------------------------------------------

void filter_alignments_lists(char report_all, size_t report_n_best, size_t report_n_hits,
                  int report_best, size_t num_lists, array_list_t **mapping_lists) {
  array_list_t *mapping_list;

  for (int i = 0; i < num_lists; i++) {
    mapping_list = mapping_lists[i];
    filter_alignments(report_all, report_n_best, report_n_hits,report_best,mapping_list);
  }
}

//====================================================================================
// main functions: apply pair and prperare alignments
//====================================================================================

void clean_apply_pair_stage_workspace(void* workspace) {
  if (workspace) {
    apply_pair_bs_stage_workspace_t *wf = (apply_pair_bs_stage_workspace_t*)workspace;

    if (wf->mapped1) {
      free(wf->mapped1);
    }

    if (wf->mapped2) {
      free(wf->mapped2);
    }

    free(wf);
  }
}

//====================================================================================

int apply_pair(pair_server_input_t* input, batch_t *batch, apply_pair_bs_stage_workspace_t *workspace) {
  struct timeval prepair_start, prepair_end;
  float prepair_time = 0.0f;
  
  if (time_on) {
    start_timer(prepair_start);
  }




  mapping_batch_t *mapping_batch = batch->mapping_batch;
  array_list_t *fq_batch = mapping_batch->fq_batch;

  int pair_mode = input->pair_mng->pair_mode;
  size_t min_distance = input->pair_mng->min_distance;
  size_t max_distance = input->pair_mng->max_distance;
  int distance;


  size_t num_items1, num_items2, num_reads = array_list_size(fq_batch);

  int flag1, flag2;
  array_list_t *list1, *list2;

  size_t end1, start2;
  short int chr1, chr2, strand1, strand2;

  size_t mapped1_counter = 0, mapped2_counter = 0;

  if (workspace->allocated_mapped1 == 0) {
    workspace->allocated_mapped1 = 100;
  }

  if (workspace->allocated_mapped2 == 0) {
    workspace->allocated_mapped2 = 100;
  }

  size_t allocated_mapped1 = workspace->allocated_mapped1;
  size_t allocated_mapped2 = workspace->allocated_mapped2;

  if (workspace->mapped1 == NULL) {
    workspace->mapped1 = malloc(allocated_mapped1 * sizeof(size_t));
  }

  if (workspace->mapped2 == NULL) {
    workspace->mapped2 = malloc(allocated_mapped2 * sizeof(size_t));
  }

  size_t *mapped1 = workspace->mapped1;
  size_t *mapped2 = workspace->mapped2;

  int pair_found;

  alignment_t *alig;
  cal_t *cal;

  for (size_t i = 0; i < num_reads; i += 2) {

	  //In order to differentiate list mode (mapping_list and mapping_list2, or vice versa)
	  for (int option = 0; option < 2; option++) {




	if (option == 0)
	{

		list1 = mapping_batch->mapping_lists[i];
		list2 = mapping_batch->mapping_lists2[i + 1];

		flag1 = array_list_get_flag(list1);
		flag2 = array_list_get_flag(list2);

		num_items1 = 0;

		if (list1 != NULL) {
			num_items1 = array_list_size(list1);
		}

		num_items2 = 0;

		if (list2 != NULL) {
			num_items2 = array_list_size(list2);
		}


	}

	else if (option == 1)
    {

    	list1 = mapping_batch->mapping_lists2[i];
    	list2 = mapping_batch->mapping_lists[i + 1];

    	flag1 = array_list_get_flag(list1);
    	flag2 = array_list_get_flag(list2);

    	num_items1 = 0;

    	if (list1 != NULL) {
    	    num_items1 = array_list_size(list1);
    	}

    	num_items2 = 0;

    	if (list2 != NULL) {
    	    num_items2 = array_list_size(list2);
    	}
    }

    if (num_items1 > 0 && num_items2 > 0)
    		//&& (num_items1 > 1 || num_items2 > 1))
    {
      // initializes memory and counters
      mapped1_counter = 0;

      if (allocated_mapped1 < num_items1) {
        workspace->mapped1 = realloc(workspace->mapped1, num_items1 * sizeof(size_t));
        workspace->allocated_mapped1 = num_items1;

        mapped1 = workspace->mapped1;
        allocated_mapped1 = workspace->allocated_mapped1;
      }

      memset(mapped1, 0, num_items1 * sizeof(size_t));
      mapped2_counter = 0;

      if (allocated_mapped2 < num_items2) {
        workspace->mapped2 = realloc(workspace->mapped2, num_items2 * sizeof(size_t));
        workspace->allocated_mapped2 = num_items2;

        mapped2 = workspace->mapped2;
        allocated_mapped2 = workspace->allocated_mapped2;
      }

      memset(mapped2, 0, num_items2 * sizeof(size_t));
      pair_found = 0;





      // search for pairs properly aligned
      for (size_t j1 = 0; j1 < num_items1; j1++) {
        if (flag1 == ALIGNMENTS_FOUND) {
          alig = (alignment_t *) array_list_get(j1, list1);
          chr1 = alig->chromosome;
          strand1 = alig->seq_strand;
          end1 = alig->position + strlen(alig->sequence);
        } else if (flag1 == MULTIPLE_ANCHORS || flag1 == SINGLE_ANCHORS) {
          cal = (cal_t *) array_list_get(j1, list1);
          chr1 = cal->chromosome_id - 1;
          strand1 = cal->strand;
          end1 = cal->end;
        }
        else if (flag1 == NOT_ANCHORS)
        {
      	  chr1 = -1;
      	  strand1 = -1;
      	  end1 = 0;
      	  break;

        }
      	 else
      	 {
      		 printf("Error in pair_server.c, apply_sw function (pair2)\n");
      		 abort();
        }

        for (size_t j2 = 0; j2 < num_items2; j2++) {
          if (mapped2[j2] == 1) {
            continue;
          }

          if (flag2 == ALIGNMENTS_FOUND) {
            alig = (alignment_t *) array_list_get(j2, list2);
            chr2 = alig->chromosome;
            strand2 = alig->seq_strand;
            start2 = alig->position;
          } else if (flag2 == MULTIPLE_ANCHORS || flag2 == SINGLE_ANCHORS) {
            cal = (cal_t *) array_list_get(j2, list2);
            chr2 = cal->chromosome_id - 1;
            strand2 = cal->strand;
            start2 = cal->start;
          }
          else if (flag2 == NOT_ANCHORS)
          {
        	  chr2 = -1;
        	  strand2 = -1;
        	  start2 = 0;
        	  break;

          }
        	 else
        	 {
        		 printf("Error in pair_server.c, apply_sw function (pair2)\n");
        		 abort();
          }

          // computes distance between alignments,
          // is a valid distance ?
          distance = (start2 > end1 ? start2 - end1 : end1 - start2);

          if ((chr1 == chr2)                                           &&
              (distance >= min_distance) && (distance <= max_distance) &&
             ((strand1 != strand2 && pair_mode == PAIRED_END_MODE)     ||
              (strand1 == strand2 && pair_mode == MATE_PAIR_MODE))) {
            mapped1[j1] = 1;
            mapped2[j2] = 1;
            mapped1_counter++;
            mapped2_counter++;

            pair_found = 1;
            break;
          }
        } // end for j2..num_items2
      } // end for j1..num_item1




      if (pair_found) {
        // removing no valid items
        if (mapped1_counter != num_items1) {
        	if (option == 0)
        		mapping_batch->mapping_lists[i] = create_new_list(mapped1, mapped1_counter, list1);
        	else
        		mapping_batch->mapping_lists2[i] = create_new_list(mapped1, mapped1_counter, list1);
        }

        if (mapped2_counter != num_items2) {
        	if (option == 0)
        		mapping_batch->mapping_lists2[i + 1] = create_new_list(mapped2, mapped2_counter, list2);
        	else
        		mapping_batch->mapping_lists[i + 1] = create_new_list(mapped2, mapped2_counter, list2);
        }
        //Set pair_found variable
        if (option == 0)
        {
        	mapping_batch->mapping_lists[i]->paired = 1;
        	mapping_batch->mapping_lists2[i+1]->paired = 1;
        }
        else
        {
        	mapping_batch->mapping_lists2[i]->paired = 1;
        	mapping_batch->mapping_lists[i+1]->paired = 1;
        }


      }      
    }
   }
  }

  if (time_on) {
    stop_timer(prepair_start, prepair_end, prepair_time);
    timing_add(prepair_time, PRE_PAIR_TIME, timing);
  }

  // go to the next stage

if ((batch->mapping_batch->num_targets  > 0) ||
	(batch->mapping_batch->num_targets2 > 0)) {
    return BS_SW_STAGE;
  }

  return BS_POST_PAIR_STAGE;


}

//====================================================================================
// main functions: apply pair and prepare alignments
//====================================================================================

void clean_prepare_alignments_bs_stage_workspace(void *workspace) {
  if (workspace) {
    prepare_alignments_bs_workspace_t *wf = (prepare_alignments_bs_workspace_t*)workspace;

    free(wf);
  }
}

//====================================================================================

int prepare_alignments_bs(pair_server_input_t *input, batch_t *batch, prepare_alignments_bs_workspace_t *workspace) {
  struct timeval post_pair_time_start, post_pair_time_end;
  float post_pair_time = 0.0f;
  alignment_t *pair1,*pair2;
  array_list_t *list1, *list2;
  size_t num_items1, num_items2;

  int option = 0;   //In order to differentiate list mode (mapping_list and mapping_list2, or vice versa)

  if (time_on) {
    start_timer(post_pair_time_start);
  }

  if (input->pair_mng->pair_mode == SINGLE_END_MODE) {
    // Filter alignments by best-alignments, n-hits, all and unpaired reads
    filter_alignments_lists(input->report_optarg->all, 
			    input->report_optarg->n_best, 
			    input->report_optarg->n_hits,
			    input->report_optarg->best,
			    array_list_size(batch->mapping_batch->fq_batch),
			    batch->mapping_batch->mapping_lists);

    filter_alignments_lists(input->report_optarg->all, 
			    input->report_optarg->n_best, 
			    input->report_optarg->n_hits,
			    input->report_optarg->best,
			    array_list_size(batch->mapping_batch->fq_batch),
			    batch->mapping_batch->mapping_lists2);
  }


  if (time_on) {
    stop_timer(post_pair_time_start, post_pair_time_end, post_pair_time);
    timing_add(post_pair_time, POST_PAIR_TIME, timing);
  }

  if (batch->write_mcontext) {
    // Return for the analysis without the post-process
	  return BS_STATUS_STAGE;
  } else {
    // Return for the analysis with the postprocess
	  return CONSUMER_STAGE;
  }
}

//------------------------------------------------------------------------------------

void pair_server_input_init(pair_mng_t *pair_mng, report_optarg_t *report_optarg, 
			    list_t* pair_list, list_t *sw_list,
			    list_t *write_list, pair_server_input_t* input) {

  input->report_optarg = report_optarg;

  input->pair_mng = pair_mng;
  input->pair_list = pair_list;
  input->sw_list = sw_list;
  input->write_list = write_list;
}
