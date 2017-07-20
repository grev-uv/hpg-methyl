#include "cal_seeker.h"

//====================================================================================
// apply_caling bs
//====================================================================================

array_list_t *filter_cals(size_t num_cals, size_t read_length, array_list_t *list,
                    caling_bs_stage_workspace_t *workspace) {
  cal_t *cal;
  int min_seeds, max_seeds;
  array_list_t *cal_list;

  //filter-incoherent CALs
  int founds[num_cals], found = 0;

  for (size_t j = 0; j < num_cals; j++) {
    founds[j] = 0;
    cal = array_list_get(j, list);
   
    if (cal->sr_list->size > 0) {
      int start = 0;
      int first = 1;
      size_t genome_start = 0;
      
      for (linked_list_item_t *list_item = cal->sr_list->first; 
           list_item != NULL; 
           list_item = list_item->next) {
	      seed_region_t *s = list_item->item;
	 
        if (start > s->read_start || s->read_start >= s->read_end) {
	        found++;
	        founds[j] = 1;
	      }
	
	      if (!first && 
	         ((s->genome_start < genome_start) || 
	          (s->genome_start - genome_start) > 2 * read_length)) {
	        found++;
	        founds[j] = 1;
	      }
	
	      first = 0;
	      start = s->read_end + 1;
	      genome_start = s->genome_end + 1;
      }
    } else {
      found++;
      founds[j] = 1;
    }
  }
  
  if (found) {
    min_seeds = 100000;
    max_seeds = 0;
    
    cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

    for (size_t j = 0; j < num_cals; j++) {
      if (!founds[j]) {
	      cal = array_list_get(j, list);
	      cal->num_seeds = cal->sr_list->size;

	      if (cal->num_seeds > max_seeds) {
          max_seeds = cal->num_seeds;
        }

	      if (cal->num_seeds < min_seeds) {
          min_seeds = cal->num_seeds;
        }

	      array_list_insert(cal, cal_list);
	      array_list_set(j, NULL, list);
      }
    }

    array_list_free(list, (void *) cal_free);
    num_cals = array_list_size(cal_list);
    list = cal_list;
  }
  
  num_cals = array_list_size(list);
  
  //Ricardo max_cals 100, put a variable
  int max = 100;		

  if (num_cals > max) {
    for(int j = num_cals - 1; j >= max; j--) {
      cal_free(array_list_remove_at(j, list));
    }
  }
 
  return list;
}

//------------------------------------------------------------------------------------

void clean_apply_caling_bs_stage_workspace(void* workspace) {
  if (workspace) {
    caling_bs_stage_workspace_t *wf = (caling_bs_stage_workspace_t*)workspace;
    
    free(wf);
  }
}

//------------------------------------------------------------------------------------

int apply_caling_bs(cal_seeker_input_t* input, batch_t *batch, caling_bs_stage_workspace_t *workspace) {
  LOG_DEBUG("========= APPLY CALING BS =========\n");
  
  struct timeval start, end;
  double time;

  if (time_on) { 
    start_timer(start); 
  }

  mapping_batch_t *mapping_batch = batch->mapping_batch;
  size_t num_cals;
  size_t total_reads = 0;
  size_t num_targets, target_pos = 0, target_pos2 = 0;
  fastq_read_t *read;
  size_t target_index;
  array_list_t *list;

  size_t *targets = mapping_batch->targets;

  num_targets = mapping_batch->num_targets;
  total_reads += num_targets;

  mapping_batch->extra_stage_do = 1;

  extern pthread_mutex_t mutex_sp;
  extern size_t TOTAL_READS_SEEDING;

  pthread_mutex_lock(&mutex_sp);
  TOTAL_READS_SEEDING += num_targets;
  pthread_mutex_unlock(&mutex_sp);
  
  float Ncg, Ngc;
  float margen = mapping_batch->margin;

  for (size_t i = 0; i < num_targets; i++) {
	  target_index = targets[i];

    Ncg = mapping_batch->histogram_sw[target_index];
    Ngc = 1 - Ncg;
    
    //==========================================================================================
    // GA reads
    LOG_DEBUG("searching CALs for GA reads\n");

    read = array_list_get(target_index, mapping_batch->GA_rev_fq_batch);
    list = mapping_batch->mapping_lists[target_index];

    if (Ngc <= margen) {
      num_cals = array_list_size(list);

      //Ricardo - sin filtrar
      list = filter_cals(num_cals, read->length, list, workspace);

      // and update targets
      mapping_batch->mapping_lists[target_index] = list;

      if (array_list_size(list) > 0) {
	      mapping_batch->targets[target_pos++] = target_index;
      }
    } else {
      array_list_clear(list, (void *) NULL);
    }
    
    LOG_DEBUG("searching CALs for CT reads\n");

    read = array_list_get(target_index, mapping_batch->CT_rev_fq_batch);
    list = mapping_batch->mapping_lists2[target_index];

    if (Ncg <= margen) {
      num_cals = array_list_size(list);

      //Ricardo - filtrado cals
      list = filter_cals(num_cals, read->length, list, workspace);

      // and update targets
      mapping_batch->mapping_lists2[target_index] = list;
      
      if (array_list_size(list) > 0) {
	      mapping_batch->targets2[target_pos2++] = target_index;
      }
    } else {
      array_list_clear(list, (void *) NULL);
    }
  } 

  // updating number of targets for the next stage
  mapping_batch->num_targets = target_pos;
  mapping_batch->num_targets2 = target_pos2;

  if (time_on) { 
    stop_timer(start, end, time); 
    timing_add(time, CAL_SEEKER, timing); 
  }

  LOG_DEBUG("========= END OF APPLYING CALING BS =========\n");

  if (batch->pair_input->pair_mng->pair_mode != SINGLE_END_MODE) {
    return BS_PRE_PAIR_STAGE;
  } else if ((batch->mapping_batch->num_targets  > 0) ||
	           (batch->mapping_batch->num_targets2 > 0)) {
    return BS_SW_STAGE;
  }
  
  return BS_POST_PAIR_STAGE;
}

//------------------------------------------------------------------------------------
// cal_seeker_input functions: init
//------------------------------------------------------------------------------------

void cal_seeker_input_init(list_t *regions_list, cal_optarg_t *cal_optarg, 
			   list_t* write_list, unsigned int write_size, 
			   list_t *sw_list, list_t *pair_list, 
			   genome_t *genome, bwt_optarg_t *bwt_optarg, 
			   bwt_index_t *index, metaexons_t *metaexons, 
			   cal_seeker_input_t *input) {
  input->regions_list = regions_list;
  input->cal_optarg = cal_optarg;
  input->batch_size = write_size;
  input->sw_list = sw_list;
  input->pair_list = pair_list;
  input->write_list = write_list;
  input->genome = genome;
  input->bwt_optarg = bwt_optarg;
  input->index = index;
  input->metaexons = metaexons;
}
