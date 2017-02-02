#include "cal_seeker.h"

//====================================================================================
// apply_caling bs
//====================================================================================

array_list_t *filter_cals(size_t num_cals, size_t read_length, array_list_t *list) {
  cal_t *cal;
  int min_seeds, max_seeds;
  array_list_t *cal_list;
  size_t select_cals;

  //filter-incoherent CALs
  int founds[num_cals], found = 0;

  for (size_t j = 0; j < num_cals; j++) {
    founds[j] = 0;
    cal = array_list_get(j, list);

    LOG_DEBUG_F("\tcal %i of %i: sr_list size = %i (cal->num_seeds = %i) %i:%lu-%lu\n", 
		            j, num_cals, cal->sr_list->size, cal->num_seeds,
		            cal->chromosome_id, cal->start, cal->end);
   
    if (cal->sr_list->size > 0) {
      int start = 0;
      int first = 1;
      size_t genome_start = 0;
      
      for (linked_list_item_t *list_item = cal->sr_list->first; 
           list_item != NULL; 
           list_item = list_item->next) {
	      seed_region_t *s = list_item->item;
	
	      LOG_DEBUG_F("\t\t:: star %lu > %lu s->read_start\n", start, s->read_start);
	      LOG_DEBUG_F("\t\t:: read_star %lu > read_end %lu \n", s->read_start, s->read_end);
	      
        if (start > s->read_start || s->read_start >= s->read_end) {
	        LOG_DEBUG("\t\t\t:: remove\n");
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
  
  //Ricardo max_cals  100, put a variable
  int max = 100;		

  if (num_cals > max) {
    select_cals = num_cals - max;

    for(int j = num_cals - 1; j >= max; j--) {
      cal_free(array_list_remove_at(j, list));
    }
  }
 
  return list;
}

//------------------------------------------------------------------------------------

int apply_caling_bs(cal_seeker_input_t* input, batch_t *batch) {
  LOG_DEBUG("========= APPLY CALING BS =========\n");
  
  struct timeval start, end;
  double time;

  if (time_on) { 
    start_timer(start); 
  }

  metaexons_t *metaexons = input->metaexons;
  bwt_optarg_t *bwt_optarg = input->bwt_optarg;
  bwt_index_t *bwt_index = input->index;
  bwt_index_t *bwt_index2 = input->index2;
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  array_list_t *allocate_cals;
  size_t num_cals, total_cals = 0;
  size_t num_batches = 0, num_reads_unmapped = 0, num_without_cals = 0;
  size_t max_seeds, total_reads = 0;
  size_t num_targets, target_pos = 0, target_pos2 = 0;
  fastq_read_t *read, *read2;
  genome_t *genome = input->genome;
  size_t *targets_aux, target_index;
  int seed_size = input->cal_optarg->seed_size;
  array_list_t *list;
  region_t *bwt_region_back, *bwt_region_forw;
  linked_list_t *linked_list;
  seed_region_t *seed_region_start, *seed_region_end, *seed_region;
  int gap_nt, anchor_nt;
  bwt_anchor_t *bwt_anchor_back, *bwt_anchor_forw;

  size_t *targets = mapping_batch->targets;
  size_t *targets2 = mapping_batch->targets2;

  num_targets = mapping_batch->num_targets;
  total_reads += num_targets;

  mapping_batch->extra_stage_do = 1;

  extern pthread_mutex_t mutex_sp;
  extern size_t TOTAL_READS_SEEDING, TOTAL_READS_SEEDING2;

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
    read2 = array_list_get(target_index, mapping_batch->GA_fq_batch);
    list = mapping_batch->mapping_lists[target_index];

    if (Ngc <= margen) {
      max_seeds = (read->length / 15) * 2 + 10;
      num_cals = array_list_size(list);

      //Ricardo - sin filtrar
      list = filter_cals(num_cals, read->length, list);

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
    read2 = array_list_get(target_index, mapping_batch->CT_fq_batch);
    list = mapping_batch->mapping_lists2[target_index];

    if (Ncg <= margen) {
      max_seeds = (read->length / 15) * 2 + 10;
      num_cals = array_list_size(list);

      //Ricardo - filtrado cals
      list = filter_cals(num_cals, read->length, list);

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
