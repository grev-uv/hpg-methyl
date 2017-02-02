#include "bwt_server.h"

#define SIMPLE_SW 1
#define SP_SW 2

#define SW_NORMAL 0 
#define SW_FINAL 1 

#define MIN_ANCHOR 25
#define CAL_FACTOR 10
#define MIN_CAL_DEFINED 30
#define MAX_CAL_DEFINED 200


//Maximum number of anchors allowed for a read
#define MAX_NUM_ANCHOR 512 //1024;

//#define MAX_BWT_REGIONS 100	//100
#define MAX_BWT_ANCHOR_DISTANCE 500000
#define MAX_BIG_ANCHOR_SIZE 20 //16

//====================================================================================
// bwt_server_input functions: init
//====================================================================================

void bwt_server_input_init(list_t* read_list_p, unsigned int batch_size, bwt_optarg_t *bwt_optarg_p, 
			   bwt_index_t *bwt_index_p, list_t* write_list_p, unsigned int write_size, 
			   list_t* unmapped_read_list_p, metaexons_t *metaexons, sw_optarg_t *sw_optarg,
			   genome_t *genome, bwt_server_input_t* input_p) {
  input_p->read_list_p = read_list_p;
  input_p->batch_size = batch_size;
  input_p->bwt_optarg_p = bwt_optarg_p;
  input_p->write_list_p = write_list_p;
  input_p->write_size = write_size;
  input_p->bwt_index_p = bwt_index_p;
  input_p->unmapped_read_list_p = unmapped_read_list_p;
  input_p->metaexons = metaexons;
  input_p->sw_optarg = sw_optarg;
  input_p->genome = genome;
}

//====================================================================================
// apply_bwt
//====================================================================================

cal_t *convert_bwt_anchor_to_CAL(bwt_anchor_t *bwt_anchor, size_t read_start, size_t read_end) {
  linked_list_t *linked_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
  seed_region_t *seed_region = seed_region_new(read_start, read_end,
					       bwt_anchor->start, bwt_anchor->end, 0);

  linked_list_insert_first(seed_region, linked_list);

  cal_t *cal = cal_new(bwt_anchor->chromosome + 1, bwt_anchor->strand,
		       bwt_anchor->start, bwt_anchor->end,
		       1, linked_list,
		       linked_list_new(COLLECTION_MODE_ASYNCHRONIZED));

  return cal;
}

//------------------------------------------------------------------------------------

size_t bwt_search_pair_anchors(array_list_t *list, unsigned int read_length, size_t min_cal_size, double umbral_cal_length_factor, int min_read_discard) {
  bwt_anchor_t *bwt_anchor;
  bwt_anchor_t *max_anchor = NULL;
  bwt_anchor_t *bwt_anchor_back, *bwt_anchor_forw;

	int max_double_anchor = 0, max_anchor_length = 0;
  int anchor_length_tmp, anchor_back, anchor_forw;
  int strand, type;
  int found_anchor = 0, found_double_anchor = 0;
	int seed_size, gap_read, gap_genome;
  int MAX_BWT_REGIONS = array_list_size(list);

  array_list_t *anchor_list_tmp, *forward_anchor_list, *backward_anchor_list;
  cal_t *cal;

  if (MAX_BWT_REGIONS > MAX_NUM_ANCHOR) {
	  MAX_BWT_REGIONS = MAX_NUM_ANCHOR;
  }

  array_list_t *backward_anchor_list_0 = array_list_new(MAX_BWT_REGIONS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *forward_anchor_list_0 = array_list_new(MAX_BWT_REGIONS, 1.25f , COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *backward_anchor_list_1 = array_list_new(MAX_BWT_REGIONS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *forward_anchor_list_1 = array_list_new(MAX_BWT_REGIONS, 1.25f , COLLECTION_MODE_ASYNCHRONIZED);

  array_list_t *big_anchor_list = array_list_new(MAX_BWT_REGIONS, 1.25f , COLLECTION_MODE_ASYNCHRONIZED);


  bool discard = false;

  for (int i = 0; i < array_list_size(list); i++) {
    discard = false;
		bwt_anchor = array_list_get(i, list);
    anchor_length_tmp = bwt_anchor->end - bwt_anchor->start + 1;

    if ((read_length >= min_read_discard)  && (anchor_length_tmp <=  min_cal_size)) {
    	discard = true;
		}

    if ((array_list_size(forward_anchor_list_1) >= MAX_BWT_REGIONS)  || 
				(array_list_size(forward_anchor_list_0) >= MAX_BWT_REGIONS)  ||
    		(array_list_size(backward_anchor_list_1) >= MAX_BWT_REGIONS) ||
				(array_list_size(backward_anchor_list_0) >= MAX_BWT_REGIONS)) {
    	discard = true;
    	printf("Discarding Anchor - Increase MAX_NUM_ANCHOR variable in bwt_server_cpu.c\n");
    }

		// - Who: Ricardo
    if (!discard)	{
    	if (bwt_anchor->strand == 1) {
        if (bwt_anchor->type == FORWARD_ANCHOR) {
    			array_list_insert(bwt_anchor, forward_anchor_list_1);
    		} else {
    			array_list_insert(bwt_anchor, backward_anchor_list_1);
    		}
    	} else {
        if (bwt_anchor->type == FORWARD_ANCHOR) {
    			array_list_insert(bwt_anchor, forward_anchor_list_0);
    		} else {
    			array_list_insert(bwt_anchor, backward_anchor_list_0);
    		}
    	}

    	anchor_length_tmp = bwt_anchor->end - bwt_anchor->start + 1;
    	found_anchor = 1;
			strand = bwt_anchor->strand;
			type = bwt_anchor->type;

			if (read_length - anchor_length_tmp < MAX_BIG_ANCHOR_SIZE) { //16 --> 20
    		array_list_insert(bwt_anchor, big_anchor_list);
    	}
    } else {
			// Borrar bwt_anchor
    	bwt_anchor_free(bwt_anchor);
    }

  }

  array_list_clear(list, NULL);

  if (array_list_size(big_anchor_list) > 0) {
    for (int i = array_list_size(big_anchor_list) - 1; i >= 0; i--) {
      bwt_anchor = array_list_remove_at(i, big_anchor_list);
      size_t seed_size = bwt_anchor->end - bwt_anchor->start;

      if (bwt_anchor->type == FORWARD_ANCHOR) {
				// Now it starts in "seq_start" not in 0
				// - Who: Ricardo
    	  cal = convert_bwt_anchor_to_CAL(bwt_anchor, bwt_anchor->seq_start, bwt_anchor->seq_start + seed_size);
	    } else {
				// Now it starts in "seq_end" not in read_length-1 
				// - Who: Ricardo
				cal = convert_bwt_anchor_to_CAL(bwt_anchor, bwt_anchor->seq_end - seed_size, bwt_anchor->seq_end); 
     }

      array_list_insert(cal, list);
    }

    array_list_set_flag(SINGLE_ANCHORS, list);

		//Big anchor is almost a matching, so it only will be one if exists
    goto exit; 
  }
 
  if (!found_double_anchor && found_anchor) {
		const int nstrands = 2;
		const int nchromosomes = 30; //TODO: Parameter

		linked_list_t ***new_cals_list = (linked_list_t ***)malloc(sizeof(linked_list_t **)*nstrands);

		for (unsigned int ii = 0; ii < nstrands; ii++) {
	    new_cals_list[ii] = (linked_list_t **)malloc(sizeof(linked_list_t *)*nchromosomes);

	    for (unsigned int jj = 0; jj < nchromosomes; jj++) {
	       new_cals_list[ii][jj] = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	    }
		}

		// We pass a unique id for identify anchor
		// - Who: Ricardo
		int start_id = 0;
		int suma;
	
		if ((suma = array_list_size(forward_anchor_list_1)) > 0) {
			insert_seeds_and_merge_anchor(forward_anchor_list_1, new_cals_list,  read_length, start_id);
			start_id += suma;
		}

		if ((suma = array_list_size(backward_anchor_list_1)) > 0) {
			insert_seeds_and_merge_anchor(backward_anchor_list_1, new_cals_list,  read_length, start_id);
			start_id += suma;
		}

		if ((suma = array_list_size(forward_anchor_list_0)) > 0) {
			insert_seeds_and_merge_anchor(forward_anchor_list_0, new_cals_list,  read_length, start_id);
			start_id += suma;
		}
		
		if ((suma = array_list_size(backward_anchor_list_0)) > 0) {
			insert_seeds_and_merge_anchor(backward_anchor_list_0, new_cals_list,  read_length, start_id);
			start_id += suma;
		}

	 	// Store CALs in Array List for return results
	  size_t start_cal, end_cal;
	  size_t seq_start, seq_end;
	  seed_region_t *s, *s_first, *s_last, *s_aux, *seed_region;
	  cal_t *cal;
	  short_cal_t *short_cal, *short_cal_aux;
	  linked_list_iterator_t itr, itr2, itr3;
	  linked_list_item_t *list_item_cal, *list_item_aux;

	  int best_cal_num_seeds = 0;
	  int best_cal_length = 0;
	  int actual_cal_length = 0;
	  int cals_added = 0;

	  // Cal is accepted like candidate if they are the largest cal mapped or if it has been mapped more than umbral_cal_length
	  int umbral_cal_length = read_length / umbral_cal_length_factor;		//2

	  const int max_intron_size = 500000; //TODO: Parameter
	 
	  for (unsigned int jj = 0; jj < nchromosomes; jj++) {
	    for (unsigned int ii = 0; ii < nstrands; ii++) {
	      linked_list_iterator_init(new_cals_list[ii][jj], &itr);
	      list_item_cal = linked_list_iterator_list_item_curr(&itr);

	      while ((list_item_cal != NULL )) {
					short_cal = (short_cal_t *)list_item_cal->item;

		  		int actual_cal_length = 0;
		  		linked_list_iterator_init(short_cal->sr_list, &itr3);
		  		s = linked_list_iterator_curr(&itr3);
 	      	
					while (s != NULL) {
 	    	  	actual_cal_length += (s->genome_end - s->genome_start);
		  	  	s =  	     	      linked_list_iterator_next(&itr3);
 	      	}

 	      	if ((actual_cal_length >= best_cal_length) || (actual_cal_length >= umbral_cal_length)) {
 	    	 		linked_list_t *list_aux = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);

 	    	  	while (s = (seed_region_t *)linked_list_remove_last(short_cal->sr_list)) {
 	    		  	//TODO: Change all parameters to seed_region_t
 	    		  	append_seed_region_linked_list(list_aux,
						   				s->read_start, s->read_end,
						   				s->genome_start, s->genome_end,
						   				s->id);

 	    		  	seed_region_free(s);
 	    	  	}
						
						//Delete previous added cals
 	    	  	if ((actual_cal_length > best_cal_length)) {
 	    		  	for (int kk=0; kk<array_list_size(list);kk++) {
 	    			  	cal = array_list_get(kk, list);
 	    			  	
								if ( (cal->end - cal->start + 1 < umbral_cal_length)) {
 	    				  	cal = (cal_t *) array_list_remove_at(kk, list);
 	    				  	cal_free(cal);
 	    			  	}
 	    		  	}
 	    	  	}

 	    	  	cal = cal_new(jj+1, ii, short_cal->start, short_cal->end, short_cal->num_seeds, 
					 								list_aux, short_cal->sr_duplicate_list);

 	    	  	array_list_insert(cal, list);

 	    	    if (actual_cal_length > best_cal_length) {
 	    		  	best_cal_length = actual_cal_length;
						}

 	    	 		short_cal->sr_duplicate_list = NULL;
					} else {
						// Borrar
						while (s = (seed_region_t *)linked_list_remove_last(short_cal->sr_list)) {
							seed_region_free(s);
						}
					}

					linked_list_iterator_next(&itr);
					list_item_cal = linked_list_iterator_list_item_curr(&itr);
				}
			}
		}

		//Ricardo - now can be multiple anchors, not single one
		array_list_set_flag(MULTIPLE_ANCHORS, list);			

		for (unsigned int i = 0; i < nstrands; i++) {
			for (unsigned int j = 0; j < nchromosomes; j++) {
				linked_list_free(new_cals_list[i][j], (void *)short_cal_free);
			}
			free(new_cals_list[i]);
		}

		free(new_cals_list);
	}

 exit:
  array_list_free(forward_anchor_list_1, (void *)bwt_anchor_free);
  array_list_free(backward_anchor_list_1,  (void *)bwt_anchor_free);
  array_list_free(forward_anchor_list_0,  (void *)bwt_anchor_free);
  array_list_free(backward_anchor_list_0,  (void *)bwt_anchor_free);
  array_list_free(big_anchor_list,  (void *)bwt_anchor_free);

  return array_list_size(list);
}

//====================================================================================
// BS - Burrows-Wheeler Transform
//====================================================================================

int apply_bwt_bs(bwt_server_input_t* input, batch_t *batch) {
  LOG_DEBUG("========= APPLY BWT BS START =========\n");

	mapping_batch_t *mapping_batch = batch->mapping_batch;
  struct timeval start, end;
  double time, total_time;
  
  // Copy the batch reads
  size_t num_reads = array_list_size(mapping_batch->fq_batch);

  // Intialize the new mappings and indices
	mapping_batch->mapping_lists2 = (array_list_t **) calloc(num_reads, sizeof(array_list_t*));

  for (size_t i = 0; i < num_reads; i++) {
    mapping_batch->mapping_lists2[i] = array_list_new(500, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  }

  mapping_batch->CT_fq_batch     = array_list_new(num_reads + 2, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  mapping_batch->CT_rev_fq_batch = array_list_new(num_reads + 2, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  mapping_batch->GA_fq_batch     = array_list_new(num_reads + 2, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  mapping_batch->GA_rev_fq_batch = array_list_new(num_reads + 2, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  // Copy and transform the reads simultaneously
  cpy_transform_array_bs(mapping_batch->fq_batch, 
			 mapping_batch->CT_fq_batch, mapping_batch->CT_rev_fq_batch, 
			 mapping_batch->GA_fq_batch, mapping_batch->GA_rev_fq_batch);

	// Make the four searches
  alignment_t *alignment;
  size_t header_len;

  size_t num_threads = input->bwt_optarg_p->num_threads;
  num_reads = array_list_size(mapping_batch->fq_batch);

  fastq_read_t* fq_read;
  mapping_batch->num_targets = 0;

  float Ncg, Ngc;
  float margen = mapping_batch->margin;
	histogram_input_t hist_input;

	struct timeval histTimeStart, histTimeEnd;
	float histTime = 0.0f, histIterTime = 0.0f;

  for (size_t i = 0; i < num_reads; i++) {
    // Obtain histogram of each read to filter the number of searches realized
    LOG_DEBUG_F("========= OBTAIN HISTOGRAM OF READ %lu =========\n", i);

    fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->fq_batch);
		if (time_on) {
			start_timer(histTimeStart);
		}

		// Obtain the histogram
		histogram_input_init(&hist_input, fq_read);
		histogram_apply(&hist_input);
		mapping_batch->histogram_sw[i] = hist_input.out_ncg;

		if (time_on) {
			stop_timer(histTimeStart, histTimeEnd, histIterTime);
			start_timer(start);
			histTime += histIterTime;
		}

    Ncg = hist_input.out_ncg;
    Ngc = hist_input.out_ngc;
    
    array_list_set_flag(1, mapping_batch->mapping_lists[i]);
    array_list_set_flag(1, mapping_batch->mapping_lists2[i]);

    int readStart = 0;
    int readEnd = 0; //399
    interval_t* interval = NULL;

    // If not defined min_cal_size, use one adapted to the read length
    int MIN_SINGLE_ANCHOR = 0;

    if (input->bwt_optarg_p->max_inner_gap <= 0) {
    	MIN_SINGLE_ANCHOR = fq_read->length / CAL_FACTOR;

    	if (MIN_SINGLE_ANCHOR < MIN_CAL_DEFINED) {
    		MIN_SINGLE_ANCHOR = MIN_CAL_DEFINED;
    	} else {
    		if (MIN_SINGLE_ANCHOR > MAX_CAL_DEFINED) {
    			MIN_SINGLE_ANCHOR = MAX_CAL_DEFINED;
    		}
			}
    } else {
    	MIN_SINGLE_ANCHOR = input->bwt_optarg_p->max_inner_gap;
    }

    //Min_cal configure, at the moment same parameters than MIN_SINGLE_ANCHOR
    int MIN_CAL_SIZE;

    if (input->bwt_optarg_p->min_cal_size <= 0) {
			MIN_CAL_SIZE = fq_read->length / CAL_FACTOR;
			
			if (MIN_CAL_SIZE < MIN_CAL_DEFINED) {
				MIN_CAL_SIZE = MIN_CAL_DEFINED;
			} else {
				if (MIN_CAL_SIZE > MAX_CAL_DEFINED) {
					MIN_CAL_SIZE = MAX_CAL_DEFINED;
				}
			}
		} else {
			MIN_CAL_SIZE = input->bwt_optarg_p->min_cal_size;
		}

    int intervalUmbral = MIN_SINGLE_ANCHOR;

    if (Ngc <= margen) {
      // First search the reverse of the G->A transformation
      fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->GA_rev_fq_batch);

      readStart = 0;
      readEnd = fq_read->length - 1; //399

			// BWT with interior intervals
      // - Who: Ricardo
      do {
				interval = bwt_map_inexact_read_bs_new(fq_read,
							input->bwt_optarg_p, input->bwt_index2_p,
							mapping_batch->mapping_lists[i], 1, readStart, readEnd);
				
				// TODO: Purge
				readStart = interval->start + 1;	//Continues after the error
				readEnd = interval->end - 1;		//Continues after the error

      	free(interval);
      }
      while ((mapping_batch->mapping_lists[i]->flag == BWT_BS_MATCHING_SEED) && 
						((readEnd - readStart) > intervalUmbral));

      if (array_list_get_flag(mapping_batch->mapping_lists[i]) != BWT_BS_MATCHING_EXCEEDED) {
    	  readStart = 0;
    	  readEnd = fq_read->length - 1; //399

				do {
					// Next search the direct of the G->A transformation
					fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->GA_fq_batch);
					interval = bwt_map_inexact_read_bs_new(fq_read,
        			  input->bwt_optarg_p, input->bwt_index_p,
							  mapping_batch->mapping_lists[i], 0, readStart, readEnd);

					readStart = interval->start + 1;	//Continues after the error
					readEnd = interval->end - 1;		//Continues after the error
					free(interval);
				}
        while ((mapping_batch->mapping_lists[i]->flag == BWT_BS_MATCHING_SEED) && 
							((readEnd - readStart) > intervalUmbral));
      }
    }

    if (Ncg <= margen) {
      // First search the reverse of the C->T transformation
      fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->CT_rev_fq_batch);

      readStart = 0;
      readEnd = fq_read->length - 1; //399

      do {
    	  interval = bwt_map_inexact_read_bs_new(fq_read,
			      input->bwt_optarg_p, input->bwt_index_p,
			      mapping_batch->mapping_lists2[i], 1, readStart, readEnd);

    	  readStart = interval->start + 1;	//Continues after the error
    	  readEnd = interval->end - 1;		//Continues after the error
    	  free(interval);  
      } while ((mapping_batch->mapping_lists2[i]->flag == BWT_BS_MATCHING_SEED) && 
							((readEnd - readStart) > intervalUmbral));

      if (array_list_get_flag(mapping_batch->mapping_lists2[i]) != BWT_BS_MATCHING_EXCEEDED) {
    	  readStart = 0;
    	  readEnd = fq_read->length - 1; //399

    	  do {
    		  // Next search the direct of the C->T transformation
    		  fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->CT_fq_batch);
    		  interval = bwt_map_inexact_read_bs_new(fq_read,
    				 		 							input->bwt_optarg_p, input->bwt_index2_p,
					  									mapping_batch->mapping_lists2[i], 0, readStart, readEnd);
    		  
    		  readStart = interval->start + 1;	//Continues after the error
    		  readEnd = interval->end - 1;		//Continues after the error
    		  free(interval);
				} while ((mapping_batch->mapping_lists2[i]->flag == BWT_BS_MATCHING_SEED) && 
								((readEnd - readStart) > intervalUmbral));
      }
    }

		// Post-process alignment
    if ((array_list_get_flag(mapping_batch->mapping_lists[i]) == BWT_BS_MATCHING_EXACT  && 
	 			 mapping_batch->mapping_lists[i]->size)										  										|| 
				(array_list_get_flag(mapping_batch->mapping_lists2[i]) == BWT_BS_MATCHING_EXACT && 
	 			 mapping_batch->mapping_lists2[i]->size)) {
      if (array_list_get_flag(mapping_batch->mapping_lists[i]) == BWT_BS_MATCHING_SEED) {
				array_list_clear( mapping_batch->mapping_lists[i], (void *)bwt_anchor_free);
      } else if (array_list_get_flag(mapping_batch->mapping_lists2[i]) == BWT_BS_MATCHING_SEED) {
				array_list_clear(mapping_batch->mapping_lists2[i], (void *)bwt_anchor_free);
      }

      array_list_set_flag(ALIGNMENTS_FOUND, mapping_batch->mapping_lists[i]);
      array_list_set_flag(ALIGNMENTS_FOUND, mapping_batch->mapping_lists2[i]);
    } else {
      if (array_list_get_flag(mapping_batch->mapping_lists[i]) != BWT_BS_MATCHING_EXCEEDED &&
	  			array_list_get_flag(mapping_batch->mapping_lists2[i]) != BWT_BS_MATCHING_EXCEEDED) {
				if (mapping_batch->mapping_lists[i]->size) {
	  			bwt_search_pair_anchors(mapping_batch->mapping_lists[i], fq_read->length, 
											MIN_CAL_SIZE, input->bwt_optarg_p->umbral_cal_length_factor, 
											input->bwt_optarg_p->min_read_discard); 
				}

				if (!mapping_batch->mapping_lists[i]->size) {
	  			array_list_clear( mapping_batch->mapping_lists[i], (void *)bwt_anchor_free);
	  			array_list_set_flag(NOT_ANCHORS, mapping_batch->mapping_lists[i]);	  
				}

				if (mapping_batch->mapping_lists2[i]->size) {
					bwt_search_pair_anchors(mapping_batch->mapping_lists2[i], fq_read->length, MIN_CAL_SIZE, 
											input->bwt_optarg_p->umbral_cal_length_factor, 
											input->bwt_optarg_p->min_read_discard); 
				}

				if (!mapping_batch->mapping_lists2[i]->size) {
	  			array_list_clear( mapping_batch->mapping_lists2[i], (void *)bwt_anchor_free);
	  			array_list_set_flag(NOT_ANCHORS, mapping_batch->mapping_lists2[i]);	  
				}

				if ((mapping_batch->mapping_lists[i]->size) || (mapping_batch->mapping_lists2[i]->size)) {
					mapping_batch->targets[(mapping_batch->num_targets)++] = i;
				}
      } else { 
		    array_list_clear( mapping_batch->mapping_lists[i], (void *)bwt_anchor_free);
   		  array_list_clear( mapping_batch->mapping_lists2[i], (void *)bwt_anchor_free);
				array_list_set_flag(ALIGNMENTS_EXCEEDED, mapping_batch->mapping_lists[i]);
				array_list_set_flag(ALIGNMENTS_EXCEEDED, mapping_batch->mapping_lists2[i]);
      }
    }

		if (time_on) {
			stop_timer(start, end, time);
			total_time += time;
		}
  }

  size_t num_mapped_reads = array_list_size(mapping_batch->fq_batch) - mapping_batch->num_targets;
  mapping_batch->num_to_do = num_mapped_reads;

  if (time_on) {
		timing_add(time, BWT_SERVER, timing);
		timing_add(histTime, BWT_SERVER_HISTOGRAM, timing);
	}

  LOG_DEBUG("========= APPLY BWT BS END =========\n");

  if (batch->mapping_batch->num_targets > 0) {
    for (int i = 0; i < batch->mapping_batch->num_targets; i++) {
      batch->mapping_batch->bwt_mappings[batch->mapping_batch->targets[i]] = 1;
    }
		
    return BS_CAL_STAGE;
  }

  return BS_POST_PAIR_STAGE;
}
