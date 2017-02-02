#include "sw_server.h"


//------------------------------------------------------------------------------------

void seed_region_free2(seed_region_t *seed_region) {
  if (seed_region) {
    if (seed_region->info) {
      cigar_code_free(seed_region->info);
    }

    free(seed_region);
  }
}

//------------------------------------------------------------------------------------

void cal_free2(cal_t *cal) {
  if (cal) {
    if (cal->sr_list) {
      linked_list_free(cal->sr_list, (void *) seed_region_free2);
    }

    if (cal->sr_duplicate_list) {
      linked_list_free(cal->sr_duplicate_list, (void *) seed_region_free2);
    }

    if (cal->candidates_seeds_start) {
      array_list_free(cal->candidates_seeds_start, (void *)seed_region_free2);
    }

    if (cal->candidates_seeds_end) {
      array_list_free(cal->candidates_seeds_end, (void *)seed_region_free2);
    }

    if (cal->info) {
      cigar_code_free(cal->info);
    }

    free(cal);
  }
}

//====================================================================================
//  Input structure for Smith-Waterman server
//====================================================================================

void sw_optarg_init(float gap_open, float gap_extend, 
		    float match, float mismatch, sw_optarg_t *sw_optarg) {
  sw_optarg->gap_open = gap_open;
  sw_optarg->gap_extend = gap_extend;

  sw_optarg->subst_matrix['A']['A'] = match;
  sw_optarg->subst_matrix['C']['A'] = mismatch;
  sw_optarg->subst_matrix['T']['A'] = mismatch;
  sw_optarg->subst_matrix['G']['A'] = mismatch;
  sw_optarg->subst_matrix['N']['A'] = mismatch;

  sw_optarg->subst_matrix['A']['C'] = mismatch;
  sw_optarg->subst_matrix['C']['C'] = match;
  sw_optarg->subst_matrix['T']['C'] = mismatch;
  sw_optarg->subst_matrix['G']['C'] = mismatch;
  sw_optarg->subst_matrix['N']['C'] = mismatch;

  sw_optarg->subst_matrix['A']['T'] = mismatch;
  sw_optarg->subst_matrix['C']['T'] = mismatch;
  sw_optarg->subst_matrix['T']['T'] = match;
  sw_optarg->subst_matrix['G']['T'] = mismatch;
  sw_optarg->subst_matrix['N']['T'] = mismatch;

  sw_optarg->subst_matrix['A']['G'] = mismatch;
  sw_optarg->subst_matrix['C']['G'] = mismatch;
  sw_optarg->subst_matrix['T']['G'] = mismatch;
  sw_optarg->subst_matrix['G']['G'] = match;
  sw_optarg->subst_matrix['N']['G'] = mismatch;

  sw_optarg->subst_matrix['A']['N'] = mismatch;
  sw_optarg->subst_matrix['C']['N'] = mismatch;
  sw_optarg->subst_matrix['T']['N'] = mismatch;
  sw_optarg->subst_matrix['G']['N'] = mismatch;
  sw_optarg->subst_matrix['N']['N'] = match;
}

//------------------------------------------------------------------------------------

void sw_server_input_init(list_t* sw_list, list_t* alignment_list, unsigned int write_size, 
			  float match, float mismatch, float gap_open, float gap_extend, 
			  float min_score, unsigned int flank_length, genome_t* genome, 
			  size_t max_intron_size, int min_intron_size, 
			  size_t seed_max_distance, bwt_optarg_t* bwt_optarg_p, 
			  cal_optarg_t *cal_optarg_p, bwt_index_t *bwt_index_p,
			  metaexons_t *metaexons, linked_list_t *buffer, 
			  linked_list_t *buffer_hc, FILE *f_sa, FILE *f_hc,
			  int pair_mode, sw_server_input_t* input) { 
  input->sw_list_p = sw_list;
  input->alignment_list_p = alignment_list;
  input->write_size = write_size;
  input->genome_p = genome;
  input->max_intron_size = max_intron_size;
  input->min_intron_size = min_intron_size;
  input->seed_max_distance = seed_max_distance;
  input->bwt_optarg_p =  bwt_optarg_p; 

  // Smith-Waterman parameters
  sw_optarg_init(gap_open, gap_extend, 
		 match, mismatch, &input->sw_optarg);
  input->match = match;
  input->mismatch = mismatch;
  input->gap_open = gap_open;
  input->gap_extend = gap_extend;
  input->min_score = min_score;

  // CAL
  input->flank_length = flank_length;

  input->cal_optarg_p = cal_optarg_p;
  input->bwt_index_p = bwt_index_p;
  input->metaexons = metaexons;

  input->buffer = buffer;
  input->buffer_hc = buffer_hc;

  input->f_sa = f_sa;
  input->f_hc = f_hc;

  input->pair_mode = pair_mode;
}

//====================================================================================
// main sw function
//====================================================================================

int apply_sw_bs(sw_server_input_t* input, batch_t *batch) {
  LOG_DEBUG("========= APPLY SW BS START =========\n");

  // Should this be an option?
  int sw_3_nucleotides = 0;
  
  struct timeval sw_time_start, sw_time_end;
  float sw_time = 0.0f;

  if (time_on) {
    start_timer(sw_time_start);
  }

  if (sw_3_nucleotides) {
    apply_sw_bs_3nt(input, batch);
  } else {
    apply_sw_bs_4nt(input, batch);
  }

  if (time_on) {
    stop_timer(sw_time_start, sw_time_end, sw_time);
    timing_add(sw_time, SW_STAGE_TIME, timing);
  }
  
  LOG_DEBUG("========= END OF APPLY SW BS =========\n");
  return BS_POST_PAIR_STAGE;
}

//--------------------------------------------------------------------------------------

void fill_matrix(subst_matrix_t subst_matrix, float match, float mismatch, int type, float factor_match, float factor_mismatch) {
  subst_matrix['A']['A'] = match;
  subst_matrix['C']['A'] = mismatch;
  subst_matrix['G']['A'] = mismatch;
  subst_matrix['T']['A'] = mismatch;
  subst_matrix['N']['A'] = mismatch;

  subst_matrix['A']['C'] = mismatch;
  subst_matrix['C']['C'] = match;
  subst_matrix['G']['C'] = mismatch;
  subst_matrix['T']['C'] = mismatch;
  subst_matrix['N']['C'] = mismatch;

  subst_matrix['A']['T'] = mismatch;
  subst_matrix['C']['T'] = mismatch;
  subst_matrix['G']['T'] = mismatch;
  subst_matrix['T']['T'] = match;
  subst_matrix['N']['T'] = mismatch;

  subst_matrix['A']['G'] = mismatch;
  subst_matrix['C']['G'] = mismatch;
  subst_matrix['G']['G'] = match;
  subst_matrix['T']['G'] = mismatch;
  subst_matrix['N']['G'] = mismatch;

  subst_matrix['A']['N'] = mismatch;
  subst_matrix['C']['N'] = mismatch;
  subst_matrix['T']['N'] = mismatch;
  subst_matrix['G']['N'] = mismatch;
  subst_matrix['N']['N'] = match;

  float x = 500;
  float y = -400;
  float z = 2.5;
  float w = -2;

  if (type == 1) {
    subst_matrix['C']['A'] = y;
    subst_matrix['C']['C'] = x;
    subst_matrix['C']['G'] = y;
    subst_matrix['C']['T'] = y;

    subst_matrix['T']['A'] = w;
    subst_matrix['T']['C'] = z;
    subst_matrix['T']['G'] = w;
    subst_matrix['T']['T'] = z;

  } else {
    subst_matrix['G']['T'] = y;
    subst_matrix['G']['G'] = x;
    subst_matrix['G']['C'] = y;
    subst_matrix['G']['A'] = y;

    subst_matrix['A']['T'] = w;
    subst_matrix['A']['G'] = z;
    subst_matrix['A']['C'] = w;
    subst_matrix['A']['A'] = z;
  }
}

//--------------------------------------------------------------------------------------

void apply_sw_bs_3nt(sw_server_input_t* input, batch_t *batch) {
  int tid = omp_get_thread_num();

  mapping_batch_t *mapping_batch = batch->mapping_batch;
  cal_t *cal = NULL;
  array_list_t *cal_list = NULL, *mapping_list = NULL;
    
  array_list_t *fq_batch = mapping_batch->fq_batch;
  fastq_read_t *fq_read;
    
  // added by PP for bisulfite
  array_list_t *CT_fq_batch = mapping_batch->CT_fq_batch;
  array_list_t *GA_fq_batch = mapping_batch->GA_fq_batch;
  array_list_t *CT_rev_fq_batch = mapping_batch->CT_rev_fq_batch;
  array_list_t *GA_rev_fq_batch = mapping_batch->GA_rev_fq_batch;
  fastq_read_t *fq_read2;
  // end added by PP for bisulfite
    
  size_t start, end;
  size_t start2, end2;

  // added by PP for bisulfite
  genome_t *genome1 = input->genome1_p;
  genome_t *genome2 = input->genome2_p;
  // end added by PP for bisulfite
    
  size_t flank_length = input->flank_length;
    
  // SIMD support for Smith-Waterman
  float score, min_score = input->min_score;  
  sw_output_t *sw_output;
  size_t read_index, num_cals;
    
  size_t num_targets = mapping_batch->num_targets;
  size_t new_num_targets = 0;

  // added by PP for bisulfite
  size_t num_targets2 = mapping_batch->num_targets2;
  size_t new_num_targets2 = 0;
  // added by PP for bisulfite
    
  // added by PP for bisulfite
  size_t sw_total1 = mapping_batch->num_to_do;
  size_t sw_total2 = mapping_batch->num_to_do2;
  size_t sw_total = sw_total1 + sw_total2;
  // end added by PP for bisulfite
    
  // set to zero
  mapping_batch->num_to_do = 0;

  // added by PP for bisulfite
  mapping_batch->num_to_do2 = 0;
  int g[sw_total];
  // end added by PP for bisulfite
    
  sw_optarg_t *sw_optarg = &input->sw_optarg;
  sw_multi_output_t *output = sw_multi_output_new(sw_total);
  char *q[sw_total], *r[sw_total];
  uint8_t strands[sw_total], chromosomes[sw_total];
  size_t starts[sw_total];
  size_t sw_count = 0, read_indices[sw_total], sw_count2 = 0;
  int read_len, ref_len, max_ref_len;
    
  
  // initialize query and reference sequences to Smith-Waterman
  for (size_t i = 0; i < num_targets; i++) {
    read_index = mapping_batch->targets[i];
      
    // to use with the three nucleotides searches
    fq_read  = (fastq_read_t *) array_list_get(read_index, GA_fq_batch);
    fq_read2 = (fastq_read_t *) array_list_get(read_index, GA_rev_fq_batch);
  
    cal_list = mapping_batch->mapping_lists[read_index];
    num_cals = array_list_size(cal_list);
      
    read_len = fq_read->length;
    
    // processing each CAL from this read
    for (size_t j = 0; j < num_cals; j++) {
      // get cal and read index
      cal = array_list_get(j, cal_list);
      read_indices[sw_count] = read_index;
	
      if (flank_length >= cal->start) {
	      start = 0;
      } else {
	      start = cal->start - flank_length;
      }
	
      end = cal->end + flank_length;

      if (end >= genome1->chr_size[cal->chromosome_id - 1]) {
	      end = genome1->chr_size[cal->chromosome_id - 1] - 1;
      }
	
      ref_len = end - start + 2;
	
      // query sequence, revcomp if necessary
      q[sw_count] = (char *) calloc((read_len + 1), sizeof(char));
	
      // to use with the three nucleotides searches
      if (cal->strand == 0) {
	      memcpy(q[sw_count], fq_read->sequence, read_len);
      } else {
	      memcpy(q[sw_count], fq_read2->sequence, read_len);
      }
	
      // reference sequence
      r[sw_count] = calloc(1, end - start + 2);
	
      // to use with the three nucleotides searches
      if (cal->strand == 0) {
	      genome_read_sequence_by_chr_index(r[sw_count], 0, cal->chromosome_id - 1, &start, 
                                          &end, genome1);
      } else {
	      genome_read_sequence_by_chr_index(r[sw_count], 0, cal->chromosome_id - 1, &start, 
                                          &end, genome2);
      }

      // save some stuff, we'll use them after...
      strands[sw_count] = cal->strand;
      chromosomes[sw_count] = cal->chromosome_id;
      starts[sw_count] = start;

      // increase counter
      sw_count++;
    }
      
    // free cal_list
    array_list_clear(cal_list, (void *) cal_free2);
  }
  
  sw_count2 = sw_count;
    
  for (size_t i = 0; i < num_targets2; i++) {
    read_index = mapping_batch->targets2[i];
      
    // to use with the three nucleotides searches
    fq_read  = (fastq_read_t *) array_list_get(read_index, CT_fq_batch);
    fq_read2 = (fastq_read_t *) array_list_get(read_index, CT_rev_fq_batch);  
    cal_list = mapping_batch->mapping_lists2[read_index];
    num_cals = array_list_size(cal_list);     
    read_len = fq_read->length;

    // processing each CAL from this read
    for(size_t j = 0; j < num_cals; j++) {
	
      // get cal and read index
      cal = array_list_get(j, cal_list);
      read_indices[sw_count] = read_index;
	
      if (flank_length >= cal->start) {
	      start = 0;
      } else {
	      start = cal->start - flank_length;
      }
	
      end = cal->end + flank_length;

      if (end >= genome1->chr_size[cal->chromosome_id - 1]) {
	      end = genome1->chr_size[cal->chromosome_id - 1] - 1;
      }
	
      ref_len = end - start + 2;

      // query sequence, revcomp if necessary
      q[sw_count] = (char *) calloc((read_len + 1), sizeof(char));
	
      // to use with the three nucleotides searches
      if (cal->strand == 0) {
	      memcpy(q[sw_count], fq_read->sequence, read_len);
      } else {
	      memcpy(q[sw_count], fq_read2->sequence, read_len);
      }

      // reference sequence
      r[sw_count] = calloc(1, end - start + 2);
	
      // to use with the three nucleotides searches
      if (cal->strand == 0) {
	      genome_read_sequence_by_chr_index(r[sw_count], 0, cal->chromosome_id - 1, &start, 
                                          &end, genome2);
      } else {
	      genome_read_sequence_by_chr_index(r[sw_count], 0, cal->chromosome_id - 1, &start, 
                                          &end, genome1);
      }

      // save some stuff, we'll use them after...
      strands[sw_count] = cal->strand;
      chromosomes[sw_count] = cal->chromosome_id;
      starts[sw_count] = start;
	
      // increase counter
      sw_count++;
    }
      
    // free cal_list
    array_list_clear(cal_list, (void *) cal_free2);
  }
    
  // run Smith-Waterman
  smith_waterman_mqmr(q, r, sw_count, sw_optarg, 1, output);
    
  for (size_t i = 0; i < sw_count; i++) {
    LOG_DEBUG_F("cal: start = %lu, strand = %i\n", starts[i], strands[i]);
    LOG_DEBUG_F("\tquery : %s\n", q[i]); 
    LOG_DEBUG_F("\tref.  : %s\n", r[i]); 
    LOG_DEBUG_F("\tquery map: %s (start: %i)\n", output->query_map_p[i], 
                output->query_start_p[i]);
    LOG_DEBUG_F("\tref. map : %s (start: %i)\n", output->ref_map_p[i],
                output->ref_start_p[i]);
    LOG_DEBUG("\n");
  }
    
  double norm_score;
  
  for (size_t i = 0; i < sw_count2; i++) { 
    read_index = read_indices[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, GA_fq_batch);
    fq_read2 = (fastq_read_t *) array_list_get(read_index, GA_rev_fq_batch);
      
    read_len = fq_read->length;
    norm_score = NORM_SCORE(output->score_p[i], read_len, input->match);
      
    if (norm_score >= min_score) {
      // valid mappings, 
      //insert in the list for further processing
      mapping_list = mapping_batch->mapping_lists[read_index];
      array_list_set_flag(0, mapping_list);
	
      if (array_list_size(mapping_list) == 0) {
	      mapping_batch->targets[new_num_targets++] = read_index;
      }
	
      sw_output = sw_output_new(strands[i],
				chromosomes[i],
				starts[i],
				strlen(r[i]),
				strlen(output->query_map_p[i]),
				output->query_start_p[i],
				output->ref_start_p[i],
				output->score_p[i],
				norm_score,
				output->query_map_p[i],
				output->ref_map_p[i]);
      array_list_insert(sw_output, mapping_list);
	
      mapping_batch->num_to_do++;
    }
      
    // free query and reference
    free(q[i]);
    free(r[i]);
  }

  mapping_batch->num_targets = new_num_targets;
    
  for (size_t i = sw_count2; i < sw_count; i++) {  
    read_index = read_indices[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, CT_fq_batch);
    fq_read2 = (fastq_read_t *) array_list_get(read_index, CT_rev_fq_batch);
      
    read_len = fq_read->length;
    norm_score = NORM_SCORE(output->score_p[i], read_len, input->match);
      
    if (norm_score >= min_score) {
      // valid mappings, 
      //insert in the list for further processing
      mapping_list = mapping_batch->mapping_lists2[read_index];
      array_list_set_flag(0, mapping_list);
	
      if (array_list_size(mapping_list) == 0) {
	      mapping_batch->targets2[new_num_targets2++] = read_index;
      }
	
      sw_output = sw_output_new(strands[i],
				chromosomes[i],
				starts[i],
				strlen(r[i]),
				strlen(output->query_map_p[i]),
				output->query_start_p[i],
				output->ref_start_p[i],
				output->score_p[i],
				norm_score,
				output->query_map_p[i],
				output->ref_map_p[i]);
      array_list_insert(sw_output, mapping_list);
	
      mapping_batch->num_to_do2++;
    }
      
    // free query and reference
    free(q[i]);
    free(r[i]);
  }

  mapping_batch->num_targets2 = new_num_targets2;
    
  // free
  sw_multi_output_free(output);

  LOG_DEBUG("------------ END OF APPLY SW BS 3NT -----------------------\n");    
}
//--------------------------------------------------------------------------------------

float cigar_code_get_sw_score(float match, float missm, 
			      float gap_open, float gap_extend, 
			      cigar_code_t *p) {
  float score = 0.0f;

  if (p) {    
    cigar_op_t *cigar_op;
    size_t num_ops = array_list_size(p->ops);
    int gaps = 0, M = 0;

    for (int i = 0; i < num_ops; i++) {
      cigar_op = array_list_get(i, p->ops);

      if (cigar_op->name == 'M') {
      	M += cigar_op->number;
      } else if (cigar_op->name == 'I' || cigar_op->name == 'D') {
        gaps += cigar_op->number;
        score -= (gap_open + (cigar_op->number - 1) * gap_extend);
      }
    }

    int num_mismatches = p->distance;
    int num_matches = M - num_mismatches;
    score += (num_matches * match + num_mismatches * missm);

    if (score < 0) {
      score = 0;
    }
  }
  return score;
}

//--------------------------------------------------------------------------------------

void apply_sw_bs_4nt(sw_server_input_t* input, batch_t *batch) {
  LOG_DEBUG("------------ APPLY SW BS 4NT START -----------------------\n");
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  genome_t *genome  = input->genome_p;
  genome_t *genome1 = input->genome1_p;
  genome_t *genome2 = input->genome2_p;
  sw_optarg_t *sw_optarg = &input->sw_optarg;

  ///////////////////////////////
  // Esta sección debe sacarse fuera de este fichero y ponerlo en el inicio de la ejecución
  ///////////////////////////////
  sw_optarg_t sw_optarg1;
  sw_optarg_t sw_optarg2;
 
  float match = sw_optarg->subst_matrix['A']['A'];
  float missm = sw_optarg->subst_matrix['C']['A'];
  float gap_open = sw_optarg->gap_open;
  float gap_extend = sw_optarg->gap_extend;

  sw_optarg1.gap_open   = gap_open;
  sw_optarg1.gap_extend = gap_extend;
  sw_optarg2.gap_open   = gap_open;
  sw_optarg2.gap_extend = gap_extend;
  
  fill_matrix(sw_optarg1.subst_matrix, match, missm, 0, 8, 2);
  fill_matrix(sw_optarg2.subst_matrix, match, missm, 1, 8, 2);
  ///////////////////////////////

  LOG_DEBUG("++++++++++ FILL GAPS 0     ++++++++++\n");
  fill_gaps_bs(mapping_batch, sw_optarg, genome1, genome2, 20, 5, 0, &sw_optarg1, &sw_optarg2);
  merge_seed_regions_bs(mapping_batch, 0);
  fill_end_gaps_bs(mapping_batch, sw_optarg, genome1, genome2, 20, 5, 0);
  LOG_DEBUG("++++++++++  END LIST 0     ++++++++++\n");

  LOG_DEBUG("++++++++++ FILL GAPS 1     ++++++++++\n");
  fill_gaps_bs(mapping_batch, sw_optarg, genome2, genome1, 20, 5, 1, &sw_optarg2, &sw_optarg1);
  merge_seed_regions_bs(mapping_batch, 1);
  fill_end_gaps_bs(mapping_batch, sw_optarg, genome2, genome1, 20, 5, 1);
  LOG_DEBUG("++++++++++  END LIST 1     ++++++++++\n");

  // now we can create the alignments
  fastq_read_t *read;
  array_list_t *fq_batch = mapping_batch->fq_batch;
  
  char *match_seq, *match_qual;
  size_t read_index, read_len, match_len, match_start;
  
  cal_t *cal;
  array_list_t *cal_list = NULL;
  size_t num_cals;
  
  seed_region_t *s;
  cigar_code_t *cigar_code;
  cigar_op_t *first_op;

  float score, norm_score, min_score = input->min_score;

  alignment_t *alignment;
  array_list_t *alignment_list;

  char *p, *optional_fields;
  int optional_fields_length, AS;

  array_list_t **mapping_lists;
  size_t num_targets;
  size_t *targets;

  for (int bs_id = 0; bs_id < 2; bs_id++) {
    if (bs_id == 0) {
      mapping_lists = mapping_batch->mapping_lists;
      num_targets = mapping_batch->num_targets;
      targets = mapping_batch->targets;
    } else {
      mapping_lists = mapping_batch->mapping_lists2;
      num_targets = mapping_batch->num_targets2;
      targets = mapping_batch->targets2;
    }

    for (size_t i = 0; i < num_targets; i++) {
      read_index = targets[i];
      read = (fastq_read_t *) array_list_get(read_index, fq_batch);
      
      cal_list = mapping_lists[read_index];
      num_cals = array_list_size(cal_list);
      
      if (num_cals <= 0) continue;
    
      read_len = read->length;
    
      alignment_list = array_list_new(num_cals, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

      // processing each CAL from this read
      for (size_t j = 0; j < num_cals; j++) {
        // get cal and read index
        cal = array_list_get(j, cal_list);

        if (cal->sr_list->size == 0) {
          continue;
        }

        if (cal->info == NULL) {
          s = (seed_region_t *) linked_list_get_first(cal->sr_list);
          cal->info = s->info;
          s->info = NULL;
        }
        
        cigar_code = (cigar_code_t *) cal->info;
        cigar_code_update(cigar_code);
        score = cigar_code_get_sw_score(match, missm, gap_open, gap_extend, cigar_code);
        norm_score = score / (match * read_len);

        LOG_DEBUG_F("score = %0.2f\n", norm_score);

        // filter by SW score
        if (norm_score > min_score) {
          // update cigar and sequence and quality strings
          LOG_DEBUG_F("\tcigar code = %s\n", new_cigar_code_string(cigar_code));
          match_start = 0;
          match_len = cigar_code_nt_length(cigar_code); 
          first_op = cigar_code_get_first_op(cigar_code);
          match_start = (first_op && first_op->name == 'H' ? first_op->number : 0);
          
          match_seq = (char *) malloc((match_len + 1)* sizeof(char));
          memcpy(match_seq, &read->sequence[match_start], match_len);
          match_seq[match_len] = 0;
          
          match_qual = (char *) malloc((match_len + 1)* sizeof(char));
          memcpy(match_qual, &read->quality[match_start], match_len);
          match_qual[match_len] = 0;
          
          // set optional fields
          optional_fields_length = 100;
          optional_fields = (char *) calloc(optional_fields_length, sizeof(char));
          
          p = optional_fields;
          AS = (int) norm_score * 100;
        
          sprintf(p, "ASi");
          p += 3;
          memcpy(p, &AS, sizeof(int));
          p += sizeof(int);
          
          sprintf(p, "NHi");
          p += 3;
          memcpy(p, &num_cals, sizeof(int));
          p += sizeof(int);
          
          sprintf(p, "NMi");
          p += 3;
          memcpy(p, &cigar_code->distance, sizeof(int));
          p += sizeof(int);
        
          // create an alignment and insert it into the list
          alignment = alignment_new();

          size_t header_len = strlen(read->id);
          char *head_id = (char *) malloc(header_len + 1);
          
          get_to_first_blank(read->id, header_len, head_id);
        
          alignment_init_single_end(head_id, match_seq, match_qual, 
                  cal->strand, cal->chromosome_id - 1, cal->start - 1,
                  strdup(new_cigar_code_string(cigar_code)), 
                  cigar_code_get_num_ops(cigar_code), 
                  norm_score * 254, 1, (num_cals > 1),
                  optional_fields_length, optional_fields, alignment);
          
          array_list_insert(alignment, alignment_list);

          LOG_DEBUG_F("creating alignment (bs_id = %i)...\n", bs_id);
        }
      }
      
      // free the cal list, and update the mapping list with the alignment list
      array_list_free(cal_list, (void *) cal_free2);
      mapping_lists[read_index] = alignment_list;
    }
  }

  LOG_DEBUG("------------ END OF APPLY SW BS 4NT -----------------------\n");
}

//--------------------------------------------------------------------------------------

sw_output_t *sw_output_new(int strand, size_t chrom, size_t ref_start, size_t ref_len,
			   size_t mref_len, size_t mquery_start, size_t mref_start,
			   float score, float norm_score, char* mquery, char* mref) {

  sw_output_t *p = (sw_output_t *) calloc(1, sizeof(sw_output_t));

  p->strand = strand;
  p->chromosome = chrom;
  p->ref_start = ref_start;
  p->ref_len = ref_len;
  p->mref_len = mref_len;
  p->mquery_start = mquery_start;
  p->mref_start = mref_start;
  p->score = score;
  p->norm_score = norm_score;
  p->mquery = strdup(mquery);
  p->mref = strdup(mref);

  return p;
}
