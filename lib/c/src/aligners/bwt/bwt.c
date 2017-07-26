#include "bwt.h"


//------------------------------------------------------------------------------

alignment_t* add_optional_fields(alignment_t *alignment, size_t n_mappings);

//------------------------------------------------------------------------------

void *__bwt_generate_anchor_list(size_t k_start, size_t l_start, int len_calc, bwt_optarg_t *bwt_optarg, 
				 bwt_index_t *index, int type, array_list_t *anchor_list, int type_anchor,
				 int dsp, size_t seq_start, size_t seq_end);

//------------------------------------------------------------------------------

void append_seed_region_linked_list(linked_list_t* sr_list,
				    size_t read_start, size_t read_end, 
				    size_t genome_start, size_t genome_end, 
				    int seed_id);

//------------------------------------------------------------------------------

void seed_region_select_linked_list(linked_list_t* sr_list, linked_list_t* sr_duplicate_list, 
				    size_t read_start, size_t read_end,
				    size_t genome_start, size_t genome_end,
				    int seed_id, unsigned char *seeds_ids_array);

//------------------------------------------------------------------------------

void my_new_cp_list_append_linked_list_v2(linked_list_t* list_p, bwt_anchor_t *bwt_anchor, size_t max_cal_distance, int max_seeds, int id);
//------------------------------------------------------------------------------
// Paratemers for the candidate alignment localizations (CALs)
//------------------------------------------------------------------------------



cal_optarg_t *cal_optarg_new(const size_t min_cal_size, 
			     const size_t max_cal_distance, 
			     const size_t num_seeds,
			     const size_t min_num_seeds_in_cal,
			     const size_t seed_size,
			     const size_t min_seed_size,
			     const size_t num_errors){
			       
  cal_optarg_t *cal_optarg_p = (cal_optarg_t *)malloc(sizeof(cal_optarg_t));
  cal_optarg_p->min_cal_size = min_cal_size;
  cal_optarg_p->max_cal_distance = max_cal_distance;
  cal_optarg_p->num_seeds = num_seeds;
  cal_optarg_p->min_num_seeds_in_cal = min_num_seeds_in_cal;
  cal_optarg_p->min_seed_size = min_seed_size;
  cal_optarg_p->seed_size = seed_size;
  cal_optarg_p->num_errors = num_errors;

  return cal_optarg_p;
}
    
void cal_optarg_free(cal_optarg_t *optarg) {
  free(optarg);
}

//------------------------------------------------------------------------------

cal_t *cal_new(const size_t chromosome_id, 
	       const short int strand,
	       const size_t start, 
	       const size_t end,
	       const size_t num_seeds,
	       const linked_list_t *sr_list,
	       const linked_list_t *sr_duplicate_list) {
		 
  cal_t *cal = (cal_t *)malloc(sizeof(cal_t));  

  cal->chromosome_id = chromosome_id;
  cal->end = end;
  cal->start = start;
  cal->strand = strand;
  cal->num_seeds = num_seeds;
  cal->sr_list = sr_list;
  cal->sr_duplicate_list = sr_duplicate_list;
  cal->read_area = 0;
  cal->info = NULL;
  cal->fill_gaps = 0;
  cal->num_targets;
  cal->candidates_seeds_start = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  cal->candidates_seeds_end = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  //TODO: Good CAL Selector??Mmm...
  //First, seed size sumation  

  //Second, read distance between the first element of the list and the last
  seed_region_t *s_first = linked_list_get_first(cal->sr_list);
  seed_region_t *s_last = linked_list_get_last(cal->sr_list);
  
  if (s_first && s_last) {
    cal->read_area += (s_last->read_end - s_first->read_start);
  }

  return cal;
}


void cal_free(cal_t *cal) {
  if (cal) {
    if (cal->sr_list) linked_list_free(cal->sr_list, (void *) seed_region_free);
    if (cal->sr_duplicate_list) linked_list_free(cal->sr_duplicate_list, (void *) seed_region_free);

    if (cal->candidates_seeds_start) array_list_free(cal->candidates_seeds_start, (void *)seed_region_free);
    if (cal->candidates_seeds_end) array_list_free(cal->candidates_seeds_end, (void *)seed_region_free);

    free(cal);
  }
}


//------------------------------------------------------------------------------

seed_region_t *seed_region_new(size_t read_start, size_t read_end, size_t genome_start, size_t genome_end, int id) {
  seed_region_t *seed_region = (seed_region_t *)malloc(sizeof(seed_region_t));
  seed_region->read_start = read_start;
  seed_region->read_end = read_end;
  seed_region->genome_start = genome_start;
  seed_region->genome_end = genome_end;
  seed_region->id = id;
  seed_region->info = NULL;
  seed_region->fusion_left = 0;
  seed_region->fusion_right = 0;

  return seed_region;
}

void seed_region_free(seed_region_t *seed_region) {
  free(seed_region);
}

//------------------------------------------------------------------------------
//Ricardo -- modified
short_cal_t *short_cal_new(const size_t start, 
	                   const size_t end,
			   const size_t seq_start,
			   const size_t seq_end,
			   const size_t seq_len,
			   const int max_seeds,
			   const int id) {

  short_cal_t *short_cal = (short_cal_t *)malloc(sizeof(short_cal_t));  
  
  short_cal->end = end;
  short_cal->start = start;
  short_cal->seq_len = seq_len;
  short_cal->num_seeds = 1;
  short_cal->seq_start = seq_start;	//Ricardo
  short_cal->seq_end = seq_end; //Ricardo
  short_cal->seeds_ids_array = (unsigned char *)calloc(max_seeds, sizeof(unsigned char));

  if (short_cal->seeds_ids_array == NULL) LOG_FATAL("NO MORE MEMORY\n");
  if (id >= max_seeds) LOG_FATAL_F("STORAGE SEED ID OVERFLOW: %i > %i\n", id, max_seeds);

  short_cal->seeds_ids_array[id] = 1;

  short_cal->sr_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
  seed_region_t *seed_region = seed_region_new(seq_start, seq_end, start, end, id);
  linked_list_insert(seed_region, short_cal->sr_list);
  short_cal->sr_duplicate_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);

  return short_cal;
}


void short_cal_free(short_cal_t *short_cal){
  if (short_cal) {
    if (short_cal->seeds_ids_array) free(short_cal->seeds_ids_array);
    if (short_cal->sr_list) linked_list_free(short_cal->sr_list, (void *) seed_region_free);
    if (short_cal->sr_duplicate_list) linked_list_free(short_cal->sr_duplicate_list, (void *) seed_region_free);

    free(short_cal);
  }
}

//------------------------------------------------------------------------------

bwt_anchor_t *bwt_anchor_new(int strand, int chromosome, size_t start, size_t end, int type, size_t seq_start, size_t seq_end) {
  bwt_anchor_t *bwt_anchor = (bwt_anchor_t *)malloc(sizeof(bwt_anchor_t));

  bwt_anchor->strand = strand;
  bwt_anchor->chromosome = chromosome;
  bwt_anchor->start = start;
  bwt_anchor->end = end;
  bwt_anchor->type = type;

  bwt_anchor->seq_start = seq_start;
  bwt_anchor->seq_end = seq_end;
  bwt_anchor->seq_len = seq_end - seq_start + 1;

  return bwt_anchor;
}

void bwt_anchor_free(bwt_anchor_t *bwt_anchor) {
  free(bwt_anchor);
}


//------------------------------------------------------------------------------

bwt_optarg_t *bwt_optarg_new(const size_t num_errors, const size_t num_threads, const int filter_read_mappings, 
    const int filter_seed_mappings, const size_t min_cal_size, const double umbral_cal_length_factor,
		const int min_read_discard, const int max_inner_gap) {
  bwt_optarg_t *bwt_optarg = (bwt_optarg_t *) calloc(1, sizeof(bwt_optarg_t));
  
  bwt_optarg->num_errors = num_errors;
  bwt_optarg->num_threads = num_threads;
  bwt_optarg->filter_read_mappings = filter_read_mappings;
  bwt_optarg->filter_seed_mappings = filter_seed_mappings;
  bwt_optarg->min_cal_size = min_cal_size;		//Ricardo
  bwt_optarg->umbral_cal_length_factor = umbral_cal_length_factor; //Ricardo
  bwt_optarg->min_read_discard = min_read_discard; //Ricardo
  bwt_optarg->max_inner_gap = max_inner_gap; //Ricardo

  return bwt_optarg;
}

//-----------------------------------------------------------------------------

void bwt_optarg_free(bwt_optarg_t *optarg) {
  if (optarg != NULL) {
    free(optarg);
  }
}

//-----------------------------------------------------------------------------

void bwt_cigar_cpy(alignment_t *mapping, char *quality) {
  
  unsigned int quality_type;
  size_t quality_len;

  quality_len = strlen(quality);
  quality_type = atoi(mapping->quality);
  free(mapping->quality);

  mapping->quality = (char *)malloc(sizeof(char)*(quality_len + 1));

  if (mapping->seq_strand == 0) {
       memcpy(mapping->quality, quality, quality_len);
  } else {
       reverse_str(quality,
		   mapping->quality, quality_len);
  }
}

//-----------------------------------------------------------------------------

unsigned int alignmentcmp(alignment_t *alignment_1, alignment_t *alignment_2) {
  
  size_t cigar1_len = strlen(alignment_1->cigar);
  size_t cigar2_len = strlen(alignment_2->cigar);

  if (alignment_1->map_quality > alignment_2->map_quality) { 
    return 1;
  } else if (alignment_1->map_quality < alignment_2->map_quality) { 
    return 2; 
  } else {
    if (alignment_1->num_cigar_operations == 1 || alignment_2->num_cigar_operations == 1) {
      if (alignment_1->cigar[cigar1_len - 1] == '=' && alignment_2->cigar[cigar1_len - 1] == '=') { 
        return 0;
      } if (alignment_1->cigar[cigar1_len - 1] == '=') { 
        return 1; 
      } else if (alignment_2->cigar[cigar2_len - 1] == '=') { 
        return 2; 
      }
      
      if (alignment_1->cigar[cigar1_len - 1] == 'M' && alignment_2->cigar[cigar1_len - 1] == 'M') { 
        return 0; 
      }
      
      if (alignment_1->cigar[cigar1_len - 1] == 'M') { 
        return 1; 
      } else if (alignment_2->cigar[cigar2_len - 1] == 'M') { 
        return 2; 
      }  
    } else {
      if (alignment_1->num_cigar_operations < alignment_2->num_cigar_operations) {
	      return 1;
      } else if (alignment_1->num_cigar_operations > alignment_2->num_cigar_operations) {
	      return 2;
      } else { 
        return 0; 
      }
    }
  }

  return 0;
}

//-----------------------------------------------------------------------------

char* reverse_str(char *src, char *dsp, size_t length) {
  size_t j = length - 1;
  size_t i = 0;

  memcpy(dsp, src, length);
  
  for (i = 0; i < length; i++) {
    dsp[i] = src[j--];
  }

  dsp[i] = '\0';
  return dsp;
}

//-----------------------------------------------------------------------------

bwt_index_t *bwt_index_new(const char *dirname) {

  bwt_index_t *index = (bwt_index_t*) calloc(1, sizeof(bwt_index_t));

  index->dirname = strdup(dirname);

  // be careful: added by JTPP to support 3-nt genomes (useful for bs)
  index->nucleotides = readNucleotide(dirname, "Nucleotide");
  bwt_init_replace_table(index->nucleotides, index->table, index->rev_table);

  readUIntVector(&index->h_C, dirname, "C");
  readUIntVector(&index->h_C1, dirname, "C1");
  readCompMatrix(&index->h_O, dirname, "O");
  readUIntCompVector(&index->S, dirname, "Scomp");

  reverseStrandC(&index->h_rC, &index->h_C,
  		 &index->h_rC1, &index->h_C1);

  reverseStrandO(&index->h_rO, &index->h_O);

  // 1-error handling 
  readCompMatrix(&index->h_Oi, dirname, "Oi");
  readUIntCompVector(&index->Si, dirname, "Scompi");

  reverseStrandO(&index->h_rOi, &index->h_Oi);

  // load karyotype
  char path[strlen(dirname) + 512];
  //  sprintf(path, "%s/index.txt", dirname);
  //load_exome_file(&index->karyotype, path);

  load_exome_file(&index->karyotype, dirname);

  return index;
}

//-----------------------------------------------------------------------------

void bwt_index_free(bwt_index_t *index) {
  if (index == NULL) return;

  free(index->dirname);
  // be careful: added by JTPP to support 3-nt genomes (useful for bs)
  free(index->nucleotides);

  free(index->h_C.vector);
  free(index->h_C1.vector);
  free(index->h_rC.vector);
  free(index->h_rC1.vector);

  freeCompMatrix(&index->h_O);
  free(index->S.vector);

  // 1-error handling
  freeCompMatrix(&index->h_Oi);
  free(index->Si.vector);

  free(index);
}

//-----------------------------------------------------------------------------

void bwt_generate_index_files(char *ref_file, char *output_dir, 
			      unsigned int s_ratio) {

  byte_vector X, B, Bi;
  vector C, C1;
  comp_vector S, Si, Scomp, Scompi;
  comp_vector R, Ri, Rcomp, Rcompi;
  comp_matrix O, Oi;

  exome ex;
  initReplaceTable_bs(NULL);

  // Calculating BWT
  calculateBWT(&B, &S, &X, 0, &ex, ref_file);

  save_exome_file(&ex, output_dir);

  saveCharVector(&X, output_dir, "X");
  free(X.vector);

  printUIntVector(S.vector, S.n);
  printUIntVector(B.vector, B.n);

  // Calculating prefix-trie matrices C and O
  calculateC(&C, &C1, &B, 0);
  calculateO(&O, &B);

  printUIntVector(C.vector, C.n);
  printUIntVector(C1.vector, C1.n);
  printCompMatrix(O);

  saveCharVector(&B, output_dir, "B");
  free(B.vector);
  saveUIntVector(&C, output_dir, "C");
  free(C.vector);
  saveUIntVector(&C1, output_dir, "C1");
  free(C1.vector);
  saveCompMatrix(&O, output_dir, "O");
  freeCompMatrix(&O);

  // Calculating R
  calculateR(&S, &R);

  printUIntVector(R.vector, R.n);

  // Calculating Scomp Rcomp
  calculateSRcomp(&S, &Scomp, s_ratio);
  printUIntVector(Scomp.vector, Scomp.n);
  calculateSRcomp(&R, &Rcomp, s_ratio);
  printUIntVector(Rcomp.vector, Rcomp.n);


  saveUIntCompVector(&S, output_dir, "S");
  free(S.vector);
  saveUIntCompVector(&R, output_dir, "R");
  free(R.vector);
  saveUIntCompVector(&Scomp, output_dir, "Scomp");
  free(Scomp.vector);
  saveUIntCompVector(&Rcomp, output_dir, "Rcomp");
  free(Rcomp.vector);

  //Calculating BWT of reverse reference
  calculateBWT(&Bi, &Si, &X, 1, NULL, ref_file);

  saveCharVector(&X, output_dir, "Xi");
  free(X.vector);

  printUIntVector(Bi.vector, Bi.n);
  printUIntVector(Si.vector, Si.n);

  //Calculating inverted prefix-trie matrix Oi
  calculateO(&Oi, &Bi);

  printCompMatrix(Oi);

  saveCharVector(&Bi, output_dir, "Bi");
  free(Bi.vector);

  saveCompMatrix(&Oi, output_dir, "Oi");
  freeCompMatrix(&Oi);

  //Calculating Ri
  calculateR(&Si, &Ri);

  printUIntVector(Ri.vector, Ri.n);

  // Calculating Scompi Rcompi
  calculateSRcomp(&Si, &Scompi, s_ratio);
  printUIntVector(Scompi.vector, Scompi.n);
  calculateSRcomp(&Ri, &Rcompi, s_ratio);
  printUIntVector(Rcompi.vector, Rcompi.n);

  saveUIntCompVector(&Si, output_dir, "Si");
  free(Si.vector);
  saveUIntCompVector(&Ri, output_dir, "Ri");
  free(Ri.vector);
  saveUIntCompVector(&Scompi, output_dir, "Scompi");
  free(Scompi.vector);
  saveUIntCompVector(&Rcompi, output_dir, "Rcompi");
  free(Rcompi.vector);
}

//-----------------------------------------------------------------------------

void *__bwt_generate_anchor_list(size_t k_start, size_t l_start, int len_calc, bwt_optarg_t *bwt_optarg, 
				 bwt_index_t *index, int type, array_list_t *anchor_list, int type_anchor,
				 int dsp, size_t seq_start, size_t seq_end) {
  size_t idx, key, direction;
  size_t start_mapping;
  const int MAX_BWT_ANCHORS = 100;
  bwt_anchor_t *bwt_anchor;

  if (type) { 
    //(+)
    if (type_anchor == BACKWARD_ANCHOR) {
      direction = BACKWARD_ANCHOR;
    } else {
      direction = FORWARD_ANCHOR; 
    }
  } else { 
    //(-)
    if (type_anchor == BACKWARD_ANCHOR) {
      direction = BACKWARD_ANCHOR; 
    } else {
      direction = FORWARD_ANCHOR;
    }
  }

  if (l_start - k_start < MAX_BWT_ANCHORS) {
    for (size_t j = k_start; j <= l_start; j++) {
      if (index->S.ratio == 1) {
        key = (direction) ? index->Si.siz - index->Si.vector[j] - len_calc - 1 : index->S.vector[j];
      } else {
        key = (direction) ? index->Si.siz - getScompValue(j, &index->Si, &index->h_C,
              &index->h_Oi) - len_calc - 1: getScompValue(j, &index->S, &index->h_C, &index->h_O);
      }

      idx = binsearch(index->karyotype.offset, index->karyotype.size, key);

      if (key + len_calc <= index->karyotype.offset[idx]) {
        start_mapping = index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]);	
        bwt_anchor = bwt_anchor_new(!type, idx - 1, start_mapping + 1, start_mapping + len_calc, type_anchor, seq_start, seq_end);
        LOG_DEBUG_F("anchor: %i:%lu-%lu\n", idx, bwt_anchor->start, bwt_anchor->end);
        array_list_insert(bwt_anchor, anchor_list);
      }
    }
  }
}

//-----------------------------------------------------------------------------

alignment_t* add_optional_fields(alignment_t *alignment, size_t n_mappings) {
  bam_int_t distance;
  bam_int_t AS = strlen(alignment->sequence);
  bam_int_t NH = n_mappings;

  bam_tag_t* as_tag = bam_tag_init(AS_TAG_NAME, BAM_TAG_TYPE_INT, 0, 0);
  bam_tag_t* nh_tag = bam_tag_init(NH_TAG_NAME, BAM_TAG_TYPE_INT, 0, 0);
  bam_tag_t* nm_tag = bam_tag_init(NM_TAG_NAME, BAM_TAG_TYPE_INT, 0, 0);

  if (alignment->cigar[strlen(alignment->cigar) - 1] == '=') {
    distance = 0;
  } else {
    distance = 1;
  }

  bam_tag_set_scalar(as_tag, &AS);
  bam_tag_set_scalar(nh_tag, &NH);
  bam_tag_set_scalar(nm_tag, &distance);

  array_list_insert(as_tag, alignment->optional_tags);
  array_list_insert(nh_tag, alignment->optional_tags);
  array_list_insert(nm_tag, alignment->optional_tags);
  
  return alignment;
}

//-----------------------------------------------------------------------------
//Ricardo - new function in order to use bwt_anchor instead of region
void insert_seeds_and_merge_anchor(array_list_t *mapping_list, linked_list_t ***cals_list,  size_t max_cal_distance, int id_start) {
  for (int m = array_list_size(mapping_list) - 1; m >= 0; m--) {
	bwt_anchor_t *bwt_anchor = array_list_remove_at(m, mapping_list);
    int  chromosome_id = bwt_anchor->chromosome;
    int strand = bwt_anchor->strand;


    my_new_cp_list_append_linked_list_v2(cals_list[strand][chromosome_id], bwt_anchor, max_cal_distance, 1000, id_start + m);		//m

    bwt_anchor_free(bwt_anchor);
  }
}

//-----------------------------------------------------------------------------

void print_se_region(short_cal_t *short_cal){
  printf("(%lu-%lu)->", short_cal->start, short_cal->end);
}


void seed_region_select_linked_list(linked_list_t* sr_list, linked_list_t* sr_duplicate_list, 
				    size_t read_start, size_t read_end,
				    size_t genome_start, size_t genome_end,
				    int seed_id, unsigned char *seeds_ids_array) {
  //printf("\tInsert [Seed:=%lu-%lu](%i): ", read_start, read_end, seed_id);
  seed_region_t *item;
  if (!seeds_ids_array[seed_id]) { 
    //printf(" Not in!!\n");
    seeds_ids_array[seed_id]++;
    item = seed_region_new(read_start, read_end, genome_start, genome_end, seed_id);
    linked_list_insert(item, sr_list);
  } else {
    //printf(" In!!\n");
    if (seeds_ids_array[seed_id] == 1) { 
      linked_list_iterator_t* itr = linked_list_iterator_new(sr_list);
      item = (seed_region_t *)linked_list_iterator_curr(itr);
      while (item != NULL) {
	if (item->id == seed_id) {
	  item = linked_list_iterator_remove(itr);
	  linked_list_insert(item, sr_duplicate_list);
	  //printf("\tRemove [Seed:=%lu-%lu](%i)\n", item->read_start, item->read_end, item->id);
	  break;
	}
	linked_list_iterator_next(itr);
	item = linked_list_iterator_curr(itr);
      }
      linked_list_iterator_free(itr);
    }
    item = seed_region_new(read_start, read_end, genome_start, genome_end, seed_id);
    linked_list_insert(item, sr_duplicate_list);
    seeds_ids_array[seed_id]++;
  }


  //printf("\t sr_size = %i, sr_duplicate = %i\n", linked_list_size(sr_list), linked_list_size(sr_duplicate_list));

  /*
  unsigned char actualization = 0;
  seed_region_t *item, *item_aux, *new_item, *item_free;
  linked_list_iterator_t* itr = linked_list_iterator_new(list_p);

  printf("\tInsert [Seed:=%lu-%lu]\n", read_start, read_end);
  if (linked_list_size(list_p) <= 0) {
    new_item = seed_region_new(read_start, read_end, genome_start, genome_end);
    linked_list_insert(new_item, list_p);
  } else {
    item = (seed_region_t *)linked_list_iterator_curr(itr);
    while (item != NULL) {
      //printf("\t compare with %lu\n", item->start);
      if (read_start <= item->read_start) {
	new_item = seed_region_new(read_start, read_end, genome_start, genome_end);
	linked_list_iterator_insert(new_item, itr);
	break;
      }      
      //continue loop...
      linked_list_iterator_next(itr);
      item = linked_list_iterator_curr(itr);      
    }// end while

    if (item == NULL) {
      new_item = seed_region_new(read_start, read_end, genome_start, genome_end);
      //printf("\tInsert at END\n");
      linked_list_insert_last(new_item, list_p);
    }
    //printf("Insert OK! and now actualization\n");
  }

  linked_list_iterator_free(itr);
  */
}


//-----------------------------------------------------------------------------

void append_seed_region_linked_list(linked_list_t* sr_list,
				    size_t read_start, size_t read_end, 
				    size_t genome_start, size_t genome_end, 
				    int seed_id) {
  unsigned char actualization = 0;
  seed_region_t *item, *item_aux, *new_item, *item_free;
  linked_list_iterator_t* itr = linked_list_iterator_new(sr_list);
  
  if (linked_list_size(sr_list) <= 0) {
    new_item = seed_region_new(read_start, read_end, genome_start, genome_end, seed_id);
    linked_list_insert(new_item, sr_list);
    //printf("Call linked_list_insert\n");
  } else {
    item = (seed_region_t *)linked_list_iterator_curr(itr);
    while (item != NULL) {
      //printf("\t compare with %lu\n", item->start);
      if (genome_start < item->genome_start) {
	if (genome_end + 1 < item->genome_start) {
	  /*********************************************
	   *    Case 1: New item insert before item.   *
           *                                           *
           *        new item     item                  *
           *       |-------| |--------|                *
           ********************************************/
	  //printf("\t Insert now before %lu\n", item->start);
	  new_item = seed_region_new(read_start, read_end, genome_start, genome_end, seed_id);
	  linked_list_iterator_insert(new_item, itr);
	  //printf("Call linked_list_iterator_insert\n");
	  linked_list_iterator_prev(itr);
	} else {
	  /********************************************
           *  Case 2: Actualization item start        *
           *           new item                       *
           *          |-------|   item                *
           *                   |--------|             *                            
           ********************************************/
	  /*printf("\tFusion Case 2! [%lu|%i - %i|%lu], [%lu|%i - %i|%lu]\n", 
		 genome_start, read_start, genome_end, read_end, 
		 item->genome_start, item->read_start, item->genome_end, item->read_end);*/
	  item->read_start = read_start;
	  item->genome_start = genome_start;
	  if (genome_end > item->genome_end) {
	    /**************************************************
             *  Case 3: Actualization item start and item end *
             *          new item                              *
             *         |------------|                         *
             *              item                              *    
             *           |--------|                           *                                    
             **************************************************/
	    //printf("\tFusion!\n");
	    item->read_end = read_end;
	    item->genome_end = genome_end;
	    actualization = 1;
	  }
	}
	break;
      } else {
	if (genome_end <= item->genome_end) {
	  /**************************************************                                       
           *  Case 4: The new item don't insert in the list *                             
           *              item                              * 
           *         |-------------|                        * 
           *             new item                           * 
           *            |--------|                          * 
           **************************************************/
	  //printf("\tFusion!\n");
	  break;
	} else if (item->genome_end + 1 >= genome_start) {
	  /********************************************                                              
           *  Case 5: Actualization item end          *
           *            item                          *                                              
           *          |-------| new item              *                                            
           *                 |--------|               *                                              
           ********************************************/
	  //printf("\tFusion!\n");
	  /*printf("\tFusion Case 5! New Item:[%lu|%i - %i|%lu], Item:[%lu|%i - %i|%lu]\n", 
		 genome_start, read_start, read_end, genome_end,  
		 item->genome_start, item->read_start, item->read_end, item->genome_end);
	  */
	  item->read_end = read_end;
	  item->genome_end = genome_end;
	  actualization = 1;
	  break;
	}
      } // end else

      //continue loop...
      linked_list_iterator_next(itr);
      item = linked_list_iterator_curr(itr);
      
    } // end while

    if (item == NULL) {
      /******************************************************* 
       * Case 5: Insert new item at the end of the list      * 
       *                 item    new item                    * 
       *              |-------| |--------|                   *    
       *******************************************************/
      new_item = seed_region_new(read_start, read_end, genome_start, genome_end, seed_id);
      //printf("\tInsert at END\n");
      linked_list_insert_last(new_item, sr_list);
      //printf("Call linked_list_insert_last\n");
    }
    //printf("Insert OK! and now actualization\n");
    if (actualization == 1) {
      //printf("\tActualization RIGHT items (Next). Current item [%d-%d]\n", item->start, item->end);
      //printf("List before actualization\n");
      //linked_list_print(list_p, print_se_region);
      //printf("\n\n");
 
      linked_list_iterator_next(itr);
      item_aux = linked_list_iterator_curr(itr);

      while (item_aux != NULL) {
	//printf("\t\tFusion right items. item->end=%d < item_aux->start=%d?\n", item->end, item_aux->start);
	//printf("\tIterator are in [%lu-%lu] and compare with item [%lu-%lu]\n", item_aux->start, item_aux->end, item->start, item->end);
	if (item->genome_end + 1 < item_aux->genome_start) {
	  //printf("\t\tSTOP Actualization\n");
	  break;
	} else {
	  //printf("\t\tCONTINUE Actualization. item->end=%d < item_aux->end=%d?\n", item->end, item_aux->end);
	  if (item->genome_end + 1 < item_aux->genome_end) {
	    //printf("\t\tActualization end value %d\n", item_aux->end);
	    item->read_end = item_aux->read_end;
	    item->genome_end = item_aux->genome_end;
	  }
          //printf("\t\tDelete item %d-%d\n", item_aux->start, item_aux->end);
	  item_free = linked_list_iterator_remove(itr);
	  if (item_free) { seed_region_free(item_free); }
	  //printf("\t\tDelete OK!\n");
	  item_aux = linked_list_iterator_curr(itr);
	}                                                                                             
      }
      /*printf("List after actualization\n");
      linked_list_print(list_p, print_se_region);
      printf("\n\n");*/
    }
  }//end first else
  //printf("End insert and actualization\n");
  linked_list_iterator_free(itr);

  /*printf("Status Seed Region list %lu [%i-%i][%i-%i]:\n", linked_list_size(sr_list), 
	 ((seed_region_t *)sr_list->first->item)->read_start, ((seed_region_t *)sr_list->first->item)->read_end,
	 ((seed_region_t *)sr_list->last->item)->read_start,  ((seed_region_t *)sr_list->last->item)->read_end);
  */
  /*for (linked_list_item_t *list_item = sr_list->first; list_item != NULL; list_item = list_item->next) {
    printf("Status Seed Region list %lu:\n", linked_list_size(sr_list));
    seed_region_t *s = list_item->item;
    printf("[%i|%i - %i|%i]  ", s->genome_start, s->read_start, s->read_end, s->genome_end);
  }
  printf("\n");
  */
}

//-----------------------------------------------------------------------------
//Versión nueva Ricardo -- v2
//No se puede mapear más que el tamaño de la read (max_cal_distance = tamaño read, aunque podría ser distinta ya que se pone como parámetro)

void my_new_cp_list_append_linked_list_v2(linked_list_t* list_p, bwt_anchor_t *bwt_anchor, size_t max_cal_distance, int max_seeds, int id) {
  unsigned char actualization = 0;
  short_cal_t *item, *item_aux, *new_item_p, *item_free;

  //int strand = bwt_anchor->strand;
  size_t start = bwt_anchor->start;
  size_t end = bwt_anchor->end;
  size_t seq_start = bwt_anchor->seq_start;
  size_t seq_end = bwt_anchor->seq_end;
  size_t seq_len = bwt_anchor->seq_len;
  linked_list_iterator_t* itr = linked_list_iterator_new(list_p);
  //printf("LINKED LIST INSERT %lu|%i-%i|%lu:\n", start, seq_start, seq_end, end);

 //printf("\nmax_cals_distance =%i\n", max_cal_distance);
 //printf("start: %i, end: %i, seq_start: %i, seq_end: %i, seq_len: %i\n", start, end, seq_start, seq_end, seq_len);



  if (linked_list_size(list_p) <= 0) {
    new_item_p = short_cal_new(start, end, seq_start, seq_end, seq_len, max_seeds, id);
    linked_list_insert(new_item_p, list_p);
  } else {
    item = (short_cal_t *)linked_list_iterator_curr(itr);
    //printf("item_start: %i, item_end: %i, iterm_seq_start: %i, item_seq_end: %i, item_seq_len: %i\n", item->start, item->end, item->seq_start, item->seq_end, item->seq_len);
    while (item != NULL) {
      if (start < item->start) {
    	 if (start + max_cal_distance - seq_start < item->end) {
	  /*********************************************
	   *    Case 1: New item insert before item.   *
           *       | -- max cal--|                                    *
           *        new item         item                  *
           *       |-------|        |--------|             *
           ********************************************/
      /*********************************************
      *    Case 1: New item insert before item.        *
    	   *          | -- max cal--|                                 *
    	   *           new item                        *
           *          |-------|       item                 *
           *                 |----------|                *
           *********************************************/
      // printf("Caso1\n");
	  //printf("\t Insert now before %lu\n", item->start);
	  new_item_p = short_cal_new(start, end, seq_start, seq_end, seq_len, max_seeds, id);
	  linked_list_iterator_insert(new_item_p, itr);
	  linked_list_iterator_prev(itr);
	} else {
	  /********************************************
           *  Case 2: Actualization item start        *
           *          | ---- max cal------|                                 *
           *           new item                       *
           *          |-------------|  			*
           *          			 item                *
           *                 |--------|             *
           ********************************************/
	  //printf("\tFusion!\n");
		//printf("Caso2\n");
	  item->start = start;
	  item->seq_start = seq_start;
	  item->num_seeds++;
	  //item->seq_start = seq_start;
	  seed_region_select_linked_list(item->sr_list, item->sr_duplicate_list,
					 seq_start, seq_end, start, end, id,
					 item->seeds_ids_array);
	  if (end > item->end) {
	    /**************************************************
             *  Case 3: Actualization item start and item end *
             *         | ---- max cal------|                                 *
             *          new item                              *
             *         |------------|                         *
             *              item                              *
             *           |--------|                           *
             **************************************************/
	    //printf("\tFusion!\n");
	//	  printf("Caso3\n");
	    item->end = end;
	    item->seq_end = seq_end;
	    actualization = 1;
	  }
	}
	break;
      } else {
	if (end <= item->end) {
	  /**************************************************
           *  Case 4: The new item don't insert in the list *
           *              item                              *
           *         |-------------|                        *
           *             new item                           *
           *            |--------|                          *
           **************************************************/
	  //printf("\tFusion!\n");
		//printf("Caso4\n");
	  item->num_seeds++;
	  seed_region_select_linked_list(item->sr_list, item->sr_duplicate_list, seq_start,
					 seq_end, start, end, id, item->seeds_ids_array);
	  break;
	} else if (item->start + max_cal_distance - item->seq_start >= end){

	  /********************************************
           *  Case 5: Actualization item end          *
           *          | ---- max cal------|                                 *
           *            item                          *
           *          |-------| new item              *
           *                 |--------|               *
           ********************************************/
	  //printf("\tFusion!\n");
		//printf("Caso5\n");
	  item->end = end;
	  item->seq_end = seq_end;
	  seed_region_select_linked_list(item->sr_list, item->sr_duplicate_list, seq_start,
					 seq_end, start, end, id, item->seeds_ids_array);
	  actualization = 1;
	  item->num_seeds++;
	  break;
	}
      } // end else

      //continue loop...
      linked_list_iterator_next(itr);
      item = linked_list_iterator_curr(itr);

    } // end while

    if (item == NULL) {
      /*******************************************************
       * Case 6: Insert new item at the end of the list      *
       *   *          | --- max cal----|                                 *
       *                 item            new item                    *
       *              |-------|          |--------|                   *
       *******************************************************/
    	//printf("Caso6\n");
      new_item_p = short_cal_new(start, end, seq_start, seq_end, seq_len, max_seeds, id);
      //printf("\tInsert at END\n");
      linked_list_insert_last(new_item_p, list_p);
    }
    //printf("Insert OK! and now actualization\n");
 //RICARDO
    if (actualization == 1) {
  //  	printf("Caso7\n");
      //printf("\tActualization RIGHT items (Next). Current item [%d-%d]\n", item->start, item->end);
      //printf("List before actualization\n");
      //linked_list_print(list_p, print_se_region);
      //printf("\n\n");

      linked_list_iterator_next(itr);
      item_aux = linked_list_iterator_curr(itr);

      while (item_aux != NULL) {
	//printf("\t\tFusion right items. item->end=%d < item_aux->start=%d?\n", item->end, item_aux->start);
	//printf("\tIterator are in [%lu-%lu] and compare with item [%lu-%lu]\n", item_aux->start, item_aux->end, item->start, item->end);
	//if (item->end + max_cal_distance < item_aux->start) {
    if (item->start + max_cal_distance - item->seq_start < item_aux->end){
	  //printf("\t\tSTOP Actualization\n");
	  break;
	} else {
	  //printf("\t\tCONTINUE Actualization. item->end=%d < item_aux->end=%d?\n", item->end, item_aux->end);
	  if (item->end < item_aux->end) {
	    //printf("\t\tActualization end value %d\n", item_aux->end);
	    item->end = item_aux->end;
	    //item->seq_end = item_aux->seq_end;
	    seed_region_t *seed_region_aux;
	    while (seed_region_aux = linked_list_remove_first(item_aux->sr_list)) {
	      seed_region_select_linked_list(item->sr_list, item->sr_duplicate_list, seed_region_aux->read_start,
					     seed_region_aux->read_end, seed_region_aux->genome_start,
					     seed_region_aux->genome_end, seed_region_aux->id, item->seeds_ids_array);
	      seed_region_free(seed_region_aux);
	    }

	  }
          //printf("\t\tDelete item %d-%d\n", item_aux->start, item_aux->end);
	  item_free = linked_list_iterator_remove(itr);
	  if (item_free) { short_cal_free(item_free); }
	  //printf("\t\tDelete OK!\n");
	  item_aux = linked_list_iterator_curr(itr);
	}
      }
      /*printf("List after actualization\n");
      linked_list_print(list_p, print_se_region);
      printf("\n\n");*/
    } //Fin if actualia
  }//end first else
  //printf("End insert and actualization\n");
  linked_list_iterator_free(itr);
}
