#include "bwt.h"

//-----------------------------------------------------------------------------
// exact functions
//-----------------------------------------------------------------------------

size_t bwt_map_exact_seq(char *seq, 
			 bwt_optarg_t *bwt_optarg, 
			 bwt_index_t *index, 
			 array_list_t *mapping_list);

size_t bwt_map_exact_read(fastq_read_t *read, 
			  bwt_optarg_t *bwt_optarg, 
			  bwt_index_t *index, 
			  array_list_t *mapping_list);

size_t bwt_map_exact_seqs(char **seqs, 
			  size_t num_reads,
			  bwt_optarg_t *bwt_optarg, 
			  bwt_index_t *index, 
			  char *out_status,
			  array_list_t *mapping_list);

size_t bwt_map_exact_batch(fastq_batch_t *batch,
			   bwt_optarg_t *bwt_optarg, 
			   bwt_index_t *index, 
			   fastq_batch_t *unmapped_batch,
			   array_list_t *mapping_list);

//------------------------------------------------------------------------------

size_t bwt_map_exact_seed(char *seq, size_t seq_len, 
			  size_t start, size_t end,
			  bwt_optarg_t *bwt_optarg,
			  bwt_index_t *index,
			  array_list_t *mapping_list,
			  int id);

//-----------------------------------------------------------------------------
// inexact functions
//-----------------------------------------------------------------------------

size_t bwt_map_inexact_seq(char *seq, 
			   bwt_optarg_t *bwt_optarg, 
			   bwt_index_t *index, 
			   array_list_t *mapping_list);


size_t bwt_map_inexact_seqs(char **seqs, 
			    size_t num_reads,
			    bwt_optarg_t *bwt_optarg, 
			    bwt_index_t *index, 
			    char *out_status,
			    array_list_t *mapping_list);

//------------------------------------------------------------------------------

char *bwt_error_type(char error_kind);

//------------------------------------------------------------------------------

size_t __bwt_map_inexact_read(fastq_read_t *read, 
			      bwt_optarg_t *bwt_optarg, 
			      bwt_index_t *index, 
			      array_list_t *mapping_list);

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
  
  //printf("List Size %i:\n", cal->sr_list->size);
  if (s_first && s_last) {
    cal->read_area += (s_last->read_end - s_first->read_start);
    //printf("\t %i - %i = %i\n", s_last->read_end, s_first->read_start, cal->read_area);
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

void cal_print(cal_t *cal) {
  printf(" CAL (%c)[%lu:%lu-%lu]:\n", cal->strand==0?'+':'-', 
	 cal->chromosome_id, cal->start, cal->end);
  printf("\t SEEDS LIST: ");
  if (cal->sr_list == NULL || cal->sr_list->size == 0) {
    printf(" NULL\n");
  } else {
    for (linked_list_item_t *item = cal->sr_list->first; 
	 item != NULL; item = item->next) {
      seed_region_t *seed = item->item;
      printf(" [%lu|%lu - %lu|%lu] ", seed->genome_start, seed->read_start, 
	     seed->read_end, seed->genome_end);
    }
    printf("\n");
  }

  if (cal->sr_duplicate_list != NULL && cal->sr_duplicate_list->size > 0) {
    printf("\t SEEDS DUPLICATE LIST: ");
    for (linked_list_item_t *item = cal->sr_duplicate_list->first; 
	 item != NULL; item = item->next) {
      seed_region_t *seed = item->item;
      printf(" [%lu|%lu - %lu|%lu] ", seed->genome_start, seed->read_start, 
	     seed->read_end, seed->genome_end);
    }
    printf("\n");
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

  int jj= 0;
  jj++;
  jj++;
  jj++;
  jj++;
  jj++;
  jj++;
  jj++;
  jj++;
  jj++;
  jj++;
  jj++;
  jj++;
  jj++;
  jj++;
  jj++;
  jj++;
  seed_region->read_start = read_start;
   seed_region->read_end = read_end;
   seed_region->genome_start = genome_start;
   seed_region->genome_end = genome_end;
   seed_region->id = id;
   seed_region->info = NULL;
   seed_region->fusion_left = 0;
   seed_region->fusion_right = 0;
  jj++;
  jj++;
  jj++;
  jj++;
  jj++;

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
  //  printf("\tInsert New [Seed:=%lu-%lu]\n", seq_start, seq_end);
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

region_t *region_bwt_new(const size_t chromosome_id, 
			 const short int strand,
			 const size_t start, 
			 const size_t end,
			 const size_t seq_start,
			 const size_t seq_end,
			 const size_t seq_len,
			 const int id) {
		 
  region_t *region = (region_t *) malloc(sizeof(region_t));  

  region->chromosome_id = chromosome_id;
  region->end = end;
  region->start = start;
  region->strand = strand;
  region->seq_start = seq_start;
  region->seq_end = seq_end;
  region->seq_len = seq_len;
  region->id = id;

  return region;
}


void region_bwt_free(region_t *region) {
  free(region);
}

//------------------------------------------------------------------------------

read_cals_t *read_cals_new(fastq_read_t *read) {

    read_cals_t *read_cals = (read_cals_t *) malloc(sizeof(read_cals_t));

    read_cals->cal_list = array_list_new(100000, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
    read_cals->read = read;

    return read_cals;
}

void read_cals_free(read_cals_t *read_cals){
  fastq_read_free(read_cals->read);
  array_list_free(read_cals->cal_list, (void *)cal_free);
  free(read_cals);
}

//------------------------------------------------------------------------------
//Options settings Burrows-Wheeler Transform 
//------------------------------------------------------------------------------

bwt_optarg_t *bwt_optarg_new(const size_t num_errors,
			     const size_t num_threads,
 			     const int filter_read_mappings, 
			     const int filter_seed_mappings,
				 const size_t min_cal_size,
				 const double umbral_cal_length_factor,
				 const int min_read_discard,
				 const int max_inner_gap) {
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

void bwt_optarg_free(bwt_optarg_t *optarg) {
  if (optarg != NULL) {
    free(optarg);
  }
}

//-----------------------------------------------------------------------------
// Index for the Burrows-Wheeler transform
//-----------------------------------------------------------------------------

char *bwt_error_type(char error_kind){

  char *msg_eror;
  switch(error_kind){
    case DELETION:
      msg_eror = "Deletion";
      break;
    case INSERTION:
      msg_eror = "Insertion";
      break;
    case MISMATCH:
      msg_eror = "Mismatch";
      break;
    default:
      msg_eror = "Unknown";
      break;
  }
  return msg_eror;
}

//-----------------------------------------------------------------------------
char *strcpy_capitalize(char *dest, const char *src, size_t n) {
  size_t i;
  
  for (i = 0 ; i < n && src[i] != '\0' ; i++) {
    if (src[i] > 96) {
      dest[i] = src[i] - 32;
    } else {
      dest[i] = src[i];
    }
  }

  for ( ; i < n ; i++) {
    dest[i] = '\0';
  }

  return dest;
}

//-----------------------------------------------------------------------------

void bwt_cigar_cpy(alignment_t *mapping, char *quality) {
  
  unsigned int quality_type;
  size_t quality_len;

  quality_len = strlen(quality);
  quality_type = atoi(mapping->quality);
  //printf("Quality len from batch: %i\n", quality_len);
  free(mapping->quality);
  mapping->quality = (char *)malloc(sizeof(char)*(quality_len + 1));
  //printf("Read:%s\n", mapping->sequence);
/*
  if (quality_type == START_HARD_CLIPPING){
    if (mapping->seq_strand == 0) {
      memcpy(mapping->quality, 
	     quality + 1, 
	     quality_len - 1);
    } else {
      reverse_str(quality + 1,
		  mapping->quality, quality_len - 1);
    }
    mapping->quality[quality_len - 1] = '\0';	    
    //printf("HARD START : %s\n", mapping->quality);
  }else if(quality_type == END_HARD_CLIPPING){
    if (mapping->seq_strand == 0) {
      memcpy(mapping->quality, 
	     quality,
	     quality_len - 1);
    } else {
      reverse_str(quality,
		  mapping->quality, quality_len - 1);
    }
    mapping->quality[quality_len - 1] = '\0';
    //printf("HARD END : %s\n", mapping->quality);
  }else{*/
    //printf("ELSE....\n");
  if (mapping->seq_strand == 0) {
       memcpy(mapping->quality, quality, quality_len);
  } else {
       reverse_str(quality,
		   mapping->quality, quality_len);
  }
  //mapping->quality[quality_len] = '\0';
  //printf("(%i)NORMAL : %s\n", mapping->seq_strand, mapping->quality);
//}
  //array_list_insert( mapping, mapping_list);
}
//-----------------------------------------------------------------------------

void bwt_cigar_cpy_batch(alignment_t *mapping, size_t read_i, fastq_batch_t *batch) {

  unsigned int quality_type;
  size_t quality_len;
  quality_len = batch->data_indices[read_i + 1] - batch->data_indices[read_i] - 1;
  quality_type = atoi(mapping->quality);
  //printf("Quality len from batch: %i\n", quality_len);                                                                                                                                                                                                                                                                                    
  free(mapping->quality);
  mapping->quality = (char *)malloc(sizeof(char)*(quality_len + 1));
  //printf("Read:%s\n", mapping->sequence);                                                                                                                                                                                                                                                                                                 

  if (quality_type == START_HARD_CLIPPING){
    if (mapping->seq_strand == 0) {
      memcpy(mapping->quality ,
             &(batch->quality[batch->data_indices[read_i]]) + 1,
             quality_len - 1);
    } else {
      reverse_str(&(batch->quality[batch->data_indices[read_i]]) + 1,
                  mapping->quality, quality_len - 1);
    }
    mapping->quality[quality_len - 1] = '\0';
    //printf("HARD START : %s\n", mapping->quality);                                                                                                                                                                                                                                                                                        
  }else if(quality_type == END_HARD_CLIPPING){
    if (mapping->seq_strand == 0) {
      memcpy(mapping->quality,
             &(batch->quality[batch->data_indices[read_i]]),
             quality_len - 1);
    } else {
      reverse_str(&(batch->quality[batch->data_indices[read_i]]),
                  mapping->quality, quality_len - 1);
    }
    mapping->quality[quality_len - 1] = '\0';
    //printf("HARD END : %s\n", mapping->quality);                                                                                                                                                                                                                                                                                          
  }else{
    //printf("ELSE....\n");                                                                                                                                                                                                                                                                                                                 
    if (mapping->seq_strand == 0) {
      memcpy(mapping->quality, &(batch->quality[batch->data_indices[read_i]]), quality_len);
    } else {
      reverse_str(&(batch->quality[batch->data_indices[read_i]]),
                  mapping->quality, quality_len);
    }
    //mapping->quality[quality_len] = '\0';                                                                                                                                                                                                                                                                                                 
    //printf("(%i)NORMAL : %s\n", mapping->seq_strand, mapping->quality);                                                                                                                                                                                                                                                                   
  }
  //array_list_insert( mapping, mapping_list);                                                                                                                                                                                                                                                                                              
}

//-----------------------------------------------------------------------------

unsigned int alignmentcmp(alignment_t *alignment_1, alignment_t *alignment_2) {
  
  size_t cigar1_len = strlen(alignment_1->cigar);
  size_t cigar2_len = strlen(alignment_2->cigar);

  if (alignment_1->map_quality > alignment_2->map_quality) { return 1; }
  else if (alignment_1->map_quality < alignment_2->map_quality) { return 2; }
  else {
    if (alignment_1->num_cigar_operations == 1 ||
	alignment_2->num_cigar_operations == 1) {
      if (alignment_1->cigar[cigar1_len - 1] == '=' && 
	  alignment_2->cigar[cigar1_len - 1] == '=') { return 0; }
      
      if (alignment_1->cigar[cigar1_len - 1] == '=') { return 1; }
      else if (alignment_2->cigar[cigar2_len - 1] == '=') { return 2; }
      
      if (alignment_1->cigar[cigar1_len - 1] == 'M' && 
	  alignment_2->cigar[cigar1_len - 1] == 'M') { return 0; }
      
      if (alignment_1->cigar[cigar1_len - 1] == 'M') { return 1; }
      else if (alignment_2->cigar[cigar2_len - 1] == 'M') { return 2; }
      
    }else {
      if (alignment_1->num_cigar_operations < 
	  alignment_2->num_cigar_operations) {
	return 1;
      } else if (alignment_1->num_cigar_operations > 
		 alignment_2->num_cigar_operations) {
	return 2;
      } else { return 0; }
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

void seq_reverse_complementary(char *seq, unsigned int len){
  unsigned int j = 0;
  
  char *seq_tmp = (char *)malloc(len*sizeof(char));
  memcpy(seq_tmp, seq, len);
  for (int i = len - 1; i >= 0; i--){
    if (seq_tmp[i] == 'A' || seq_tmp[i] == 'a') {
      seq[j] = 'T';
    }
    else if (seq_tmp[i] == 'C' || seq_tmp[i] == 'c'){
       seq[j] = 'G';
    }
    else if (seq_tmp[i] == 'G' || seq_tmp[i] == 'g'){
       seq[j] = 'C';
    }
    else if (seq_tmp[i] == 'T' || seq_tmp[i] == 't'){
       seq[j] = 'A';
    }else {
       seq[j] = 'N';
    }
    j++;  
  }

  free(seq_tmp);

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
  //ex.chromosome = (char *) calloc(INDEX_EXOME*IDMAX, sizeof(char));
  //ex.start = (unsigned int *) calloc(INDEX_EXOME, sizeof(unsigned int));
  //ex.end = (unsigned int *) calloc(INDEX_EXOME, sizeof(unsigned int));
  //ex.offset = (unsigned int *) calloc(INDEX_EXOME, sizeof(unsigned int));
  
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
// general functions
//-----------------------------------------------------------------------------
/*
size_t bwt_map_seq(char *seq, 
		   bwt_optarg_t *bwt_optarg, 
		   bwt_index_t *index, 
		   array_list_t *mapping_list) {
  
  if (bwt_optarg != NULL && bwt_optarg->num_errors == 0) {
    return bwt_map_exact_seq(seq, bwt_optarg, index, mapping_list);
  }

  return bwt_map_inexact_seq(seq, bwt_optarg, index, mapping_list);  
}
*/
//-----------------------------------------------------------------------------
/*
size_t bwt_map_read(fastq_read_t *read, 
		    bwt_optarg_t *bwt_optarg, 
		    bwt_index_t *index, 
		    array_list_t *mapping_list) {
  
  if (bwt_optarg != NULL && bwt_optarg->num_errors == 0) {
    return bwt_map_exact_read(read, bwt_optarg, index, mapping_list);
  }

  return bwt_map_inexact_read(read, bwt_optarg, index, mapping_list);  
}
*/
//-----------------------------------------------------------------------------

/*size_t bwt_map_seqs(char **seqs, 
		    size_t num_reads,
		    bwt_optarg_t *bwt_optarg, 
		    bwt_index_t *index, 
		    char *out_status,
		    array_list_t *mapping_list) {

  if (bwt_optarg != NULL && bwt_optarg->num_errors == 0) {
    return bwt_map_exact_seqs(seqs, num_reads, bwt_optarg, index, 
			      out_status, mapping_list);
  }
  
  return bwt_map_inexact_seqs(seqs, num_reads, bwt_optarg, index, 
			      out_status, mapping_list);  
			      }*/

//-----------------------------------------------------------------------------
/*
size_t bwt_map_batch(fastq_batch_t *batch,
		     bwt_optarg_t *bwt_optarg, 
		     bwt_index_t *index, 
		     fastq_batch_t *unmapped_batch,
		     array_list_t *mapping_list) {
  
  if (bwt_optarg != NULL && bwt_optarg->num_errors == 0) {
    return bwt_map_exact_batch(batch, bwt_optarg, index, 
			       unmapped_batch, mapping_list);
  }

  return bwt_map_inexact_batch(batch, bwt_optarg, index, 
			       unmapped_batch, mapping_list);  
}*/

//-----------------------------------------------------------------------------
// exact functions
//-----------------------------------------------------------------------------

size_t bwt_map_exact_seq(char *seq, 
			 bwt_optarg_t *bwt_optarg, 
			 bwt_index_t *index, 
			 array_list_t *mapping_list) {
  printf("\tIn function ...\n");
  unsigned int len = strlen(seq);
  int start = 0;
  int end = len - 1;
  char *code_seq = (char *) malloc(len * sizeof(char));
  result result;
  size_t l_aux, k_aux;
  unsigned int num_mappings = 0;
  char plusminus[2] = "-+";
  size_t idx, key, error, pos;
  alignment_t *alignment;
  char *cigar_p, *seq_dup, *seq_strand;
  unsigned int start_mapping;
  //char *chromosome = (char *)malloc(sizeof(char)*10);
  seq_strand = strdup(seq);
  // be careful, now we use the table from the index
  //replaceBases(seq, code_seq, len);
  bwt_encode_Bases(code_seq, seq, len, &index->table);
  struct timeval t_start, t_end;
  char *optional_fields /*= (char *)malloc(sizeof(char)*256)*/;
  //printf("---> EXACT Search Read (%d): %s\n", len, seq);

  
  for (short int type = 1; type >= 0; type--) {
    result.k = 0;
    result.l = index->h_O.siz - 2;
    result.start = start;
    result.end = end;
    
    if (type == 1) {
      result.pos = end;
      //printf("Before call BWT\n");
      start_timer(t_start);
      BWExactSearchBackward(code_seq, &index->h_C, &index->h_C1, &index->h_O, &result);
      stop_timer(t_start, t_end, time_bwt);
      //printf("After call BWT\n");
    } else {
      result.pos = start;
      start_timer(t_start);
      BWExactSearchForward(code_seq, &index->h_rC, &index->h_rC1, &index->h_rO, &result);
      stop_timer(t_start, t_end, time_bwt);
      if(l_aux - k_aux + 1 < bwt_optarg->filter_read_mappings) {
	seq_reverse_complementary(seq_strand, len);
      }
      
    }
      
    start_timer(t_start);
    k_aux = result.k;
    l_aux = result.l;
    //printf("\tk=%d - l=%d\n", k_aux, l_aux);      
    if (l_aux - k_aux + 1 >= bwt_optarg->filter_read_mappings) {
      l_aux = k_aux + 1;
    }
    
    
    for (size_t j = k_aux; j <= l_aux; j++) {
      
      if (index->S.ratio == 1) {
	key = index->S.vector[j];
      } else {
	key = getScompValue(j, &index->S, &index->h_C, &index->h_O);
      }
      
      idx = binsearch(index->karyotype.offset, index->karyotype.size, key);
      //chromosome = index->karyotype.chromosome + (idx-1) * IDMAX;
      
      if(key + len <= index->karyotype.offset[idx]) {
	start_mapping = index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]);
	
	cigar_p = (char *)malloc(sizeof(char)*len);
	sprintf(cigar_p, "%d=\0", len);
	
	// save all into one alignment structure and insert to the list
	alignment = alignment_new();
	alignment_init_single_end(NULL, strdup(seq_strand), NULL, !type, 
				  idx - 1,				      
				  start_mapping, 
				  cigar_p, 1, 254, 1, (num_mappings > 0), 0, optional_fields, alignment);
	  
	if(!array_list_insert((void*) alignment, mapping_list)){
	  printf("Error to insert item into array list\n");
	}
	
	num_mappings++;
      }
    }
    
    stop_timer(t_start, t_end, time_search);
  }
  
  /*
  printf("Time BWT Search %.3fs\n", time_bwt/1000000);
  printf("Time S Search %.3fs\n", time_search/1000000);
  */
  free(seq_strand);
  free(code_seq);

  //printf("\tOut function\n");
  return num_mappings;
}

//-----------------------------------------------------------------------------

size_t bwt_map_exact_read(fastq_read_t *read, 
			  bwt_optarg_t *bwt_optarg, 
			  bwt_index_t *index, 
			  array_list_t *mapping_list) {
  
  size_t padding = array_list_size(mapping_list);
  
  size_t num_mappings = bwt_map_exact_seq(read->sequence, 
					  bwt_optarg, index, 
					  mapping_list);
  alignment_t *mapping;
  
  unsigned int header_len, quality_len;
  
  for (int i = 0; i < num_mappings ; i++) {
    //printf("access to pos:%d-padding:%d\n", i + padding, padding);
    mapping = (alignment_t *) array_list_get(i + padding, mapping_list);
    if (mapping != NULL) {
      header_len = strlen(read->id) + 1;
      mapping->query_name = (char *)malloc(sizeof(char)*header_len);
      get_to_first_blank(read->id, header_len, mapping->query_name);
      //printf("--->header id %s to %s\n",read->id, header_id);
      //free(read->id);
      
      quality_len = strlen(read->quality) + 1;
      mapping->quality = (char *)malloc(sizeof(char)*quality_len);
      if (mapping->seq_strand == 0) {
	memcpy(mapping->quality, read->quality, quality_len);
      }else {
	reverse_str(read->quality,
		    mapping->quality, quality_len);
      }
      mapping->quality[quality_len - 1] = '\0';
      //mapping->quality = read->quality;
    } else {
      printf("Error to extract item\n");
    }
    //printf("Store Ok!\n"); 
    //printf("\tSEQ:%s-ID:%s\n", mapping->sequence, mapping->query_name);
  }
  
  return num_mappings;
}

//-----------------------------------------------------------------------------

size_t bwt_map_exact_seqs(char **seqs, 
			  size_t num_reads,
			  bwt_optarg_t *bwt_optarg, 
			  bwt_index_t *index, 
			  char *out_status,
			  array_list_t *mapping_list) {
  
  unsigned int num_threads = bwt_optarg->num_threads;
  array_list_t **individual_mapping_list_p = (array_list_t **)malloc(sizeof(array_list_t *)*num_threads);
  const size_t chunk = MAX(1, num_reads/( num_threads*10));
  size_t num_mappings = 0; 

  short int th_id;
  alignment_t *mapping;
  
  for(int i = 0; i < num_threads; i++){
   individual_mapping_list_p[i]  = array_list_new(100000/num_threads, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
  }
  
  omp_set_num_threads(num_threads);
  
  //printf("Initialization ok!Process %d reads x %d threads\n", num_reads, num_threads);
  
  #pragma omp parallel for private(th_id, num_mappings) schedule(dynamic, chunk)
  for(int i = 0; i < num_reads; i++){
    th_id = omp_get_thread_num();
    //printf("I'm thread %d and process read %d\n", th_id,  i);
    num_mappings = bwt_map_exact_seq(seqs[i], bwt_optarg, index, individual_mapping_list_p[th_id]);
    if(num_mappings){
      out_status[i] = 1;
    }else{
      out_status[i] = 0;
    }
    num_mappings = 0;
  }

  
  //printf("Process reads ok! Total mappings %d\n\n", num_mappings_tot);
  
  //Join results
  for(int i = 0; i < num_threads; i++ ){
    num_mappings = array_list_size(individual_mapping_list_p[i]);
   
    for(int j = 0; j < num_mappings; j++){
      mapping = (alignment_t *)array_list_get(j, individual_mapping_list_p[i]);
      array_list_insert( mapping, mapping_list);
    }
  }

  free(individual_mapping_list_p);

  return array_list_size(mapping_list);
}


//-----------------------------------------------------------------------------

size_t bwt_map_exact_batch(fastq_batch_t *batch,
			   bwt_optarg_t *bwt_optarg, 
			   bwt_index_t *index, 
			   fastq_batch_t *unmapped_batch,
			   array_list_t *mapping_list) {
  
  //size_t length_header, length_seq;
  //fastq_read_t *fq_read;
  //char header[1024], seq[1024], quality[1024]; 
  size_t num_mappings, num_mappings_tot; //num_unmapped;
  size_t read_pos;
  size_t num_threads = bwt_optarg->num_threads;
  //unsigned int th_id;
  size_t num_reads = batch->num_reads;
  //size_t batch_individual_size = batch->data_size / num_threads;
  size_t chunk = MAX(1, num_reads/(num_threads*10));

  array_list_t **individual_mapping_list_p   = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  //array_list_t **individual_unmapping_list_p = (array_list_t **)malloc(sizeof(array_list_t *)*num_threads);
  
  alignment_t *mapping;
  
  //struct timeval start_time, end_time;
  //double timer, parallel_t, sequential_t;
  unsigned int j, header_len, quality_len;
  size_t read_id = 0;
  for (unsigned int i = 0; i < num_reads; i++) {
    individual_mapping_list_p[i] = array_list_new(bwt_optarg->filter_read_mappings, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
    //individual_unmapping_list_p[i] = array_list_new(100000/num_threads, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
  }  
  
  //printf("%d Threads %d chunk\n", num_threads, chunk);
  //num_mappings = 0;
  omp_set_num_threads(num_threads);
  
  //  start_timer(start_time);
  //printf("Process %d reads\n", num_reads);
  #pragma omp parallel for schedule(dynamic, chunk)
  for (unsigned int i = 0; i < num_reads; i++) {
 
    bwt_map_exact_seq(&(batch->seq[batch->data_indices[i]]), 
		      bwt_optarg, index, 
		      individual_mapping_list_p[i]);
  }

  //printf("Out mappings = %d\n", num_mappings_tot);

  read_pos = 0;
  unmapped_batch->header_indices[read_pos] = 0;
  unmapped_batch->data_indices[read_pos] = 0;
  for (unsigned int i = 0; i < num_reads; i++) {
    num_mappings = array_list_size(individual_mapping_list_p[i]);
    if(num_mappings){
      //printf("@DNA\n%s\n+\n%s\n", &(batch->seq[batch->data_indices[i]]), &(batch->seq[batch->data_indices[i]]));
      for (unsigned int j = 0; j < num_mappings; j++) {
	mapping = (alignment_t *) array_list_get(j, individual_mapping_list_p[i]);
	if(mapping != NULL){
	  header_len = batch->header_indices[i + 1] - batch->header_indices[i];
	  mapping->query_name = (char *)malloc(sizeof(char)*header_len);
	  get_to_first_blank(&(batch->header[batch->header_indices[i]]), header_len, mapping->query_name);
	  //printf("--->header id %s to %s\n",read->id, header_id);
	  
	  quality_len = batch->data_indices[i + 1] - batch->data_indices[i];
	  mapping->quality = (char *)malloc(sizeof(char)*quality_len);
	  memcpy(mapping->quality, &(batch->quality[batch->data_indices[i]]), quality_len);
	  //array_list_insert( mapping, mapping_list);
	  array_list_insert( mapping, mapping_list);
	}else{
	  printf("Error to extract item\n");
	}
      }
    }else{
      //fq_read = (fastq_read_t *)array_list_get(j, individual_unmapping_list_p[i]);
      read_pos++;
      header_len = batch->header_indices[i + 1] - batch->header_indices[i];
      quality_len = batch->data_indices[i + 1] - batch->data_indices[i];
      
      memcpy(&(unmapped_batch->header[unmapped_batch->header_indices[read_pos - 1]]), &(batch->header[batch->header_indices[i]]), header_len);
      memcpy(&(unmapped_batch->seq[unmapped_batch->data_indices[read_pos - 1]]),      &(batch->seq[batch->data_indices[i]]),      quality_len);
      memcpy(&(unmapped_batch->quality[unmapped_batch->data_indices[read_pos - 1]]),  &(batch->quality[batch->data_indices[i]]),  quality_len);
      
      unmapped_batch->data_indices[read_pos]   = unmapped_batch->data_indices[read_pos - 1]   + quality_len;
      unmapped_batch->header_indices[read_pos] = unmapped_batch->header_indices[read_pos - 1] + header_len;
      
    }
  }
  //printf("\tReads unmapped %d\n", read_pos);  
  unmapped_batch->num_reads = read_pos;

  for (unsigned int i = 0; i < num_reads; i++) {
    array_list_free(individual_mapping_list_p[i], NULL);
    //individual_unmapping_list_p[i] = array_list_new(100000/num_threads, 1.25f, COLLECTION_MODE_SYNCHRONIZED);
  }  

  free(individual_mapping_list_p);
  //free(individual_unmapping_list_p);
  
  return array_list_size(mapping_list);

}

//-----------------------------------------------------------------------------

size_t bwt_map_exact_seed(char *seq, size_t seq_len,
			  size_t seq_start, size_t seq_end,
			  bwt_optarg_t *bwt_optarg, 
			  bwt_index_t *index,
			  array_list_t *mapping_list,
			  int id/*, int *limit_exceeded*/) {
  
  //printf("Process New Seeds\n");

  region_t *region;
  size_t start = 0;
  size_t end = seq_end - seq_start;
  result result;
  char *code_seq = &seq[seq_start];
  size_t len = seq_end - seq_start;
  size_t num_mappings = 0;
  char plusminus[2] = "-+";
  size_t idx, key, direction, error, pos;
  results_list *r_list;
  //  result *r;
  size_t l_aux, k_aux;
  alignment_t *alignment;
  //  size_t len = strlen(seq);
  //  size_t start = 0;
  //  size_t end = len - 1;
  size_t start_mapping;
  size_t aux_seq_start, aux_seq_end;
  array_list_t *mappings = array_list_new(bwt_optarg->filter_read_mappings, 1.25f, 
						     COLLECTION_MODE_ASYNCHRONIZED);
  int discard_seed = 0;
  int actual_mappings = 0;
  struct timeval t_start, t_end;
  //*limit_exceeded = 0;
  for (short int type = 1; type >= 0; type--) {
    result.k = 0;
    result.l = index->h_O.siz - 2;
    result.start = start;
    result.end = end;
    if (type == 1) {
      // strand +
      aux_seq_start = seq_start;
      aux_seq_end = seq_end;
      result.pos = end;	
      start_timer(t_start);
      BWExactSearchBackward(code_seq, &index->h_C, &index->h_C1, &index->h_O, &result);
      //BWExactSearchBackward(code_seq, start, end, &index->h_C, &index->h_C1, &index->h_O, result_p);
      stop_timer(t_start, t_end, time_bwt_seed);
    } else {
      // strand -
      aux_seq_start = seq_len - seq_end - 1;
      aux_seq_end = seq_len - seq_start - 1;
      //printf("Translate coords %i-%i(id %i)\n", aux_seq_start, aux_seq_end, id);
      result.pos = start;      
      start_timer(t_start);
      BWExactSearchForward(code_seq, &index->h_rC, &index->h_rC1, &index->h_rO, &result);
      //BWExactSearchForward(code_seq, start, end, &index->h_rC, &index->h_rC1, &index->h_rO, result_p);
      stop_timer(t_start, t_end, time_bwt_seed);
    }

    //start_timer(t_start);

    k_aux = result.k;
    l_aux = result.l;
    actual_mappings += (result.l - result.k + 1);

    if (actual_mappings > bwt_optarg->filter_seed_mappings) {
      //discard_seed = 1;
      //printf("Limit exceded LIMIT = %i, MAP = %i\n", bwt_optarg->filter_seed_mappings, actual_mappings);
      continue;
      k_aux = result.k;
      l_aux = result.k + 10;      
      //limit_exceeded = 1;
      break;
    } else {
	//printf("\tk=%d - l=%d\n", r->k, r->l);      
      k_aux = result.k;
      l_aux = result.l;
    }
      
    for (size_t j = k_aux; j <= l_aux; j++) {            
      key = getScompValue(j, &index->S, &index->h_C, &index->h_O);
      idx = binsearch(index->karyotype.offset, index->karyotype.size, key);     
      
      if (key + len <= index->karyotype.offset[idx]) {
	start_mapping = index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]) + 1;
	// save all into one alignment structure and insert to the list
	//printf("\t\t Region%i[%i:%lu-%lu]\n", !type, idx, start_mapping, start_mapping + len);
	region = region_bwt_new(idx, !type, start_mapping, start_mapping + len, aux_seq_start, aux_seq_end, seq_len, id);
	
	if (!array_list_insert((void*) region, mapping_list)){
	  printf("Error to insert item into array list\n");
	}
	
	num_mappings++;
      }
    }
    //stop_timer(t_start, t_end, time_search_seed);
  }
  //  free(result_p);
  /*
  if (discard_seed) {
    array_list_clear(mappings, region_bwt_free);
  } else {
    for (int i = num_mappings - 1; i >= 0; i--) {
      region = array_list_remove_at(i, mappings);
      array_list_insert(region, mapping_list);
    }
  }
  */
  array_list_free(mappings, NULL);

  return num_mappings;  
}

size_t bwt_map_exact_seed_by_region(char *seq, size_t seq_len,
				    size_t seq_start, size_t seq_end,
				    bwt_optarg_t *bwt_optarg, 
				    bwt_index_t *index,
				    array_list_t *mapping_list,
				    int id, 
				    int strand_target,
				    int chromosome_target,
				    size_t start_target, 
				    size_t end_target) {  
  region_t *region;
  size_t start = 0;
  size_t end = seq_end - seq_start;
  result result;
  char *code_seq = &seq[seq_start];

  size_t len = seq_end - seq_start;
  size_t num_mappings = 0;
  char plusminus[2] = "-+";
  size_t idx, key, direction, error, pos;
  results_list *r_list;
  size_t l_aux, k_aux;
  alignment_t *alignment;

  size_t start_mapping;
  size_t aux_seq_start, aux_seq_end;
  array_list_t *mappings = array_list_new(bwt_optarg->filter_read_mappings, 1.25f, 
						     COLLECTION_MODE_ASYNCHRONIZED);
  int discard_seed = 0;
  int actual_mappings = 0;
  struct timeval t_start, t_end;

  for (short int type = 1; type >= 0; type--) {
    if (!type != strand_target) { continue; }

    result.k = 0;
    result.l = index->h_O.siz - 2;
    result.start = start;
    result.end = end;
    if (type == 1) {
      // strand +
      aux_seq_start = seq_start;
      aux_seq_end = seq_end;
      result.pos = end;	
      start_timer(t_start);
      BWExactSearchBackward(code_seq, &index->h_C, &index->h_C1, &index->h_O, &result);
      //BWExactSearchBackward(code_seq, start, end, &index->h_C, &index->h_C1, &index->h_O, result_p);
      stop_timer(t_start, t_end, time_bwt_seed);
    } else {
      // strand -
      aux_seq_start = seq_len - seq_end - 1;
      aux_seq_end = seq_len - seq_start - 1;

      result.pos = start;      
      start_timer(t_start);
      BWExactSearchForward(code_seq, &index->h_rC, &index->h_rC1, &index->h_rO, &result);
      //BWExactSearchForward(code_seq, start, end, &index->h_rC, &index->h_rC1, &index->h_rO, result_p);
      stop_timer(t_start, t_end, time_bwt_seed);
    }

    //start_timer(t_start);

    k_aux = result.k;
    l_aux = result.l;
    actual_mappings += (result.l - result.k + 1);

    if (actual_mappings > 500) {
      k_aux = result.k;
      l_aux = result.k + 10;
    } else {
      k_aux = result.k;
      l_aux = result.l;
    }
      
    for (size_t j = k_aux; j <= l_aux; j++) {            
      key = getScompValue(j, &index->S, &index->h_C, &index->h_O);
      idx = binsearch(index->karyotype.offset, index->karyotype.size, key);     
      if (idx != chromosome_target) { continue; }

      if (key + len <= index->karyotype.offset[idx]) {
	start_mapping = index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]) + 1;
	if (start_mapping < start_target || start_mapping > end_target) { continue; }
	// save all into one alignment structure and insert to the list
	printf(" \t::: %lu:%lu :::\n", idx, start_mapping);
	region = region_bwt_new(idx, !type, start_mapping, start_mapping + len, aux_seq_start, aux_seq_end, seq_len, id);
	
	if (!array_list_insert((void*) region, mapping_list)){
	  printf("Error to insert item into array list\n");
	}
	
	num_mappings++;
      }
    }
    //stop_timer(t_start, t_end, time_search_seed);
  }
  //  free(result_p);
  array_list_free(mappings, NULL);

  return num_mappings;  
}


//-----------------------------------------------------------------------------
// inexact functions
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//Ricardo - Modified in order to allow internal bwt intervals (seq_start, seq_end)
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

size_t bwt_map_inexact_read(fastq_read_t *read, 
			    bwt_optarg_t *bwt_optarg, 
			    bwt_index_t *index, 
			    array_list_t *mapping_list) {

  return __bwt_map_inexact_read(read, 
			       bwt_optarg, 
			       index, 
			       mapping_list);  
}

size_t __bwt_map_inexact_read(fastq_read_t *read, 
			      bwt_optarg_t *bwt_optarg, 
			      bwt_index_t *index, 
			      array_list_t *mapping_list) {
   
     char *seq = read->sequence;
     alignment_t *alignment;
     size_t len = strlen(seq);

     if (len < 5) {
	  char aux[len + 2];
	  sprintf(aux, "%luX", len);

	  char *quality_clipping = (char *)malloc(sizeof(char)*50);
	  sprintf(quality_clipping, "%i", NONE_HARD_CLIPPING);

	  alignment = alignment_new();
	  alignment_init_single_end(NULL,
				    strdup(seq),
				    quality_clipping,
				    0,
				    -1,
				    -1,
				    strdup(aux), 1, 0, 0, 0, 0, NULL, alignment);
	  array_list_insert((void*) alignment, mapping_list);
	  return 1;
     }
 
     char *seq_dup, *seq_strand;
     size_t start = 0;
     size_t end = len - 1;
     size_t len_calc = len;

     char *code_seq = (char *) malloc(len * sizeof(char));

     // be careful, now we use the table from the index
     //replaceBases(seq, code_seq, len);
     bwt_encode_Bases(code_seq, seq, len, &index->table);

     // calculate vectors k and l
     size_t *k0 = (size_t *) malloc(len * sizeof(size_t));
     size_t *l0 = (size_t *) malloc(len * sizeof(size_t));
     size_t *k1 = (size_t *) malloc(len * sizeof(size_t));
     size_t *l1 = (size_t *) malloc(len * sizeof(size_t));
     size_t *ki0 = (size_t *) malloc(len * sizeof(size_t));
     size_t *li0 = (size_t *) malloc(len * sizeof(size_t));
     size_t *ki1 = (size_t *) malloc(len * sizeof(size_t));
     size_t *li1 = (size_t *) malloc(len * sizeof(size_t));


     size_t last_k0, last_l0;
     size_t last_ki0, last_li0;

     size_t last_k1, last_l1;
     size_t last_ki1, last_li1;

     int back0_nt, forw0_nt, back1_nt, forw1_nt;

     BWExactSearchVectorBackward(code_seq, start, end, 0, index->h_O.siz - 2,
				 k1, l1, &index->h_C, &index->h_C1, &index->h_O,
				 &last_k1, &last_l1, &back1_nt);
     
     BWExactSearchVectorForward(code_seq, start, end, 0, index->h_Oi.siz - 2,
				ki1, li1, &index->h_C, &index->h_C1, &index->h_Oi,
				&last_ki1, &last_li1, &forw1_nt);
  
     BWExactSearchVectorForward(code_seq, start, end, 0, index->h_rO.siz - 2,
				k0, l0, &index->h_rC, &index->h_rC1, &index->h_rO,
				&last_k0, &last_l0, &forw0_nt);
  
     BWExactSearchVectorBackward(code_seq, start, end, 0, index->h_rOi.siz - 2,
				 ki0, li0, &index->h_rC, &index->h_rC1, &index->h_rOi,
				 &last_ki0, &last_li0, &back0_nt);
  

     //compare the vectors k and l to get mappings in the genome
     size_t num_mappings = 0;
     char plusminus[2] = "-+";
     size_t idx, key, direction;
     char error;
     int pos;
     //results_list *r_list;
     results_list r_list;
     result *r;
     char error_debug = 0;
     char *cigar_dup;
     char cigar[1024];
     size_t cigar_len, num_cigar_ops;
     array_list_t *mapping_list_filter;
     alignment_t *best_alignment, *aux_alignment;
     size_t best_pos, array_size;
     //size_t i, j, z;
     size_t *allocate_pos_alignments;
     size_t k_start, l_start;

     size_t start_mapping;
     size_t tot_alignments = 0;
     const int MAX_BWT_ALIGNMENTS = 10;
     int filter_exceeded = 0;
     //seq_dup = (char *)malloc(sizeof(char)*(len + 1));
     seq_strand = strdup(seq);
     error = MISMATCH;


     char *quality_clipping = (char *) malloc(sizeof(char) * 50);
     seq_dup = (char *) malloc(sizeof(char) * (len + 1));

     new_results_list(&r_list, bwt_optarg->filter_read_mappings);

     array_list_t *tmp_mapping_list = array_list_new(bwt_optarg->filter_read_mappings, 1.25f, 
						     COLLECTION_MODE_ASYNCHRONIZED);

     for (int type = 1; type >= 0; type--) {
       //printf("Process Strand(%i)\n", !type);
	  r_list.num_results = 0;
	  r_list.read_index = 0;

	  //    printf("*** bwt.c: calling BWSearch1 with type = %d...\n", type);
	  if (type == 1) {
	       BWSearch1(code_seq, start, end, k1, l1, ki1, li1, 
			 &index->h_C, &index->h_C1, &index->h_O, &index->h_Oi, &r_list);
	  } else {
	       BWSearch1(code_seq, start, end, ki0, li0, k0, l0, 
			 &index->h_rC, &index->h_rC1, &index->h_rOi, &index->h_rO, &r_list);      
	       if (r_list.num_results) {
		    seq_reverse_complementary(seq_strand, len);
	       }
	  }

	  //printf("*** bwt.c: calling BWSearch1 with type = %d (num_results = %d). Done !!\n", type, r_list.num_results);

	  for (size_t ii = 0; ii < r_list.num_results; ii++) {
	    r = &r_list.list[ii];

	       //printf("Errors number %d\n", r->num_mismatches);
	       //if(r == NULL){printf("ERROR: Result list position null");exit(-1);}
	       if(!r->num_mismatches)
		    error = 0;
	       else
		    error = r->err_kind[0];

	       pos = r->err_pos[0];
	       //pos = r->position[0];
	  

	       if (type) {
		    direction = r->dir;
	       } else {
		    direction = !r->dir;
	       }

	       len_calc = len;
	       if (error == DELETION) {
		    len_calc--;
	       } else if (error == INSERTION) {
		    len_calc++;
	       }

	       // generating cigar
	       sprintf(quality_clipping, "%i", NONE_HARD_CLIPPING);
	       if (error == 0) {
		 sprintf(cigar, "%lu=\0", len);
		 num_cigar_ops = 1;
		 memcpy(seq_dup, seq_strand, len);
		 seq_dup[len] = '\0';
		 
	       } else if (error == MISMATCH) {
		 if (pos == 0) {
		   //Positive strand
		   if(type) { 
		     sprintf(cigar, "1S%luM\0", len-1); 
		     start_mapping++;
		   }
		   else { 
		     sprintf(cigar, "%luM1S\0", len-1); 
		   }
		   num_cigar_ops = 2;
		 } else if (pos == len - 1) {
		   //Positive strand
		   if(type) { 
		     sprintf(cigar, "%luM1S\0", len - 1); 
		   }
		   else{ 
		     sprintf(cigar, "1S%luM\0", len-1); 
		     start_mapping++;
		   }
		   num_cigar_ops = 2;
		 } else {
		   sprintf(cigar, "%luM\0", len);
		   num_cigar_ops = 1;
		 }
		 memcpy(seq_dup, seq_strand, len);
		 seq_dup[len] = '\0';
		 //printf("MISMATCH\n");
		 
	       } else if (error == INSERTION) {
		 //printf("INSERTION\n");
		 if (pos == 0) {
		   if(type) {
		     sprintf(cigar, "1M1D%luM\0", len - 1); 
		   }
		   else{ 
		     sprintf(cigar, "%luM1D1M\0", len - 1); 
		   }	      
		 } else if (pos == len - 1) {
		   if(type) { 
		     sprintf(cigar, "%luM1D1M\0", len - 1); 
		   }
		   else{ 
		     sprintf(cigar, "1M1D%luM\0", len - 1); 
		   }
		 } else {
		   
		   if(type) {
		     if(r->dir)
		       sprintf(cigar, "%iM1D%luM\0", pos, len - pos);
		     else
		       sprintf(cigar, "%iM1D%luM\0", pos + 1, len - pos - 1);
		   } else { 
		     if(r->dir)
		       sprintf(cigar, "%luM1D%dM\0", len - pos, pos);
		     else
		       sprintf(cigar, "%luM1D%dM\0", len - pos - 1, pos + 1);
		   }
		 }
		 num_cigar_ops = 3;
		 memcpy(seq_dup, seq_strand, len);
		 seq_dup[len] = '\0';
	       } else if (error == DELETION) {	     
		 //printf("DELETION\n");
		 if (pos == 0) {
		   if(type) { sprintf(cigar, "1I%luM\0", len -1); }
		   else{ 
		     sprintf(cigar, "%luM1I\0", len -1); 
		     //		   start_mapping++;
		   }
		   
		   num_cigar_ops = 2;		
		 } else if (pos == len - 1) {
		   if(type) { 
		     sprintf(cigar, "%luM1I\0", len -1); 
		     //		   start_mapping++;
		   }
		   else{ 
		     sprintf(cigar, "1I%luM\0", len -1); 
		   }
		   num_cigar_ops = 2;
		 } else {
		   if(type) { 
		     sprintf(cigar, "%dM1I%luM\0", pos, len - pos - 1); 
		   }
		   else{ 
		     sprintf(cigar, "%luM1I%dM\0", len - pos - 1, pos); 
		   }
		   num_cigar_ops = 3;
		 }
		 memcpy(seq_dup, seq_strand , len );
		 seq_dup[len] = '\0';
		 
	       }else{
		 printf("NUM MAPPINGS %lu -> POS %d -> ERROR %d -> (%lu):%s", num_mappings, pos, error, len, seq);
		 continue;
		 //exit(-1);
		 //error_debug = 1;
	       }
	       k_start = r->k;
	       l_start = r->l;
	       //}

	       //printf("%lu-%lu :: num_results %lu\n", k_start, l_start, r_list.num_results);
	       tot_alignments += (l_start - k_start);
	       if (tot_alignments >  bwt_optarg->filter_read_mappings) {
		 filter_exceeded = 1;
		 LOG_DEBUG_F("Filter exceeded: num. read mappings: %i (total: %i > filter %i)\n", 
			     l_start - k_start, tot_alignments, bwt_optarg->filter_read_mappings);
		 break;
	       }

	       for (size_t j = k_start; j <= l_start; j++) {
		    if (index->S.ratio == 1) {
			 key = (direction)
			      ? index->Si.siz - index->Si.vector[j] - len_calc - 1
			      : index->S.vector[j];
		    } else {
			 key = (direction)
			      ? index->Si.siz - getScompValue(j, &index->Si, &index->h_C,
							      &index->h_Oi) - len_calc - 1
			      : getScompValue(j, &index->S, &index->h_C, &index->h_O);
		    }
		    idx = binsearch(index->karyotype.offset, index->karyotype.size, key);
		    if(key + len_calc <= index->karyotype.offset[idx]) {
		         start_mapping = index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]);
			 // save all into one alignment structure and insert to the list
			 //printf("*****Alignments %i:%lu\n", idx, start_mapping);
			 alignment = alignment_new();
			 alignment_init_single_end(NULL, strdup(seq_dup), strdup(quality_clipping), !type, 
						   idx - 1, //index->karyotype.chromosome + (idx-1) * IDMAX,
						   start_mapping + 1, //index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]), 
						   strdup(cigar), num_cigar_ops, 254, 1, (num_mappings > 0), 0, NULL, alignment);

			 array_list_insert((void*) alignment, tmp_mapping_list);
		    }
	       }//end for k and l
	  }//end for 
	  
	  if (filter_exceeded) {
	    array_list_clear(tmp_mapping_list, (void *)alignment_free);
	    array_list_set_flag(2, mapping_list);
	    break;
	  }

	  //
	  //
	  //Search for equal BWT mappings and set the mappings that will be delete
	  int n_mappings = array_list_size(tmp_mapping_list);
	  alignment_t *alig_1, *alig_2;
	  unsigned int *delete_mark = (unsigned int *)calloc(n_mappings, sizeof(unsigned int));
	  const int max_distance = 10;
	  //printf("------------------Num mappings %i---------------\n", n_mappings);
     
	  for (int a1 = n_mappings - 1; a1 >= 1; a1--) {
	    if (!delete_mark[a1]) {
	      alig_1 = array_list_get(a1, tmp_mapping_list);
	      //printf("Alig1(%i): %i chromosome, %i seq_strand, %s cigar:\n", a1, alig_1->chromosome, alig_1->seq_strand, alig_1->cigar);
	      for (int a2 = a1 - 1; a2 >= 0; a2--) {
		alig_2 = array_list_get(a2, tmp_mapping_list);
		size_t dist = abs(alig_1->position - alig_2->position);
		//printf("\t Alig2(%i): %i chromosome, %i seq_strand, %i dist, %i delete mark, %s cigar\n", a2, alig_2->chromosome, alig_2->seq_strand, dist, delete_mark[a2],alig_2->cigar );
		if (alig_1->chromosome == alig_2->chromosome &&
		    dist < max_distance && 
		    !delete_mark[a2]) {
		  //Same chromosome && same position
		  if (alig_1->num_cigar_operations < alig_2->num_cigar_operations) {
		    //printf("\tSet read %i\n", a2);
		    delete_mark[a2] = 1;
		  } else {
		    //printf("\tSet read %i\n", a2);
		    delete_mark[a1] = 1;
		  }
		}
	      }
	    }
	  }
	  //Delete all set mappings
	  int primary_delete = 0;
	  for (int m = n_mappings - 1; m >= 0; m--) {
	    alig_1 = array_list_remove_at(m, tmp_mapping_list);
	    if (delete_mark[m]) {
	      if (!is_secondary_alignment(alig_1)) { primary_delete = 1; }
	      alignment_free(alig_1);
	    } else {
	      array_list_insert((void*) alig_1, mapping_list);
	      set_secondary_alignment(num_mappings > 0, alig_1);
	      num_mappings++;
	    }
	  }
	  if (primary_delete) { 
	    alig_1 = array_list_get(0, mapping_list); 
	    set_secondary_alignment(0, alig_1);
	  }
	  free(delete_mark);	  	  
     } // end for type 

     if (array_list_size(mapping_list) == 0) {
       if (array_list_get_flag(mapping_list) == 1) {
	 //====================================================================================
	 //printf("BWT(+): FORWARD(k-l)%i:%lu-%lu BACKWARD(k-l)%i:%lu-%lu\n", forw1_nt, last_ki1, last_li1, back1_nt, last_k1, last_l1);
	 //printf("BWT(-): FORWARD(k-l)%i:%lu-%lu BACKWARD(k-l)%i:%lu-%lu\n", forw0_nt, last_ki0, last_li0, back0_nt, last_k0, last_l0);
	 array_list_t *forward_anchor_list, *backward_anchor_list;       
	 int type = 1;//(+)

	 //printf("BACKWARD (+)\n");
	 __bwt_generate_anchor_list(last_k1, last_l1, back1_nt, bwt_optarg, 
				    index, type, mapping_list, BACKWARD_ANCHOR, 0, start, end);
	 //printf("FORWARD (+)\n");
	 __bwt_generate_anchor_list(last_ki1, last_li1, forw1_nt, bwt_optarg, 
				    index, type,  mapping_list, FORWARD_ANCHOR, 0, start, end);
	 type = 0;//(-)
	 //printf("BACKWARD (-)\n");
	 /* __bwt_generate_anchor_list(last_k0, last_l0, back0_nt, bwt_optarg, 
				    index, type,  mapping_list, BACKWARD_ANCHOR, back0_nt - forw0_nt);
	 //printf("FORWARD (-)\n");
	 __bwt_generate_anchor_list(last_ki0, last_li0, forw0_nt, bwt_optarg, 
	 index, type,  mapping_list, FORWARD_ANCHOR, back0_nt - forw0_nt);*/
	 
	 __bwt_generate_anchor_list(last_k0, last_l0, forw0_nt, bwt_optarg, 
				    index, type,  mapping_list, BACKWARD_ANCHOR, back0_nt - forw0_nt, start, end);
	 //printf("FORWARD (-)\n");
	 __bwt_generate_anchor_list(last_ki0, last_li0, back0_nt, bwt_optarg, 
				    index, type,  mapping_list, FORWARD_ANCHOR, back0_nt - forw0_nt, start, end);
       
       }
     } else {
       int header_len;
       size_t num_mapping = array_list_size(mapping_list);
       //printf("BWT:Report Total mappings %i\n", num_mappings);
       for (int i = 0; i < num_mapping; i++) {
	 alignment = array_list_get(i, mapping_list);
	 header_len = strlen(read->id);
	 alignment->query_name = (char *) malloc(sizeof(char) * (header_len + 1));
	 get_to_first_blank(read->id, header_len, alignment->query_name);
	 bwt_cigar_cpy(alignment, read->quality);
	 //alignment->quality = strdup(&(batch->quality[batch->data_indices[i]]));                                                                                     
	 //************************* OPTIONAL FIELDS ***************************//
	 alignment = add_optional_fields(alignment, num_mappings);
       }
       array_list_set_flag(0, mapping_list);
       //printf("########## EXACT! #########\n");
     }

     free(r_list.list);
     free(code_seq);
     free(seq_strand);
     free(k0);
     free(l0);
     free(k1);
     free(l1);
     free(ki0);
     free(li0);
     free(ki1);
     free(li1);

     free(seq_dup);
     free(quality_clipping);
     array_list_free(tmp_mapping_list, NULL);

     return array_list_size(mapping_list);

}

//-----------------------------------------------------------------------------

alignment_t* add_optional_fields(alignment_t *alignment, size_t n_mappings) {
  char *p, *optional_fields;
  size_t optional_fields_length = 100;
  int distance;
  int AS = 254;
  int cigar_len;
  
  optional_fields = (char *)calloc(optional_fields_length, sizeof(char));
  p = optional_fields;
  
  sprintf(p, "ASi");
  p += 3;
  memcpy(p, &AS, sizeof(int));
  p += sizeof(int);
  
  sprintf(p, "NHi");
  p += 3;
  memcpy(p, &n_mappings, sizeof(int));
  p += sizeof(int);
  
  sprintf(p, "NMi");
  p += 3;
  cigar_len = strlen(alignment->cigar);
  
  if (alignment->cigar[cigar_len - 1] == '=') {
    distance = 0;
  } else {
    distance = 1;
  }
  
  memcpy(p, &distance, sizeof(int));
  p += sizeof(int);
  alignment->optional_fields_length = p - optional_fields;
  alignment->optional_fields = optional_fields;
  
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
// 11 -> 9 seeds (first pass), and 22 (second pass)
inline size_t seedingOK(char *code_seq, size_t seq_len, size_t num_seeds,
		      size_t seed_size, size_t min_seed_size,
		      bwt_optarg_t *bwt_optarg, bwt_index_t *index,
		      array_list_t *mapping_list) {

  size_t num_mappings = 0;
  size_t offset, offset_inc, offset_end = seq_len - min_seed_size;
  int seed_id = 0;

  if (seed_size * num_seeds > seq_len) {
    offset_inc = seed_size - (((seed_size * num_seeds) - seq_len) / num_seeds);
  } else {
    offset_inc = seq_len / num_seeds;
  }

  for (offset = 0; offset < offset_end; offset += offset_inc) {
    //printf("\nseed [%i - %i]\n", offset, offset + seed_size - 1);
    num_mappings += bwt_map_exact_seed(code_seq, seq_len, offset, offset + seed_size - 1,
				       bwt_optarg, index, mapping_list, seed_id);
    seed_id++;
  }

  return num_mappings;
}

//-----------------------------------------------------------------------------

inline size_t seeding(char *code_seq, size_t seq_len, size_t num_seeds,
		      size_t seed_size, size_t min_seed_size,
		      bwt_optarg_t *bwt_optarg, bwt_index_t *index,
		      array_list_t *mapping_list) {

  size_t n_seeds, total_mappings = 0, num_mappings = 0;
  //  size_t offset, offset_inc, offset_end = seq_len - min_seed_size;
  size_t start, end;
  size_t offset = 0, offset_inc, offset_end = seq_len - min_seed_size;
  int seed_id = 0;

  n_seeds = num_seeds;
  offset_inc = ceil(1.0f * seq_len / (num_seeds + 1));
  if (offset_inc <= 0) offset_inc = 1;
  /*
  if (seed_size * num_seeds > seq_len) {
    size_t max_seeds = seq_len - min_seed_size;
    if (num_seeds >= max_seeds) {
      n_seeds = max_seeds;
      offset_inc = 1;
    } else {
      n_seeds = num_seeds;
      offset_inc = seq_len / num_seeds;
    }
  } else {
    n_seeds = num_seeds;
    offset_inc = seq_len / num_seeds;
  }
  */

  start = 0;
  for (size_t i = 0; i < n_seeds; i++) {
    end = start + seed_size;
    if (end >= seq_len) end = seq_len;
    num_mappings = bwt_map_exact_seed(code_seq, seq_len, start, end - 1,
				      bwt_optarg, index, mapping_list, seed_id);
    seed_id++;
    total_mappings += num_mappings;
    //    LOG_DEBUG_F("\tseed %i\t[%i - %i], length read = %i, num. mappings = %i\n", 
    //		i + 1, start, end, seq_len, num_mappings);
    start += offset_inc;
    if (start > offset_end) {
      if (offset_inc == 1) break;
      start = offset_inc / 2;
    }
    /*
    if (start > offset_end) {
      offset++;
      start = offset;
    }
    */
  }

  //  LOG_DEBUG_F("\t\ttotal mappings = %i\n", total_mappings);

  return total_mappings;
}

//-----------------------------------------------------------------------------

size_t bwt_map_exact_seeds_seq_by_num(char *seq, size_t num_seeds, 
				      size_t seed_size, size_t min_seed_size,
				      bwt_optarg_t *bwt_optarg, bwt_index_t *index,
				      array_list_t *mapping_list) {
  size_t seq_len = strlen(seq);
  size_t num_mappings = 0;
  if (seed_size <= 10) { seed_size = 16; }
  char *code_seq = (char *) calloc(seq_len, sizeof(char));

  // be careful, now we use the table from the index
  //replaceBases(seq, code_seq, strlen(seq));
  bwt_encode_Bases(code_seq, seq, seq_len, &index->table);

  num_mappings = seeding(code_seq, seq_len, num_seeds, seed_size, min_seed_size,
			 bwt_optarg, index, mapping_list);
  //  printf("\tfirst, num_mappings = %d\n", num_mappings);
  /*
  if (num_mappings < 10) {
    num_mappings += seeding(code_seq, seq_len, max_num_seeds, seed_size - 4, min_seed_size - 4,
			    bwt_optarg, index, mapping_list);
    //    printf("\tsecond -4, num_mappings (include first) = %d\n", num_mappings);
    //  } else if (num_mappings >= 10000) {
  } else if (num_mappings >= bwt_optarg->filter_read_mappings) {
    array_list_clear(mapping_list, (void *) region_bwt_free);
    num_mappings = seeding(code_seq, seq_len, max_num_seeds, seed_size + 2, min_seed_size + 2,
			   bwt_optarg, index, mapping_list); 
    //    printf("\tthird +2, num_mappings = %d\n", num_mappings);
 }
  */
  free(code_seq);
  return num_mappings;
}

//-----------------------------------------------------------------------------
// CAL functions
//-----------------------------------------------------------------------------

int print_short_cal(void *item, void *dummy){
        short_cal_t *coordenate_p = (short_cal_t *)item;
	printf("[%lu-%lu]-> ",  coordenate_p->start, coordenate_p->end);

        return 0;
}

int cal_location_compare(short_cal_t *a, short_cal_t *b) {
	int res = 1;

	if(a->start == b->start){
		if(a->end == b->end){
			res = 0;
		}else if(a->end < b->end){
			res = -1;
		}else{
			res = 1;
		}
	}else if(a->start < b->start){
		res = -1;
	}

	return res;

}

short_cal_t* cal_location_dup(short_cal_t *a) {
	short_cal_t *exon_location_p = (short_cal_t *)malloc(sizeof(short_cal_t));
	exon_location_p->end = a->end;
	exon_location_p->start = a->start;

	return exon_location_p;
}

void print_se_region(short_cal_t *short_cal){
  printf("(%lu-%lu)->", short_cal->start, short_cal->end);
}


void print_se_region_cp(short_cal_t *short_cal, void* dummy){
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

void my_cp_list_append_linked_list(linked_list_t* list_p, region_t *region, size_t max_cal_distance, int max_seeds) {  
  unsigned char actualization = 0;
  short_cal_t *item, *item_aux, *new_item_p, *item_free;
  
  //int strand = region->strand;
  size_t start = region->start;
  size_t end = region->end;
  size_t seq_start = region->seq_start;
  size_t seq_end = region->seq_end;
  size_t seq_len = region->seq_len;
  linked_list_iterator_t* itr = linked_list_iterator_new(list_p);
  //printf("LINKED LIST INSERT %lu|%i-%i|%lu:\n", start, seq_start, seq_end, end);

  if (linked_list_size(list_p) <= 0) {
    new_item_p = short_cal_new(start, end, seq_start, seq_end, seq_len, max_seeds, region->id);
    linked_list_insert(new_item_p, list_p);
  } else {
    item = (short_cal_t *)linked_list_iterator_curr(itr);
    while (item != NULL) {
      if (start < item->start) {
	if (end + max_cal_distance < item->start) {
	  /*********************************************
	   *    Case 1: New item insert before item.   *
           *                                           *
           *        new item     item                  *
           *       |-------| |--------|                *
           ********************************************/

	  //printf("\t Insert now before %lu\n", item->start);
	  new_item_p = short_cal_new(start, end, seq_start, seq_end, seq_len, max_seeds, region->id);
	  linked_list_iterator_insert(new_item_p, itr);
	  linked_list_iterator_prev(itr);
	} else {
	  /********************************************
           *  Case 2: Actualization item start        *
           *           new item                       *
           *          |-------|   item                *
           *                   |--------|             *                            
           ********************************************/
	  //printf("\tFusion!\n");
	  item->start = start;
	  item->num_seeds++;
	  //item->seq_start = seq_start;
	  seed_region_select_linked_list(item->sr_list, item->sr_duplicate_list, 
					 seq_start, seq_end, start, end, region->id,
					 item->seeds_ids_array);
	  if (end > item->end) {
	    /**************************************************
             *  Case 3: Actualization item start and item end *
             *          new item                              *
             *         |------------|                         *
             *              item                              *    
             *           |--------|                           *                                    
             **************************************************/
	    //printf("\tFusion!\n");
	    item->end = end;
	    //item->seq_end = seq_end;
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
	  item->num_seeds++;
	  seed_region_select_linked_list(item->sr_list, item->sr_duplicate_list, seq_start,
					 seq_end, start, end, region->id, item->seeds_ids_array);
	  break;
	} else if (item->end + max_cal_distance >= start) {
	  /********************************************                                              
           *  Case 5: Actualization item end          *
           *            item                          *                                              
           *          |-------| new item              *                                            
           *                 |--------|               *                                              
           ********************************************/
	  //printf("\tFusion!\n");
	  item->end = end;
	  //item->seq_end = seq_end;
	  seed_region_select_linked_list(item->sr_list, item->sr_duplicate_list, seq_start,
					 seq_end, start, end, region->id, item->seeds_ids_array);
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
       * Case 5: Insert new item at the end of the list      * 
       *                 item    new item                    * 
       *              |-------| |--------|                   *    
       *******************************************************/
      new_item_p = short_cal_new(start, end, seq_start, seq_end, seq_len, max_seeds, region->id);
      //printf("\tInsert at END\n");
      linked_list_insert_last(new_item_p, list_p);
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
	if (item->end + max_cal_distance < item_aux->start) {
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
    }
  }//end first else
  //printf("End insert and actualization\n");
  linked_list_iterator_free(itr);
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
