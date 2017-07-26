#include "bwt.h"
#include <stdbool.h>



void __bwt_map_inexact_read_bs_new(fastq_read_t *read,
				 bwt_optarg_t *bwt_optarg,
				 bwt_index_t *index,
				 array_list_t *mapping_list,
				 int type,
				 size_t  start_read, size_t  end_read,
         bwt_stage_bs_workspace_t *workspace,
         interval_t *interval);

//-----------------------------------------------------------------------------

void bwt_init_replace_table(const char *str, int *table, int *rev_table) {
  if (str == NULL) {
    nA = 4;
    AA = 0; CC = 1; GG = 2; TT = 3;

    table['a'] = AA; table['A'] = AA;
    table['c'] = CC; table['C'] = CC;
    table['t'] = TT; table['T'] = TT;
    table['g'] = GG; table['G'] = GG;
    table['n'] = AA; table['N'] = AA;

    rev_table[AA] = 'A';
    rev_table[CC] = 'C';
    rev_table[GG] = 'G';
    rev_table[TT] = 'T';
  } else {
    nA = strlen(str);

    for (int i = 0; i < nA; i++) {
      rev_table[i] = toupper(str[i]);

      table[toupper(str[i])] = i;
      table[tolower(str[i])] = i;

      if      (toupper(str[i]) == 'A') AA = i;
      else if (toupper(str[i]) == 'C') CC = i;
      else if (toupper(str[i]) == 'G') GG = i;
      else if (toupper(str[i]) == 'T') TT = i;
    }
  }
}

//-----------------------------------------------------------------------------

char* bwt_encode_Bases(char* dest, char* src, unsigned int length, int *table) {
  size_t i;

  for (i=0; i<length; i++) {
    dest[i] = table[(int)src[i]];
  }

  return dest;
}

//-----------------------------------------------------------------------------

void bwt_generate_index_files_bs(char *ref_file, char *output_dir, 
				 unsigned int s_ratio, char *bases) {

  byte_vector X, B, Bi;
  vector C, C1;
  comp_vector S, Si, Scomp, Scompi;
  comp_vector R, Ri, Rcomp, Rcompi;
  comp_matrix O, Oi;

  exome ex;

  initReplaceTable_bs(bases);
  saveNucleotide(bases, output_dir, "Nucleotide");

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

void bwt_map_inexact_read_bs_new(fastq_read_t *read,
			       bwt_optarg_t *bwt_optarg,
			       bwt_index_t *index,
			       array_list_t *mapping_list,
			       int type, size_t  start_read, size_t end_read,
             bwt_stage_bs_workspace_t *workspace,
             interval_t *interval) {

  __bwt_map_inexact_read_bs_new(read,
				   bwt_optarg, 
				   index, 
				   mapping_list,
				   type,
				   start_read, end_read,
           workspace,
           interval);
}

//-----------------------------------------------------------------------------

//New inexact functions, used for apply bwt several times into the reads (when two errors take place, this function can be applied to the interior interval //Ricardo
void __bwt_map_inexact_read_bs_new(fastq_read_t *read,
				 bwt_optarg_t *bwt_optarg,
				 bwt_index_t *index,
				 array_list_t *mapping_list,
				 int type, size_t start_read, size_t end_read,
         bwt_stage_bs_workspace_t *workspace,
         interval_t *interval) {
	interval->start = start_read;
	interval->end = end_read;

	alignment_t *alignment;
  char *seq = read->sequence;
  size_t total_read_len = read->length;

  size_t start = start_read;
  size_t end = end_read;
  size_t len = end_read - start_read + 1;
  size_t len_calc = len;

  // Only a entire read can be mapped, if not the read will be included in an anchor
  bool mapp = true;

  if (len < total_read_len) {
    mapp = false;
  }

  if (len < 5) {
    return interval;
  }

  char *seq_dup, *seq_strand, *quality_clipping;

  if (workspace->map_inexact_read_code_seq == NULL) {
    workspace->map_inexact_read_code_seq = malloc(total_read_len * sizeof(char));
  }

  char *code_seq = workspace->map_inexact_read_code_seq;

  //encode entire read (400)
  bwt_encode_Bases(code_seq, seq, total_read_len, &index->table); 		

  // We work with the entire read, but the kl, ll, kil and lil arrays only fill len elements 
  // (end - start) <= 400
  // Calculate vectors k and l
  if (workspace->map_inexact_read_k1 == NULL) {
    workspace->map_inexact_read_k1 = malloc(len * sizeof(size_t));
  }

  if (workspace->map_inexact_read_l1 == NULL) {
    workspace->map_inexact_read_l1 = malloc(len * sizeof(size_t));
  }

  if (workspace->map_inexact_read_ki1 == NULL) {
    workspace->map_inexact_read_ki1 = malloc(len * sizeof(size_t));
  }

  if (workspace->map_inexact_read_li1 == NULL) {
    workspace->map_inexact_read_li1 = malloc(len * sizeof(size_t));
  }

  size_t *k1 = workspace->map_inexact_read_k1;
  size_t *l1 = workspace->map_inexact_read_l1;
  size_t *ki1 = workspace->map_inexact_read_ki1;
  size_t *li1 = workspace->map_inexact_read_li1;

  size_t last_k1, last_l1;
  size_t last_ki1, last_li1;

  int back1_nt, forw1_nt;

  //To store the last index that matches
  size_t min_start = start;

  //To store the last index that matches	
  size_t max_end = end; 
  size_t aux_index;

  aux_index = BWExactSearchVectorBackward(code_seq, start, end, 0, index->h_O.siz - 2,
            k1, l1, &index->h_C, &index->h_C1, &index->h_O,
            &last_k1, &last_l1, &back1_nt);

  if (aux_index < max_end) {
    max_end = aux_index;
  }

  aux_index = BWExactSearchVectorForward(code_seq, start, end, 0, index->h_Oi.siz - 2,
          ki1, li1, &index->h_C, &index->h_C1, &index->h_Oi,
          &last_ki1, &last_li1, &forw1_nt);

  if (aux_index > min_start) {
    min_start = aux_index;
  }

  // compare the vectors k and l to get mappings in the genome
  size_t num_mappings = 0;
  size_t idx, key, direction;
  char error;
  int pos;
  results_list r_list;
  result *r;
  char cigar[1024];
  size_t num_cigar_ops;
  size_t k_start, l_start;

  size_t start_mapping;
  size_t tot_alignments = 0;
  int filter_exceeded = 0;
  seq_strand = strdup(seq);
  error = MISMATCH;

  if (workspace->map_inexact_read_seq_dup == NULL) {
    workspace->map_inexact_read_seq_dup = malloc(sizeof(char) * (len + 1));
  }

  if (workspace->map_inexact_read_quality_clipping == NULL) {
    workspace->map_inexact_read_quality_clipping = calloc(len + 1, sizeof(char));
  } else {
    memset(workspace->map_inexact_read_quality_clipping, 0, len + 1);
  }

  seq_dup = workspace->map_inexact_read_seq_dup;
  quality_clipping = workspace->map_inexact_read_quality_clipping;

  new_results_list(&r_list, bwt_optarg->filter_read_mappings);

  if (workspace->map_inexact_read_tmp_mapping_list == NULL) {
    workspace->map_inexact_read_tmp_mapping_list = array_list_new(
        bwt_optarg->filter_read_mappings, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  } else {
    array_list_clear(workspace->map_inexact_read_tmp_mapping_list, alignment_free);
  }

  array_list_t *tmp_mapping_list = workspace->map_inexact_read_tmp_mapping_list;

  r_list.num_results = 0;
  r_list.read_index = 0;

  BWSearch1(code_seq, start, end, k1, l1, ki1, li1,
      &index->h_C, &index->h_C1, &index->h_O, &index->h_Oi, &r_list);

  for (size_t ii = 0; ii < r_list.num_results; ii++) {
    r = &r_list.list[ii];

    if(r == NULL) {
      printf("ERROR: Result list position null");exit(-1);
    }

    if(!r->num_mismatches) {
      error = 0;
    } else {
      error = r->err_kind[0];
    }

    pos = r->err_pos[0];
    direction = r->dir;
    len_calc = len;

    if (error == DELETION) {
      len_calc--;
    } else if (error == INSERTION) {
      len_calc++;
    }

    // generating cigar
    sprintf(quality_clipping, "%i", NONE_HARD_CLIPPING);

    if (error == 0) {
      sprintf(cigar, "%lu=", len);
      num_cigar_ops = 1;
      memcpy(seq_dup, seq_strand, len);
      seq_dup[len] = '\0';
    } else if (error == MISMATCH) {
      if (pos == 0) {
        //Positive strand
        sprintf(cigar, "1S%luM", len-1);
        start_mapping++;
        num_cigar_ops = 2;
      } else if (pos == len - 1) {
        //Positive strand
        sprintf(cigar, "%luM1S", len - 1);
        num_cigar_ops = 2;
      } else {
        sprintf(cigar, "%luM", len);
        num_cigar_ops = 1;
      }

      memcpy(seq_dup, seq_strand, len);
      seq_dup[len] = '\0';
    } else if (error == INSERTION) {
      if (pos == 0) {
        sprintf(cigar, "1M1D%luM", len - 1);
      } else if (pos == len - 1) {
        sprintf(cigar, "%luM1D1M", len - 1);
      } else {
        if(r->dir) {
          sprintf(cigar, "%iM1D%luM", pos, len - pos);
        } else {
          sprintf(cigar, "%iM1D%luM", pos + 1, len - pos - 1);
        }
      }

      num_cigar_ops = 3;
      memcpy(seq_dup, seq_strand, len);
      seq_dup[len] = '\0';
    } else if (error == DELETION) {
      if (pos == 0) {
        sprintf(cigar, "1I%luM", len -1);
        num_cigar_ops = 2;
      } else if (pos == len - 1) {
        sprintf(cigar, "%luM1I", len -1);
        num_cigar_ops = 2;
      } else {
        sprintf(cigar, "%dM1I%luM", pos, len - pos - 1);
        num_cigar_ops = 3;
      }

      memcpy(seq_dup, seq_strand , len );
      seq_dup[len] = '\0';
    } else {
      continue;
    }

    k_start = r->k;
    l_start = r->l;
    tot_alignments += (l_start - k_start);

    // check filter_read_mappings for bisulfite case
    if (tot_alignments >  bwt_optarg->filter_read_mappings) {
      filter_exceeded = 1;
      LOG_DEBUG_F("Filter exceeded: num. read mappings: %i (total: %i > filter %i)\n",
                  l_start - k_start, tot_alignments, bwt_optarg->filter_read_mappings);
      break;
    }

    for (size_t j = k_start; j <= l_start; j++) {
      if (index->S.ratio == 1) {
        key = (direction) ? index->Si.siz - index->Si.vector[j] - len_calc - 1 : index->S.vector[j];
      } else {
        key = (direction) ? index->Si.siz - getScompValue(j, &index->Si, &index->h_C,&index->h_Oi) - len_calc - 1
                          : getScompValue(j, &index->S, &index->h_C, &index->h_O);
      }

      idx = binsearch(index->karyotype.offset, index->karyotype.size, key);

      if (key + len_calc <= index->karyotype.offset[idx]) {
        start_mapping = index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]);

        // save all into one alignment structure and insert to the list
        alignment = alignment_new();

        alignment_init_single_end(NULL, strdup(seq_dup), strdup(quality_clipping), type,
                          idx - 1, start_mapping, strdup(cigar), num_cigar_ops, 254, 1, 
                          (num_mappings > 0), 0, NULL, alignment);
        array_list_insert((void*) alignment, tmp_mapping_list);
      }
    }//end for k and l
  }//end for ii (550)

  if (filter_exceeded) {
    array_list_clear(tmp_mapping_list, (void *)alignment_free);

    if (!mapping_list->size) {
      array_list_set_flag(2, mapping_list);
    }

    goto exit;
  }

  //Search for equal BWT mappings and set the mappings that will be delete
  int n_mappings = array_list_size(tmp_mapping_list);
  alignment_t *alig_1, *alig_2;

  unsigned int *delete_mark = workspace->map_inexact_read_delete_mark;
  memset(workspace->map_inexact_read_delete_mark, 0, BWT_MAX_INEXACT_READ_DELETE_MARK);

  const int max_distance = 10;

  for (int a1 = n_mappings - 1; a1 >= 1; a1--) {
    if (!delete_mark[a1]) {
      alig_1 = array_list_get(a1, tmp_mapping_list);

      for (int a2 = a1 - 1; a2 >= 0; a2--) {
        alig_2 = array_list_get(a2, tmp_mapping_list);
        size_t dist = abs(alig_1->position - alig_2->position);

        if (alig_1->chromosome == alig_2->chromosome &&
            dist < max_distance                      &&
            !delete_mark[a2]) {
          //Same chromosome && same position
          if (alig_1->num_cigar_operations < alig_2->num_cigar_operations) {
            delete_mark[a2] = 1;
          } else {
            delete_mark[a1] = 1;
          }
        }
      }
    }
  }

  if (mapp && array_list_get_flag(mapping_list) == 1 && n_mappings) {
    array_list_clear(mapping_list, (void *)bwt_anchor_free);
  }

 
  if (mapp) {
    //Delete all set mappings
    int primary_delete = 0;
    int header_len;

    for (int m = n_mappings - 1; m >= 0; m--) {
      alig_1 = array_list_remove_at(m, tmp_mapping_list);

      if (delete_mark[m]) {
        if (!is_secondary_alignment(alig_1)) { 
          primary_delete = 1; 
        }

        alignment_free(alig_1);
      } else {
        set_secondary_alignment(num_mappings > 0, alig_1);
        num_mappings++;
        
        header_len = strlen(read->id);
        alig_1->query_name = (char *) malloc(sizeof(char) * (header_len + 1));
        get_to_first_blank(read->id, header_len, alig_1->query_name);

        bwt_cigar_cpy(alig_1, read->quality);
        //************************* OPTIONAL FIELDS ***************************//
        alig_1 = add_optional_fields(alig_1, n_mappings);
        array_list_insert((void*) alig_1, mapping_list);
      }
    }

    if (primary_delete) {
      alig_1 = array_list_get(0, mapping_list);
      set_secondary_alignment(0, alig_1);
    }
  }

  if (n_mappings == 0) {
    if (array_list_get_flag(mapping_list) == 1) {

      //1 == (+)
      int new_type = !type;

      __bwt_generate_anchor_list(last_k1, last_l1, back1_nt, bwt_optarg,
                                index, new_type, mapping_list, BACKWARD_ANCHOR, 
                                0, max_end+1, end_read);
      __bwt_generate_anchor_list(last_ki1, last_li1, forw1_nt, bwt_optarg,
                                index, new_type,  mapping_list, FORWARD_ANCHOR, 
                                0, start_read, min_start-1);
    }
  } else {
    if (mapp) {
      array_list_set_flag(0, mapping_list);
    } else {
      //it is mapped completely, but is a seed, no the entire read, so it must be added like an anchor
      //RICARDO,
      array_list_set_flag(1,mapping_list);

      int new_type = !type; //1 == (+)

      //Added like one seed if not error
      if ((start_read >= max_end) && (end_read <= min_start)) {	
        //Solo es necesario añadirla en una lista
        if (start_read > (total_read_len - end_read)) {
          __bwt_generate_anchor_list(last_k1, last_l1, back1_nt, bwt_optarg,
                                    index, new_type, mapping_list, BACKWARD_ANCHOR, 
                                    0, start_read, end_read);
        } else {
          __bwt_generate_anchor_list(last_ki1, last_li1, forw1_nt, bwt_optarg,
                                     index, new_type,  mapping_list, FORWARD_ANCHOR, 
                                     0, start_read, end_read);
        }
      } else {
        //Mapped with one error, we add two seeds without error
        __bwt_generate_anchor_list(last_k1, last_l1, back1_nt, bwt_optarg,
                                   index, new_type, mapping_list, BACKWARD_ANCHOR, 
                                   0, max_end+1, end_read);

        __bwt_generate_anchor_list(last_ki1, last_li1, forw1_nt, bwt_optarg,
                                   index, new_type,  mapping_list, FORWARD_ANCHOR, 
                                   0, start_read, min_start-1);
      }
    }
  }

exit:
  free(r_list.list);
  free(seq_strand);

  interval->start = min_start;
  interval->end = max_end;
}


//-----------------------------------------------------------------------------

char * readNucleotide(const char *directory, const char *name) {
  FILE *fp;

  char path[500];
  char tmp[5];
  //char *tmp = malloc(5 * sizeof(char));
  //char *tmp;

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".txt");

  fp  = fopen(path,  "r");
  checkFileOpen(fp, path);

  fgets(tmp, 4, fp);

  fclose(fp);

  return strdup(tmp);
  //return tmp;
}

//-----------------------------------------------------------------------------

void saveNucleotide(char *nucleotide, const char *directory, const char *name) {
  FILE *fp;

  char path[500];
  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".txt");

  fp  = fopen(path,  "w");
  checkFileOpen(fp, path);

  fputs(nucleotide, fp);

  fclose(fp);
}

//-----------------------------------------------------------------------------
