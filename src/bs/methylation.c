#include "methylation.h"


//------------------------------------------------------------------------------------

void replace(char * refs, int len, int type) {
  char c1, c2;

  switch (type) {
    case ACGT:
      // Case with no transformation required
      return;

    case AGT:
      c1 = 'C';
      c2 = 'T';
      break;

    case ACT:
      c1 = 'G';
      c2 = 'A';
      break;

    case AT:
      // Apply both transformations, and end the execution
      replace(refs, len, 1);
      replace(refs, len, 2);
      return;

    default:
      // If the type value is not recognised, nothing is done
      return;
  }

  // Transforms the refs sequence, replacing the caracter c1 with the caracter c2
  for (int j = 0; j < len; j++) {
    if (refs[j] == c1) {
      refs[j] = c2;
    }
  }

  return;
}

//------------------------------------------------------------------------------------

char complement (char c) {
  // Return the complementary base
  switch (c) {
    case 'A':
      return 'T';
    case 'T':
      return 'A';
    case 'C':
      return 'G';
    case 'G':
      return 'C';
    default:
      return c;
  }
}

//------------------------------------------------------------------------------------

void rev_comp(char *orig, char *dest, int len) {
  // Put the reverse complementary sequence of orig in dest
  for (int i = 0; i < len; i++) {
    dest[len - i - 1] = complement(orig[i]);
  }

//RICARDO -- must be dest[len] = '\0', but strndup must do it
  //dest[len - 1] = '\0';
}

//------------------------------------------------------------------------------------

void comp(char *seq, int len) {
  // obtain the complementary sequence of seq
  for (int i = 0; i < len - 1; i++) {
    seq[i] = complement(seq[i]);
  }
}


//------------------------------------------------------------------------------------

void cpy_transform_array_bs(array_list_t *src, array_list_t *dest_ct, array_list_t *dest_ct_rev, array_list_t *dest_ga, array_list_t *dest_ga_rev) {
  size_t num_reads = array_list_size(src);
  fastq_read_t *fq_read_src;
  fastq_read_t *fq_read_dest;
  fastq_read_t *fq_read_tmp;

  // Read element by element in the array, transform each one, and put in the new arrays
  for (size_t i = 0; i < num_reads; i++) {
    fq_read_src  = (fastq_read_t *) array_list_get(i, src);

    fq_read_dest = fastq_read_dup(fq_read_src);
    replace(fq_read_dest->sequence, fq_read_dest->length, AGT);
    array_list_insert(fq_read_dest, dest_ct);

    fq_read_tmp = fastq_read_dup(fq_read_src);
    rev_comp(fq_read_dest->sequence, fq_read_tmp->sequence, fq_read_dest->length);
    array_list_insert(fq_read_tmp, dest_ct_rev);

    fq_read_dest = fastq_read_dup(fq_read_src);
    replace(fq_read_dest->sequence, fq_read_dest->length, ACT);
    array_list_insert(fq_read_dest, dest_ga);

    fq_read_tmp = fastq_read_dup(fq_read_src);
    rev_comp(fq_read_dest->sequence, fq_read_tmp->sequence, fq_read_dest->length);
    array_list_insert(fq_read_tmp, dest_ga_rev);
  }
}

//------------------------------------------------------------------------------------

void revert_mappings_seqs(array_list_t **src1, array_list_t **src2, array_list_t *orig) {
  size_t num_mappings;
  alignment_t *align_tmp;
  fastq_read_t *fastq_orig;
  size_t num_reads = array_list_size(orig);

  // Go over all the sequences
  for (size_t i = 0; i < num_reads; i++) {
    fastq_orig  = (fastq_read_t *) array_list_get(i, orig);
 
    // Go over all the alignments in list 1
    num_mappings = array_list_size(src1[i]);
    
    for (size_t j = 0; j < num_mappings; j++) {
      align_tmp  = (alignment_t *) array_list_get(j, src1[i]);

      // Free existing memory and copy the original
      if (align_tmp->sequence != NULL) {
        free(align_tmp->sequence);
      }

      if (align_tmp->num_cigar_operations > 1)	//Check for Hard clipping in order to remove bases and Soft clipping in order to restore quality
           {
         	  //Find hard clipping

         	  int i=0;
         	  char aux[10000], aux2[10000];

         	  char operation;
         	  int bases;
         	  int inicio=0;
         	  int fin = strlen(fastq_orig->sequence);
         	  bool restore_quality = false;

         	  strcpy(aux, align_tmp->cigar);
         	  for (i=0; i<align_tmp->num_cigar_operations;i++)
         	  {
         		  sscanf(aux,"%d%c%s", &bases, &operation, aux2);
         		  if (operation == 'H')	//Hard Clipping
         		  {
         			  if (i==0) //Hard clipping in the start of the read
         			  {
         				inicio = bases;

         			  }
         			  else	//Hard clipping in the end of the read
         			  {
         				 fin = fin - bases;
         			  }
         		  }
         		  else  if (operation == 'S')	//Soft Clipping
         		  {
         			  restore_quality = true;
         		  }

         		  strcpy(aux,aux2);


         	  }

         	  align_tmp->sequence = strndup(fastq_orig->sequence+inicio,fin-inicio);
         	  if (restore_quality)
         	  {
         		 if (align_tmp->quality != NULL) {
         		         free(align_tmp->quality);
         		       }
         		 align_tmp->quality = strndup(fastq_orig->quality+inicio,fin-inicio);
         	  }
           }

      else{
    	        align_tmp->sequence = strdup(fastq_orig->sequence);
      }
    }

    // Go over all the alignments in list 2
    num_mappings = array_list_size(src2[i]);

    for (size_t j = 0; j < num_mappings; j++) {
      align_tmp  = (alignment_t *) array_list_get(j, src2[i]);

      // Free existing memory and copy the original
      if (align_tmp->sequence != NULL) {
        free(align_tmp->sequence);
      }

      if (align_tmp->num_cigar_operations > 1)	//Check for Hard clipping in order to remove bases and Soft clipping in order to restore quality
      {
    	  //Find hard clipping
    	  int i=0;
    	  char aux[10000], aux2[10000];

    	  char operation;
    	  int bases;
    	  int inicio=0;
    	  int fin = strlen(fastq_orig->sequence);
    	  bool restore_quality = false;

    	  strcpy(aux, align_tmp->cigar);
    	  for (i=0; i<align_tmp->num_cigar_operations;i++)
    	  {
    		  sscanf(aux,"%d%c%s", &bases, &operation, aux2);
    		  if (operation == 'H')	//Hard Clipping
    		  {
    			  if (i==0) //Hard clipping in the start of the read
    			  {
    				inicio = bases;

    			  }
    			  else	//Hard clipping in the end of the read
    			  {
    				 fin = fin - bases;
    			  }
    		  }
    		  else  if (operation == 'S')	//Soft Clipping
    		  {
    			  restore_quality = true;
    		  }

    		  strcpy(aux,aux2);


    	  }

    	  align_tmp->sequence = strndup(fastq_orig->sequence+inicio,fin-inicio);
    	  if (restore_quality)
    	  {
    	  	 if (align_tmp->quality != NULL) {
    	         free(align_tmp->quality);
    	     }
    	  	 align_tmp->quality =  strndup(fastq_orig->quality+inicio,fin-inicio);
    	  }
      }
      else
      {
    	  align_tmp->sequence = strdup(fastq_orig->sequence);
      }
    }
  }
}

//------------------------------------------------------------------------------------

char *obtain_seq(alignment_t *alig, fastq_read_t * orig) {
  int cont, pos = 0, pos_read = 0, num, operations;
  int offset = 0, current_offset = 0; 
  int seq_len = orig->length;
  char car;
  
  char *read = orig->sequence;
  char *cigar = strdup(alig->cigar);
  char *seq = (char *)calloc(seq_len, sizeof(char));
  
  for (operations = 0; operations < alig->num_cigar_operations; operations++) {
    sscanf(&cigar[offset], "%i%c%n", &num, &car, &current_offset);
    offset += current_offset;

    if (car == 'M' || car == '=' || car == 'X') {
      for (cont = 0; cont < num; cont++, pos++, pos_read++) {
        if (pos_read < seq_len && pos < seq_len) {
	        seq[pos] = read[pos_read];
        }
      }
    } else {
      if (car == 'D' || car == 'N') {
	      pos_read += num - 1;
      } else {
        if (car == 'I') {
          for (cont = 0; cont < num; cont++, pos++) {
            if (pos < seq_len) {
              seq[pos] = '-';
            }
          }
        }
        else{
        	if (car == 'H' || car == 'S') {
              	for (cont = 0; cont < num; cont++, pos++, pos_read++) {
        		if (pos < seq_len) {
        			seq[pos] = '-';
        		}
        	}
        }
        }
      }

    }
  }

  //Ricardo, sobreescribe elemento del read, con inserciones se pierde el final del genoma original, con delecciones tamaño más pequeño
  if (pos < seq_len)
	  seq[pos] = '\0';
	  //  seq[seq_len - 1] = '\0';

  free(cigar);
  return seq;
}

//------------------------------------------------------------------------------------

void metil_file_init(metil_file_t *metil_file, char *dir, genome_t *genome) {
  char *name_tmp = malloc(128 * sizeof(char));
  int file_error;

  sprintf(name_tmp, "%s/CpG_context.txt", dir);
  metil_file->filenameCpG = strdup(name_tmp);
  
  sprintf(name_tmp, "%s/CHG_context.txt", dir);
  metil_file->filenameCHG = strdup(name_tmp);

  sprintf(name_tmp, "%s/CHH_context.txt", dir);
  metil_file->filenameCHH = strdup(name_tmp);

  sprintf(name_tmp, "%s/MUT_context.txt", dir);
  metil_file->filenameMUT = strdup(name_tmp);

  sprintf(name_tmp, "%s/Statistics.txt", dir);
  metil_file->filenameSTAT = strdup(name_tmp);


  metil_file->genome = genome;

  metil_file->CpG  = fopen(metil_file->filenameCpG,  "w");
  metil_file->CHG  = fopen(metil_file->filenameCHG,  "w");
  metil_file->CHH  = fopen(metil_file->filenameCHH,  "w");
  metil_file->MUT  = fopen(metil_file->filenameMUT,  "w");
  metil_file->STAT = fopen(metil_file->filenameSTAT, "w");

  metil_file->skip_cpg = 0;
  metil_file->skip_chh = 0;
  metil_file->skip_chg = 0;
  metil_file->skip_mut = 0;

  FILE *a;

  a = metil_file->CpG;
  file_error = fprintf(a, "File for Cytosines in CpG context\n\n");

  if (file_error < 0) {
    printf("Error al escribir\n");
    exit(-1);
  }

  a = metil_file->CHG;
  file_error = fprintf(a, "File for Cytosines in CHG context\n\n");

  if (file_error < 0) {
    printf("Error al escribir\n");
    exit(-1);
  }

  a = metil_file->CHH;
  file_error = fprintf(a, "File for Cytosines in CHH context\n\n");

  if (file_error < 0) {
    printf("Error al escribir\n");
    exit(-1);
  }

  a = metil_file->MUT;
  file_error = fprintf(a, "File for Cytosines mutated\n\n");

  if (file_error < 0) {
    printf("Error al escribir\n");
    exit(-1);
  }

  a = metil_file->STAT;
  file_error = fprintf(a, "File for Methylation Statistics\n\n");

  if (file_error < 0) {
    printf("Error al escribir\n");
    exit(-1);
  }

  free(name_tmp);

  metil_file->CpG_methyl   = 0;
  metil_file->CpG_unmethyl = 0;
  metil_file->CHG_methyl   = 0;
  metil_file->CHG_unmethyl = 0;
  metil_file->CHH_methyl   = 0;
  metil_file->CHH_unmethyl = 0;
  metil_file->MUT_methyl   = 0;
  metil_file->num_bases    = 0;

  metil_file->methyl_reads = calloc(genome->num_chromosomes, sizeof(uint32_t));
}

//------------------------------------------------------------------------------------

void metil_file_free(metil_file_t *metil_file) {
  free(metil_file->filenameCpG);
  free(metil_file->filenameCHG);
  free(metil_file->filenameCHH);
  free(metil_file->filenameMUT);
  free(metil_file->filenameSTAT);
  free(metil_file->methyl_reads);

  if (metil_file->CpG  != NULL) fclose(metil_file->CpG);
  if (metil_file->CHG  != NULL) fclose(metil_file->CHG);
  if (metil_file->CHH  != NULL) fclose(metil_file->CHH);
  if (metil_file->MUT  != NULL) fclose(metil_file->MUT);
  if (metil_file->STAT != NULL) fclose(metil_file->STAT);

  free(metil_file);
}

//====================================================================================

void remove_duplicates(size_t reads, array_list_t **list, array_list_t **list2) {
  size_t num_items, num_items2;
  alignment_t *alig, *alig2;

  for (size_t i = 0; i < reads; i++) {
    num_items = array_list_size(list[i]);

    for (size_t j = 0; j < num_items; j++) {
      alig = (alignment_t *) array_list_get(j, list[i]);

      if (alig != NULL && alig->is_seq_mapped) {
        num_items2 = array_list_size(list2[i]);

        for (size_t k = 0; k < num_items2; k++) {
          alig2 = (alignment_t *) array_list_get(k, list2[i]);

          if (alig->position      == alig2->position
              && alig->chromosome == alig2->chromosome
              && alig->seq_strand == alig2->seq_strand
              && alig->num_cigar_operations == alig2->num_cigar_operations) {
            alig2 = (alignment_t *) array_list_remove_at(k, list2[i]);
            alignment_free(alig2);
            k--;
            num_items2 = array_list_size(list2[i]);
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------------

void clean_methylation_stage_workspace(void *workspace) {
  if (workspace) {
    methylation_stage_workspace_t* wf = workspace;

    if (wf->add_status_seq) {
      free(wf->add_status_seq);
    }

    if (wf->add_status_gen) {
      free(wf->add_status_gen);
    }

    if (wf->add_status_seq_dup) {
      free(wf->add_status_seq_dup);
    }

    free(wf);
  }
}

//------------------------------------------------------------------------------------

int methylation_status_report(sw_server_input_t* input, batch_t *batch, methylation_stage_workspace_t *workspace) {

  LOG_DEBUG("========= METHYLATION STATUS REPORT START =========\n");

  mapping_batch_t *mapping_batch = (mapping_batch_t *) batch->mapping_batch;
  array_list_t **mapping_lists;
  size_t num_items;
  size_t num_reads = array_list_size(mapping_batch->fq_batch);
  genome_t *genome = input->genome_p;

  struct timeval methylation_time_start, methylation_time_end;
  float methylation_time = 0.0f;

  if (time_on) {
    start_timer(methylation_time_start);
  }

  // Inicializar listas para guardar datos de c's metiladas/no metiladas
  //bs_context_t *bs_context = bs_context_new(num_reads, genome->num_chromosomes);
  bs_context_t *bs_context = bs_context_new(10000, genome->num_chromosomes);

  mapping_batch->bs_context = bs_context;

  remove_duplicates(num_reads, mapping_batch->mapping_lists, mapping_batch->mapping_lists2);
  
  for (int k = 0; k < 2; k++) {
    mapping_lists = (k == 0) ? mapping_batch->mapping_lists : mapping_batch->mapping_lists2;
    
    for (size_t i = 0; i < num_reads; i++) {
      num_items = array_list_size(mapping_lists[i]);
      
      // Mapped or not mapped ?
      if (num_items != 0) {
	      add_metilation_status(mapping_lists[i], bs_context, genome, mapping_batch->fq_batch, i, k, workspace);
      }
    }
  }

  if (time_on) {
    stop_timer(methylation_time_start, methylation_time_end, methylation_time);
    timing_add(methylation_time, METHYLATION_REP_TIME, timing);
  }

  return CONSUMER_STAGE;
}

//------------------------------------------------------------------------------------

void add_metilation_status(array_list_t *array_list, bs_context_t *bs_context, 
    genome_t * genome, array_list_t * orig_seq, size_t index, int conversion,
    methylation_stage_workspace_t *workspace) {
  size_t num_items = array_list_size(array_list);
  size_t len, end, start;

  // Number of methylated entries, used for the ZM tag
  size_t num_methyl = 0;

  alignment_t *alig;
  fastq_read_t *orig;

  char *seq, *gen;
  
  int new_strand;
  int write_file = 1;

  orig = (fastq_read_t *) array_list_get(index, orig_seq);

  for (size_t j = 0; j < num_items; j++) {
    alig = (alignment_t *) array_list_get(j, array_list);

    if (alig != NULL && alig->is_seq_mapped) {
      seq = obtain_seq(alig, orig);

      if (alig->seq_strand == 1) {
        if (workspace->add_status_seq_dup == NULL) {
          workspace->add_status_seq_dup = strdup(seq);
        } else {
        	if (strlen(workspace->add_status_seq_dup)== strlen(seq))
        	{
        		strcpy(workspace->add_status_seq_dup, seq);
        	}
        	else
        	{
        		free(workspace->add_status_seq_dup);
        		workspace->add_status_seq_dup = strdup(seq);
        	}
        }

        char *seq_dup = workspace->add_status_seq_dup;
        rev_comp(seq_dup, seq, orig->length);
      }

      // Increase the counter number of bases
      len = orig->length;

      if (workspace->add_status_gen == NULL) {
        workspace->add_status_gen = calloc(len + 6, sizeof(char));
      }
      else if (strlen(workspace->add_status_gen) != (len + 6)) {
    	    free(workspace->add_status_gen);
            workspace->add_status_gen = calloc(len + 6, sizeof(char));
          }

      gen = workspace->add_status_gen;

      start = alig->position + 1;
      end = start + len + 4;

      if (end >= genome->chr_size[alig->chromosome]) {
        end = genome->chr_size[alig->chromosome] - 1;
      }

      genome_read_sequence_by_chr_index(gen, alig->seq_strand, alig->chromosome, &start, &end, genome);

      // Initialize the auxiliary tag data for the alignment
      bam_tag_t *tag_zm = bam_tag_init(ZM_TAG_NAME, BAM_TAG_TYPE_INT, 0, 0);
      bam_tag_t *tag_xm = bam_tag_init(XM_TAG_NAME, BAM_TAG_TYPE_STRING, 0, len);
      bam_tag_t *tag_xg = bam_tag_init(XG_TAG_NAME, BAM_TAG_TYPE_STRING, 0, 2);
      bam_tag_t *tag_xr = bam_tag_init(XR_TAG_NAME, BAM_TAG_TYPE_STRING, 0, 2);

      // Set the alignment conversion data
      if (conversion) {
        bam_tag_str_insert(tag_xg, XG_CONVERSION_CT, 2);
        bam_tag_str_insert(tag_xr, XR_CONVERSION_CT, 2);
      } else {
        bam_tag_str_insert(tag_xg, XG_CONVERSION_GA, 2);
        bam_tag_str_insert(tag_xr, XR_CONVERSION_GA, 2);
      }

      for (size_t i = 0; i < len; i++) {
        if ((conversion == 1 && alig->seq_strand == 0) || (conversion == 0 && alig->seq_strand == 1)) {
          // Methylated/unmethylated cytosines are located in the same strand as the alignment
          if (gen[i] == 'C') {
            if (gen[i + 1] == 'G') {
              if (seq[i] == 'C') {
                if (write_file == 1) {
                  postprocess_bs(alig->query_name, '+', alig->chromosome, start + i, XM_METHYLATED_CPG, alig->seq_strand, 0, bs_context->context_CpG);

                  // Update the alignment methylation tag
                  bam_tag_str_append(tag_xm, XM_METHYLATED_CPG);
                  ++num_methyl;
                }

                bs_context->CpG_methyl++;
              } else if (seq[i] == 'T') {
                if (write_file == 1) {
                  postprocess_bs(alig->query_name, '-', alig->chromosome, start + i, XM_UNMETHYLATED_CPG, alig->seq_strand, 0,bs_context->context_CpG);

                  // Update the alignment methylation tag
                  bam_tag_str_append(tag_xm, XM_UNMETHYLATED_CPG);
                }

                bs_context->CpG_unmethyl++;
              } else {
                if (write_file == 1) {
                  postprocess_bs(alig->query_name, '.', alig->chromosome, start + i, XM_METHYLATED_MUT, alig->seq_strand, 3,bs_context->context_MUT);

                  // Update the alignment methylation tag
                  bam_tag_str_append(tag_xm, XM_METHYLATED_MUT);
                  ++num_methyl;
                }
                
                bs_context->MUT_methyl++;
              }
            } else {
              if (gen[i + 2] == 'G') {
                if (seq[i] == 'C') {
                  if (write_file == 1) {
                    postprocess_bs(alig->query_name, '+', alig->chromosome, start + i, XM_METHYLATED_CHG, alig->seq_strand, 1,bs_context->context_CHG);

                    // Update the alignment methylation tag
                    bam_tag_str_append(tag_xm, XM_METHYLATED_CHG);
                    ++num_methyl;
                  }

                  bs_context->CHG_methyl++;
                } else if (seq[i] == 'T') {
                  if (write_file == 1) {
                    postprocess_bs(alig->query_name, '-', alig->chromosome, start + i, XM_UNMETHYLATED_CHG, alig->seq_strand, 1,bs_context->context_CHG);

                    // Update the alignment methylation tag
                    bam_tag_str_append(tag_xm, XM_UNMETHYLATED_CHG);
                  }

                  bs_context->CHG_unmethyl++;
                } else {
                  if (write_file == 1) {
                    postprocess_bs(alig->query_name, '.', alig->chromosome, start + i, XM_METHYLATED_MUT, alig->seq_strand, 3,bs_context->context_MUT);

                    // Update the alignment methylation tag
                    bam_tag_str_append(tag_xm, XM_METHYLATED_MUT);
                    ++num_methyl;
                  }

                  bs_context->MUT_methyl++;
                }
              } else {
                if (seq[i] == 'C') {
                  if (write_file == 1) {
                    postprocess_bs(alig->query_name, '+', alig->chromosome, start + i, XM_METHYLATED_CHH, alig->seq_strand, 2,bs_context->context_CHH);

                    // Update the alignment methylation tag
                    bam_tag_str_append(tag_xm, XM_METHYLATED_CHH);
                    ++num_methyl;
                  }

                  bs_context->CHH_methyl++;
                } else if (seq[i] == 'T') {
                  if (write_file == 1) {
                    postprocess_bs(alig->query_name, '-', alig->chromosome, start + i, XM_UNMETHYLATED_CHH, alig->seq_strand, 2,bs_context->context_CHH);

                    // Update the alignment methylation tag
                    bam_tag_str_append(tag_xm, XM_UNMETHYLATED_CHH);
                  }

                  bs_context->CHH_unmethyl++;
                } else {
                  if (write_file == 1) {
                    postprocess_bs(alig->query_name, '.', alig->chromosome, start + i, XM_METHYLATED_MUT, alig->seq_strand, 3,bs_context->context_MUT);

                    // Update the alignment methylation tag
                    bam_tag_str_append(tag_xm, XM_METHYLATED_MUT);
                    ++num_methyl;
                  }
                      
                  bs_context->MUT_methyl++;
                }
              }
            }
          } else {
            // Skip one methylation entry if the current
            // nucleotide is not a cytosine
            bam_tag_str_append(tag_xm, XM_NON_RELEVANT);
          }
        } else {
          // Methylated/unmethylated cytosines are located in the other strand
          if (alig->seq_strand == 0) {
            new_strand = 1;
          } else {
            new_strand = 0;
          }
          
          if (gen[i+2] == 'G') {
            if (gen[i + 1] == 'C') {
              if (seq[i] == 'G') {
                if (write_file == 1) {
                  postprocess_bs(alig->query_name, '+', alig->chromosome, start + i, XM_METHYLATED_CPG, new_strand, 0,bs_context->context_CpG);

                  // Update the alignment methylation tag
                  bam_tag_str_append(tag_xm, XM_METHYLATED_CPG);
                  ++num_methyl;
                }

                bs_context->CpG_methyl++;
              } else if (seq[i] == 'A') {
                if (write_file == 1) {
                  postprocess_bs(alig->query_name, '-', alig->chromosome, start + i, XM_UNMETHYLATED_CPG, new_strand, 0,bs_context->context_CpG);

                  // Update the alignment methylation tag
                  bam_tag_str_append(tag_xm, XM_UNMETHYLATED_CPG);
                }

                bs_context->CpG_unmethyl++;
              } else {
                if (write_file == 1) {
                  postprocess_bs(alig->query_name, '.', alig->chromosome, start + i, XM_METHYLATED_MUT, alig->seq_strand, 3,bs_context->context_MUT);

                  // Update the alignment methylation tag
                  bam_tag_str_append(tag_xm, XM_METHYLATED_MUT);
                  ++num_methyl;
                }
                
                bs_context->MUT_methyl++;
              }
            } else {
              if (gen[i] == 'C') {
                if (seq[i] == 'G') {
                  if (write_file == 1) {
                    postprocess_bs(alig->query_name, '+', alig->chromosome, start + i, XM_METHYLATED_CHG, new_strand, 1,bs_context->context_CHG);

                    // Update the alignment methylation tag
                    bam_tag_str_append(tag_xm, XM_METHYLATED_CHG);
                    ++num_methyl;
                  }
                      
                  bs_context->CHG_methyl++;
                } else if (seq[i] == 'A') {
                  if (write_file == 1) {
                    postprocess_bs(alig->query_name, '-', alig->chromosome, start + i, XM_UNMETHYLATED_CHG, new_strand, 1,bs_context->context_CHG);

                    // Update the alignment methylation tag
                    bam_tag_str_append(tag_xm, XM_UNMETHYLATED_CHG);
                  }

                  bs_context->CHG_unmethyl++;
                } else {
                  if (write_file == 1) {
                    postprocess_bs(alig->query_name, '.', alig->chromosome, start + i, XM_METHYLATED_MUT, alig->seq_strand, 3,bs_context->context_MUT);

                    // Update the alignment methylation tag
                    bam_tag_str_append(tag_xm, XM_METHYLATED_MUT);
                    ++num_methyl;
                  }

                  bs_context->MUT_methyl++;
                }
              } else {
                if (seq[i] == 'G') {
                  if (write_file == 1) {
                    postprocess_bs(alig->query_name, '+', alig->chromosome, start + i, XM_METHYLATED_CHH, new_strand, 2,bs_context->context_CHH);

                    // Update the alignment methylation tag
                    bam_tag_str_append(tag_xm, XM_METHYLATED_CHH);
                    ++num_methyl;
                  }

                  bs_context->CHH_methyl++;
                } else if (seq[i] == 'A') {
                  if (write_file == 1) {
                    postprocess_bs(alig->query_name, '-', alig->chromosome, start + i, XM_UNMETHYLATED_CHH, new_strand, 2,bs_context->context_CHH);

                    // Update the alignment methylation tag
                    bam_tag_str_append(tag_xm, XM_UNMETHYLATED_CHH);
                  }

                  bs_context->CHH_unmethyl++;
                } else {
                  if (write_file == 1) {
                    postprocess_bs(alig->query_name, '.', alig->chromosome, start + i, XM_METHYLATED_MUT, alig->seq_strand, 3,bs_context->context_MUT);

                    // Update the alignment methylation tag
                    bam_tag_str_append(tag_xm, XM_METHYLATED_MUT);
                    ++num_methyl;
                  }

                  bs_context->MUT_methyl++;
                }
              }
            }
          } else {
            // Skip one methylation entry if the current
            // nucleotide is not a cytosine
            bam_tag_str_append(tag_xm, XM_NON_RELEVANT);
          }
        }
      }

      // Set the number of methylated C's to the ZM tag
      bam_tag_set_scalar(tag_zm, &num_methyl);

      if (num_methyl != 0) {
        bs_context->methyl_reads[alig->chromosome]++;
      }
      
      num_methyl = 0;

      // Add the methylation context tags to the current alignment
      array_list_insert(tag_zm, alig->optional_tags); 
      array_list_insert(tag_xg, alig->optional_tags);
      array_list_insert(tag_xr, alig->optional_tags);
      array_list_insert(tag_xm, alig->optional_tags);
    }
  }
}

//------------------------------------------------------------------------------------

void write_bs_context(metil_file_t *metil_file, bs_context_t *bs_context, size_t num_chromosomes) {
  array_list_t *context_CpG = bs_context->context_CpG;
  array_list_t *context_CHG = bs_context->context_CHG;
  array_list_t *context_CHH = bs_context->context_CHH;
  array_list_t *context_MUT = bs_context->context_MUT;

  size_t num_items;
  int file_error;
  metil_data_t *metil_data;

  FILE * CpG = metil_file->CpG;
  FILE * CHG = metil_file->CHG;
  FILE * CHH = metil_file->CHH;
  FILE * MUT = metil_file->MUT;

  metil_file->CpG_methyl   += bs_context->CpG_methyl;
  metil_file->CpG_unmethyl += bs_context->CpG_unmethyl;
  metil_file->CHG_methyl   += bs_context->CHG_methyl;
  metil_file->CHG_unmethyl += bs_context->CHG_unmethyl;
  metil_file->CHH_methyl   += bs_context->CHH_methyl;
  metil_file->CHH_unmethyl += bs_context->CHH_unmethyl;
  metil_file->MUT_methyl   += bs_context->MUT_methyl;
  metil_file->num_bases    += bs_context->num_bases;

  // Accumulate per-read methylation statistics
  for (size_t i = 0; i < num_chromosomes; ++i) {
    metil_file->methyl_reads[i] += bs_context->methyl_reads[i];
  }

  if (CpG == NULL && !metil_file->skip_cpg) {
    CpG = fopen(metil_file->filenameCpG, "a");
  }

  if (CHG == NULL && !metil_file->skip_chg) {
    CHG = fopen(metil_file->filenameCHG, "a");
  }

  if (CHH == NULL && !metil_file->skip_chh) {
    CHH = fopen(metil_file->filenameCHH, "a");
  }

  if (MUT == NULL && !metil_file->skip_mut) {
    MUT = fopen(metil_file->filenameMUT, "a");
  }

  if (context_CpG != NULL && CpG) {
    num_items = array_list_size(context_CpG);

    for (int i = num_items - 1; i >= 0; i--) {
      metil_data = (metil_data_t *)array_list_get(i, context_CpG);
      file_error = fprintf(CpG, "%s\t%c\t%i\t%lu\t%c\t%i\n",
                          metil_data->query_name, metil_data->status,
                          metil_data->chromosome, metil_data->start,
                          metil_data->context, metil_data->strand);
      metil_data_free(metil_data);
      
      if (file_error < 0) {
        printf("Error writting the CpG context file: %s\n\n", strerror(errno));
        printf("Written %i of %lu CpG methylated bases, skipping the remaining bases.\nCheck the BAM file to retrieve all the information.\n\n", i, num_items);

        fclose(metil_file->CpG);
        metil_file->skip_cpg = 1;

        break;
      }
    }
  }

  if (context_CHG != NULL && CHG) {
    num_items = array_list_size(context_CHG);

    for (int i = num_items - 1; i >= 0; i--) {
      metil_data = (metil_data_t *)array_list_get(i, context_CHG);
      file_error = fprintf(CHG, "%s\t%c\t%i\t%lu\t%c\t%i\n",
                          metil_data->query_name, metil_data->status,
                          metil_data->chromosome, metil_data->start,
                          metil_data->context, metil_data->strand);
      metil_data_free(metil_data);
      
      if (file_error < 0) {
        printf("Error writting the CHG context file: %s\n\n", strerror(errno));
        printf("Written %i of %lu CHG methylated bases, skipping the remaining bases.\nCheck the BAM file to retrieve all the information.\n\n", i, num_items);
        
        fclose(metil_file->CHG);
        metil_file->skip_chg = 1;

        break;
      }
    }
  }

  if (context_CHH != NULL && CHH) {
    num_items = array_list_size(context_CHH);

    for (int i = num_items - 1; i >= 0; i--) {
      metil_data = (metil_data_t *)array_list_get(i, context_CHH);
      file_error = fprintf(CHH, "%s\t%c\t%i\t%lu\t%c\t%i\n",
                          metil_data->query_name, metil_data->status,
                          metil_data->chromosome, metil_data->start,
                          metil_data->context, metil_data->strand);
      metil_data_free(metil_data);
      
      if (file_error < 0) {
        printf("Error writting the CHH context file: %s\n\n", strerror(errno));
        printf("Written %i of %lu CHH methylated bases, skipping the remaining bases.\nCheck the BAM file to retrieve all the information.\n\n", i, num_items);
        
        fclose(metil_file->CHH);
        metil_file->skip_chh = 1;

        break;
      }
    }
  }
 
  if (context_MUT != NULL && MUT) {
    num_items = array_list_size(context_MUT);
    for (int i = num_items - 1; i >= 0; i--) {
      metil_data = (metil_data_t *)array_list_get(i, context_MUT);
      file_error = fprintf(MUT, "%s\t%c\t%i\t%lu\t%c\t%i\n",
                          metil_data->query_name, metil_data->status,
                          metil_data->chromosome, metil_data->start,
                          metil_data->context, metil_data->strand);
      metil_data_free(metil_data);
      
      if (file_error < 0) {
        printf("Error writting the MUT context file: %s\n\n", strerror(errno));
        printf("Written %i of %lu MUT methylated bases, skipping the remaining bases.\nCheck the BAM file to retrieve all the information.\n\n", i, num_items);
        
        fclose(metil_file->MUT);
        metil_file->skip_mut = 1;
        
        break;
      }
    }
  }

  if (context_CpG) array_list_free(context_CpG, NULL);
  if (context_CHG) array_list_free(context_CHG, NULL);
  if (context_CHH) array_list_free(context_CHH, NULL);
  if (context_MUT) array_list_free(context_MUT, NULL);
  if (bs_context)  bs_context_free(bs_context);
}

//------------------------------------------------------------------------------------

void metil_data_init(metil_data_t *metil_data, char *query, char status, int chromosome, size_t start, char context, int strand, int zone) {
  if (metil_data == NULL)
    metil_data = (metil_data_t *)malloc(sizeof(metil_data_t));

  metil_data->query_name = strdup(query);
  metil_data->status = status;
  metil_data->chromosome = chromosome;
  metil_data->start = start;
  metil_data->context = context;
  metil_data->strand = strand;
  metil_data->zone = zone;
}

//------------------------------------------------------------------------------------

void metil_data_free(metil_data_t *metil_data) {
  if (metil_data != NULL) {
    free(metil_data->query_name);
    free(metil_data);
  }
}

//------------------------------------------------------------------------------------

void postprocess_bs(char *query_name, char status, size_t chromosome, size_t start, char context, int strand, int region,
		    array_list_t *list) {

  metil_data_t *metil_data = (metil_data_t *) malloc(sizeof(metil_data_t));
  memset(metil_data, 0, sizeof(metil_data_t));

  metil_data_init(metil_data, query_name, status, chromosome, start, context, strand, region);
  array_list_insert(metil_data, list);
}

//------------------------------------------------------------------------------------

int encode_context(char* filename, char* directory, int max_num_chromosomes) {
  printf("Init Genome Compresion\n");

  FILE *f1, *f2, *f3, *f4;
  unsigned long long value, value2;
  size_t *size = (size_t*)malloc(max_num_chromosomes*sizeof(size_t)); //Tamaño máximo de cromosoma 24 para humanos
  for (int i=0;i<max_num_chromosomes;i++)
	  size[i] = 0;



  size_t size2;
  int chromosome, i, cont, cont2;
  char *line1, *line2;
  char *tmp;
  size_t contador = 100000000;

  int elem = sizeof(unsigned long long) << 2;

  size2 = strlen(directory);
  tmp = (char *)malloc((size2 + 40) * sizeof(char));
  line1 = (char *)malloc(512 * sizeof(char));
  line2 = (char *)malloc(512 * sizeof(char));

  f1 = fopen (filename, "r");

  if (f1 == NULL) {
    perror("No se puede abrir el fichero de entrada");
    return -1;
  }

  sprintf(tmp, "%s/Genome_context_CT.bin", directory);
  f2 = fopen (tmp, "wb");

  if (f2 == NULL) {
    perror("No se puede abrir el fichero de contexto CT");
    return -1;
  }

  sprintf(tmp, "%s/Genome_context_GA.bin", directory);
  f3 = fopen (tmp, "wb");

  if (f3 == NULL) {
    perror("No se puede abrir el fichero de contexto GA");
    return -1;
  }

  sprintf(tmp, "%s/Genome_context_size.txt", directory);
  f4 = fopen (tmp, "w");

  if (f4 == NULL) {
    perror("No se puede abrir el fichero de tamaños");
    return -1;
  }
  
  chromosome = 0;
  size[chromosome] = 0;
  cont  = 0;
  cont2 = 0;

  // Descartar la primera linea
  do {
    fgets(line1, 512, f1);
  } while(line1[0] == '>');
  
  
  while (fgets(line2, 512, f1) != NULL && contador) {
	  if (line2[0] == '\n')
	  {
		  char* cad = fgets(line2, 512, f1);
		  if (cad == NULL)
			  break;
	  }

    contador--;

    if (line1[0] == '>') {
      size[chromosome]++;
      printf("chromosomes = %2i, size = %10lu\n", chromosome, size[chromosome]);
      chromosome++;

      for (; cont > 0 && cont <= elem; cont++) {
	      value = value << 2;
      }

      for (; cont2 > 0 && cont2 <= elem; cont2++) {
	      value2 = value2 << 2;
      }

      fwrite(&value,  sizeof(unsigned long long), 1, f2);
      fwrite(&value2, sizeof(unsigned long long), 1, f3);

      cont  = 0;
      cont2 = 0;
    } else {
      size2 = strlen(line1);

      // Value for C->T conversion
      for (i = 0; i < (int)size2 - 2; i++, cont++) {
        if (line1[i] == 'C') {
          if (line1[i + 1] == 'G') {
            value += 1;
          } else {
            if (line1[i + 2] == 'G') {
              value += 2;
            } else {
              value += 3;
            }
          }
        }
        if (cont == elem) {
          cont = 0;
          size[chromosome]++;
          fwrite(&value, sizeof(unsigned long long), 1, f2);
        } else {
          value = value << 2;
        }
      }
      
      if (line1[size2 - 2] == 'C') {
        if (line1[size2 - 1] == 'G') {
          value += 1;
        } else {
          if (line2[0] == 'G') {
            value += 2;
          } else {
            value += 3;
          }
        }
      }

      if (cont == elem) {
        cont = 0;
        size[chromosome]++;
        fwrite(&value, sizeof(unsigned long long), 1, f2);
      } else {
	      value = value << 2;
      }

      cont++;

      if (line1[size2 - 1] == 'C') {
        if (line2[0] == 'G') {
          value += 1;
        } else {
          if (line2[1] == 'G') {
            value += 2;
          } else {
            value += 3;
          }
        }
      }

      if (cont == elem) {
        cont = 0;
        size[chromosome]++;
        fwrite(&value, sizeof(unsigned long long), 1, f2);
      } else {
	      value = value << 2;
      }
      cont++;
      // End value for C->T conversion

      // Value for G->A conversion
      for (i = 2; i < size2; i++, cont2++) {
        if (line1[i] == 'G') {
          if (line1[i - 1] == 'C') {
            value2 += 1;
          } else {
            if (line1[i - 2] == 'C') {
              value2 += 2;
            } else {
              value2 += 3;
            }
          }
        }
        if (cont2 == elem) {
          cont2 = 0;
          fwrite(&value2, sizeof(unsigned long long), 1, f3);
        } else {
          value2 = value2 << 2;
        }
      }
      
      if (line2[0] == 'G') {
        if (line1[size2 - 1] == 'C') {
          value2 += 1;
        } else {
          if (line1[size2 - 2] == 'C') {
            value2 += 2;
          } else {
            value2 += 3;
          }
        }
      }

      if (cont2 == elem) {
        cont2 = 0;
        fwrite(&value2, sizeof(unsigned long long), 1, f3);
      } else {
	      value2 = value2 << 2;
      }

      cont2++;

      if (line2[1] == 'G') {
        if (line2[0] == 'C') {
          value2 += 1;
        } else {
          if (line1[size2 - 1] == 'C') {
            value2 += 2;
          } else {
            value2 += 3;
          }
        }
      }

      if (cont2 == elem) {
        cont2 = 0;
        fwrite(&value2, sizeof(unsigned long long), 1, f3);
      } else {
	      value2 = value2 << 2;
      }

      cont2++;
      // End value for G->A conversion
    }

    free(line1);
    line1 = strdup(line2);
  }
  printf("chromosomes = %2i, size = %10lu\n", chromosome, size[chromosome]);

  fprintf(f4, "%i\n", chromosome);

  for (i = 0; i < chromosome; i++) {
    fprintf(f4, "%lu\n", size[i]);
  }

  free(size);
  free(tmp);
  free(line1);
  free(line2);
  fclose(f1);
  fclose(f2);
  fclose(f3);
  fclose(f4);

  printf("End Genome Compresion\n");
  return 0;
}

//------------------------------------------------------------------------------------

int load_encode_context(char* directory, unsigned long long **valuesCT, unsigned long long **valuesGA) {
  FILE *f1, *f2, *f3;
  size_t size, size2 = strlen(directory);
  char *tmp = (char *)malloc((size2 + 50) * sizeof(char));
  int ret, chromosome, i;

  sprintf(tmp, "%s/Genome_context_size.txt", directory);
  f1 = fopen (tmp, "r");

  if (f1 == NULL) {
    perror("No se puede abrir el fichero de tamaños");
    ret = -1;
    goto end1;
  }

  sprintf(tmp, "%s/Genome_context_CT.bin", directory);
  f2 = fopen (tmp, "rb");
  
  if (f2 == NULL) {
    perror("No se puede abrir el fichero de contexto CT");
    ret = -1;
    goto end2;
  }

  sprintf(tmp, "%s/Genome_context_GA.bin", directory);
  f3 = fopen (tmp, "rb");
  
  if (f3 == NULL) {
    perror("No se puede abrir el fichero de contexto GA");
    ret = -1;
    goto end3;
  }

  fscanf(f1, "%i\n", &chromosome);

  for (i = 0; i < chromosome; i++) {
    fscanf(f1, "%lu\n", &size);

    valuesCT[i] = (unsigned long long *)calloc(size, sizeof(unsigned long long));
    fread (valuesCT[i], sizeof(unsigned long long), size, f2);

    valuesGA[i] = (unsigned long long *)calloc(size, sizeof(unsigned long long));
    fread (valuesGA[i], sizeof(unsigned long long), size, f3);
  }

  ret = chromosome;

  fclose(f3);
end3:
  fclose(f2);
end2:
  fclose(f1);
end1:
  free(tmp);

  return ret;
}
