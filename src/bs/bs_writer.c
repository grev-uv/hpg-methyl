#include "bs_writer.h"
#include "methylation.h"


//------------------------------------------------------------------------------------

int bs_writer(void *data) {
  struct timeval start, end;
  double time;
  
  if (time_on) { 
    start_timer(start); 
  }
  
  batch_t *batch = (batch_t *) data;
  fastq_read_t *fq_read;
  size_t num_items;
  
  mapping_batch_t *mapping_batch = (mapping_batch_t *) batch->mapping_batch;
  bs_context_t *bs_context = (mapping_batch->bs_context);

  // Set the sequences of the mapping to the original //RICARDO adapted to avoid cigar problem with hard clipping
  revert_mappings_seqs(mapping_batch->mapping_lists, mapping_batch->mapping_lists2, mapping_batch->fq_batch);
  
  batch_writer_input_t *writer_input = batch->writer_input;
  bam_file_t *bam_file = writer_input->bam_file;     
  size_t num_reads = array_list_size(mapping_batch->fq_batch);
  size_t num_mapped_reads = 0;
  size_t total_mappings = 0;
  
  
  writer_input->total_batches++;
  
  array_list_t **mapping_lists;
  int *found = (int *) calloc(num_reads, sizeof(int));
  alignment_t ** best_alignments = ((alignment_t **) calloc (num_reads, sizeof(alignment_t *)));
  metil_file_t *metil_file = writer_input->metil_file;

  size_t num_chromosomes = batch->bwt_input->genome->num_chromosomes;



  pair_server_input_t *input = batch->pair_input;

  size_t end1, start2;
  short int chr1, chr2, strand1, strand2;
  int pair_mode = input->pair_mng->pair_mode;
  size_t min_distance = input->pair_mng->min_distance;
  size_t max_distance = input->pair_mng->max_distance;
  int distance;

  if (input->pair_mng->pair_mode == SINGLE_END_MODE) {

  // Process mapping_lists and mapping_lists2
  for (int k = 0; k < 2; k++) {
    mapping_lists = (k == 0) ? mapping_batch->mapping_lists : mapping_batch->mapping_lists2;
    
    for (size_t i = 0; i < num_reads; i++) {
      num_items = array_list_size(mapping_lists[i]);
      total_mappings += num_items;
      fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->fq_batch);




      // Mapped or not mapped ?	 
      if (num_items == 0) {
        if (mapping_lists[i]) {
          array_list_free(mapping_lists[i], NULL);
        }
      } else {
        found[i] = 1;
        LOG_DEBUG_F("printing alignment (bs_id = %i)...\n", k);
        write_mapped_read(mapping_lists[i], bam_file);
      }
    }
  }
  
  for (size_t i = 0; i < num_reads; i++) {
    if (found[i]) {
      num_mapped_reads++;
    } else {
      total_mappings++;
      fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->fq_batch);
      write_unmapped_read(fq_read, bam_file);
    }
  }

  }
  else //Paired-end-mode
  {
	  array_list_t *list1, *list2;
	  size_t num_items1, num_items2;
	  alignment_t *alig1, *alig2;
	  //Ricardo
	  for (size_t i = 0; i < num_reads; i += 2) {

	 	  //In order to differentiate list mode (mapping_list and mapping_list2, or vice versa)
	 	for (int option = 0; option < 2; option++) {


	 	if (option == 0)
	 	{

	 		list1 = mapping_batch->mapping_lists[i];
	 		list2 = mapping_batch->mapping_lists2[i + 1];


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


	     	num_items1 = 0;

	     	if (list1 != NULL) {
	     	    num_items1 = array_list_size(list1);
	     	}

	     	num_items2 = 0;

	     	if (list2 != NULL) {
	     	    num_items2 = array_list_size(list2);
	     	}
	     }

	 	 total_mappings += num_items1 + num_items2;

	     if (num_items1 > 0 && num_items2 > 0)
	     {
	    	 found[i] = 1;
	    	 found[i+1] = 1;


	    	 if (list1->paired)
	    	 {	//read and mate map in the same area
	    		if (num_items1 == 1 && num_items2 == 1)	//only one paired option
	    		{
	    			 alig1 = (alignment_t *) array_list_get(0, list1);
	    			 alig2 = (alignment_t *) array_list_get(0, list2);
	    			 alignment_set_paired_end_data(alig1, alig2);


	    		}
	    		else
	    		{	//looking for paired read

	    			 for (size_t j1 = 0; j1 < num_items1; j1++) {
	    			       alig1 = (alignment_t *) array_list_get(j1, list1);
	    			       chr1 = alig1->chromosome;
	    			       strand1 = alig1->seq_strand;
	    			       end1 = alig1->position + strlen(alig1->sequence);


	    			       for (size_t j2 = 0; j2 < num_items2; j2++) {


	    			           alig2 = (alignment_t *) array_list_get(j2, list2);
	    			           chr2 = alig2->chromosome;
	    			           strand2 = alig2->seq_strand;
	    			           start2 = alig2->position;


	    			          // computes distance between alignments,
	    			          // is a valid distance ?
	    			          distance = (start2 > end1 ? start2 - end1 : end1 - start2);

	    			          if ((chr1 == chr2)                                           &&
	    			              (distance >= min_distance) && (distance <= max_distance) &&
	    			             ((strand1 != strand2 && pair_mode == PAIRED_END_MODE)     ||
	    			              (strand1 == strand2 && pair_mode == MATE_PAIR_MODE))) {



	    			            alignment_set_paired_end_data(alig1, alig2);

	    			            break;
	    			          }
	    			        } // end for j2..num_items2
	    			      } // end for j1..num_item1


	    		}

	    	 }
	    	 else
	    	 {	//read and mate not map in the same area (but both map)

	    		 //Looking for the best quality mapping
	    		 int best_list1_index = 0;
	    		 int best_list2_index = 0;
	    		 int max_quality = -1;

	    		 for (size_t j1 = 0; j1 < num_items1; j1++) {
	    			  alig1 = (alignment_t *) array_list_get(j1, list1);
	    			  if (alig1->map_quality > max_quality)
	    			  {
	    				  max_quality = alig1->map_quality;
	    				  best_list1_index = j1;
	    			  }
	    		 }

	    		 max_quality = -1;
	    		 for (size_t j2 = 0; j2 < num_items2; j2++) {
	    			   alig2 = (alignment_t *) array_list_get(j2, list2);
	    			   if (alig2->map_quality > max_quality)
	    			   {
	    			   	  max_quality = alig2->map_quality;
	    			   	  best_list2_index = j2;
	    			   }
	    		}

	    		for (size_t j1 = 0; j1 < num_items1; j1++) {
	    			  alig1 = (alignment_t *) array_list_get(j1, list1);
	    			  alig2 = (alignment_t *) array_list_get(best_list2_index, list2);
	    			  alignment_set_paired_end_data_best(alig1, alig2,1);
	    		}

	    		for (size_t j2 = 0; j2 < num_items2; j2++) {
	    			  alig2 = (alignment_t *) array_list_get(j2, list2);
	    			  alig1 = (alignment_t *) array_list_get(best_list1_index, list1);
	    			  alignment_set_paired_end_data_best(alig2, alig1, 2);
	    		}



	    	 }



	     }
	     else  //only one (read or mate) map
	     {
	    	 	 int best_list_index = 0;
	    	 	 int max_quality = -1;

	    		 if (num_items1 > 0) {	//pair1 map, pair2 not map
	    			 found[i] = 1;
	    			 for (size_t j1 = 0; j1 < num_items1; j1++) {
	    				 alig1 = (alignment_t *) array_list_get(j1, list1);
	    				 alignment_set_paired_end_data_pair2none(alig1, 1);

	    				 if (alig1->map_quality > max_quality)
	    				 {
	    				 	    max_quality = alig1->map_quality;
	    				 	    best_list_index = j1;
	    				 }
	    			 }

	    			 if (best_alignments[i] != NULL)
	    			 {
	    				 if (max_quality > best_alignments[i]->map_quality)
	    					 best_alignments[i] = (alignment_t *) array_list_get(best_list_index, list1);
	    			 }
	    			 else
	    				 best_alignments[i] = (alignment_t *) array_list_get(best_list_index, list1);

	    		 }


	    		 if (num_items2 > 0) {	//pair2 map, pair 1 not map
	    			 found[i+1] = 1;
	    			 for (size_t j2 = 0; j2 < num_items2; j2++) {
	    				 alig2 = (alignment_t *) array_list_get(j2, list2);
	    				 alignment_set_paired_end_data_pair2none(alig2, 2);

	    				 if (alig2->map_quality > max_quality)
	    				 {
	    					    max_quality = alig2->map_quality;
	    					    best_list_index = j2;
	    				 }
	    			 }

	    			 if (best_alignments[i+1] != NULL)
	    			 {
	    				 if (max_quality > best_alignments[i+1]->map_quality)
	    					 best_alignments[i+1] = (alignment_t *) array_list_get(best_list_index, list2);
	    			 }
	    			 else
	    				 best_alignments[i+1] = (alignment_t *) array_list_get(best_list_index, list2);


	    		 }

	    }

	   }//end_options


	  }//end num_reads



	  for (size_t i = 0; i < num_reads; i+=2) {

		  //First write unmapped
		 if ((!found[i]) && (!found[i+1])) //none of them (read and mate) map
		 {
			  total_mappings++;
			  fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->fq_batch);
			  write_unmapped_read_paired_none(fq_read, bam_file, 1);

			  total_mappings++;
			  fq_read = (fastq_read_t *) array_list_get(i+1, mapping_batch->fq_batch);
			  write_unmapped_read_paired_none(fq_read, bam_file, 2);
		 }
		 else{
			 if ((found[i]) && (!found[i+1])) //pair1 map, pair2 not map
			 {
				 total_mappings++;
				 fq_read = (fastq_read_t *) array_list_get(i+1, mapping_batch->fq_batch);
				 write_unmapped_read_paired_map(fq_read, bam_file, 2, best_alignments[i]);

			 }
			 else if ((!found[i]) && (found[i+1])) //pair1 not map, pair2 map
			 {
				 total_mappings++;
				 fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->fq_batch);
				 write_unmapped_read_paired_map(fq_read, bam_file, 1, best_alignments[i+1]);

			 }

		 }

		 //write mapped

		  if (found[i]) {
	        num_mapped_reads++;

	        if ((mapping_batch->mapping_lists[i] != NULL) && (array_list_size(mapping_batch->mapping_lists[i])>0)){
	        	 write_mapped_read(mapping_batch->mapping_lists[i], bam_file);
	        }
	        else
	        	 array_list_free(mapping_batch->mapping_lists[i], NULL);

	        if ((mapping_batch->mapping_lists2[i] != NULL) && (array_list_size(mapping_batch->mapping_lists2[i])>0)){
	       	     write_mapped_read(mapping_batch->mapping_lists2[i], bam_file);
	       	}
	        else
	       	     array_list_free(mapping_batch->mapping_lists2[i], NULL);
		  }
		  else
		  {
			  array_list_free(mapping_batch->mapping_lists[i], NULL);
			  array_list_free(mapping_batch->mapping_lists2[i], NULL);
		  }


		  if (found[i+1]) {
		 	num_mapped_reads++;

		 	if ((mapping_batch->mapping_lists[i+1] != NULL) && (array_list_size(mapping_batch->mapping_lists[i+1])>0)){
		 	  	 write_mapped_read(mapping_batch->mapping_lists[i+1], bam_file);
		 	}
		 	else
		 		 array_list_free(mapping_batch->mapping_lists[i+1], NULL);

		 	if ((mapping_batch->mapping_lists2[i+1] != NULL) && (array_list_size(mapping_batch->mapping_lists2[i+1])>0)){
		 	     write_mapped_read(mapping_batch->mapping_lists2[i+1], bam_file);
		 	}
		 	else
		 	   	 array_list_free(mapping_batch->mapping_lists2[i+1], NULL);

		  }
		  else
		  {
			  	 array_list_free(mapping_batch->mapping_lists[i+1], NULL);
		 		 array_list_free(mapping_batch->mapping_lists2[i+1], NULL);
		  }


	   }

  }//end paired-end mode


  // To use with the postprocess,  write metil context files
  if (batch->write_mcontext) {
	  write_bs_context(metil_file, bs_context, num_chromosomes);
  }
  
  if (basic_st->total_reads >= writer_input->limit_print) {
    LOG_DEBUG_F("TOTAL READS PROCESS: %lu\n", basic_st->total_reads);
    LOG_DEBUG_F("\tTotal Reads Mapped: %lu(%.2f%%)\n", 
		basic_st->num_mapped_reads, 
		(float) (basic_st->num_mapped_reads*100)/(float)(basic_st->total_reads));
    writer_input->limit_print += 1000000;
  }
  //printf("HOLA\n");
  if (mapping_batch) {
    mapping_batch_free(mapping_batch);
  }
  
  if (batch) {
    batch_free(batch);
  }
  
  if (found) {
    free(found);
  }


  if (best_alignments) {
    free(best_alignments);
  }

  basic_statistics_add(num_reads, num_mapped_reads, total_mappings, basic_st);
  
  if (time_on) { 
    stop_timer(start, end, time); 
    timing_add(time, BAM_WRITER, timing); 
  }

  return 0;
}

//--------------------------------------------------------------------

void write_mapped_read(array_list_t *array_list, bam_file_t *bam_file) {
  size_t num_items = array_list_size(array_list);
  alignment_t *alig;
  bam1_t *bam1;

  for (size_t j = 0; j < num_items; j++) {
    alig = (alignment_t *) array_list_get(j, array_list);
    LOG_DEBUG("writting bam..\n");

    if (alig != NULL) {
      //Ricardo - quitar
    	//alig->quality[0] = "\0";

      bam1 = convert_to_bam(alig, 33);
      bam_fwrite(bam1, bam_file);
      bam_destroy1(bam1);
      alignment_free(alig);
    } else {
      LOG_FATAL_F("alig is NULL, num_items = %lu\n", num_items)
    }
  }

  if (array_list) { 
    array_list_free(array_list, NULL); 
  }
}

//--------------------------------------------------------------------

void write_unmapped_read(fastq_read_t *fq_read, bam_file_t *bam_file) {
  static char aux[1024];
  alignment_t *alig;
  size_t header_len;
  char *id;
  bam1_t *bam1;

  // Calculating cigar
  sprintf(aux, "%iX", fq_read->length);

  alig = alignment_new();
  header_len = strlen(fq_read->id);
  id = (char *) malloc(sizeof(char) * (header_len + 1));

  get_to_first_blank(fq_read->id, header_len, id);
  
  alignment_init_single_end(id, fq_read->sequence, fq_read->quality,
			    0, -1, -1, aux, 1, 0, 0, 0, 0, NULL, alig);

  bam1 = convert_to_bam(alig, 33);
  bam_fwrite(bam1, bam_file);
  bam_destroy1(bam1);

  alig->sequence = NULL;
  alig->quality = NULL;
  alig->cigar = NULL;

  alignment_free(alig);
}

//both (read and mate) not map
void write_unmapped_read_paired_none(fastq_read_t *fq_read, bam_file_t *bam_file, int pair_num) {
  static char aux[1024];
  alignment_t *alig;
  size_t header_len;
  char *id;
  bam1_t *bam1;

  // Calculating cigar
  sprintf(aux, "%iX", fq_read->length);

  alig = alignment_new();
  header_len = strlen(fq_read->id);
  id = (char *) malloc(sizeof(char) * (header_len + 1));

  get_to_first_blank(fq_read->id, header_len, id);

  alignment_init_paired_end_none(id, fq_read->sequence, fq_read->quality,
			    0, -1, -1, aux, 1, 0, 0, 0, 0, NULL, alig, pair_num);

  bam1 = convert_to_bam(alig, 33);
  bam_fwrite(bam1, bam_file);
  bam_destroy1(bam1);

  alig->sequence = NULL;
  alig->quality = NULL;
  alig->cigar = NULL;

  alignment_free(alig);
}


//read not map, and mate map
void write_unmapped_read_paired_map(fastq_read_t *fq_read, bam_file_t *bam_file, int pair_num, alignment_t* alig2) {
  static char aux[1024];
  alignment_t *alig;
  size_t header_len;
  char *id;
  bam1_t *bam1;

  // Calculating cigar
  sprintf(aux, "%iX", fq_read->length);

  alig = alignment_new();
  header_len = strlen(fq_read->id);
  id = (char *) malloc(sizeof(char) * (header_len + 1));

  get_to_first_blank(fq_read->id, header_len, id);

  alignment_init_paired_end_map(id, fq_read->sequence, fq_read->quality,
			    0, -1, -1, aux, 1, 0, 0, 0, 0, NULL, alig, pair_num,
				alig2->seq_strand, alig2->chromosome, alig2->position);



  bam1 = convert_to_bam(alig, 33);
  bam_fwrite(bam1, bam_file);
  bam_destroy1(bam1);

  alig->sequence = NULL;
  alig->quality = NULL;
  alig->cigar = NULL;

  alignment_free(alig);
}
