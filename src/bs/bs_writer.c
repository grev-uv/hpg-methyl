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

  // Set the sequences of the mapping to the original
  revert_mappings_seqs(mapping_batch->mapping_lists, mapping_batch->mapping_lists2, mapping_batch->fq_batch);
  
  batch_writer_input_t *writer_input = batch->writer_input;
  bam_file_t *bam_file = writer_input->bam_file;     
  size_t num_reads = array_list_size(mapping_batch->fq_batch);
  size_t num_mapped_reads = 0;
  size_t total_mappings = 0;
  
  
  writer_input->total_batches++;
  
  array_list_t **mapping_lists;
  int *found = (int *) calloc(num_reads, sizeof(int));
  metil_file_t *metil_file = writer_input->metil_file;

  size_t num_chromosomes = batch->bwt_input->genome->num_chromosomes;

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

  if (mapping_batch) {
    mapping_batch_free(mapping_batch);
  }
  
  if (batch) {
    batch_free(batch);
  }
  
  if (found) {
    free(found);
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
