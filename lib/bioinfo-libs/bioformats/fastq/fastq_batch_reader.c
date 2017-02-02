#include "fastq_batch_reader.h"

//------------------------------------------------------------------------------------
// init structure for the batch reader
//------------------------------------------------------------------------------------

void fastq_batch_reader_input_init(char *filename1, char *filename2,
				   int flags, int batch_size, list_t *list, 
				   fastq_batch_reader_input_t *input) {
  input->filename1 = filename1;
  input->filename2 = filename2;
  input->flags = flags;
  input->batch_size = batch_size;
  input->list = list;

  // internal
  input->fq_file1 = NULL;
  input->fq_file2 = NULL;
}

//------------------------------------------------------------------------------------
// this functions reads read from disk, and save them
// into a list to be processed by the next thread
//------------------------------------------------------------------------------------

void fastq_batch_reader_single(fastq_batch_reader_input_t* input);
void fastq_batch_reader_pair(fastq_batch_reader_input_t* input);
void fastq_batch_reader_aligner_pair(fastq_batch_reader_input_t* input);
//------------------------------------------------------------------------------------

void fastq_batch_reader(fastq_batch_reader_input_t* input) {
  if (input->flags == SINGLE_END_MODE) {
    fastq_batch_reader_single(input);
  } else {
    fastq_batch_reader_pair(input);
  }
}

//------------------------------------------------------------------------------------

void fastq_batch_reader_aligner(fastq_batch_reader_input_t* input) {
  if (input->flags == SINGLE_END_MODE) {
    fastq_batch_reader_single(input);
  } else {
    fastq_batch_reader_aligner_pair(input);
  }
}


//------------------------------------------------------------------------------------
      
void fastq_batch_reader_single(fastq_batch_reader_input_t* input) {
  unsigned int total_reads = 0;

  /*struct timespec ts;
  ts.tv_sec = 1;
  ts.tv_nsec = 0;*/

  char *filename = input->filename1;
  int flags = input->flags & 3;
  size_t bytes = input->batch_size;
  list_t *list = input->list;

  printf("fastq_batch_reader SINGLE mode (%i): START, for file %s\n", 
	 omp_get_thread_num(), filename);
		
  size_t num_reads = 0, num_batches = 0;
  
  array_list_t *reads; //= array_list_new(1000, 
    //		       1.25f, 
    //				       COLLECTION_MODE_ASYNCHRONIZED);

  list_item_t *item = NULL;

  fastq_file_t *file = fastq_fopen(filename);

  while (1) {
    // read reads from file
    reads = array_list_new(10000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    fastq_fread_bytes_se(reads, bytes, file);

    num_reads = array_list_size(reads);
    total_reads	+= num_reads;
    
    // if no reads, free memory and go out....
    if (num_reads == 0) {
      array_list_free(reads, (void *)fastq_read_free);
      break;
    }

    // otherwise, create a new batch object..
    // and insert this batch to the corresponding list
    item = list_item_new(num_batches, flags, reads);
    //printf("fastq_ex_batch_reader: reading batch %i done !!\n", num_batches);
    list_insert_item(item, list);
    
    
    num_batches++;
  } // end of batch loop
  
  list_decr_writers(list);
  fastq_fclose(file);
  printf("fastq_batch_reader SINGLE mode: END, %i total reads (%i batches), for file %s\n", 
	 total_reads, num_batches, filename);
}

//------------------------------------------------------------------------------------

void fastq_batch_reader_pair(fastq_batch_reader_input_t* input) {
 unsigned int total_reads = 0;
  /*struct timespec ts;
  ts.tv_sec = 1;
  ts.tv_nsec = 0;*/

  char *filename1 = input->filename1;
  char *filename2 = input->filename2;
  int flags = input->flags;
  size_t bytes = input->batch_size;
  list_t *list = input->list;

  printf("fastq_batch_reader SINGLE mode (%i): START, for file %s\n", 
	 omp_get_thread_num(), filename1);
		
  size_t num_reads = 0, num_batches = 0;
  
  array_list_t *reads;

  list_item_t *item = NULL;

  fastq_file_t *file1 = fastq_fopen(filename1);
  fastq_file_t *file2 = fastq_fopen(filename2);

  while (1) {
    // read reads from file
    reads = array_list_new(10000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    //fastq_fread_bytes_se(reads, bytes, file);
    fastq_fread_bytes_pe(reads, bytes, file1, file2, flags);
    num_reads = array_list_size(reads);
    total_reads	+= num_reads;
    
    // if no reads, free memory and go out....
    if (num_reads == 0) {
      array_list_free(reads, (void *)fastq_read_free);
      break;
    }

    // otherwise, create a new batch object..
    // and insert this batch to the corresponding list
    item = list_item_new(num_batches, flags, reads);

    list_insert_item(item, list);
  
    //printf("fastq_ex_batch_reader: reading batch %i done !!\n", num_batches);
    num_batches++;
  } // end of batch loop
  
  list_decr_writers(list);
  fastq_fclose(file1);
  printf("fastq_batch_reader SINGLE mode: END, %i total reads (%i batches), for file %s\n", 
	 total_reads, num_batches, filename1);

}
//------------------------------------------------------------------------------------


void fastq_batch_reader_aligner_pair(fastq_batch_reader_input_t* input) {
 unsigned int total_reads = 0;
  /*struct timespec ts;
  ts.tv_sec = 1;
  ts.tv_nsec = 0;*/

  char *filename1 = input->filename1;
  char *filename2 = input->filename2;
  int flags = input->flags;
  size_t bytes = input->batch_size;
  list_t *list = input->list;

  printf("fastq_batch_reader PAIR mode (%i): START for files %s, %s\n", 
	 omp_get_thread_num(), filename1, filename2);
		
  size_t num_reads = 0, num_batches = 0;
  
  array_list_t *reads;

  list_item_t *item = NULL;

  fastq_file_t *file1 = fastq_fopen(filename1);
  fastq_file_t *file2 = fastq_fopen(filename2);

  while (1) {
    // read reads from file
    reads = array_list_new(10000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    //fastq_fread_bytes_se(reads, bytes, file);
    fastq_fread_bytes_aligner_pe(reads, bytes, file1, file2);
    num_reads = array_list_size(reads);
    total_reads	+= num_reads;
    
    // if no reads, free memory and go out....
    if (num_reads == 0) {
      array_list_free(reads, (void *)fastq_read_free);
      break;
    }

    // otherwise, create a new batch object..
    // and insert this batch to the corresponding list
    item = list_item_new(num_batches, flags, reads);

    list_insert_item(item, list);
  
    //printf("fastq_ex_batch_reader: reading batch %i done !!\n", num_batches);
    num_batches++;
  } // end of batch loop
  
  list_decr_writers(list);
  fastq_fclose(file1);
  fastq_fclose(file2);
  printf("fastq_batch_reader PAIR mode: END, %i total reads (%i batches), for files %s, %s\n", 
	 total_reads, num_batches, filename1, filename2);

}
