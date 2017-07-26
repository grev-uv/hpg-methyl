#ifndef FASTQ_BATCH_READER_H
#define FASTQ_BATCH_READER_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <omp.h>

#include "containers/list.h"
#include "fastq_batch.h"
#include "fastq_file.h"
#include "fastq_gzfile.h"

//#include "timing.h"
//#include "buffers.h"
//#include "statistics.h"

//====================================================================================

#define READ_ITEM       1

//====================================================================================

#define SINGLE_END_MODE 0
#define PAIRED_END_MODE 1
#define MATE_PAIR_MODE  2

//====================================================================================

#define SINGLE_END_FLAG    4
#define PAIRED_END_FLAG    8
#define MATE_PAIR_FLAG    16
#define PAIR1_FLAG        32
#define PAIR2_FLAG        64

//====================================================================================
//  structures and prototypes
//====================================================================================

typedef struct fastq_batch_reader_input {
  char *filename1;
  char *filename2;
  int flags;
  int batch_size;
  int gzip;
  list_t *list;

  // internal for fastq
  fastq_file_t *fq_file1;
  fastq_file_t *fq_file2;

  // internal for gzip fastq
  fastq_gzfile_t *fq_gzip_file1;
  fastq_gzfile_t *fq_gzip_file2;  

} fastq_batch_reader_input_t;


//------------------------------------------------------------------------------------

void fastq_batch_reader_input_init(char *filename1, char *filename2, int flags,
				   int batch_size, list_t *list, int gzip,
				   fastq_batch_reader_input_t* input);

//-----------------------------------------------------
// main function
//-----------------------------------------------------

void fastq_batch_reader(fastq_batch_reader_input_t *input);
void fastq_batch_reader_aligner(fastq_batch_reader_input_t* input);

//-----------------------------------------------------
//-----------------------------------------------------



#ifdef __cplusplus
}
#endif

#endif	/*  FASTQ_BATCH_READER_H  */
