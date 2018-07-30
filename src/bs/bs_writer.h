#ifndef BS_WRITER_H
#define BS_WRITER_H

#include <stdio.h>

#include "commons/log.h"
#include "commons/file_utils.h"
#include "commons/workflow_scheduler.h"

#include "bioformats/fastq/fastq_batch_reader.h"
#include "bioformats/bam/bam_file.h"

#include "buffers.h"
#include "batch_writer.h"

int bs_writer(void *data);

//--------------------------------------------------------------------

void write_mapped_read(array_list_t *array_list, bam_file_t *bam_file);
void write_unmapped_read(fastq_read_t *fq_read, bam_file_t *bam_file);
void write_unmapped_read_paired_none(fastq_read_t *fq_read, bam_file_t *bam_file, int pair_num);

#endif // BS_WRITER_H
