#ifndef GFF_RAGEL_H
#define GFF_RAGEL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <containers/list.h>
#include <cprops/linked_list.h>

#include "gff_file_structure.h"
#include "gff_file.h"
#include "gff_read.h"
#include "gff_batch.h"

int gff_ragel_read(list_t *batches_list, size_t batch_size, gff_file_t *file);

#ifdef __cplusplus
}
#endif

#endif
