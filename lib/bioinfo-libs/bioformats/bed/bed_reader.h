#ifndef BED_RAGEL_H
#define BED_RAGEL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <containers/list.h>
#include <containers/cprops/linked_list.h>

#include "bed_file_structure.h"
#include "bed_file.h"
#include "bed_read.h"
#include "bed_batch.h"

int bed_ragel_read(list_t *batches_list, size_t batch_size, bed_file_t *file);

#endif
