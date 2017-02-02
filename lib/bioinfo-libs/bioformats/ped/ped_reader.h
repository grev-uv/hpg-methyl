#ifndef PED_RAGEL_H
#define PED_RAGEL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <containers/list.h>
#include <containers/cprops/linked_list.h>

#include "ped_error.h"
#include "ped_file_structure.h"
#include "ped_file.h"
#include "ped_read.h"
#include "ped_batch.h"

enum PED_Field { FAMILY_ID, INDIVIDUAL_ID, FATHER_ID, MOTHER_ID, SEX, PHENOTYPE, OTHER };

int ped_ragel_read(list_t *batches_list, size_t batch_size, ped_file_t *file);

#endif
