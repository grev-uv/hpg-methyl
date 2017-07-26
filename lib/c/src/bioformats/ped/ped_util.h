#ifndef PED_UTIL_H
#define PED_UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <bioformats/family/family.h>
#include "ped_file.h"
enum Condition get_condition_from_phenotype(char* phenotype, ped_file_t *ped_file);

#ifdef __cplusplus
}
#endif

#endif
