#ifndef PED_READ_H
#define PED_READ_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <bioformats/family/family.h>
#include <commons/log.h>

#include "ped_file_structure.h"

//====================================================================================
//  ped_read.h
//
//  ped reading functions prototypes
//====================================================================================


ped_record_t* create_ped_record();

void set_ped_record_family_id(char* family_id, ped_record_t* ped_record);

void set_ped_record_individual_id(char* individual_id, ped_record_t* ped_record);

void set_ped_record_father_id(char* father_id, ped_record_t* ped_record);

void set_ped_record_mother_id(char* mother_id, ped_record_t* ped_record);

void set_ped_record_sex(enum Sex sex, ped_record_t* ped_record);

void set_ped_record_phenotype(char* phenotype, ped_record_t* ped_record, ped_file_t* ped_file);

void set_ped_record_custom_field(char* field, ped_record_t* ped_record, ped_file_t *ped_file);

#ifdef __cplusplus
}
#endif

#endif
