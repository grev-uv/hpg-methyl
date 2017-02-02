#ifndef PED_READ_H
#define PED_READ_H

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

void set_ped_record_phenotype(float phenotype, ped_record_t* ped_record);

#endif
