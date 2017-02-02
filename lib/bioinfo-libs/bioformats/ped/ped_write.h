#ifndef PED_WRITE_H
#define PED_WRITE_H

#include <stdio.h>
#include <string.h>

#include "ped_file_structure.h"
#include "ped_file.h"
#include "ped_batch.h"

//====================================================================================
//  ped_write.h
//
//  ped writing functions prototypes
//====================================================================================

int ped_write_to_file(ped_file_t *ped_file, FILE *fd);


/* **************** Family in PED functions **********************/

void write_ped_individual(individual_t *individual, FILE *fd);


/* ************* Record management functions *********************/

void write_ped_batch(ped_batch_t *ped_batch, FILE *fd);

void write_ped_record(ped_record_t* ped_record, FILE *fd);


#endif 
