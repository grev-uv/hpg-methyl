#ifndef PED_FILE_H
#define PED_FILE_H

#include <assert.h>

#include <fcntl.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <omp.h>

#include <commons/file_utils.h>
#include <commons/log.h>
#include <containers/cprops/hashtable.h>
#include <containers/cprops/linked_list.h>

#include <bioformats/family/family.h>

#include "ped_batch.h"
#include "ped_error.h"
#include "ped_file_structure.h"
#include "ped_reader.h"
#include "ped_util.h"
#include "ped_write.h"

#define INIT_RECORD_SIZE    100


//====================================================================================
//  ped_file.h
//
//  ped_file_t structures and prototypes
//====================================================================================


//====================================================================================
//  Function prototypes
//====================================================================================

/*
 * Physical file management
 */

/**
 * Open a file stream in a given file mode (r/w/a) and initialize the associated 
 * ped_file_t structure.
 * 
 * @param filename The name of the file to open
 * @param mode Open mode (read/write/append)
 */
ped_file_t *ped_open(char *filename);

/**
 * Close the file stream associated to a ped_file_t type.
 * 
 * @param ped_file ped_file_t whose file stream is about to being closed
 * @param free_families flag that marks if the families must be also freed
 */
void ped_close(ped_file_t *ped_file, int free_families);

void ped_record_free(ped_record_t *ped_record);

/**
 * Fill the fields of the ped_file_t given as argument reading data from a file.
 * 
 * TODO breaks because a batches list isn't provided! 
 * a solution would be creating a temp batches list, then copying its contents to the ped_file_t
 * 
 * @param ped_file The ped_file_t whose fields will be set
 */
int ped_read(ped_file_t *ped_file);

int ped_read_batches(list_t *batches_list, size_t batch_size, ped_file_t *ped_file);

/**
 * Write the contents of the ped_file_t given as argument to the given path.
 * 
 * @param ped_file The ped_file_t whose information will be written
 * @param filename The path of the file to write the information to
 */
int ped_write(ped_file_t *ped_file, char *filename);

/*
 * Header and record entries management
 */

int add_family(family_t *family, ped_file_t *ped_file);

int get_num_families(ped_file_t *ped_file);

int add_ped_record(ped_record_t* record, ped_file_t *ped_file);

#endif
