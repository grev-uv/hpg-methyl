#ifndef BED_FILE_H
#define BED_FILE_H

#include <fcntl.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <commons/file_utils.h>
#include <commons/log.h>
#include <containers/linked_list.h>
#include <containers/cprops/linked_list.h>

#include "bed_file_structure.h"
#include "bed_reader.h"
#include "bed_write.h"

#define INIT_RECORD_SIZE    100


//====================================================================================
//  bed_file.h
//
//  bed_file_t structures and prototypes
//====================================================================================


//====================================================================================
//  Function prototypes
//====================================================================================

/*
 * Physical file management
 */

/**
 * Open a file stream in a given file mode (r/w/a) and initialize the associated 
 * bed_file_t structure.
 * 
 * @param filename The name of the file to open
 * @param mode Open mode (read/write/append)
 */
bed_file_t *bed_open(char *filename);

/**
 * Close the file stream associated to a bed_file_t type.
 * 
 * @param bed_file bed_file_t whose file stream is about to being closed
 * @param free_records flag that marks if the BED records must be also freed
 */
void bed_close(bed_file_t *bed_file, int free_records);

void bed_header_entry_free(bed_header_entry_t *bed_header_entry);

void bed_record_free(bed_record_t *bed_record);

/**
 * Fill the fields of the bed_file_t given as argument reading data from a file.
 * 
 * TODO breaks because a batches list isn't provided! 
 * a solution would be creating a temp batches list, then copying its contents to the bed_file_t
 * 
 * @param bed_file The bed_file_t whose fields will be set
 */
int bed_read(bed_file_t *bed_file);

int bed_read_batches(list_t *batches_list, size_t batch_size, bed_file_t *bed_file);

/**
 * Write the contents of the bed_file_t given as argument to the given path.
 * 
 * @param bed_file The bed_file_t whose information will be written
 * @param filename The path of the file to write the information to
 */
int bed_write(bed_file_t *bed_file, char *filename);

/*
 * Header and record entries management
 */

int add_bed_header_entry(bed_header_entry_t *header_entry, bed_file_t *bed_file);

int add_bed_record(bed_record_t* record, bed_file_t *bed_file);

#endif
