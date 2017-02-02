#ifndef GFF_FILE_H
#define GFF_FILE_H

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

#include <containers/cprops/linked_list.h>

#include "gff_file_structure.h"
#include "gff_reader.h"
#include "gff_write.h"

#define INIT_RECORD_SIZE    100


//====================================================================================
//  gff_file.h
//
//  gff_file_t structures and prototypes
//====================================================================================


//====================================================================================
//  Function prototypes
//====================================================================================

/*
 * Physical file management
 */

/**
 * Open a file stream in a given file mode (r/w/a) and initialize the associated 
 * gff_file_t structure.
 * 
 * @param filename The name of the file to open
 * @param mode Open mode (read/write/append)
 */
gff_file_t *gff_open(char *filename);

/**
 * Close the file stream associated to a gff_file_t type.
 * 
 * @param gff_file gff_file_t whose file stream is about to being closed
 * @param free_records flag that marks if the GFF records must be also freed
 */
void gff_close(gff_file_t *gff_file, int free_records);

void gff_header_entry_free(gff_header_entry_t *gff_header_entry);

void gff_record_free(gff_record_t *gff_record);

/**
 * Fill the fields of the gff_file_t given as argument reading data from a file.
 * 
 * TODO breaks because a batches list isn't provided! 
 * a solution would be creating a temp batches list, then copying its contents to the gff_file_t
 * 
 * @param gff_file The gff_file_t whose fields will be set
 */
int gff_read(gff_file_t *gff_file);

int gff_read_batches(list_t *batches_list, size_t batch_size, gff_file_t *gff_file);

/**
 * Write the contents of the gff_file_t given as argument to the given path.
 * 
 * @param gff_file The gff_file_t whose information will be written
 * @param filename The path of the file to write the information to
 */
int gff_write(gff_file_t *gff_file, char *filename);

/*
 * Header and record entries management
 */

int add_gff_header_entry(gff_header_entry_t *header_entry, gff_file_t *gff_file);

int add_gff_record(gff_record_t* record, gff_file_t *gff_file);

#endif
