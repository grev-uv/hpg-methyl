#ifndef PED_FILE_H
#define PED_FILE_H

#ifdef __cplusplus
extern "C" {
#endif

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
#include <cprops/hashtable.h>
#include <cprops/linked_list.h>

#include <bioformats/family/family.h>

#include "ped_batch.h"
#include "ped_error.h"
#include "ped_file_structure.h"
#include "ped_reader.h"
#include "ped_util.h"
#include "ped_write.h"


/* **************************************
 *      Physical file management        *
 * **************************************/

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
 * @param free_phenotypes flag that marks if the phenotypes must be also freed
 */
void ped_close(ped_file_t *ped_file, int free_families, int free_phenotypes);

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



/* **************************************
 * Header and record entries management *
 * **************************************/

int add_family(family_t *family, ped_file_t *ped_file);

int get_num_families(ped_file_t *ped_file);

int ped_add_individual(individual_t *individual, ped_file_t *ped_file);

int add_ped_record(ped_record_t* record, ped_file_t *ped_file);

/**
 * Given a collection of families read from a PED file and possibly with multiple generations,
 * returns them transformed into single-generation combinations: father-mother-children 
 * (like trios, quartets, and so on).
 * 
 * @param ped_file The file the families were read from
 * @param num_families[inout] The number of families after reorganizing
 * @return The families flattened into single generation ones
 */
family_t **ped_flatten_families(ped_file_t *ped_file, int *num_families);


/* **************************************
 *        Phenotypes management         *
 * **************************************/

khash_t(str)* get_phenotypes(ped_file_t *ped_file);

/**
 * Select the affected/unaffected ID for the samples in the selected phenotype field
 * Only affects to the "condition" in individual struct
 * */
void set_unaffected_phenotype(const char* id, ped_file_t *ped_file);

void set_affected_phenotype(const char* id, ped_file_t *ped_file);

/**
 * Changes the variable field. By default, "PHENO" or field number 6
 * */
void set_variable_field(const char* id, int num_field, ped_file_t *ped_file);

int set_phenotype_group(char** ids, int n , ped_file_t *ped_file);

khash_t(str)* get_phenotypes(ped_file_t *ped_file);

int get_num_variables(ped_file_t* ped_file);

#ifdef __cplusplus
}
#endif

#endif
