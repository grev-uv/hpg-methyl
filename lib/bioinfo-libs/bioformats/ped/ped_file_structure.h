#ifndef PED_FILE_STRUCTURE_H
#define PED_FILE_STRUCTURE_H

#include <sys/types.h>


#include <commons/file_utils.h>

#include <containers/cprops/hashtable.h>

#include <bioformats/family/family.h>

/**
 * Entry in the PED document body.
 */
typedef struct ped_record {
    char *family_id;
    char *individual_id;
    char *father_id;
    char *mother_id;
    enum Sex sex;
    float phenotype;
} ped_record_t;

/**
 * PED file structure. The physical file is defined by its file descriptor, its 
 * filename and the mode it has been open.
 *
 * It contains a collection of families with a unique ID, which are specified via 
 * a group of records which contain their members.
 */
typedef struct ped_file {
    char* filename;
    char* mode;
    char *data;
    size_t data_len;
    
    cp_hashtable *families;
} ped_file_t;

#endif
