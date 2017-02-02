#ifndef BED_FILE_STRUCTURE_H
#define BED_FILE_STRUCTURE_H

#include <sys/types.h>

#include <commons/file_utils.h>
#include <containers/linked_list.h>

/**
 * Entry in the BED document header.
 */
typedef struct bed_header_entry {
    char *text;
    int text_len;
} bed_header_entry_t;


/**
 * Entry in the BED document body.
 */
typedef struct bed_record {
    char *sequence;
    char *name;
    char *item_rgb;
    char *block_sizes;
    char *block_starts;
    unsigned long start;
    unsigned long end;
    unsigned long thickstart;
    unsigned long thickend;
    float score;
    int block_count;
    char strand;
    
    int sequence_len;           /**< Length of the seqname field */
    int name_len;               /**< Length of the name field */
    int item_rgb_len;           /**< Length of the item_rgb field */
    int block_sizes_len;        /**< Length of the block_sizes field */
    int block_starts_len;        /**< Length of the block_start field */
} bed_record_t;


/**
 * BED file structure. The physical file is defined by its file descriptor, its 
 * filename and the mode it has been open.
 *
 * It contains a header with several entries, and a body with several records.
 * For last version of the spec see: http://www.sanger.ac.uk/resources/software/bed/spec.html
 */
typedef struct bed_file {
    char* filename;
    char* mode;
    char *data;
    size_t data_len;
    
    linked_list_t *header_entries;
    linked_list_t *records;
} bed_file_t;

#endif
