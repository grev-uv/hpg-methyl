#ifndef GFF_FILE_STRUCTURE_H
#define GFF_FILE_STRUCTURE_H

#include <sys/types.h>

#include <commons/file_utils.h>
#include <containers/linked_list.h>

/**
 * Entry in the GFF document header.
 */
typedef struct gff_header_entry {
    char *text;
    int text_len;
} gff_header_entry_t;


/**
 * Entry in the GFF document body.
 */
typedef struct gff_record {
    char *sequence;
    char *source;
    char *feature;
    char *attribute;
    unsigned long start;
    unsigned long end;
    float score;
    char frame;
    char strand;
    
    int sequence_len;           /**< Length of the seqname field */
    int source_len;             /**< Length of the source field */
    int feature_len;            /**< Length of the feature field */
    int attribute_len;          /**< Length of the attributes field */
} gff_record_t;


/**
 * GFF file structure. The physical file is defined by its file descriptor, its 
 * filename and the mode it has been open.
 *
 * It contains a header with several entries, and a body with several records.
 * For last version of the spec see: http://www.sanger.ac.uk/resources/software/gff/spec.html
 */
typedef struct gff_file {
    char* filename;
    char* mode;
    char *data;
    size_t data_len;
    
    linked_list_t *header_entries;
    linked_list_t *records;
} gff_file_t;

#endif
