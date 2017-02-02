#ifndef ALIGNMENTS_H
#define ALIGNMENTS_H

#include <stdio.h>

#include "commons/commons.h"

#include "samtools/bam.h"

//#include "bam_commons.h"

#define MIN_ALLOCATED_SIZE_FOR_CIGAR_STRING  	5
#define HUMAN      				1 //Specie: Human
#define NCBI37      				1 //Assembly: NCBI 3.7

//========CIGAR AUTOMATA STATUS========//
#define CIGAR_MATCH_MISMATCH  	0
#define CIGAR_INSERTION  	1
#define CIGAR_DELETION  	2
#define CIGAR_BIG_DELETION    	3
#define CIGAR_SKIPPED  		4
#define CIGAR_PADDING 		5
#define CIGAR_PERFECT_MATCH   	6
//=====================================//


/* **************************************

OA *      	Structures    		*
 * *************************************/

/**
* @brief Structure for storing alignments
*
* Structure for storing alignments
*/
typedef struct alignment {
    short int optional_fields_length; 	/**< Length of optional fields. */
    short int chromosome;		/**< Chromosome. */
    int position;			/**< Mapping position. */
    int mate_position;			/**< Mapping position of the mate. */
    short int mate_chromosome;		/**< Chromosome of the mate. */
    short int template_length;		/**< Template length. */
    int map_quality;			/**< Map quality. */
    int num_cigar_operations;		/**< Number of CIGAR operations. */

    //flags
    uint8_t is_paired_end;		/**< 0: single end, 1: paired end. */ 
    uint8_t is_paired_end_mapped;	/**< 0: pair not mapped, 1: pair mapped. */
    uint8_t is_seq_mapped;		/**< 0: seq not mapped, 1: seq mapped. */ 
    uint8_t is_mate_mapped;		/**< 0: mate unmapped, 1, mate mapped. */ 
    uint8_t seq_strand;			/**< 0 (forward) or 1 (reverse). */
    uint8_t mate_strand;		/**< Strand of the mate read, same values. */
    uint8_t pair_num;			/**< Paired-end number 1 or 2. */
    uint8_t secondary_alignment;		/**< 0: not primary, 1 primary. */
    uint8_t fails_quality_check;	/**< 0: meet quality checks, 1: quality checks not meeted. */ 
    uint8_t pc_optical_duplicate;	/**< 0: not duplicate, 1: duplicate. */

    char* query_name;			/**< Query template name. */
    char* sequence;			/**< Sequence of nts. */
    char* quality;			/**< Quality of nts. */
    char* cigar;			/**< CIGAR string. */
    uint8_t* optional_fields;		/**< Optional fields. */
} alignment_t;

/**
* @brief Batch of alignments
*
* Batch of alignments
*/
typedef struct alignment_batch {
    int type; 				/**< SINGLE or MULTIPLE CHROMOSOMES (alignments). */
    int num_alignments;			/**< Number of alignments in the batch. */

    int allocated_alignments;		/**< Number of allocated alignments for memory reserve. */
    alignment_t** alignments_p;		/**< Pointers to alignments. */
} alignment_batch_t;

/* **************************************
 *      	Functions    		*
 * *************************************/

/**
*  @brief Creates a new alignment
*  @return pointer to the created alignment
*  
*  Creates and returns a new alignment
*/
alignment_t* alignment_new();

/**
*  @brief Inits an alignment with a single end mapping
*  @param query_name name of the mapped read
*  @param sequence sequence string
*  @param quality quality string
*  @param strand strand of the alignment (1: forward, 0: reverse)
*  @param chromosome chromosome of the alignment
*  @param position start position of the alignment
*  @param cigar cigar (string format)
*  @param num_cigar_operations number of cigar operations
*  @param map_quality quality of the mapping 
*  @param is_seq_mapped flag indicating if sequence is mapped
*  @param secondary_alignment flag indicating if alignment is primary
*  @param[in,out] alignment_p pointer to the alignment to init
*  @return pointer to the created qc hash list item
*  
*  Creates and returns a new qc hash list item
*/
void alignment_init_single_end(char* query_name, char* sequence, char* quality, short int strand, short int chromosome, int position, char* cigar, short int num_cigar_operations, int map_quality, short int is_seq_mapped, short int secondary_alignment, int optional_fields_length, char *optional_fields, alignment_t* alignment_p);

/**
*  @brief Inits an alignment with a single end mapping
*  @param query_name name of the mapped read
*  @param sequence1 sequence string for paired end 1
*  @param sequence2 sequence string for paired end 2
*  @param quality1 quality string for paired end 1
*  @param quality2 quality string for paired end 2
*  @param strand1 strand of the alignment (1: forward, 0: reverse) for paired end 1
*  @param strand2 strand of the alignment (1: forward, 0: reverse) for paired end 2
*  @param chromosome1 chromosome of the alignment for paired end 1
*  @param chromosome2 chromosome of the alignment for paired end 2
*  @param position1 start position of the alignment for paired end 1
*  @param position2 start position of the alignment for paired end 2
*  @param cigar1 cigar (string format) for paired end 1
*  @param cigar2 cigar (string format) for paired end 2
*  @param num_cigar_operations1 number of cigar operations for paired end 1
*  @param num_cigar_operations2 number of cigar operations for paired end 2
*  @param map_quality1 quality of the mapping for paired end 1
*  @param map_quality2 quality of the mapping for paired end 2
*  @param secondary_alignment1 flag indicating if alignment is primary for paired end 1
*  @param secondary_alignment1 flag indicating if alignment is primary for paired end 2
*  @param[in,out] alignment1_p pointer to the alignment to init for paired end 1
*  @param[in,out] alignment2_p pointer to the alignment to init for paired end 2
*  @return pointer to the created qc hash list item
*  
*  Creates and returns a new qc hash list item
*/
void alignment_init_paired_end(char* query_name, char* sequence1, char* sequence2, char* quality1, char* quality2, short int strand1, short int strand2, short int chromosome1, int position1, int position2, short int chromosome2, char* cigar1, char* cigar2, short int num_cigar_operations1, short int num_cigar_operations2, short int map_quality1, short int map_quality2, short int secondary_alignment1, short int secondary_alignment2,  alignment_t* alignment1_p, alignment_t* alignment2_p);

void alignment_update_paired_end(alignment_t* alignment1_p, alignment_t* alignment2_p);

/**
*  @brief Creates an alignment from a bam1_t structure
*  @param bam_p pointer to the bam1_t structure
*  @param base_quality base quality for quality values normalization
*  @return pointer to the created alignment
*  
*  Creates and returns an alignment from a bam1_t structure
*/
alignment_t* alignment_new_by_bam(bam1_t* bam_p, int base_quality);

/**
*  @brief Frees a given alignment
*  @param[in,out] alignment_p pointer to the alignment to free 
*  @return void
*  
*  Frees a given alignment
*/
void alignment_free(alignment_t* alignment_p);

/**
*  @brief Converts an alignment into a bam1_t structure
*  @param alignment_p pointer to the alignment
*  @param base_quality base quality for quality values normalization
*  @return pointer to the created bam1_t
*  
*  Creates and returns an alignment from a bam1_t structure
*/
bam1_t* convert_to_bam(alignment_t* alignment_p, int base_quality);

/**
*  @brief Prints the content of an alignment
*  @param alignment_p pointer to the alignment to print
*  @return void
*  
*  Prints the content of an alignment
*/
void alignment_print(alignment_t* aligment_p);

/**
*  @brief Prints the content of a bam1_t structure
*  @param bam_p pointer to the bam1_t to print
*  @param base_quality base quality for quality values normalization
*  @return void
*  
*  Prints the content of a bam1_t structure
*/
void bam_print(bam1_t* bam_p, int base_quality);

/**
*  @brief Creates a bam_header_t
*  @param specie specie
*  @param assembly number of version of the reference genome
*  @return pointer to the created bam_header_t
*  
*  Creates a bam_header_t for the specified specie and assembly
*  Headers are stored in /bam_headers 
*/
bam_header_t* bam_header_new(int specie, int assembly, char* file_path);

void bam_header_free(bam_header_t *header);

short int is_secondary_alignment(alignment_t *alignment);

void set_secondary_alignment(short int set, alignment_t *alignment);

/* **********************************************************************
 *      	Functions to manage bam1_t coded fields    		*
 * *********************************************************************/

/**
*  @brief Converts a uint8_t coded sequence to string
*  @param sequence_p sequence coded in uint8_t type
*  @param sequence_length length of the sequence
*  @return pointer to the char sequence
*  
*  Converts a uint8_t sequence to string
*/
char* convert_to_sequence_string(uint8_t* sequence_p, int sequence_length);

/**
*  @brief Converts a uint8_t coded quality to string
*  @param quality_p quality coded in uint8_t type
*  @param quality_length length of the quality
*  @param base_quality base quality for quality values normalization
*  @return pointer to the char quality
*  
*  Converts a uint8_t quality to string
*/
char* convert_to_quality_string_length(char* quality_desp, uint8_t* quality_src, int quality_length, int base_quality);

//char* convert_from_uint8_to_string(uint8_t* sequence_p, int sequence_length);

/**
*  @brief Converts a uint8_t coded quality to string
*  @param quality_p quality coded in uint8_t type
*  @return pointer to the char quality
*  
*  Converts a uint8_t quality to string
*/
char* convert_to_quality_string(uint8_t* quality_p);

/**
*  @brief Converts a uint32_t coded cigar to string
*  @param cigar_p cigar coded in uint32_t type
*  @param num_cigar_operations number of operations in the cigar
*  @return pointer to the char cigar
*  
*  Converts a uint32_t coded cigar to string
*/
char* convert_to_cigar_string(uint32_t* cigar_p, int num_cigar_operations);

/**
*  @brief Converts a cigar string to uint32_t type
*  @param data[in,out] pointer to the cigar in uint32_t type
*  @param cigar cigar in string format
*  @param num_cigar_operations number of operations in the cigar
*  @return void
*  
*  Converts a cigar string to uint8_t type
*/
void convert_to_cigar_uint32_t(uint8_t* data, char* cigar, int num_cigar_operations);

/**
*  @brief Converts a sequence string to uint8_t type
*  @param data[in,out] pointer to the sequence in uint8_t type
*  @param sequence_p sequence in string format
*  @param sequence_length length of the sequence
*  @return void
*  
*  Converts a string sequence to uint8_t coded sequence
*/
void convert_to_sequence_uint8_t(uint8_t* data, char* sequence_p, int sequence_length);

/**
*  @brief Converts a quality string to uint8_t type
*  @param data[in,out] pointer to the quality in uint8_t type
*  @param quality_p quality in string format
*  @param quality_length length of the quality
*  @param base_quality base quality for quality values normalization
*  @return void
*  
*  Converts a string quality to uint8_t coded quality
*/
void convert_to_quality_uint8_t(uint8_t* data, char* quality_p, int quality_length, int base_quality);


/* **********************************************************************
 *      	Functions to manage CIGAR                		*
 * *********************************************************************/

char select_op(unsigned char status);

/**
 */
char* generate_cigar_str(char *str_seq_p, char *str_ref_p, unsigned int start_seq, 
			 unsigned int seq_orig_len, unsigned int length, 
			 int *distance, int *number_op_tot);

#endif /* ALIGNMENTS_H */

