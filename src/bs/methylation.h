#ifndef METHYLATION_H
#define METHYLATION_H

//====================================================================================

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "commons/commons.h"
#include "commons/system_utils.h"

#include "containers/array_list.h"
#include "containers/list.h"

#include "aligners/bwt/bwt.h"
#include "aligners/bwt/genome.h"
#include "aligners/sw/smith_waterman.h"

#include "bioformats/fastq/fastq_file.h"

#include "buffers.h"
#include "pair_server.h"

#include "hash_table.h"

#include "bioformats/bam/bam_tags.h"
#include "bwt_server.h"

//====================================================================================

// values for the diferent combination of nucleotides
#define ACGT 0
#define AGT  1
#define ACT  2
#define AT   3


typedef struct metil_data {
  char*  query_name;
  char   status;
  int    chromosome;
  size_t start;
  char   context;
  int    strand;
  int    zone;
} metil_data_t;


//====================================================================================

typedef struct metil_file {
  char* filenameCpG;                     /**< Metilation file name for CpG regions. */
  char* filenameCHG;                     /**< Metilation file name for CpG regions. */
  char* filenameCHH;                     /**< Metilation file name for CpG regions. */
  char* filenameMUT;                     /**< Metilation file name for Mutations.   */
  char* filenameSTAT;                    /**< Metilation file name for statistics.  */
  char* mode;                            /**< Open mode ("r", "w").                 */ // not in use
  genome_t *genome;                      /**< Reference to the original genome.     */

  FILE *CpG;                             /**< File pointer to CpG output.           */
  FILE *CHG;                             /**< File pointer to CHG output.           */
  FILE *CHH;                             /**< File pointer to CHH output.           */
  FILE *MUT;                             /**< File pointer to Mutation output.      */
  FILE *STAT;                            /**< File pointer to Statistics output.    */

  size_t skip_cpg;        /**< Skip flag for CpG regions. If set, the file is skipped for the remaining batches */
  size_t skip_chg;        /**< Skip flag for CHG regions. If set, the file is skipped for the remaining batches */
  size_t skip_chh;        /**< Skip flag for CHH regions. If set, the file is skipped for the remaining batches */
  size_t skip_mut;        /**< Skip flag for MUT regions. If set, the file is skipped for the remaining batches */

  hash_table_t *table_isles;             /**< Structure with the hash table values  */

  // check the values in order to change the type of the counter variables
  size_t CpG_methyl;                 /**< Global Counter for methylated Cytosines in CpG context   */
  size_t CpG_unmethyl;               /**< Global Counter for unmethylated Cytosines in CpG context */
  size_t CHG_methyl;                 /**< Global Counter for methylated Cytosines in CHG context   */
  size_t CHG_unmethyl;               /**< Global Counter for unmethylated Cytosines in CHG context */
  size_t CHH_methyl;                 /**< Global Counter for methylated Cytosines in CHH context   */
  size_t CHH_unmethyl;               /**< Global Counter for unmethylated Cytosines in CHH context */
  size_t MUT_methyl;                 /**< Global Counter for mutated Cytosines                     */
  size_t CUN_methyl;				 /**< Global Counter for Unknown context (CN or CHN) - methylated */ //RICARDO
  size_t CUN_unmethyl;				 /**< Global Counter for Unknown context (CN or CHN) - unmethylated */ //RICARDO
  size_t num_bases;                  /**< Global Counter for number of bases in the batch          */

  uint32_t *methyl_reads;            /**< Array with the number of methylated reads per chromosome */
} metil_file_t;

//====================================================================================

typedef struct methylation_stage_workspace {
  char* add_status_seq;
  char* add_status_gen;
  char* add_status_seq_dup;
} methylation_stage_workspace_t;

/**
 * @brief  Initialization of the metilation structure.
 * @param  metil_file metilation structure to use in the final analysis.
 * @param  dir        output directory.
 * @param  genome     pointer to the compressed original genome.
 * 
 * Initialize the variables in the structure, and open the output files.
 */
void metil_file_init(metil_file_t *metil_file, char *dir, genome_t *genome);

/**
 * @brief  Release the memory.
 * @param  metil_file metilation structure.
 * 
 * Release the memory used by the structure, and close the files in use.
 */
void metil_file_free(metil_file_t *metil_file);

//====================================================================================

/**
 * @brief  
 * @param  
 * @param  
 * @param  
 * 
 * 
 */
void metil_data_init(metil_data_t *metil_data, char *query, char status, int chromosome, size_t start, char context, int strand, int zone);

/**
 * @brief  
 * @param  
 * 
 * 
 */
void metil_data_free(metil_data_t *metil_data);

//====================================================================================

/**
 * @brief  Transfor the sequence to the metilation case
 * @param  refs sequence of reference
 * @param  len length of the sequence
 * @param  type type of transformation
 * 
 * Replace the bases of the sequence depending the type.
 */
void replace(char *refs, int len, int type);

//====================================================================================

/**
 * @brief  Obtain the complementary.
 * @param  c base of a sequence.
 * @return complementary base.
 * 
 * Obtain the complementary base of the input.
 */
char complement (char c);

//====================================================================================

/**
 * @brief  Obtain the reverse complementary of a dna sequence.
 * @param  orig sequence to transform.
 * @param  dest sequence transformed.
 * @param  len  length of the sequence.
 * 
 * Transform the @a orig sequence of DNA, to the reverse complementary.
 */
void rev_comp(char *query, char *seq, int len);

//====================================================================================

/**
 * @brief  Transform the sequence into the complementary.
 * @param  seq sequence of nucleotides.
 * @param  len length of the sequence.
 * 
 * Complement the sequence element by element.
 */
void comp(char *seq, int len);

//====================================================================================

/**
 * @brief  Make four copies of the array of reads, and transform for the bisulfite case.
 * @param  src original data array.
 * @param  dest_ct     destination to the ct transformation.
 * @param  dest_ct_rev destination to the ct rev comp transformation.
 * @param  dest_ga     destination to the ga transformation.
 * @param  dest_ga_rev destination to the ga rev comp transformation.
 * 
 * Copy the reads and transform it before insert in the new array.
 */
void cpy_transform_array_bs(array_list_t *src, array_list_t *dest_ct, array_list_t *dest_ct_rev,
			    array_list_t *dest_ga, array_list_t *dest_ga_rev);

//====================================================================================

/**
 * @brief  Set the original read in the sequences of each mapping.
 * @param  src1  Array with the mappings of the reads transformed.
 * @param  src1  Array with the mappings of the reads transformed.
 * @param  orig Fastq batch with the original reads.
 * 
 * Go on the array of mappings, and set the sequence to the original.
 */
void revert_mappings_seqs(array_list_t **src1, array_list_t **src2, array_list_t *orig);

//====================================================================================

/**
 * @brief  Stage to determine the status of each Cytosine in the sequence alignment                                       
 * @param  batch                                                                                                           
 *                                                                                                                         
 *                                                                                                                         
 */
int methylation_status_report(sw_server_input_t* input, batch_t *batch, methylation_stage_workspace_t *workspace);
void clean_methylation_stage_workspace(void *workspace);

//====================================================================================                                     

/**
 * @brief  Add the metilation status of each Cytosine to the write list                                                    
 * @param  array_list                                                                                                      
 * @param  bs_status                                                                                                       
 *                                                                                                                         
 *                                                                                                                         
 */
//void add_metilation_status(array_list_t *array_list, array_list_t *list, bs_context_t *bs_context, genome_t * genome, array_list_t * orig_seq, size_t index, int conversion);
void add_metilation_status(array_list_t *array_list, bs_context_t *bs_context, 
        genome_t * genome, array_list_t * orig_seq, size_t index, int conversion,
        methylation_stage_workspace_t *workspace);

//====================================================================================                                     

/**
 * @brief  
 * @param  
 * @param  
 *                                                                                                                         
 *                                                                                                                         
 */
void postprocess_bs(char *query_name, char status, size_t chromosome, size_t start, char context, int strand, int region,
		    array_list_t *list);


//====================================================================================                                     

/**
 * @brief  Write the metilation status of each Cytosine
 * @param  array_list 
 * @param  metil_file 
 * 
 * 
 */
void write_bs_context(metil_file_t *metil_file, bs_context_t *bs_context, size_t num_chromosomes);

//====================================================================================

int encode_context(char* filename, char* directory, int max_num_chromosomes);

//====================================================================================

int load_encode_context(char* directory, unsigned long long **valuesCT, unsigned long long **valuesGA);

//====================================================================================

#endif // METHYLATION_H
