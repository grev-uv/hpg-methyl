#ifndef SW_SERVER_H
#define SW_SERVER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "commons/commons.h"
#include "commons/system_utils.h"
#include "containers/list.h"

#include "aligners/sw/smith_waterman.h"
#include "aligners/bwt/genome.h"
#include "bioformats/bam/bam_tags.h"

#include "timing.h"
#include "buffers.h"
#include "pair_server.h"
#include "breakpoint.h"
#include "cal_seeker.h"

//====================================================================================
//  Input structure for Smith-Waterman server
//====================================================================================

/**
 * @brief Input structure for the Smith-Waterman server function.
 *
 * This structure contains all required parameters by
 * the Smith-Waterman server function (@see sw_server).
 */
struct sw_server_input {
  float match;      /**< Penalty for match. */
  float mismatch;   /**< Penalty for mismatch. */
  float gap_open;   /**< Penalty for gap opening. */
  float gap_extend; /**< Penalty for gap extending. */
  float min_score;  /**< Minimum normalized score (0..1) for valid alignments. */
  unsigned int flank_length; /**< Length to extend the CAL region. */
  unsigned int write_size;   /**< Size of the writing batch (to disk). */
  
  sw_optarg_t sw_optarg;
  
  // for RNA 
  unsigned int max_intron_size; /**< Intron max size */
  int min_intron_size; /**< Intron max size */
  unsigned int seed_max_distance;
  int pair_mode;
  
  // to get inputs and to save outputs
  list_t* sw_list_p;    /**< Pointer to the list that contains the input sequences to align. */
  list_t* alignment_list_p; /**< Pointer to the list that contains the output aligned sequences. */
  genome_t* genome_p;   /**< Pointer to the genome structure to get the reference sequences. */
  genome_t* genome1_p;  /**< Pointer to the first genome structure to get the reference sequences for bisulfite. */
  genome_t* genome2_p;  /**< Pointer to the second genome structure to get the reference sequences for bisulfite. */
  genome_t* genome_or_p;  /**< Pointer to the original genome structure to get the reference sequences for bisulfite. */
  bwt_optarg_t* bwt_optarg_p;

  cal_optarg_t *cal_optarg_p;     /**< cal seeker configuration values */
  bwt_index_t *bwt_index_p;       /**< structure where were stored burrows wheeler transform index */
  metaexons_t *metaexons;
  linked_list_t *buffer;
  linked_list_t *buffer_hc;
  
  FILE *f_sa;
  FILE *f_hc;
  
  unsigned long long **valuesCT;
  unsigned long long **valuesGA;
};

//------------------------------------------------------------------------------------

typedef struct apply_sw_bs_stage_workspace {
  sw_optarg_t sw_optarg1;
  sw_optarg_t sw_optarg2;
  int sw_optarg_is_set;

  array_list_t *fill_gaps_sw_prepare_list1;
  array_list_t *fill_gaps_sw_prepare_list2;

  array_list_t *fill_end_gaps_sw_prepare_list1;
  array_list_t *fill_end_gaps_sw_prepare_list2;
} apply_sw_bs_stage_workspace_t;

//------------------------------------------------------------------------------------
void sw_optarg_init(float gap_open, float gap_extend, 
		    float match, float mismatch, sw_optarg_t *sw_optarg);

/**
 * @brief Initialization function for the @a sw_server_input_t structure.
 * @param sw_list_p pointer to the list that contains the input sequences to align
 * @param write_list_p pointer to the list that contains the output aligned sequences
 * @param write_size Size of the writing batch (to disk)
 * @param match Penalty for match
 * @param mismatch Penalty for mismatch
 * @param gap_open Penalty for gap opening
 * @param gap_extend Penalty for gap extending
 * @param min_score Minimum normalized score (0..1) for valid alignments
 * @param flank_length Length to extend the CAL region
 * @param genome_p pointer to the genome structure to get the reference sequences
 * @param[out] input_p pointer to the structure to initialize
 *
 * This function takes the input parameters and initializes the @a sw_server_input_t
 * structure that will be used by the Smith-Waterman server function (@see sw_server).
 */
void sw_server_input_init(list_t* sw_list_p, list_t* write_list_p, unsigned int write_size, 
			  float match, float mismatch, float gap_open, float gap_extend, 
			  float min_score, unsigned int flank_length, genome_t* genome_p,
			  size_t max_intron_size, int min_intron_size, size_t seed_max_distance, 
			  bwt_optarg_t* bwt_optarg_p,
			  cal_optarg_t *cal_optarg_p, bwt_index_t *bwt_index_p,
			  metaexons_t *metaexons, linked_list_t *buffer, 
			  linked_list_t *buffer_hc, FILE *f_sa, FILE *f_hc, 
			  int pair_mode, sw_server_input_t* input_p);

//====================================================================================
//  Smith-Waterman channel for SIMD implementation
//====================================================================================

/**
 * @brief Input structure for the Smith-Waterman server function.
 *
 * This structure contains all required parameters by
 * the Smith-Waterman server function (@see sw_server).
 * Smith-Waterman channel are re-usable by the different
 * reads to align, in order to reduce memory allocations.
 */
typedef struct sw_channel {
    size_t read_index;   /**< Index to the target read. */
    unsigned int cal_index;   /**< Index to the target CAL. */
    unsigned int header_len; /**< Read header length. */
    unsigned int read_len;  /**< Read length. */

    unsigned int allocated_ref_size;  /**< Allocated memory for the refence sequence. */
    unsigned int ref_len;            /**< Reference sequence length. */
    unsigned int start_splice;	    /**< Start splice site. */
    unsigned int end_splice;	   /**< End splice site */
    short int type; 		 /**< Type of CAL **/
    unsigned int extra_search;  /**<  Multiples intron marks found **/
    char* ref_p;                /**< Pointer to the reference sequence. */
} sw_channel_t;

//====================================================================================
//  Smith-Waterman server main function
//====================================================================================

/**
 * @brief Smith-Waterman server main function.
 * @param input_p pointer to the input structure containing the parameters required
 *    by the Smith-Waterman server (@see sw_server_input_t)
 *
 * Basically, this function removes, from an input list, batches containing
 * the sequences to align, then performs the Smith-Watermen algorithm based-on
 * SSE or AVX instructions (depending on the system architecture), and the
 * aligned sequences are packed into batches and insert these batches in an
 * output list for further writting to disk.
 */
void sw_server(sw_server_input_t* input_p);

//====================================================================================
// apply_sw
//====================================================================================

typedef struct sw_output {
  int strand;
  size_t chromosome;
  size_t ref_start;
  size_t ref_len;
  size_t mref_len;
  size_t mquery_start;
  size_t mref_start;
  float score;
  float norm_score;
  char* mquery;
  char* mref;
} sw_output_t;

sw_output_t *sw_output_new(int strand, size_t chrom, size_t ref_start, size_t ref_len,
			   size_t mref_len, size_t mquery_start, size_t mref_start,
			   float score, float norm_score, char* mquery, char* mref);

//--------------------------------------------------------------------------------------

int apply_sw_bs(sw_server_input_t* input, batch_t *batch, apply_sw_bs_stage_workspace_t *workspace);
void clean_apply_sw_stage_workspace(void* workspace);
void apply_sw_bs_4nt(sw_server_input_t* input, batch_t *batch, apply_sw_bs_stage_workspace_t *workspace);
void fill_matrix(subst_matrix_t subst_matrix, float match, float mismatch, int type, float factor_match, float factor_mismatch);

void fill_gaps_bs(mapping_batch_t *mapping_batch, sw_optarg_t *sw_optarg, 
		  genome_t *genome1, genome_t *genome2, int min_gap, int min_distance,
		  int bs_id, sw_optarg_t *sw_optarg1, sw_optarg_t *sw_optarg2,
      apply_sw_bs_stage_workspace_t *workspace);

void merge_seed_regions_bs(mapping_batch_t *mapping_batch, int bs_id);

void fill_end_gaps_bs(mapping_batch_t *mapping_batch, sw_optarg_t *sw_optarg, 
		      genome_t *genome1, genome_t *genome2, int min_H, int min_distance,
		      int bs_id, apply_sw_bs_stage_workspace_t *workspace);

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

#endif  // SW_SERVER_H
