#ifndef BWT_SERVER_H
#define BWT_SERVER_H

#include "commons/commons.h"
#include "containers/array_list.h"

#include "buffers.h"
#include "pair_server.h"
#include "bs/methylation.h"
#include "histogram.h"

//====================================================================================
//  structures and prototypes
//====================================================================================

/**
 * @brief Structure for store all parameters needed for run @a bwt_server_cpu.
 * 
 * Structure for store some configuration values and data structures like lists 
 * to insert and read data batches. 
 */
struct bwt_server_input {
  unsigned int batch_size;        /**< size of data batches*/
  unsigned int write_size;        /**< size of write batches*/
  bwt_optarg_t *bwt_optarg_p;     /**< burrows wheeler transform configuration values */
  list_t* write_list_p;           /**< list for store write batches*/
  list_t* read_list_p;            /**< list for read batches of reads */
  list_t* unmapped_read_list_p;   /**< list for store batches with the reads no mapped */
  bwt_index_t *bwt_index_p;       /**< structure where were stored burrows wheeler transform index */
  metaexons_t *metaexons;
  bwt_index_t *bwt_index2_p;      /**< structure where were stored burrows wheeler transform index */
  report_optarg_t *report_optarg;
  sw_optarg_t *sw_optarg;
  genome_t *genome;
};

//------------------------------------------------------------------------------------

/**
 * @brief  Initializer for the @a bwt_server_input_t structure.
 * @param  read_list_p list for read batches of reads
 * @param  batch_size size of data batches
 * @param  bwt_optarg_p burrows wheeler transform configuration values
 * @param  bwt_index_p structure where were stored burrows wheeler transform index
 * @param  write_list_p list for store write batches
 * @param  write_size size of write batches
 * @param  unmapped_read_list_p list for store batches with the reads no mapped
 * 
 * Initialize all @a bwt_server_input_t fields with the input parameters.
 */
void bwt_server_input_init(list_t* read_list_p, unsigned int batch_size, bwt_optarg_t *bwt_optarg_p, 
			   bwt_index_t *bwt_index_p, list_t* write_list_p, unsigned int write_size, 
			   list_t* unmapped_read_list_p, metaexons_t *metaexons, sw_optarg_t *sw_optarg,
			   genome_t *genome, bwt_server_input_t* input_p);

//====================================================================================
// apply_bwt
//====================================================================================

int apply_bwt_bs(bwt_server_input_t* input, batch_t *batch_p, bwt_stage_bs_workspace_t *workspace);
void clean_bwt_stage_bs_workspace(void *workspace);

#endif  // BW_SERVER_H
