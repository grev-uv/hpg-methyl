#ifndef CAL_SEEKER_H
#define CAL_SEEKER_H

#include "buffers.h"
#include "aligners/bwt/bwt.h"
#include "pair_server.h"

#define MAX_CALS 200
#define MAX_RNA_CALS 200

extern size_t reads_max_cals;
extern size_t reads_no_cals;
//====================================================================================
//  structures and prototypes
//====================================================================================

/**
 * @brief Structure for store all parameters needed for run @a cal_seeker_server.
 * 
 * Structure for store some configuration values and data structures like lists 
 * to insert and read data batches.
 */
struct cal_seeker_input {
  unsigned batch_size;           /**< size of data batches*/
  list_t *write_list;            /**< list for store write batches*/
  list_t *regions_list;          /**< list for read batches with all regions found for each read */
  list_t *sw_list;               /**< list to store batches with all CALs found for each read */
  list_t *pair_list;
  cal_optarg_t *cal_optarg;      /**< cal seeker configuration values */
  genome_t *genome;
  bwt_optarg_t *bwt_optarg;
  bwt_index_t *index;
  bwt_index_t *index2;
  metaexons_t *metaexons;
  genome_t *genome2;             /**< second genome for bisulfite */
};

typedef struct caling_bs_stage_workspace {
} caling_bs_stage_workspace_t;

/**
 * @brief  Initializer for the @a cal_seeker_input_t structure.
 * @param  region_list_p list for read batches with all regions found for each read
 * @param  cal_optarg_p cal seeker configuration values 
 * @param  write_list_p list for store write batches
 * @param  write_size size of write batches
 * @param  sw_list_p list to store batches with all CALs found for each read
 * @param  pair_list_p list to store batches with pairs
 *
 * 
 * Initialize all @a cal_seeker_input_t fields with the input parameters.
 */
void cal_seeker_input_init(list_t *regions_list, cal_optarg_t *cal_optarg, 
			   list_t* write_list, unsigned int write_size, 
			   list_t *sw_list, list_t *pair_list,
			   genome_t *genome, bwt_optarg_t *bwt_optarg, 
			   bwt_index_t *index, metaexons_t *metaexons, 
			   cal_seeker_input_t *input);

//--------------------------------------------------------------------------------------

/**
 * @brief CAL Seeker server make CALs from regions
 * @param input_p all configuration values and data structures needed by the server
 *
 * CAL Seeker server extract batches of reads with their regions from @a regions_list_p and
 * for each read call @a bwt_generate_cal_list for generate CALs. Finally, all reads with a lot of
 * CALs are stored in @a write_list_p because those reads aren't mapped and will be store in bam file,
 * the rest of reads are stored with their CALs in @a sw_list_p.      
 */
//void cal_seeker_server(cal_seeker_input_t* input);

//====================================================================================
// defines
//====================================================================================

#define NONE_POS   0
#define BEGIN_POS  1
#define END_POS    2

#define SINGLE_FLANK 0 //2
#define DOUBLE_FLANK 0 //4

//====================================================================================
// structures
//====================================================================================

typedef struct sw_prepare {
  int left_flank;
  int right_flank;  
  int ref_type;
  char *query;
  char *ref;
  seed_region_t *seed_region;
  cal_t *cal;
  fastq_read_t *read;
} sw_prepare_t;

//--------------------------------------------------------------------------------------
// In case the build is being compiled with all optimizations disabled (ie: debug build)
// inline functions must be forward declared to enable the compiler to find the appropriate
// simbols to link into them on all the translation units.
// - Date: 14 / 11 / 2016
// - Who: Cesar
#ifdef __GNUC__
#ifdef __NO_INLINE__
static inline sw_prepare_t *sw_prepare_new(char *query, char *ref, int left_flank, int right_flank, int ref_type) {
#else
inline sw_prepare_t *sw_prepare_new(char *query, char *ref, int left_flank, int right_flank, int ref_type) {
#endif /* __NO_INLINE__ */
#else
inline sw_prepare_t *sw_prepare_new(char *query, char *ref, int left_flank, int right_flank, int ref_type) {
#endif /* __GNUC__ */
  sw_prepare_t *p = (sw_prepare_t *) malloc(sizeof(sw_prepare_t));
  p->query = query;
  p->ref = ref;
  p->left_flank = left_flank;
  p->right_flank = right_flank;
  p->seed_region = NULL;
  p->cal = NULL;
  p->read = NULL;
  p->ref_type = ref_type;

  return p;
}

#ifdef __GNUC__
#ifdef __NO_INLINE__
static inline void sw_prepare_free(sw_prepare_t *p) {
#else
static inline void sw_prepare_free(sw_prepare_t *p) {
#endif /* __NO_INLINE__ */
#else
static inline void sw_prepare_free(sw_prepare_t *p) {
#endif /* __GNUC__ */
	if (p) {
    if (p->query) free(p->query);
    if (p->ref) free(p->ref);
    free(p);
  }
}

//====================================================================================
// apply_caling
//====================================================================================

int apply_caling_bs(cal_seeker_input_t* input, batch_t *batch, caling_bs_stage_workspace_t *workspace);
void clean_apply_caling_bs_stage_workspace(void* workspace);

//--------------------------------------------------------------------------------------

#endif  // CAL_SEEKER_H
