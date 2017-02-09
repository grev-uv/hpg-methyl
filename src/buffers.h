#ifndef BUFFERS_H
#define BUFFERS_H

#include <stdbool.h>

#include "containers/array_list.h"
#include "containers/list.h"
#include "aligners/bwt/bwt.h"
#include "bioformats/bam/alignment.h"
#include "aligners/bwt/genome.h"
#include "timing.h"
#include "statistics.h"
#include "commons/log.h"

#include "breakpoint.h"


//-----------------------------------------------

#define ALIGNMENT_TYPE       0
#define CAL_TYPE             1
#define META_ALIGNMENT_TYPE  2

//====================================================================================
//  Buffer RNA Vars
//====================================================================================

#define BITEM_NO_CALS          0
#define BITEM_SINGLE_ANCHORS   1
#define BITEM_CALS             2
#define BITEM_META_ALIGNMENTS  3

//====================================================================================


//====================================================================================
//  Workflow Vars
//====================================================================================

//================================= NORMAL WORKFLOW ==================================

//-------- DEFINE WORKFLOW COMMON VARS -----------

#define BWT_STAGE               0
#define CONSUMER_STAGE         -1


//--------  DEFINE WORKFLOW BS VARS  -----------

#define BS_HIST_MARGIN          0.80

#define BS_BWT_STAGE            0
#define BS_CAL_STAGE            1
#define BS_PRE_PAIR_STAGE       2
#define BS_SW_STAGE             3
#define BS_POST_PAIR_STAGE      4
#define BS_STATUS_STAGE         5

//--------  DEFINE WORKFLOW BS UNIFIED VARS  -----------

#define BS_UN_BWT_STAGE         0
#define BS_UN_PRE_PAIR_STAGE    1
#define BS_UN_SW_STAGE          2
#define BS_UN_POST_PAIR_STAGE   3
#define BS_UN_STATUS_STAGE      4

//--------  DEFINE WORKFLOW DNA VARS  -----------

#define CAL_STAGE               1
#define PRE_PAIR_STAGE          2
#define SW_STAGE                3
#define DNA_POST_PAIR_STAGE     4

//--------  DEFINE WORKFLOW RNA VARS   -----------

#define RNA_CAL_STAGE           1
#define RNA_STAGE               2
#define RNA_POST_PAIR_STAGE     3

//================================= LAST WORKFLOW ==================================

#define LAST_RNA_POST_PAIR_STAGE     1

//====================================================================================
//  globals
//====================================================================================

#define NOT_ANCHORS         0
#define SINGLE_ANCHORS      1
#define DOUBLE_ANCHORS      2
#define ALIGNMENTS_FOUND    3
#define ALIGNMENTS_EXCEEDED 4
#define MULTIPLE_ANCHORS    5

//------------------------------------------------------------------------------------

#define DNA_MODE           0
#define RNA_MODE           1
// added by PP
#define BS_MODE            2
#define BS_UN_MODE         3

//------------------------------------------------------------------------------------

#define DNA_FLAG           1
#define RNA_FLAG           2
#define SINGLE_END_FLAG    4
#define PAIRED_END_FLAG    8
#define MATE_PAIR_FLAG    16
#define PAIR1_FLAG        32
#define PAIR2_FLAG        64
#define WRITE_ITEM_FLAG  128
#define SW_ITEM_FLAG     256
// added by PP
#define BS_FLAG          512

//------------------------------------------------------------------------------------

#define SINGLE_END_MODE 0
#define PAIRED_END_MODE 1
#define MATE_PAIR_MODE  2

//------------------------------------------------------------------------------------

#define PAIR_1   1
#define PAIR_2   2
#define PAIR_1_2 3

//------------------------------------------------------------------------------------

#define UNKNOWN_ITEM    0
#define READ_ITEM       1
#define KL_ITEM         2
#define SEED_ITEM       3
#define SPLIT_KL_ITEM   4
#define SW_ITEM         5
#define WRITE_ITEM      6

//------------------------------------------------------------------------------------
// Define BWT matching flags

#define BWT_BS_MATCHING_EXACT     0
#define BWT_BS_MATCHING_SEED      1
#define BWT_BS_MATCHING_EXCEEDED  2

//------------------------------------------------------------------------------------

#define MISMATCH_FLAG        0
#define MATCH_FLAG           1
#define SPLICE_EXACT_FLAG    2
#define SPLICE_EXTEND_FLAG   3

//------------------------------------------------------------------------------------

#define NORMAL_MODE 0
#define SEED_MODE   1

//------------------------------------------------------------------------------------

#define FASTQ_READER         	0
#define BWT_SERVER	   	      1
#define BWT_SERVER_HISTOGRAM  2
#define CAL_SEEKER	   	      3
#define PRE_PAIR_TIME		      4
#define SW_STAGE_TIME         5
#define POST_PAIR_TIME        6
#define METHYLATION_REP_TIME  7
#define BAM_WRITER         	  8
#define BS_ALIGNMENT_TOTAL    9
#define TOTAL_TIME         	  10

//------------------------------------------------------------------------------------

#define MAX_READ_MAPPINGS 100 

//------------------------------------------------------------------------------------

#define UNKNWON_ACTION 0
#define BWT_ACTION     1
#define SEEDING_ACTION 2
#define CAL_ACTION     3
#define PAIR_ACTION    4
#define SW_ACTION      5

#define NUM_STRANDS 2

//====================================================================================
//  SPLICE JUNCTION TYPE
//====================================================================================

//--------------------------------------------//
//            Not found splice junction       //
//--------------------------------------------//

#define NOT_SPLICE	-1

//--------------------------------------------//
//      No Cannonical Splice junction         //
//--------------------------------------------//

#define UNKNOWN_SPLICE	0
                                            
//--------------------------------------------//
//        Cannonical Splice Junction          //
//--------------------------------------------//

#define GT_AG_SPLICE  	1 //+
#define CT_AC_SPLICE  	2 //-
  
//--------------------------------------------//
//      Semi-Cannonical Splice Junction       //
//--------------------------------------------//

#define AT_AC_SPLICE  	3 //+
#define GT_AT_SPLICE  	4 //-
#define GC_AG_SPLICE  	5 //+
#define CT_GC_SPLICE  	6 //-

//===============================================

//====================================================================================
//  structures and prototypes
//====================================================================================

bam_header_t *create_bam_header_by_genome(genome_t *genome);


/**
 * @brief Structure for store in files all the process data.
 * 
 * Structure for store in files all data process by each pipeline phase.
 */
typedef struct write_batch {
  unsigned char flag;           /**< type of data stored*/
  unsigned int size;            /**< actual size of the batch (in bytes)*/
  unsigned int allocated_size;  /**< maximum size of the batch (in bytes)*/
  void* buffer_p;               /**< buffer to store data*/
} write_batch_t;

//====================================================================================

typedef struct region_batch {
  array_list_t **allocate_mapping_p;
  fastq_batch_t *unmapped_batch_p;
} region_batch_t;

typedef struct report_optarg {
  int all;
  int n_best;
  int n_hits;
  int only_paired;
  int best;
} report_optarg_t;

report_optarg_t *report_optarg_new(int all, int n_best, int n_hits, int only_paired, int best);

void report_optarg_free(report_optarg_t *p);

//====================================================================================

typedef struct pair_mng {
  int pair_mode;
  size_t min_distance;
  size_t max_distance;
  int report_only_paired;
} pair_mng_t;

pair_mng_t *pair_mng_new(int pair_mode, size_t min_distance, 
			 size_t max_distance, int report_only_upaired);

void pair_mng_free(pair_mng_t *p);

//=====================================================================================

// Added by PP
//====================================================================================

typedef struct bs_context {
  size_t CpG_methyl;                 /**< Partial Counter for methylated Cytosines in CpG context   */
  size_t CpG_unmethyl;               /**< Partial Counter for unmethylated Cytosines in CpG context */
  
  size_t CHG_methyl;                 /**< Partial Counter for methylated Cytosines in CHG context   */
  size_t CHG_unmethyl;               /**< Partial Counter for unmethylated Cytosines in CpG context */
  
  size_t CHH_methyl;                 /**< Partial Counter for methylated Cytosines in CHH context   */
  size_t CHH_unmethyl;               /**< Partial Counter for unmethylated Cytosines in CpG context */
  
  size_t MUT_methyl;                 /**< Partial Counter for mutated Cytosines                     */
  size_t num_bases;                  /**< Partial Counter for number of bases in the batch          */
  
  array_list_t *context_CpG;         /**< Array with the sequences from CpG context to write */
  array_list_t *context_CHG;         /**< Array with the sequences from CHG context to write */
  array_list_t *context_CHH;         /**< Array with the sequences from CHH context to write */
  array_list_t *context_MUT;         /**< Array with the sequences from mutations to write   */

  uint32_t *methyl_reads;            /**< Array with the number of methylated reads per chromosome */
} bs_context_t;

bs_context_t *bs_context_new(size_t num_reads, size_t num_chromosomes);
void bs_context_free(bs_context_t *bs_context);

//====================================================================================

typedef struct mapping_batch {
  int action;
  size_t num_targets;
  size_t num_extra_targets;
  size_t num_allocated_targets;
  size_t num_to_do;
  unsigned char extra_stage_do;
  unsigned char was_process;

  size_t num_gaps;
  size_t num_sws;
  size_t num_ext_sws;

  unsigned char *extra_stage_id;
  array_list_t *fq_batch;
  size_t *targets;
  size_t *extra_targets;
  array_list_t **mapping_lists;
  char *status; 
  pair_mng_t *pair_mng;
  array_list_t **old_mapping_lists;
  unsigned char *bwt_mappings;

  // histogram handling for filtering
  float margin;
  float *histogram_sw;

  // bs handling
  size_t num_targets2;
  size_t num_to_do2;
  size_t *targets2;

  array_list_t **mapping_lists2;
  array_list_t *CT_fq_batch;
  array_list_t *GA_fq_batch;

  array_list_t *CT_rev_fq_batch;
  array_list_t *GA_rev_fq_batch;

  array_list_t *bs_status;
  bs_context_t *bs_context;
  //  bs_context_t bs_context;
} mapping_batch_t;

mapping_batch_t *mapping_batch_new_2(size_t num_reads, array_list_t *fq_batch, pair_mng_t *pair_mng);
mapping_batch_t *mapping_batch_new(array_list_t *fq_batch, pair_mng_t *pair_mng);
mapping_batch_t *mapping_batch_new_by_num(size_t num_reads, pair_mng_t *pair_mng);
void mapping_batch_free(mapping_batch_t *p);

//=====================================================================================
//=====================================================================================

/**
 * @brief Structure for store all process data in @a region_seeker_server.
 * 
 * Structure for store process data in @a region_seeker_server, contains an 
 * array list with all mappings found for all seeds in each read and a
 * batch with all reads unmapped.
 */
typedef struct cal_batch {
  array_list_t **allocate_mapping;  /**< array list with all mappings found for all seeds in each read*/
  fastq_batch_t *unmapped_batch;   /**< batch with all reads unmapped*/
} cal_batch_t;

//------------------------------------------------------------------------------------

/**
 * @brief Structure for store all process data in @a cal_seeker_server.
 * 
 * Structure for store process data in @a cal_seeker_server, store all cals found
 * for each read unmapped.
 */
typedef struct sw_batch {
  unsigned int num_reads;          /**< number of reads allocate in the batch*/
  array_list_t **allocate_cals_p;  /**< array list that store all cals found for each read */
  fastq_read_t **allocate_reads_p; /**< array that store for each read the header, the sequence and the quality*/
} sw_batch_t;

typedef struct bwt_server_input bwt_server_input_t;
typedef struct region_seeker_input region_seeker_input_t;
typedef struct cal_seeker_input cal_seeker_input_t;
typedef struct pair_server_input pair_server_input_t;
typedef struct sw_server_input sw_server_input_t;
typedef struct batch_writer_input batch_writer_input_t;


typedef struct batch {
  int mapping_mode;
  bwt_server_input_t *bwt_input;
  region_seeker_input_t *region_input;
  cal_seeker_input_t *cal_input;
  pair_server_input_t *pair_input;
  sw_server_input_t *sw_input;
  batch_writer_input_t *writer_input;
  mapping_batch_t *mapping_batch;
  void *optional_data;
  bool write_mcontext;
} batch_t;


batch_t *batch_new(bwt_server_input_t *bwt_input,
                   region_seeker_input_t *region_input,
                   cal_seeker_input_t *cal_input,
                   pair_server_input_t *pair_input,
                   sw_server_input_t *sw_input,
                   batch_writer_input_t *writer_input,
		               int mapping_mode,
                   mapping_batch_t *mapping_batch,
				           bool write_mcontext);

void batch_free(batch_t *b);

//======================================================================================

typedef struct buffer_item {
  fastq_read_t *read;
  array_list_t *items_list;
  void *aux_data;
} buffer_item_t;

buffer_item_t *buffer_item_complete_new(fastq_read_t *read, array_list_t *cals_list, void *aux_data);
void buffer_item_insert_new_item(fastq_read_t *fq_read,
                                 linked_list_t *items_list,
                                 void *data,
                                 int type_items,
                                 linked_list_t *buffer,
                                 linked_list_t *buffer_hc,
                                 int phase);
void buffer_item_free(buffer_item_t *buffer_item);

void insert_file_item(fastq_read_t *fq_read, array_list_t *items, FILE *f_sa);

void insert_file_item_2(fastq_read_t *fq_read, array_list_t *items, FILE *f_hc);
//======================================================================================


typedef struct meta_alignment {
  int status;
  int sp_sw;
  int num_cigars;
  int type;
  int type_cigars[10];
  int score;
  float f_score;
  int flag;
  array_list_t *cals_list;
  cigar_code_t *middle_cigars[10];
  cigar_code_t *cigar_left;
  cigar_code_t *cigar_right;
  cigar_code_t *cigar_code;
} meta_alignment_t;

typedef struct simple_alignment {
  int gap_start;
  int gap_end;
  int map_strand;
  int map_chromosome;
  size_t map_start;
  int map_distance;
  int cigar_len;
} simple_alignment_t;

typedef struct alignment_aux {
  int mapping_len;
  short int optional_fields_length;
  short int chromosome;
  int position;
  int map_quality;
  int num_cigar_operations;
  unsigned char seq_strand;	
  int cigar_len;
} alignment_aux_t;


void alignment_aux_init(alignment_t* alignment, alignment_aux_t *alignment_aux);


cal_t *convert_bwt_anchor_to_CAL(bwt_anchor_t *bwt_anchor, 
				 size_t read_start, size_t read_end);

void file_write_items(fastq_read_t *fq_read, array_list_t *items, 
		      unsigned char data_type, FILE *fd1, FILE *fd2, int mode);

#endif // BUFFERS_H
