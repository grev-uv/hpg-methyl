#ifndef BWT_H
#define BWT_H

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <pthread.h> 

#include "commons/string_utils.h"
#include "containers/array_list.h"
#include "containers/linked_list.h"
#include "bioformats/fastq/fastq_read.h"
#include "bioformats/fastq/fastq_batch.h"
#include "bioformats/bam/alignment.h"

#include "BW_io.h"
#include "BW_search.h"
#include "BW_preprocess.h"

#define NONE_HARD_CLIPPING 0

#define BACKWARD_ANCHOR 0
#define FORWARD_ANCHOR  1

#ifndef MAX
  #define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

//-----------------------------------------------------------------------------
// Paratemers for the candidate alignment localizations (CALs)
//-----------------------------------------------------------------------------

typedef struct cal_optarg {
  size_t min_cal_size;
  size_t max_cal_distance;
  size_t num_seeds;
  size_t min_num_seeds_in_cal;
  size_t seed_size;
  size_t min_seed_size;
  size_t num_errors;
} cal_optarg_t;

cal_optarg_t *cal_optarg_new(const size_t min_cal_size, 
			     const size_t max_cal_distance, 
			     const size_t num_seeds,
			     const size_t min_num_seeds_in_cal,
			     const size_t seed_size,
			     const size_t min_seed_size,
			     const size_t num_errors);

void cal_optarg_free(cal_optarg_t *optarg);

//-----------------------------------------------------------------------------

typedef struct seed_region {
  size_t read_start;
  size_t read_end;
  size_t genome_start;
  size_t genome_end;
  int id;
  int fusion_left;
  int fusion_right;
  void *info;
} seed_region_t;

seed_region_t *seed_region_new(size_t read_start, size_t read_end, 
			       size_t genome_start, size_t genome_end, int id);

void seed_region_free();

//-----------------------------------------------------------------------------

typedef struct cal {
  size_t chromosome_id;
  short int strand;
  size_t start;
  size_t end;
  size_t num_seeds;
  int read_area;
  int l_flank;
  int r_flank;
  int fill_gaps;
  int num_targets;
  linked_list_t *sr_list;
  linked_list_t *sr_duplicate_list;
  array_list_t *candidates_seeds_start;
  array_list_t *candidates_seeds_end;
  void *info;
} cal_t;

cal_t *cal_new(const size_t chromosome_id,
               const short int strand,
               const size_t start,
               const size_t end,
               const size_t num_seeds,
               const linked_list_t *sr_list,
               const linked_list_t *sr_duplicate_list);

void cal_free(cal_t *cal);


//-----------------------------------------------------------------------------
//Ricardo - Modified in order to perform BWT in interior intervals of the read
typedef struct bwt_anchor {
  int strand;
  int chromosome;
  size_t start;
  size_t end;
  int type;
  size_t seq_start;		//new for interior intervals
  size_t seq_end;		//new for interior intervals
  size_t seq_len;		//new for interior intervals
} bwt_anchor_t;

bwt_anchor_t *bwt_anchor_new(int strand, int chromosome, size_t start, size_t end, int type, size_t seq_start, size_t seq_end);
void bwt_anchor_free(bwt_anchor_t *bwt_anchor);

//-----------------------------------------------------------------------------

typedef struct region {
  size_t chromosome_id;
  short int strand;
  size_t start;
  size_t end;
  size_t seq_start;
  size_t seq_end;
  size_t seq_len;
  int id;
} region_t;


//-------------------- to store the internal bwt interval, (sequence start and end) //RICARDO
typedef struct interval {
	size_t start;
	size_t end;
} interval_t;

//-----------------------------------------------------------------------------

typedef struct short_cal {
  size_t start;
  size_t end;
  size_t seq_len;
  size_t num_seeds;
  size_t seq_start;
  size_t seq_end;
  linked_list_t *sr_list;
  linked_list_t *sr_duplicate_list;
  unsigned char *seeds_ids_array;
} short_cal_t;

short_cal_t *short_cal_new(const size_t start, 
	                   const size_t end,
			   const size_t seq_start,
			   const size_t seq_end,
			   const size_t seq_len,
			   const int max_seeds,
			   const int id);

void short_cal_free(short_cal_t *short_cal_p);

//-----------------------------------------------------------------------------

typedef struct read_cals {
  fastq_read_t *read;
  array_list_t *cal_list; // array list of cal_t structures
} read_cals_t;



//-----------------------------------------------------------------------------
// Burrows-Wheeler Transform
//-----------------------------------------------------------------------------

typedef struct bwt_optarg {
  size_t num_errors;
  size_t num_threads;
  int filter_read_mappings;
  int filter_seed_mappings;
  size_t min_cal_size;
  double umbral_cal_length_factor;
  int min_read_discard;
  int max_inner_gap;
} bwt_optarg_t;

bwt_optarg_t *bwt_optarg_new(const size_t num_errors,
			     const size_t num_threads,
			     const int filter_read_mappings, 
			     const int filter_seed_mappings,
				 const size_t min_cal_size,
				 const double umbral_cal_length_factor,
				 const int min_read_discard,
				 const int max_inner_gap);

void bwt_optarg_free(bwt_optarg_t *optarg);

//-----------------------------------------------------------------------------

typedef struct bwt_index {
  char *dirname;
  char *nucleotides;
  comp_matrix h_O, h_rO, h_Oi, h_rOi;
  vector h_C, h_rC, h_C1, h_rC1;
  byte_vector B;
  comp_vector S, Si;
  exome karyotype;
  int table[128];
  int rev_table[4];
} bwt_index_t;

bwt_index_t *bwt_index_new(const char *dirname);
//void bwt_index_new(const char *dirname, bwt_index_t **index);
void bwt_index_free(bwt_index_t *index);


void bwt_generate_index_files(char *ref_file, char *output_dir, 
			      unsigned int s_ratio);

void bwt_generate_index_files_bs(char *ref_file, char *output_dir, 
				 unsigned int s_ratio, char *bases);


//-----------------------------------------------------------------------------

typedef struct mapping {
  cal_t *cal;
  char *cigar;
} mapping_t;

//-----------------------------------------------------------------------------

typedef struct read_mappings {
  fastq_read_t *read;
  array_list_t *mapping_list; // array list of mapping_t structures
} read_mappings_t;

//-----------------------------------------------------------------------------
// general functions
//-----------------------------------------------------------------------------

alignment_t* add_optional_fields(alignment_t *alignment, size_t n_mappings);


/**
 */
char* reverse_str(char *src, char *dsp, size_t length);


size_t bwt_map_seq(char *seq, 
		   bwt_optarg_t *bwt_optarg, 
		   bwt_index_t *index, 
		   array_list_t *mapping_list);

size_t bwt_map_read(fastq_read_t *read, 
		    bwt_optarg_t *bwt_optarg, 
		    bwt_index_t *index, 
		    array_list_t *mapping_list);

size_t bwt_map_seqs(char **seqs, 
		    size_t num_reads,
		    bwt_optarg_t *bwt_optarg, 
		    bwt_index_t *index, 
		    char *out_status,
		    array_list_t *mapping_list);

size_t bwt_map_reads(fastq_read_t **reads, 
		     bwt_optarg_t *bwt_optarg, 
		     bwt_index_t *index, 
		     char *out_status,
		     array_list_t *mapping_list);

size_t bwt_map_batch(fastq_batch_t *batch,
		     bwt_optarg_t *bwt_optarg, 
		     bwt_index_t *index, 
		     fastq_batch_t *unmapped_batch,
		     array_list_t *mapping_list);

size_t bwt_map_inexact_batch(fastq_batch_t *batch,
			     bwt_optarg_t *bwt_optarg, 
			     bwt_index_t *index, 
			     fastq_batch_t *unmapped_batch,
			     array_list_t *mapping_list);

interval_t* bwt_map_inexact_read_bs_new(fastq_read_t *read,
			       bwt_optarg_t *bwt_optarg,
			       bwt_index_t *index,
			       array_list_t *mapping_list,
			       int type,
				   size_t  start_read, size_t  end_read);

//-----------------------------------------------------------------------------
// seed functions
//-----------------------------------------------------------------------------

size_t bwt_map_exact_seeds_between_coords(int start_position, int end_position, 
					  char *seq, int seed_size, int min_seed_size,
					  bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
					  array_list_t *mapping_list, int extra_seed,
					  int *last_seed_id);


//-----------------------------------------------------------------------------
// cal functions
//-----------------------------------------------------------------------------



size_t bwt_generate_cals_between_coords(int strand_target, int chromosome_target,
					size_t start_target, size_t end_target, 
					int start_position, int end_posistion, 
					char *seq, int seed_size, int min_seed_size,
					bwt_optarg_t *bwt_optarg, bwt_index_t *index, 
					array_list_t *init_list, array_list_t *cal_list);

size_t bwt_find_cals_from_reads(fastq_read_t **reads, 
				bwt_optarg_t *bwt_optarg, 
				bwt_index_t *index, 
				cal_optarg_t *cal_optarg, 
				char *out_status,
				array_list_t *cal_list);

size_t bwt_map_inexact_array_list(array_list_t *reads,
				  bwt_optarg_t *bwt_optarg, 
				  bwt_index_t *index,
				  array_list_t **lists,
				  size_t *num_unmapped, 
				  size_t *unmapped_indices);

void bwt_map_inexact_array_list_by_filter(array_list_t *reads,
					  bwt_optarg_t *bwt_optarg, 
					  bwt_index_t *index,
					  array_list_t **lists,
					  size_t *num_unmapped, 
					  size_t *unmapped_indices);


size_t bwt_map_forward_inexact_seq(char *seq, 
				   bwt_optarg_t *bwt_optarg, 
				   bwt_index_t *index, 
				   array_list_t *mapping_list);
  
size_t bwt_generate_cal_list_rna_linked_list(array_list_t *mapping_list,
					     cal_optarg_t *cal_optarg,
					     array_list_t *cal_list,
					     size_t read_length,
					     size_t nchromosomes);


size_t bwt_generate_cal_rna_list_linked_list(array_list_t *mapping_list,
                                             cal_optarg_t *cal_optarg,
                                             size_t *min_seeds, int *max_seeds,
                                             size_t nchromosomes,
                                             array_list_t *cal_list,
                                             size_t read_length);

void append_seed_region_linked_list(linked_list_t* sr_list,
				    size_t read_start, size_t read_end, 
				    size_t genome_start, size_t genome_end, 
				    int seed_id);


//-----------------------------------------------------------------------------
void initReplaceTable_bs(const char *str);

char * readNucleotide(const char *directory, const char *name);
//-----------------------------------------------------------------------------

void insert_seeds_and_merge_anchor(array_list_t *mapping_list, linked_list_t ***cals_list,  size_t max_cal_distance, int id_start);



#endif // BWT_H
