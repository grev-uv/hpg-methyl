#ifndef OPTIONS_H
#define OPTIONS_H

/*
 *  options.h
 *
 *  Created on: Aug 29, 2012
 *      Author: imedina
 */

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "argtable/argtable2.h"
#include "config/libconfig.h"
#include "commons/log.h"
#include "commons/system_utils.h"
#include "commons/file_utils.h"

//============================ DEFAULT VALUES ============================
#define DEFAULT_GPU_THREADS		         32
#define DEFAULT_CPU_THREADS		         1
#define DEFAULT_CAL_SEEKER_ERRORS	     0
#define DEFAULT_MIN_SEED_PADDING_LEFT	 5
#define DEFAULT_MIN_SEED_PADDING_RIGHT	 5
#define DEFAULT_WRITE_BATCH_SIZE	     500000
#define DEFAULT_NUM_CAL_SEEKERS		     1
#define DEFAULT_REGION_THREADS		     1
#define DEFAULT_NUM_SW_THREADS		     1

// BEGIN: Mariano (4/12/2014)
#define DEFAULT_NUM_SEEDS		         10
// END: Mariano (4/12/2014)

// BEGIN: Ricardo
#define DEFAULT_UMBRAL_CAL_LENGTH_FACTOR 4
#define DEFAULT_MIN_READ_DISCARD         100
// END: Ricardo

#define DEFAULT_MIN_NUM_SEEDS_IN_CAL	 -1
#define DEFAULT_MAX_INTRON_LENGTH	     800000
#define DEFAULT_MIN_INTRON_LENGTH	     40
#define DEFAULT_SW_MIN_SCORE		     0.8
#define DEFAULT_SW_MATCH		         5
#define DEFAULT_SW_MISMATCH		         -4
#define DEFAULT_SW_GAP_OPEN		         10
#define DEFAULT_SW_GAP_EXTEND		     0.5
#define DEFAULT_PAIR_MODE	             0
#define DEFAULT_PAIR_MIN_DISTANCE	     50
#define DEFAULT_PAIR_MAX_DISTANCE	     800
#define MINIMUM_BATCH_SIZE               10000
#define DEFAULT_FILTER_READ_MAPPINGS     50
#define DEFAULT_FILTER_SEED_MAPPINGS     500
// New variable for default uses
#define DEFAULT_NUCLEOTIDES              "ACGT"
#define DEFAULT_FILTER_READ_MAPPINGS_BS  100
#define DEFAULT_FILTER_SEED_MAPPINGS_BS  500
//=====================================================================
#define NUM_OPTIONS			             50

typedef struct options {
  char mode[64];
  unsigned char bwt_set;
  unsigned char reg_set;
  unsigned char cal_set;
  unsigned char sw_set;
  int min_intron_length;
  int num_gpu_threads;
  int num_cpu_threads;
  int min_cal_size; 
  int num_seeds; 
  int min_num_seeds_in_cal; 
  int seeds_max_distance;
  int min_seed_padding_right;
  int min_seed_padding_left;
  int batch_size;
  int write_size;
  int num_cal_seekers;
  int min_seed_size;
  int seed_size;
  int max_intron_length;
  int flank_length;
  int timming;
  int statistics;
  int rna_seq; 
  int help;
  int cal_seeker_errors;
  int pair_mode;
  int pair_min_distance;
  int pair_max_distance;
  int report_all;
  int report_n_best;
  int report_n_hits;
  int report_best;
  int report_only_paired;
  int gpu_process;
  int log_level;
  int index_ratio;
  int filter_read_mappings;
  int filter_seed_mappings;
  int bs_index;
  double min_score;
  double match;
  double mismatch;
  double gap_open;
  double gap_extend;
  char* prefix_name;
  char* in_filename;
  char* in_filename2;
  char* bwt_dirname;
  char* genome_filename;
  char* output_name;
  char* header_filename;
  char* transcriptome_filename;
  char* intron_filename;
  bool write_mcontext;
  //RICARDO
  double umbral_cal_length_factor;
  int min_read_discard;
  int max_inner_gap;
  //RICARDO
  // new variables for bisulphite case
} options_t;


options_t *options_new(void);

void options_free(options_t *options);

void options_display(options_t *options);

void usage_cli();

void validate_options(options_t *options, char *mode);

/**
 * @brief Initializes an global_options_t structure mandatory members.
 * @return A new global_options_t structure.
 *
 * Initializes the only mandatory member of a global_options_t, which is the output directory.
 */
void** argtable_options_new(void);


/**
 * @brief Free memory associated to a global_options_data_t structure.
 * @param options_data the structure to be freed
 *
 * Free memory associated to a global_options_data_t structure, including its text buffers.
 */
void argtable_options_free(void **argtable_options);


/**
 * @brief Initializes an global_options_data_t structure mandatory members.
 * @return A new global_options_data_t structure.
 *
 * Initializes the only mandatory member of a global_options_data_t, which is the output directory.
 */
options_t *read_CLI_options(void **argtable_options, options_t *options);

options_t *parse_options(int argc, char **argv);


void usage(void **argtable);

#endif
