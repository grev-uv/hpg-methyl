#include "bs/bs_aligner.h"
#include "build-index/index_builder.h"


double emboss_matrix_t = 0.0f, emboss_tracking_t = 0.0f;
double sse_matrix_t = 0.0f, sse_tracking_t = 0.0f;
double sse1_matrix_t = 0.0f, sse1_tracking_t = 0.0f;
double avx_matrix_t = 0.0f, avx_tracking_t = 0.0f;
double avx1_matrix_t = 0.0f, avx1_tracking_t = 0.0f;

//--------------------------------------------------------------------
// constants
//--------------------------------------------------------------------

#define OPTIONS 			30
#define MIN_ARGC  			5
#define NUM_SECTIONS_STATISTICS 	5
#define NUM_SECTIONS_STATISTICS_SB	21

//--------------------------------------------------------------------
// global variables for timing and capacity meausures
//--------------------------------------------------------------------

char time_on = 0;
char statistics_on = 0;
timing_t* timing = NULL;
basic_statistics_t *basic_st;
cal_st_t cal_st;
double kl_time;

pthread_mutex_t mutex_sp;

size_t bwt_correct = 0;
size_t bwt_error = 0;
size_t seeding_reads = 0;
pthread_mutex_t bwt_mutex;

// timing
size_t TOTAL_READS_SEEDING;

struct timeval time_start_alig, time_end_alig;
struct timeval total_time_start, total_time_end;
float total_time = 0.0f;
double time_alig = 0.0;

size_t w2_r = 0;
size_t w3_r = 0;
size_t w2_3_r = 0;

size_t tot_reads_in = 0;
size_t tot_reads_out = 0;

FILE *fd_log;
size_t junction_id;

//--------------------------------------------------------------------
// main parameters support
//--------------------------------------------------------------------
int main(int argc, char* argv[]) {
  pthread_mutex_init(&mutex_sp, NULL);

  // If compiling under GCC, enable the application to check the
  // CPU features at runtime  
#ifdef __GNUC__
  __builtin_cpu_init();
#endif

  TOTAL_READS_SEEDING = 0;
  
  basic_st = basic_statistics_new();

  // Init logs, after parsing the command-line
  // logs will be re-set according to the command-line
  log_level = LOG_FATAL_LEVEL;
  log_verbose = 1;
  log_file = NULL;

  if (argc <= 1) {
    LOG_FATAL("Missing command.\nValid commands are:\n \
    \tbs: to map BS sequences\n \
    \tbuild-index: to create the genome index.\n \
    Use -h or --help to display hpg-aligner options.");
  }

  if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
    usage_cli();
  }

  char *command = argv[1];  

  // We need to consume command: {dna | rna | bs | bs-un | build-index}
  argc -= 1;
  argv += 1;

  if (strcmp(command, "bs") != 0) {
    LOG_FATAL("Command unknown.\nValid commands are:\n \
    \tbs: to map BS sequences\n \
    \tbuild-index: to create the genome index.\n \
    Use -h or --help to display hpg-aligner options.");
  }

  // Parsing options
  options_t *options = parse_options(argc, argv);

  // Now, we can set logs according to the command-line
  init_log_custom(options->log_level, 1, "hpg-aligner.log", "w");

  validate_options(options, command);
  LOG_DEBUG_F("Command Mode: %s\n", command);
  
  if (!strcmp(command, "build-index")) {
    if (options->bs_index == 0) { // Regular index generation
      run_index_builder(options->genome_filename, options->bwt_dirname, options->index_ratio);
      LOG_DEBUG("Done !!\n");
      exit(0);
    } else {
      // Bisulphite index generation
      printf("\nBisulphite index generation\n");
      
      /******************************************************************************
       * 										                                    *
       * Generates the genome transform from the input and builds the index		    *
       * 										                                    *
       * The genome transformed are stored in the directory give by the user,		*
       * and the index are stored in subfolders				       		            *
       * 										                                    *
       ******************************************************************************/

      run_index_builder(options->genome_filename, options->bwt_dirname, options->index_ratio);

      // Generate binary code for original genome
      char bs_dir1[256];
      sprintf(bs_dir1, "%s/AGT_index", options->bwt_dirname);
      create_directory(bs_dir1);

      LOG_DEBUG("Generation of AGT index\n");

      char genome_1[256];
      sprintf(genome_1, "%s/AGT_genome.fa", options->bwt_dirname);

      // Substitute all cytosines by thymines in the genome file
      char gen1[256];
      sprintf(gen1, "sed 's/C/T/g' %s > %s",options->genome_filename, genome_1);
      system(gen1);

      // Generate the new index for the AGT genome file
      run_index_builder_bs(genome_1, bs_dir1, options->index_ratio, "AGT");
      LOG_DEBUG("AGT index Done !!\n");

      LOG_DEBUG("Generation of ACT index\n");
      char bs_dir2[256];
      sprintf(bs_dir2, "%s/ACT_index", options->bwt_dirname);

      create_directory(bs_dir2);

      char genome_2[256];
      sprintf(genome_2, "%s/ACT_genome.fa", options->bwt_dirname);

      // Substitute all guanines by adenines in the genome file
      char gen2[256];
      sprintf(gen2, "sed 's/G/A/g' %s > %s",options->genome_filename, genome_2);
      system(gen2);

      // Generate the new index for the ACT genome file
      run_index_builder_bs(genome_2, bs_dir2, options->index_ratio, "ACT");
      LOG_DEBUG("ACT index Done !!\n");
      
      LOG_DEBUG("Generation of CT & GA context\n");
      encode_context(options->genome_filename, options->bwt_dirname);
      LOG_DEBUG("CT & GA context Done !!\n");

      exit(0);
    }
  }

  time_on =  (unsigned int) options->timming;
  statistics_on =  (unsigned int) options->statistics;

  // Create the timing related structures
  // (check buffers.h and timing.h)
  if (time_on) {
    char* timingLabels[] = {
      "FASTQ_READER",
      "BWT_SERVER",
      "BWT_SERVER_HISTOGRAM",
      "CAL_SEEKER",
      "PRE_PAIR_TIME",
      "SW_STAGE_TIME",
      "POST_PAIR_TIME",
      "METHYLATION_REP_TIME",
      "BAM_WRITER",
      "BS_ALIGNMENT_TOTAL",
      "TOTAL_TIME"
    };

    timing = timing_new(timingLabels, TOTAL_TIME + 1);
    start_timer(total_time_start);
  }

  struct timeval time_genome_s;
  double time_genome = 0.0;

  genome_t *genome, *genome1, *genome2;
  bwt_index_t *bwt_index1, *bwt_index2;

  start_timer(time_genome_s);

  // Load BWT indices
  char bs_dir1[256];
  sprintf(bs_dir1, "%s/AGT_index", options->bwt_dirname);

  char bs_dir2[256];
  sprintf(bs_dir2, "%s/ACT_index", options->bwt_dirname);

  // Genome parameters
  LOG_DEBUG("Reading genomes...");

  genome  = genome_new("dna_compression.bin", options->bwt_dirname);
  genome1 = genome_new("dna_compression.bin", bs_dir1);
  genome2 = genome_new("dna_compression.bin", bs_dir2);

  LOG_DEBUG("Done !!");
  
  // BWT index
  LOG_DEBUG("Loading AGT index...");
  bwt_index1 = bwt_index_new(bs_dir1);
  LOG_DEBUG("Loading AGT index done !!");

  LOG_DEBUG("Loading ACT index...");
  bwt_index2 = bwt_index_new(bs_dir2);
  LOG_DEBUG("Loading ACT index done !!");

  //BWT parameters
  bwt_optarg_t *bwt_optarg = bwt_optarg_new(1, 0,
					    options->filter_read_mappings, 
					    options->filter_seed_mappings,
              options->min_cal_size,
              options->umbral_cal_length_factor,
              options->min_read_discard,
              options->max_inner_gap);
  
  // CAL parameters
  cal_optarg_t *cal_optarg = cal_optarg_new(options->min_cal_size, options->seeds_max_distance, 
					    options->num_seeds, options->min_num_seeds_in_cal,
					    options->seed_size, options->min_seed_size, 
					    options->cal_seeker_errors);
  
  // Paired mode parameters
  pair_mng_t *pair_mng = pair_mng_new(options->pair_mode, options->pair_min_distance, 
				      options->pair_max_distance, options->report_only_paired);
  
  // Report parameters
  report_optarg_t *report_optarg = report_optarg_new(options->report_all,
						     options->report_n_best,
						     options->report_n_hits, 
						     options->report_only_paired,
						     options->report_best);  


  // Start main process execution.
  if (!strcmp(command, "bs")) {
    // BS version
    run_bs_aligner(genome1, genome2, genome, bwt_index1, bwt_index2,
		   bwt_optarg, cal_optarg, pair_mng, report_optarg, options);
  }
  
  LOG_DEBUG("main done !!");

  // Free memory
  genome_free(genome);

  bwt_index_free(bwt_index1);
  genome_free(genome1);

  bwt_index_free(bwt_index2);
  genome_free(genome2);

  bwt_optarg_free(bwt_optarg);
  cal_optarg_free(cal_optarg);
  pair_mng_free(pair_mng);
  report_optarg_free(report_optarg);

  if (time_on) {
    stop_timer(total_time_start, total_time_end, total_time);
    timing_add(total_time, TOTAL_TIME, timing);

    timing_display(timing);
    timing_free(timing);
  }

  basic_statistics_display(basic_st, 0, time_alig / 1000000, time_genome / 1000000);
  free(basic_st);
  options_free(options);

  // New addition
  stop_log();
  return 0;
}
