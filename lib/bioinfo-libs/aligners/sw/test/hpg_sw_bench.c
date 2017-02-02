#include "hpg_sw_bench.h"

double *sse_matrix_t, *sse_tracking_t;
double *emboss_matrix_t, *emboss_tracking_t, *emboss_total_t;

/*
double *emboss_total_time, *emboss_matrix_time, *emboss_tracking_time;

double emboss_matrix_t = 0.0f, emboss_tracking_t = 0.0f;
double sse_matrix_t = 0.0f, sse_tracking_t = 0.0f;
double sse1_matrix_t = 0.0f, sse1_tracking_t = 0.0f;

double emboss_t = 0.0f;
double sse_t = 0.0f;
*/

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

int main(int argc, char *argv[]) {

  int c;
  char *q_filename = NULL, *r_filename = NULL;
  char *out_dirname = NULL, *sse_filename = NULL, *emboss_filename = NULL, *subst_matrix_name = "./dnafull";
  int num_threads = 4, batch_size = 2048;
  float gap_open = 10.0f;
  float gap_extend = 0.5f;
     
  while (1) {;
    static struct option long_options[] =
      {
	{"query-file",             required_argument, 0, 'q'},
	{"ref-file",               required_argument, 0, 'r'},
	{"output-dir",             required_argument, 0, 'd'},
	{"output-sse-file",        required_argument, 0, 'o'},
	{"output-emboss-file",     required_argument, 0, 'v'},
	{"gap-open-penalty",       required_argument, 0, 'p'},
	{"gap-extend-penalty",     required_argument, 0, 'e'},
	{"num-threads",            required_argument, 0, 'n'},
	{"reads-per-batch",        required_argument, 0, 'b'},
	{"help",                   no_argument,       0, 'h'},
	{0, 0, 0, 0}
      };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    
    c = getopt_long (argc, argv, "q:r:d:o:v:p:e:n:b:h",
		     long_options, &option_index);
    
    /* Detect the end of the options. */
    if (c == -1)
      break;
    
    switch (c) {
    case 0:
      /* If this option set a flag, do nothing else now. */
      if (long_options[option_index].flag != 0)
	break;
      printf ("option %s", long_options[option_index].name);
      if (optarg)
	printf (" with arg %s", optarg);
      printf ("\n");
      break;
      
    case 'q':
      q_filename = optarg;
      break;
      
    case 'r':
      r_filename = optarg;
      break;
      
    case 'd':
      out_dirname = optarg;
      break;

    case 'o':
      sse_filename = optarg;
      break;

    case 'v':
      emboss_filename = optarg;
      break;
      
    case 'p':
      sscanf(optarg, "%f", &gap_open);
      break;
      
    case 'e':
      sscanf(optarg, "%f", &gap_extend);
      break;
      
    case 'n':
      sscanf(optarg, "%i", &num_threads);
      break;

    case 'b':
      sscanf(optarg, "%i", &batch_size);
      break;

    case '?':
    case 'h':
    default:
      display_usage(NULL);
    }
  }

  // sanity cheking
  if (q_filename == NULL) {
    display_usage("\nERROR: Query file name is missing (option --query-file or -q)\n");
  } else if (!file_exists(q_filename)) {
    printf("\nERROR: Query file %s does not exist\n", q_filename);
    exit(-1);
  }

  printf("query-file = %s\n", q_filename);

  if (r_filename == NULL) {
    display_usage("\nERROR: Reference file name is missing (option --ref-file or -r)\n");
  } else if (!file_exists(r_filename)) {
    printf("\nERROR: Reference file %s does not exist\n", r_filename);
    exit(-1);
  }
  printf("ref-file = %s\n", r_filename);

  if (out_dirname == NULL) {
    display_usage("\nERROR: Output directory name is missing (option --output-dir or -d)\n");
  } else if (!file_exists(out_dirname)) {
    printf("\nERROR: Output directory %s does not exist\n", out_dirname);
    exit(-1);
  }
  printf("output-dir = %s\n", out_dirname);

  if (sse_filename == NULL) {
    display_usage("\nERROR: Ouptut SSE file name is missing (option --output-sse-file or -o)\n");
  }
  printf("output-sse-file = %s\n", sse_filename);

  if (emboss_filename == NULL) {
    display_usage("\nERROR: Ouptut EMBOSS file name is missing (option --output-emboss-file or -v)\n");
  }
  printf("output-emboss-file = %s\n", emboss_filename);

  if (gap_open < 0.0 || gap_open > 100.0) {
    display_usage("\nERROR: Penalty for the gap openning is not valid, it must be from 0.0 to 100.0 (option --gap-open-penalty or -p)\n");
  }
  printf("gap-open-penalty = %f\n", gap_open);  

  if (gap_extend < 0.0 || gap_extend > 10.0) {
    display_usage("\nERROR: Penalty for the gap extending is not valid, it must be from 0.0 to 10.0 (option --gap-extend-penalty or -e)\n");
  }
  printf("gap-extend-penalty = %f\n", gap_extend);  

  if (num_threads <= 0) {
    display_usage("\nERROR: Number of threads is not valid, it must be greater than 0 (option --num-threads or -n)\n");
  }
  printf("num-threads = %i\n", num_threads);

  if (batch_size <= 0) {
    display_usage("\nERROR: Number of reads per batch is not valid, it must be greater than 0 (option --reads-per-batch or -b)\n");
  }
  printf("reads-per-batch = %i\n", batch_size);

  char sse_out_filename[2048];
  sprintf(sse_out_filename, "%s/%s", out_dirname, sse_filename);

  printf("out SSE filename (full path) = %s\n", sse_out_filename);

  char emboss_out_filename[2048];
  sprintf(emboss_out_filename, "%s/%s", out_dirname, emboss_filename);

  printf("out EMBOSS filename (full path) = %s\n", emboss_out_filename);
  printf("\n");

  float match = 5.0f;
  float mismatch = 4.0f;


  sse_matrix_t = (double *) calloc(num_threads, sizeof(double));
  sse_tracking_t = (double *) calloc(num_threads, sizeof(double));

  emboss_matrix_t = (double *) calloc(num_threads, sizeof(double));
  emboss_tracking_t = (double *) calloc(num_threads, sizeof(double));
  emboss_total_t = (double *) calloc(num_threads, sizeof(double));

  // SSE
  run_sse(q_filename, r_filename,
	  gap_open, gap_extend, "./dnafull", batch_size, 
	  num_threads, sse_out_filename);
  
  printf("SSE + OpenMP\n");
  double max_sse = 0.0f;
  for (int i = 0; i < num_threads; i++) {
    if (sse_matrix_t[i] + sse_tracking_t[i] > max_sse)
      max_sse = sse_matrix_t[i] + sse_tracking_t[i];
  }
  printf("\tScore matrix creation time:\n");
  for (int i = 0; i < num_threads; i++) {
    printf("\t\tThread %2i:\t%0.5f s\n", i, sse_matrix_t[i]); 
  }
  printf("\tBacktracking time:\n");
  for (int i = 0; i < num_threads; i++) {
    printf("\t\tThread %2i:\t%0.5f s\n",  i, sse_tracking_t[i]); 
  }
  printf("\tTotal time:\n");
  printf("\t\tMax. time:\t%0.5f s\n", max_sse);

  // EMBOSS
  run_emboss(q_filename, r_filename,
	     gap_open, gap_extend, match, mismatch,
	     batch_size, num_threads, emboss_out_filename);
					   
  printf("EMBOSS + OpenMP\n");
  double max_emboss = 0.0f;
  for (int i = 0; i < num_threads; i++) {
    if (emboss_matrix_t[i] + emboss_tracking_t[i] > max_emboss)
      max_emboss = emboss_matrix_t[i] + emboss_tracking_t[i];
  }
  printf("\tScore matrix creation time:\n");
  for (int i = 0; i < num_threads; i++) {
    printf("\t\tThread %2i:\t%0.5f s\n", i, emboss_matrix_t[i]); 
  }
  printf("\tBacktracking time:\n");
  for (int i = 0; i < num_threads; i++) {
    printf("\t\tThread %2i:\t%0.5f s\n",  i,  emboss_tracking_t[i]); 
  }
  printf("\tTotal time:\n");
  printf("\t\tMax. time:\t%0.5f s\n", max_emboss);

  // speed-up
  printf("\nSpeed-up (EMBOS vs SSE)\n");
  printf("\tTotal             : %0.2f\n", max_emboss / max_sse);
  printf("\n");

  free(sse_matrix_t);
  free(sse_tracking_t);

  free(emboss_matrix_t);
  free(emboss_tracking_t);
  free(emboss_total_t);
}

//-------------------------------------------------------------------------

void run_emboss(char *q_filename, char *r_filename, 
		float gap_open, float gap_extend,
		float match, float mismatch,
		int batch_size, int num_threads, 
		char *out_filename) {  
  
  const int max_length = 2048;

  float score[batch_size];

  FILE *q_file = fopen(q_filename, "r");
  FILE *r_file = fopen(r_filename, "r");
  FILE *out_file = fopen(out_filename, "w");

  char *q[batch_size], *r[batch_size];
  char *m[batch_size], *n[batch_size];
  int start1[batch_size], start2[batch_size];

  for (int i = 0; i < batch_size; i++) {
    q[i] = (char *) calloc(max_length, sizeof(char));
    r[i] = (char *) calloc(max_length, sizeof(char));

    m[i] = (char *) calloc(4048, sizeof(char));
    n[i] = (char *) calloc(4048, sizeof(char));
  }
  
  unsigned int tid, count = 0, batches = 0, num_queries;

  omp_set_num_threads(num_threads);

  while (1) {
    num_queries = 0;

    // read queries
    for (int i = 0; i < batch_size; i++) {
      if (fgets(q[i], max_length, q_file) == NULL) { break; }
      str_trim(q[i]);
      num_queries++;
      count++;
    }

    // exit if no queries
    if (num_queries == 0) break;

    // read references
    for (int i = 0; i < num_queries; i++) {
      if (fgets(r[i], max_length, r_file) == NULL) { break; }
      str_trim(r[i]);

    }

    if (num_threads == 1) {
      tid = 0;
      for (int i = 0; i < num_queries; i++) {
	score[i] = sw_emboss(r[i], q[i], gap_open, gap_extend, m[i], n[i],
			     &start1[i], &start2[i], 
			     &emboss_matrix_t[tid], &emboss_tracking_t[tid],
			     &emboss_total_t[tid]);

      }
    } else {
      // call smith-waterman
      #pragma omp parallel for
      for (int i = 0; i < num_queries; i++) {
	tid = omp_get_thread_num();
	score[i] = sw_emboss(r[i], q[i], gap_open, gap_extend, m[i], n[i],
			     &start1[i], &start2[i], 
			     &emboss_matrix_t[tid], &emboss_tracking_t[tid],
			     &emboss_total_t[tid]);

      }      
    }
    for (int i = 0; i < num_queries; i++) {
	fprintf(out_file, "%s\n%s\nScore = %0.2f, query start at %i, ref. start at %i\n\n", 
		n[i], m[i], score[i],
		start2[i], start1[i]);
    }


    batches++;
  }

  /*
  for (int i = 0; i < batch_size; i++) {
    free(q[i]);
    free(r[i]);
    free(m[i]);
    free(n[i]);
  }
  */
  printf("\nEMBOSS done (%i reads)\n", count);

  fclose(q_file);
  fclose(r_file);
  fclose(out_file);
}

//-------------------------------------------------------------------------

void run_sse(char *q_filename, char *r_filename, 
	     float gap_open, float gap_extend, char *matrix_filename,
	     int batch_size, int num_threads,
	     char *out_filename) {  

  const int max_length = 2048;

  sw_optarg_t *optarg_p = sw_optarg_new(gap_open, gap_extend, matrix_filename);
  sw_multi_output_t *output_p;


  FILE *q_file = fopen(q_filename, "r");
  FILE *r_file = fopen(r_filename, "r");
  FILE *out_file = fopen(out_filename, "w");

  char *q[batch_size], *r[batch_size];
  for (int i = 0; i < batch_size; i++) {
    q[i] = (char *) calloc(max_length, sizeof(char));
    r[i] = (char *) calloc(max_length, sizeof(char));
  }
  
  int count = 0, batches = 0, num_queries;

  while (1) {
    num_queries = 0;

    // read queries
    for (int i = 0; i < batch_size; i++) {
      if (fgets(q[i], max_length, q_file) == NULL) { break; }
      str_trim(q[i]);
      num_queries++;
      count++;
    }

    // exit if no queries
    if (num_queries == 0) break;

    // read references
    for (int i = 0; i < num_queries; i++) {
      if (fgets(r[i], max_length, r_file) == NULL) { break; }
      str_trim(r[i]);
    }

    //printf("num queries = %i, count = %i\n", num_queries, count);
    output_p = sw_multi_output_new(num_queries);

    // call smith-waterman
   smith_waterman_mqmr(q, r, num_queries, optarg_p, num_threads, output_p);

    // save results
    sw_multi_output_save(num_queries, output_p, out_file);
    //break;

    sw_multi_output_free(output_p);

    batches++;
    //printf("%i batches\n", batches);
    //    if (batches == 7) break;
  }
  printf("\nSSE done (%i reads)\n", count);

  // free memory and close files
  sw_optarg_free(optarg_p);

  for (int i = 0; i < batch_size; i++) {
    free(q[i]);
    free(r[i]);
  }

  fclose(q_file);
  fclose(r_file);
  fclose(out_file);
}


//-------------------------------------------------------------------------

void display_usage(char *msg) {
  if (msg != NULL) printf("%s", msg);
  printf("Usage: ./hpg-sw -q query_filename -r ref_filename -d output_dirname -o output_sse_filename -v output_emboss_filename -p gap_open_penalty -e gap_extend_penalty -s substitution_score_matrix -n number_of_threads -b number_of_reads_per_batch\n");

  printf("\t--query-file, -q, input file name containing the queries to align\n");
  printf("\t--ref-file, -r, input file name containing the references to align\n");
  printf("\t--output-dir, -d, output directory where the results will be saved\n");
  printf("\t--output-sse-file, -o, output SSE file name where the alignments will be saved\n");
  printf("\t--output-emboss-file, -v, output EMBOSS file name where the alignments will be saved\n");
  printf("\t--gap-open-penalty, -p, penalty for the gap openning: from 0.0 to 100.0\n");

  printf("\t--gap-extend-penalty, -e, penalty for the gap extending: from 0.0 to 10.0\n");
  //  printf("\t--substitution-matrix, -s, substitution score matrix name, for DNA: dnafull, for proteins: blosum50, blosum62, blosum80\n");
  printf("\t--num-threads, -n, number of threads\n");
  printf("\t--reads-per-batch, -b, number of reads per batch\n");

  exit(-1);
}

//-------------------------------------------------------------------------

int file_exists(const char *filename) {
  FILE *file = NULL;
  if (file = fopen(filename, "r")) {
    fclose(file);
    return 1;
  }
  return 0;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
