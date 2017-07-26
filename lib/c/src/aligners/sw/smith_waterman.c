#include "smith_waterman.h"

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

#ifdef TIMING
extern double *sse_matrix_t, *sse_tracking_t;
#endif // TIMING

//------------------------------------------------------------------------------------

void init_subst_score_matrix(char *filename, subst_matrix_t matrix) {
  FILE *file = fopen(filename, "r");

  if (file == NULL) {
    printf("Error: substitution score matrix file (%s) not found\n", filename);
    exit(-1);
  }

  char *header[256],  *token[256];

  char *header_line = (char*) calloc(1, 4096);
  char *token_line = (char*) calloc(1, 4096);

  char *res = fgets(header_line, 4096, file);
  str_trim(header_line);


  res = NULL;

  // init matrix to -1000.0f
  for (int i = 0; i<128; i++) {
    for (int j = 0; j<128; j++) {
      matrix[i][j] = -1000.0f;
    }
  }
  /*
  matrix['A']['A'] = 5.0f; matrix['C']['A'] = -4.0f; matrix['T']['A'] = -4.0f; matrix['G']['A'] = -4.0f;
  matrix['A']['C'] = -4.0f; matrix['C']['C'] = 5.0f; matrix['T']['C'] = -4.0f; matrix['G']['C'] = -4.0f;
  matrix['A']['G'] = -4.0f; matrix['C']['T'] = -4.0f; matrix['T']['T'] = 5.0f; matrix['G']['T'] = -4.0f;
  matrix['A']['T'] = -4.0f; matrix['C']['G'] = -4.0f; matrix['T']['G'] = -4.0f; matrix['G']['G'] = 5.0f;
  */

  // read header row
  unsigned int num_columns = 0;

  header[num_columns] = strtok(header_line, "\t");
  //  res = strtok(header_line, "\t");
  while (header[num_columns]!= NULL) {
    num_columns++;
    header[num_columns] = strtok(NULL, "\t");
  }
  
  // read the remain rows and update matrix
  unsigned int col = 0;
  while (fgets(token_line, 4096, file) != NULL) {
    str_trim(token_line);
    col = 0;
    token[col] = strtok(token_line, "\t");
    while (token[col]!= NULL) {
      col++;
      token[col] = strtok(NULL, "\t");
    }
    if (col != num_columns) {
      printf("Error: substitution score matrix invalid format\n");
      exit(-1);
    }

    for(col = 1; col < num_columns; col++) {
      matrix[(unsigned char)header[col][0]][(unsigned char)token[0][0]] = atof(token[col]);
    }
  }

  free(header_line);
  free(token_line);

  fclose(file);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

sw_optarg_t* sw_optarg_new(float gap_open, float gap_extend, char *subst_matrix_name) {
    sw_optarg_t *optarg_p = calloc(1, sizeof(sw_optarg_t));

    optarg_p->gap_open = gap_open;
    optarg_p->gap_extend = gap_extend;
    
    init_subst_score_matrix(subst_matrix_name, optarg_p->subst_matrix);

    return optarg_p;
}

//------------------------------------------------------------------------------------

void sw_optarg_free(sw_optarg_t* optarg_p) {
    free(optarg_p);
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

sw_multi_output_t* sw_multi_output_new(unsigned int num_queries) {
    sw_multi_output_t *output_p = calloc(1, sizeof(sw_multi_output_t));
  
    output_p->num_queries = num_queries;
    output_p->query_map_p = (char**) calloc(num_queries, sizeof(char*));
    output_p->ref_map_p = (char**) calloc(num_queries, sizeof(char*));

    output_p->query_map_len_p = (unsigned int*) calloc(num_queries, sizeof(unsigned int));
    output_p->ref_map_len_p = (unsigned int*) calloc(num_queries, sizeof(unsigned int));
    
    output_p->query_start_p = (unsigned int*) calloc(num_queries, sizeof(unsigned int));
    output_p->ref_start_p = (unsigned int*) calloc(num_queries, sizeof(unsigned int));
    
    output_p->score_p = (float*) calloc(num_queries, sizeof(float));
    
    return output_p;
}

//------------------------------------------------------------------------------------

void sw_multi_output_free(sw_multi_output_t* output_p) {
    if (output_p == NULL) return;

    for (unsigned int i = 0; i < output_p->num_queries; i++) {
      if (output_p->query_map_p[i] != NULL) free(output_p->query_map_p[i]);
      if (output_p->ref_map_p[i] != NULL) free(output_p->ref_map_p[i]);
    }
    free(output_p->query_map_p);
    free(output_p->ref_map_p);

    if (output_p->query_map_len_p != NULL) free(output_p->query_map_len_p);
    if (output_p->ref_map_len_p != NULL) free(output_p->ref_map_len_p);
    
    if (output_p->query_start_p != NULL) free(output_p->query_start_p);
    if (output_p->ref_start_p != NULL) free(output_p->ref_start_p);
    
    if (output_p->score_p != NULL) free(output_p->score_p);
    
    free(output_p);
}

//------------------------------------------------------------------------------------
 
void sw_multi_output_save(int num_alignments, sw_multi_output_t* output_p, FILE *file_p) {
   unsigned int len, identity, gaps;

    if (file_p == NULL) {
      file_p = stdout;
    }
    
    for(int i = 0; i < num_alignments; i++) {
        gaps = 0;
        identity = 0;
        len = strlen(output_p->query_map_p[i]);

	fprintf(file_p, "Query: %s\tStart at %i\n", output_p->query_map_p[i], output_p->query_start_p[i]);
	fprintf(file_p,"       ");
	for(int j = 0; j < len; j++) {
	  if (output_p->query_map_p[i][j] == '-' || output_p->ref_map_p[i][j] == '-') {
	    gaps++;
	  }
	  if (output_p->query_map_p[i][j] == output_p->ref_map_p[i][j]) {
	    fprintf(file_p, "|");
	    identity++;
	  } else {
	    fprintf(file_p, "x");
	  }
	}
	fprintf(file_p, "\n");
	fprintf(file_p, "Ref. : %s\tStart at %i\n", output_p->ref_map_p[i], output_p->ref_start_p[i]);

	fprintf(file_p, "Score: %.2f\tLength: %i\tIdentity: %.2f\tGaps: %.2f\n", 
		output_p->score_p[i], len, identity * 100.0f / len, gaps * 100.0f / len);
	       
	fprintf(file_p, "\n");
    }
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void smith_waterman_mqmr(char **query_p, char **ref_p, unsigned int num_queries, 
			 sw_optarg_t *optarg_p, unsigned int num_threads, 
			 sw_multi_output_t *output_p) {
  
  const unsigned int simd_depth = 4;

  if (output_p == NULL) {
    printf("Error: output buffer is null.\n");
    exit(-1);
  }
  if (output_p->num_queries < num_queries) {
    printf("Error: num. sequencias (%i) and output buffer length (%i) mismatch !\n", 
	   num_queries, output_p->num_queries);
    exit(-1);
  }
  
  if (num_threads == 1) {

#ifdef TIMING
    double partial_t;
#endif // TIMING
    char *q_aux = NULL, *r_aux = NULL;
    int depth, aux_size = 0, H_size = 0, F_size = 0, max_q_len = 0, max_r_len = 0;
    char *q[simd_depth], *r[simd_depth];
    int len, index, q_lens[simd_depth], r_lens[simd_depth], alig_lens[simd_depth];
    float *H = NULL, *F = NULL;
    int *C = NULL;

    float gap_open = optarg_p->gap_open, gap_extend = optarg_p->gap_extend;
    float *score_p = output_p->score_p;

    depth = 0;
    for (unsigned int i = 0; i < num_queries; i++) {
      //printf("smith_waterman.c: query: #%i\n", i);
      len = strlen(query_p[i]);
      if (len > max_q_len) max_q_len = len;
      q_lens[depth] = len;
      q[depth] = query_p[i];

      len = strlen(ref_p[i]);
      if (len > max_r_len) max_r_len = len;
      r_lens[depth] = len;
      r[depth] = ref_p[i];

      depth++;
      if (depth == simd_depth) {

	index = i - depth + 1;

	reallocate_memory(max_q_len, max_r_len, simd_depth, &H_size, &H, &C, &F_size, &F, &aux_size, &q_aux, &r_aux);
	//printf("-----> max_q_len = %i, max_r_eln = %i, simd_depth = %i\n", max_q_len, max_r_len, simd_depth);
	//printf("-----> H_size = %i, F_size = %i, aux_size = %i\n", H_size, F_size, aux_size);
	//printf("num_queries = %i, gap_open = %0.2f, gap_extend = %0.2f\n", depth, gap_open, gap_extend);

	// generating score matrix
#ifdef TIMING
	partial_t = sw_tic();
#endif // TIMING
	sse_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len, 
		   optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, &score_p[index]);
#ifdef TIMING
	sse_matrix_t[0] += sw_toc(partial_t);
#endif // TIMING

	// tracebacking
#ifdef TIMING
	partial_t = sw_tic();
#endif // TIMING
	simd_traceback(simd_depth, depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
		       gap_open, gap_extend, H, C, &score_p[index],
		       &output_p->query_map_p[index], (int *)&output_p->query_start_p[index],
		       &output_p->ref_map_p[index], (int *)&output_p->ref_start_p[index], alig_lens,
		       q_aux, r_aux);
#ifdef TIMING
	sse_tracking_t[0] += sw_toc(partial_t);
#endif // TIMING
	
	depth = 0;
      }
    }

    //    printf("depth = %i\n", depth);

    if (depth > 0) {

      float max_score[simd_depth];
      
      for (unsigned int i = depth; i < simd_depth; i++) {
	q[i] = q[0];
	q_lens[i] = q_lens[0];

	r[i] = r[0];
	r_lens[i] = r_lens[0];
      }

      reallocate_memory(max_q_len, max_r_len, simd_depth, &H_size, &H, &C, &F_size, &F, &aux_size, &q_aux, &r_aux);
      
      index = num_queries - depth;

      // generating score matrix
#ifdef TIMING
      partial_t = sw_tic();
#endif // TIMING
      sse_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len, 
		 optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, max_score);
#ifdef TIMING
      sse_matrix_t[0] += sw_toc(partial_t);
#endif // TIMING

      // tracebacking
#ifdef TIMING
      partial_t = sw_tic();
#endif // TIMING
      simd_traceback(simd_depth, depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
		     gap_open, gap_extend, H, C, max_score,
		     &output_p->query_map_p[index], (int *)&output_p->query_start_p[index],
		     &output_p->ref_map_p[index], (int *)&output_p->ref_start_p[index], alig_lens,
		     q_aux, r_aux);
#ifdef TIMING
      sse_tracking_t[0] += sw_toc(partial_t);
#endif // TIMING

      for (unsigned int i = 0; i < depth; i++) {
	score_p[index +  i] = max_score[i];
      }
    }
    //    printf("Free 1\n");
    // free memory
    if (H != NULL) _mm_free(H);
    if (C != NULL) _mm_free(C);
    if (F != NULL) _mm_free(F);
    if (q_aux != NULL) free(q_aux);
    if (r_aux != NULL) free(r_aux);
    
  } else { // multi-thread

      unsigned int num_packs = num_queries / simd_depth;
      if (num_queries % simd_depth) {
        num_packs++;
      }
      unsigned int packs_per_thread = num_packs / num_threads;
      if (num_packs % num_threads) {
        packs_per_thread++;
      }

      //printf("num_packs = %i, packs_per_thread = %i\n", num_packs, packs_per_thread);
  
      #pragma omp parallel num_threads(num_threads) shared(num_packs, packs_per_thread, optarg_p, output_p)
      {
#ifdef TIMING
	double partial_t = 0.0;
#endif // TIMING
        unsigned int tid = omp_get_thread_num();
	
	unsigned int first_index = tid * packs_per_thread * simd_depth;
	unsigned int last_index = first_index + (packs_per_thread * simd_depth);
	if (last_index > num_queries) last_index = num_queries;

	//printf("tid = %i, from %i to %i\n", tid, first_index, last_index);

	char *q_aux = NULL, *r_aux = NULL;
	int depth, aux_size = 0, H_size = 0, F_size = 0, max_q_len = 0, max_r_len = 0;
	char *q[simd_depth], *r[simd_depth];
	int len, index, q_lens[simd_depth], r_lens[simd_depth], alig_lens[simd_depth];
	float *H = NULL, *F = NULL;
	int *C = NULL;
	
	float gap_open = optarg_p->gap_open, gap_extend = optarg_p->gap_extend;
	float *score_p = output_p->score_p;
	
	depth = 0;
	for (unsigned int i = first_index; i < last_index; i++) {
	  
	  len = strlen(query_p[i]);
	  if (len > max_q_len) max_q_len = len;
	  q_lens[depth] = len;
	  q[depth] = query_p[i];
	  
	  len = strlen(ref_p[i]);
	  if (len > max_r_len) max_r_len = len;
	  r_lens[depth] = len;
	  r[depth] = ref_p[i];
	  
	  depth++;
	  if (depth == simd_depth) {
	    
	    index = i - depth + 1;
	    
	    reallocate_memory(max_q_len, max_r_len, simd_depth, &H_size, &H, &C, &F_size, &F, &aux_size, &q_aux, &r_aux);
	    
	    // generating score matrix
#ifdef TIMING
	    partial_t = sw_tic();
#endif // TIMING
	    sse_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len, 
		       optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, &score_p[index]);
#ifdef TIMING
	    sse_matrix_t[tid] += sw_toc(partial_t);
#endif // TIMING
	    
	    // tracebacking
#ifdef TIMING
	    partial_t = sw_tic();
#endif // TIMING
	    simd_traceback(simd_depth, depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
			   gap_open, gap_extend, H, C, &score_p[index],
			   &output_p->query_map_p[index], (int *)&output_p->query_start_p[index],
			   &output_p->ref_map_p[index], (int *)&output_p->ref_start_p[index], alig_lens,
			   q_aux, r_aux);
#ifdef TIMING
	    sse_tracking_t[tid] += sw_toc(partial_t);
#endif // TIMING
	    
	    depth = 0;
	  }
	}
	
	//	printf("depth = %i\n", depth);
	
	if (depth > 0) {
	  
	  float max_score[simd_depth];
	  
	  for (unsigned int i = depth; i < simd_depth; i++) {
	    q[i] = q[0];
	    q_lens[i] = q_lens[0];
	    
	    r[i] = r[0];
	    r_lens[i] = r_lens[0];
	  }
	  
	  reallocate_memory(max_q_len, max_r_len, simd_depth, &H_size, &H, &C, &F_size, &F, &aux_size, &q_aux, &r_aux);
	  
	  index = num_queries - depth;
	  
	  // generating score matrix
#ifdef TIMING
	  partial_t = sw_tic();
#endif // TIMING
	  sse_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len, 
		     optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, max_score);
#ifdef TIMING
	  sse_matrix_t[tid] += sw_toc(partial_t);
#endif // TIMING
	  
	  // tracebacking
#ifdef TIMING
	  partial_t = sw_tic();
#endif // TIMING
	  simd_traceback(simd_depth, depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
			 gap_open, gap_extend, H, C, max_score,
			 &output_p->query_map_p[index], (int *)&output_p->query_start_p[index],
			 &output_p->ref_map_p[index], (int *)&output_p->ref_start_p[index], alig_lens,
			 q_aux, r_aux);
#ifdef TIMING
	  sse_tracking_t[tid] += sw_toc(partial_t);
#endif // TIMING
	  
	  for (unsigned int i = 0; i < depth; i++) {
	    score_p[index +  i] = max_score[i];
	  }
	}
	
	// free memory
	if (H != NULL) _mm_free(H);
	if (C != NULL) _mm_free(C);
	if (F != NULL) _mm_free(F);
	if (q_aux != NULL) free(q_aux);
	if (r_aux != NULL) free(r_aux);

      } // end #pragma omp parallel
  }
}

//------------------------------------------------------------------------------------

void smith_waterman_mqsr(char **query_p, char *ref_p, unsigned int num_queries, 
			 sw_optarg_t *optarg_p, unsigned int num_threads, 
			 sw_multi_output_t *output_p) {
  
  const unsigned int simd_depth = 4;

  if (output_p == NULL) {
    printf("Error: output buffer is null.\n");
    exit(-1);
  }
  if (output_p->num_queries < num_queries) {
    printf("Error: num. sequencias (%i) and output buffer length (%i) mismatch !\n", 
	   num_queries, output_p->num_queries);
    exit(-1);
  }

  int ref_len = strlen(ref_p);
  
  if (num_threads == 1) {



    char *q_aux = NULL, *r_aux = NULL;
    int depth, aux_size = 0, H_size = 0, F_size = 0, max_q_len = 0, max_r_len = ref_len;
    char *q[simd_depth], *r[simd_depth];
    int len, index, q_lens[simd_depth], r_lens[simd_depth], alig_lens[simd_depth];
    float *H = NULL, *F = NULL;
    int *C = NULL;

    float gap_open = optarg_p->gap_open, gap_extend = optarg_p->gap_extend;
    float *score_p = output_p->score_p;

    depth = 0;
    for (unsigned int i = 0; i < num_queries; i++) {

      len = strlen(query_p[i]);
      if (len > max_q_len) max_q_len = len;
      q_lens[depth] = len;
      q[depth] = query_p[i];

      r_lens[depth] = ref_len;
      r[depth] = ref_p;

      depth++;
      if (depth == simd_depth) {

	index = i - depth + 1;

	reallocate_memory(max_q_len, max_r_len, simd_depth, &H_size, &H, &C, &F_size, &F, &aux_size, &q_aux, &r_aux);
	//printf("-----> max_q_len = %i, max_r_eln = %i, simd_depth = %i\n", max_q_len, max_r_len, simd_depth);
	//printf("-----> H_size = %i, F_size = %i, aux_size = %i\n", H_size, F_size, aux_size);
	//printf("num_queries = %i, gap_open = %0.2f, gap_extend = %0.2f\n", depth, gap_open, gap_extend);

	// generating score matrix
#ifdef TIMING
	partial_t = sw_tic();
#endif // TIMING
	sse_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len, 
		   optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, &score_p[index]);
#ifdef TIMING
	sse_matrix_t[0] += sw_toc(partial_t);
#endif // TIMING

	// tracebacking
#ifdef TIMING
	partial_t = sw_tic();
#endif // TIMING
	simd_traceback(simd_depth, depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
		       gap_open, gap_extend, H, C, &score_p[index],
		       &output_p->query_map_p[index], (int *)&output_p->query_start_p[index],
		       &output_p->ref_map_p[index], (int *)&output_p->ref_start_p[index], alig_lens,
		       q_aux, r_aux);
#ifdef TIMING
	sse_tracking_t[0] += sw_toc(partial_t);
#endif // TIMING
	
	depth = 0;
      }
    }

    //    printf("depth = %i\n", depth);

    if (depth > 0) {

      float max_score[simd_depth];
      
      for (unsigned int i = depth; i < simd_depth; i++) {
	q[i] = q[0];
	q_lens[i] = q_lens[0];

	r[i] = r[0];
	r_lens[i] = r_lens[0];
      }

      reallocate_memory(max_q_len, max_r_len, simd_depth, &H_size, &H, &C, &F_size, &F, &aux_size, &q_aux, &r_aux);
      
      index = num_queries - depth;

      // generating score matrix
#ifdef TIMING
      partial_t = sw_tic();
#endif // TIMING
      sse_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len, 
		 optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, max_score);
#ifdef TIMING
      sse_matrix_t[0] += sw_toc(partial_t);
#endif // TIMING

      // tracebacking
#ifdef TIMING
      partial_t = sw_tic();
#endif // TIMING
      simd_traceback(simd_depth, depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
		     gap_open, gap_extend, H, C, max_score,
		     &output_p->query_map_p[index], (int *)&output_p->query_start_p[index],
		     &output_p->ref_map_p[index], (int *)&output_p->ref_start_p[index], alig_lens,
		     q_aux, r_aux);
#ifdef TIMING
      sse_tracking_t[0] += sw_toc(partial_t);
#endif // TIMING

      for (unsigned int i = 0; i < depth; i++) {
	score_p[index +  i] = max_score[i];
      }
    }

    // free memory
    if (H != NULL) _mm_free(H);
    if (C != NULL) _mm_free(C);
    if (F != NULL) _mm_free(F);
    if (q_aux != NULL) free(q_aux);
    if (r_aux != NULL) free(r_aux);
    
  } else { // multi-thread

      unsigned int num_packs = num_queries / simd_depth;
      if (num_queries % simd_depth) {
        num_packs++;
      }
      unsigned int packs_per_thread = num_packs / num_threads;
      if (num_packs % num_threads) {
        packs_per_thread++;
      }

      //printf("num_packs = %i, packs_per_thread = %i\n", num_packs, packs_per_thread);
  
      #pragma omp parallel num_threads(num_threads) shared(num_packs, packs_per_thread, optarg_p, output_p)
      {
        unsigned int tid = omp_get_thread_num();
	
	unsigned int first_index = tid * packs_per_thread * simd_depth;
	unsigned int last_index = first_index + (packs_per_thread * simd_depth);
	if (last_index > num_queries) last_index = num_queries;

	//printf("tid = %i, from %i to %i\n", tid, first_index, last_index);

	char *q_aux = NULL, *r_aux = NULL;
	int depth, aux_size = 0, H_size = 0, F_size = 0, max_q_len = 0, max_r_len = ref_len;
	char *q[simd_depth], *r[simd_depth];
	int len, index, q_lens[simd_depth], r_lens[simd_depth], alig_lens[simd_depth];
	float *H = NULL, *F = NULL;
	int *C = NULL;
	
	float gap_open = optarg_p->gap_open, gap_extend = optarg_p->gap_extend;
	float *score_p = output_p->score_p;
	
	depth = 0;
	for (unsigned int i = first_index; i < last_index; i++) {
	  
	  len = strlen(query_p[i]);
	  if (len > max_q_len) max_q_len = len;
	  q_lens[depth] = len;
	  q[depth] = query_p[i];
	  
	  r_lens[depth] = ref_len;
	  r[depth] = ref_p;
	  
	  depth++;
	  if (depth == simd_depth) {
	    
	    index = i - depth + 1;
	    
	    reallocate_memory(max_q_len, max_r_len, simd_depth, &H_size, &H, &C, &F_size, &F, &aux_size, &q_aux, &r_aux);
	    
	    // generating score matrix
#ifdef TIMING
	    partial_t = sw_tic();
#endif // TIMING
	    sse_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len, 
		       optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, &score_p[index]);
#ifdef TIMING
	    sse_matrix_t[tid] += sw_toc(partial_t);
#endif // TIMING
	    
	    // tracebacking
#ifdef TIMING
	    partial_t = sw_tic();
#endif // TIMING
	    simd_traceback(simd_depth, depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
			   gap_open, gap_extend, H, C, &score_p[index],
			   &output_p->query_map_p[index], (int *)&output_p->query_start_p[index],
			   &output_p->ref_map_p[index], (int *)&output_p->ref_start_p[index], alig_lens,
			   q_aux, r_aux);
#ifdef TIMING
	    sse_tracking_t[tid] += sw_toc(partial_t);
#endif // TIMING
	    
	    depth = 0;
	  }
	}
	
	//	printf("depth = %i\n", depth);
	
	if (depth > 0) {
	  
	  float max_score[simd_depth];
	  
	  for (unsigned int i = depth; i < simd_depth; i++) {
	    q[i] = q[0];
	    q_lens[i] = q_lens[0];
	    
	    r[i] = r[0];
	    r_lens[i] = r_lens[0];
	  }
	  
	  reallocate_memory(max_q_len, max_r_len, simd_depth, &H_size, &H, &C, &F_size, &F, &aux_size, &q_aux, &r_aux);
	  
	  index = num_queries - depth;
	  
	  // generating score matrix
#ifdef TIMING
	  partial_t = sw_tic();
#endif // TIMING
	  sse_matrix(depth, q, q_lens, max_q_len, r, r_lens, max_r_len, 
		     optarg_p->subst_matrix, gap_open, gap_extend, H, F, C, max_score);
#ifdef TIMING
	  sse_matrix_t[tid] += sw_toc(partial_t);
#endif // TIMING
	  
	  // tracebacking
#ifdef TIMING
	  partial_t = sw_tic();
#endif // TIMING
	  simd_traceback(simd_depth, depth, q, q_lens, max_q_len, r, r_lens, max_r_len,
			 gap_open, gap_extend, H, C, max_score,
			 &output_p->query_map_p[index], (int *)&output_p->query_start_p[index],
			 &output_p->ref_map_p[index], (int *)&output_p->ref_start_p[index], alig_lens,
			 q_aux, r_aux);
#ifdef TIMING
	  sse_tracking_t[tid] += sw_toc(partial_t);
#endif // TIMING
	  
	  for (unsigned int i = 0; i < depth; i++) {
	    score_p[index +  i] = max_score[i];
	  }
	}

	// free memory
	if (H != NULL) _mm_free(H);
	if (C != NULL) _mm_free(C);
	if (F != NULL) _mm_free(F);
	if (q_aux != NULL) free(q_aux);
	if (r_aux != NULL) free(r_aux);

      } // end #pragma omp parallel
  }
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void reallocate_memory(int max_q_len, int max_r_len, int simd_depth, 
		       int *H_size, float **H, int **C, int *F_size, float **F, 
		       int *aux_size, char **q_aux, char **r_aux) {
 
  int size_h, size_c, size_f;
  unsigned int matrix_size = max_q_len * max_r_len;
  
  size_h = simd_depth * matrix_size * sizeof(float);
  size_c = simd_depth * matrix_size * sizeof(int);
  if (matrix_size > *H_size) {
    if (*H_size > 0) {
      _mm_free(*H);
      _mm_free(*C);
    }

    *H = (float *) _mm_malloc(size_h, 16);
    *C = (int *) _mm_malloc(size_c, 16);
    //printf("new H %x, C %x\n", *H, *C);
    *H_size = matrix_size;
  }
  //  memset(*H, 0, size_h);
  //  memset(*C, 0, size_c);

  size_f = simd_depth * max_q_len * sizeof(float);
  if (max_q_len > *F_size) {
    if (*F_size > 0) _mm_free(*F);
    *F = (float *) _mm_malloc(size_f, 16);
    //printf("new F %x\n", *F);
    *F_size = max_q_len;
  }
  //  memset(*F, 0, size_f);

  int max_size = (max_r_len > max_q_len ? max_r_len : max_q_len);
  if (max_size > *aux_size) {
    if (*aux_size > 0) {
      free(*q_aux);
      free(*r_aux);
    }
    *q_aux = (char *) calloc(max_size * 2, sizeof(char));
    *r_aux = (char *) calloc(max_size * 2, sizeof(char));
    //printf("new q_aux %x, r_aux %x\n", *q_aux, *r_aux);
    *aux_size = max_size;
  }
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
