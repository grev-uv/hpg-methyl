#ifndef MACROS_H
#define MACROS_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <string.h>
#include <math.h>

//------------------------------------------------------------------------
//------------------------------------------------------------------------
#ifndef MAX
  #define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef MIN
  #define MIN(a,b) ((a) > (b) ? (b) : (a))
#endif

#ifndef NORM_SCORE
  #define NORM_SCORE(s,l,m) ((s) / ((l) * (m)))
#endif

#define ZERO 0
#define LEFT 1
#define DOWN 2
#define DIAGONAL 3

#define U_FEPS 1.192e-6F         // 1.0F + E_FEPS != 1.0F
#define E_FPEQ(a,b,e) (((b - e) < a) && (a < (b + e)))

#ifndef FLT_MAX
   #define FLT_MAX 1000000.0f
#endif

//------------------------------------------------------------------------
// utility functions
//------------------------------------------------------------------------

void str_trim(char *line);

//------------------------------------------------------------------------

double sw_tic();
double sw_toc(double t1);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

void simd_traceback(int depth, int num_seqs, 
		    char **q, int *q_len, int max_q_len,
		    char **r, int *r_len, int max_r_len,
		    float gap_open, float gap_extend,
		    float *H, int *C, 
		    float *max_score,
		    char **q_alig, int *q_start,
		    char **r_alig, int *r_start, 
		    int *len_alig,
		    char *q_aux, char *r_aux);

//------------------------------------------------------------------------

void simd_find_position(int depth, int index, char *q, int q_len, char *r, int r_len,
			float *H, int cols, int rows, float score, 
			int *q_pos, int *r_pos);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

void save_float_matrix(float *matrix, 
		       int cols, int rows,
		       char *q, int q_len, 
		       char *r, int r_len, 
		       int index, int depth,
		       char *filename);

void save_int_matrix(int *matrix, 
		     char *q, int q_len, 
		     char *r, int r_len, 
		     int index, int depth,
		     char *filename);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif // MACROS_H
