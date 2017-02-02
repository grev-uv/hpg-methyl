#ifndef SSE_H
#define SSE_H

#include <immintrin.h>

#include "macros.h"

//------------------------------------------------------------------------
// SSE functions
//------------------------------------------------------------------------

void sse_matrix(int num_seqs, 
		char **q, int *q_len, int max_q_len,
		char **r, int *r_len, int max_r_len,
		float profile[128][128], float gap_open, float gap_extend,
		float *H, float *F, int *C, float *max_score);


//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif // SSE_H
