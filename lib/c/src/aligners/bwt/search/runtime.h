#ifndef _SEARCH_RUNTIME_
#define _SEARCH_RUNTIME_

#if defined CSALIB_SEARCH
#include "../csalib/csa.h"
#else
#include "csafm.h"
#endif

#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct {

#if defined CSALIB_SEARCH
	CSA csa;
#else
	vector C, C1;
	comp_matrix O;
	comp_vector S, R;
#endif

} bwt_index;

#if defined CSALIB_SEARCH

#define BWiteration(k_in,l_in, k_out, l_out, b, index)\
	do {\
		(k_out) = k_in;\
		(l_out) = l_in;\
		(index)->csa.searchsub((b), &((index)->csa), &(k_out), &(l_out));\
	} while (0);
//printf("k-> %lu, l-> %lu\n", (k_out), (l_out));

#else

#define BWiteration(k_in,l_in, k_out, l_out, b, index)\
	do {\
		(k_out) = (index)->C1.vector[(b)] + get_O((b), (k_in) , &((index)->O));\
		(l_out) = (index)->C.vector[(b)] + get_O((b), (l_in)+1, &((index)->O));\
	} while (0);
#endif
//printf("k-> %lu, l-> %lu, C -> %u, C1 -> %u\n", (k_out), (l_out), get_O((b), (k_in), &((index)->O)), get_O((b), (l_in)+1, &((index)->O)), (index)->C.vector[(b)], (index)->C1.vector[(b)]);

#if defined CSALIB_SEARCH

#define size_SA(index) ((index)->csa.n+1)
#define get_SA(m, index) (index)->csa.lookup(&((index)->csa), (m))
#define get_ISA(m, index) (index)->csa.inverse(&((index)->csa), (m))

static inline void load_bwt_index(bwt_index *index_rev, bwt_index *index, const char *directory, int direction, bool inverse_sa, bwt_config_t config) {

	char *fname[2];

	fname[0] = (char *) malloc(500 * sizeof(char));
	fname[1] = (char *) malloc(500 * sizeof(char));
	fname[0][0]='\0';
	fname[1][0]='\0';

	if (direction) {
		strcat(fname[0], directory);
		strcat(fname[0], "/backward.bwd");
		strcat(fname[1], directory);
		strcat(fname[1], "/backward.idx");
	} else {
		strcat(fname[0], directory);
		strcat(fname[0], "/forward.bwd");
		strcat(fname[1], directory);
		strcat(fname[1], "/forward.idx");
	}

	csa_read(&(index->csa), 2, fname);

	free(fname[0]);
	free(fname[1]);

}

static inline void free_bwt_index(bwt_index *index_rev, bwt_index *index, bool inverse_sa) {
	
}

#else

#define size_SA(index) ((index)->S.siz)
#define get_SA(m, index) getScompValue((m), &((index)->S), &((index)->C), &((index)->O))
#define get_ISA(m, index) getRcompValue((m), &((index)->R), &((index)->C), &((index)->O))

static inline void load_bwt_index(bwt_index *index_rev, bwt_index *index, const char *directory, int direction, bool inverse_sa, bwt_config_t config) {

	if (direction) {
	  read_vector(&(index->C),      directory, "C");
	  read_vector(&(index->C1),     directory, "C1");
	  read_comp_matrix(&(index->O), directory, "O");
	  read_comp_vector(&(index->S), directory, "S");
	  if (inverse_sa) read_comp_vector(&(index->R), directory, "R");
	} else {
	  read_vector(&(index->C),      directory, "C");
	  read_vector(&(index->C1),     directory, "C1");
	  read_comp_matrix(&(index->O), directory, "Oi");
	  read_comp_vector(&(index->S), directory, "Si");
	  if (inverse_sa) read_comp_vector(&(index->R), directory, "Ri");
	}
	
	if (index_rev != NULL) {
	  reverse_strand_C(&(index_rev->C), &(index->C), &(index_rev->C1), &(index->C1), config);
	  reverse_strand_O(&(index_rev->O), &(index->O), config);
	  
	  index_rev->S.vector = index->S.vector;
	  index_rev->S.siz = index->S.siz;
	  index_rev->S.n = index->S.n;
	  index_rev->S.ratio = index->S.ratio;

	  if (inverse_sa) {
	    index_rev->R.vector = index->R.vector;
	    index_rev->R.siz = index->R.siz;
	    index_rev->R.n = index->R.n;
	    index_rev->R.ratio = index->R.ratio;
	  }

	}
	
}

static inline void free_bwt_index(bwt_index *index_rev, bwt_index *index, bool inverse_sa) {

	free(index->C.vector);
	free(index->C1.vector);

	if (index_rev != NULL) {
		free(index_rev->C.vector);
		free(index_rev->C1.vector);
		free_comp_matrix(&(index_rev->O), &(index->O));
	} else {
		free_comp_matrix(NULL, &(index->O));
	}

	free(index->S.vector);
	if (inverse_sa) free(index->R.vector);

}

#endif

#endif
