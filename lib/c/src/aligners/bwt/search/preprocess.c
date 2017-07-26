#include "preprocess.h"

void calculate_and_save_B(ref_vector *X, const char *directory, const char *name) {
  direct_bwt(X->vector, X->n, directory, name, true);
}

void calculate_S_and_R(comp_vector *S, comp_vector *R, ref_vector *B, vector *C, comp_matrix *O, SA_TYPE ratio) {
 
  SA_TYPE i,j;
  uint8_t b_aux;

  S->siz = B->n + 1;
  R->siz = B->n + 1;

  S->n = (S->siz / ratio);
  R->n = (S->siz / ratio);
  if (S->siz % ratio) {
    S->n++; R->n++;
  }

  S->ratio = ratio;
  R->ratio = ratio;

  S->vector = (SA_TYPE *) malloc(S->n * sizeof(SA_TYPE));
  R->vector = (SA_TYPE *) malloc(R->n * sizeof(SA_TYPE));

  if (B->dollar % S->ratio == 0) S->vector[B->dollar / S->ratio] = 0;	
  R->vector[0] = B->dollar;

  for(j = S->siz-1, i = 0; j > 0; j--) {
	
    if (i % S->ratio == 0) S->vector[i / S->ratio] = j;
    if (j % R->ratio == 0) R->vector[j / S->ratio] = i;
	  
    if (i > B->dollar) b_aux = B->vector[i-1];
    else               b_aux = B->vector[i];

    i = C->vector[b_aux] + get_O(b_aux, i+1/*0 is -1*/, O);

  }

}

void calculate_C(vector *C, vector *C1, ref_vector *B, uint8_t nA) {

  C->n  = nA;
  C1->n = nA;

  C->vector  = (SA_TYPE *) malloc(C->n  * sizeof(SA_TYPE));
  C1->vector = (SA_TYPE *) malloc(C1->n * sizeof(SA_TYPE));

  check_malloc(C->vector,  "calculateC");
  check_malloc(C1->vector, "calculateC");

  for (SA_TYPE i = 0; i < C->n; i++)
    C->vector[i] = 0;

  SA_TYPE dollar = 0;

  for (SA_TYPE i = 0; i<=B->n; i++) {
    if (i == B->dollar) {dollar = 1; continue;}
    if (B->vector[i - dollar] + 1 == nA) continue;
    C->vector[B->vector[i - dollar] + 1]++;
  }

  for (SA_TYPE i = 1; i < C->n; i++)
    C->vector[i] += C->vector[i-1];

  for (SA_TYPE i = 0; i < C->n; i++)
    C1->vector[i] = C->vector[i] + 1;

}

void calculate_O(comp_matrix *O, ref_vector *B, uint8_t nA) {

#if defined FM_COMP_32 || FM_COMP_64

	O->siz = B->n+2;          // Position 0 is -1 and the dollar, so I add two elements

	O->n_desp = nA;
	O->m_desp = (O->siz / FM_COMP_VALUE);
	if (O->siz % FM_COMP_VALUE) O->m_desp++;

	O->n_count = nA;
	O->m_count = O->m_desp;

	O->desp  = (SA_TYPE **) malloc(O->n_desp * sizeof(SA_TYPE *));
	check_malloc(O->desp, "calculateO");
	O->count = (FM_COMP_TYPE **) malloc(O->n_count * sizeof(FM_COMP_TYPE *));
	check_malloc(O->count, "calculateO");

#pragma omp parallel for
	for (SA_TYPE i=0; i<O->n_count; i++) {

		O->desp[i]  = (SA_TYPE *) malloc(O->m_desp * sizeof(SA_TYPE));
		check_malloc(O->desp[i], "calculateO");
		O->count[i] = (FM_COMP_TYPE *) malloc(O->m_count * sizeof(FM_COMP_TYPE));
		check_malloc(O->count[i], "calculateO");

		SA_TYPE k=0;
		SA_TYPE pos=0;
		SA_TYPE bit;

		O->desp[i][pos]  = 0;
		O->count[i][pos] = 0;

		/* printf("%d", 0); */

		SA_TYPE dollar = 0;

		for (SA_TYPE j=0; j<=B->n; j++) {

			bit = (j + 1) % FM_COMP_VALUE; //First column is -1 index

			if (!bit) {
				pos++;
				O->desp[i][pos]  = k;
				O->count[i][pos] = 0;          //Initialize to 0 bit-vector
			}

			if (j == B->dollar) { dollar = 1; continue;	}

			if (B->vector[j - dollar] == i) {
				k++;
				O->count[i][pos] |= ((FM_COMP_TYPE) 1) << bit;
			}

		}

	}

#else

	O->siz    = B->n+1;
	O->n_desp = nA;
	O->m_desp = O->siz;       // Position 0 is -1, so I add one element

	O->desp  = (SA_TYPE **) malloc(O->n_desp * sizeof(SA_TYPE *));
	check_malloc(O->desp, "calculateO");

#pragma omp parallel for
	for (SA_TYPE i=0; i<O->n_desp; i++) {

		O->desp[i]  = (SA_TYPE *) malloc(O->m_desp * sizeof(SA_TYPE));
		check_malloc(O->desp[i], "calculateO");

		SA_TYPE k=0;

		O->desp[i][0] = 0;     // First column is -1 index

		SA_TYPE dollar = 0;

		for (SA_TYPE j=0; j<=B->n; j++) {

			if (j == B->dollar) {
				dollar = 1;
			} else if (((SA_TYPE)B->vector[j - dollar]) == i) {
				k++;
			}

			O->desp[i][j + 1] = k; // First column is -1 index

		}

	}

#endif

}

void compress_S_or_R(comp_vector *SR, comp_vector *SRcomp, SA_TYPE ratio) {

	if (SR->ratio != 1) {
		fprintf(stderr, "compress_SR: The S or R vectors must be uncompressed\n");
		exit(1);
	}

	SRcomp->n = (SR->siz / ratio);
	if (SR->siz % ratio) SRcomp->n++;

	SRcomp->siz = SR->siz;
	SRcomp->ratio = ratio;

	SRcomp->vector = (SA_TYPE *) malloc(SRcomp->n * sizeof(SA_TYPE));
	check_malloc(SRcomp->vector, "calculateO");

#pragma omp parallel for
	for (SA_TYPE i=0; i < SRcomp->siz; i+=ratio)
		SRcomp->vector[i/ratio] = SR->vector[i];

}
