#ifndef EMBOSS_H
#define EMBOSS_H

#include "macros.h"

//------------------------------------------------------------------------
// EMBOSS functions
//------------------------------------------------------------------------

float sw_emboss(char *seq_a, char *seq_b, float gapopen, float gapextend, 
		char *m, char *n, int *start1, int *start2,
		double *matrix_time, double *tracking_time, double *total_time);

float AlignPathCalcSW(const char *a, const char *b, int lena, int lenb,
                      float gapopen, float gapextend, float* path,
                      int* compass);

void AlignWalkSWMatrix(const float* path, const int* compass,
		       float gapopen, float gapextend,
		       const char*  a, const char* b,
		       char* m, char* n,
		       int lena, int lenb,
		       int *start1, int *start2);
void revstr(char *str);

#endif // EMBOSS_H
