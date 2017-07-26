#ifndef P_ADJUST_H
#define P_ADJUST_H

#include "data/array_utils.h"
#include "data/R_utils.h"

int FDR_correction(double *pvalues, size_t n);

int Bonferroni_correction(double *pvalues, size_t n);

int BH_correction(double *pvalues, size_t n);

int Hochberg_correction(double *pvalues, size_t n);

int BY_correction(double *pvalues, size_t n);

int Holm_correction(double *pvalues, size_t n);

#endif  /* P_ADJUST_H */
