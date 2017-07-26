#ifndef FISHER_H
#define FISHER_H

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include <immintrin.h>

enum Fisher_mode { LESS = 1, GREATER, TWO_SIDED };

double fisher_test(int a, int b, int c, int d, enum Fisher_mode, double *factorial_logarithms);

double *fisher_test_vectorized(int *a, int *b, int *c, int *d, int len, enum Fisher_mode mode, double *factorial_logarithms);

double *fisher_test_omp(int *a, int *b, int *c, int *d, int len, enum Fisher_mode mode, double *factorial_logarithms);


double *init_logarithm_array(int nmax);

#endif
