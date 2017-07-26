#ifndef LS_FIT_H
#define LS_FIT_H

#include <sys/types.h>

#include <gsl/gsl_fit.h>


int linear_regression(const double *x, const size_t xstride, const double *y, const size_t ystride, const size_t n, 
                      double *c0, double *c1, double *cov00, double *cov01, double *cov11, double *sumsq);


#endif