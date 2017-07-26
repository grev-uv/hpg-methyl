#include "fitting.h"

int linear_regression(const double *x, const size_t xstride, const double *y, const size_t ystride, const size_t n, 
                      double *c0, double *c1, double *cov00, double *cov01, double *cov11, double *sumsq) {
    return gsl_fit_linear(x, xstride, y, ystride, n, c0, c1, cov00, cov01, cov11, sumsq);
}
