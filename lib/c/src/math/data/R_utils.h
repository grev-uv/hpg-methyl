#ifndef R_UTILS_H
#define R_UTILS_H

#include <stdlib.h>
#include <math.h>


double *cummin(const double *data, size_t n, double *res);

double *cummax(const double *data, size_t n, double *res);

double *cumsum(const double *data, size_t n, double *res);

double *pmin(const double *data, size_t n, double min, double *res);

double *pmax(const double *data, size_t n, double max, double *res);


#endif	/* R_UTILS_H */
