#include "R_utils.h"


double *cummin(const double *data, size_t n, double *res) {
    double min = INFINITY;
    size_t i;
    for(i = 0 ; i < n; i++ ) {
        if (isnan(data[i]) || isnan(min)) {
            min = min + data[i];  /* propagate NA and NaN */
        }else {
            min = (min < data[i]) ? min : data[i];
            res[i] = min;
        }
    }
    return res;
}


double *cummax(const double *data, size_t n, double *res) {
    double max = INFINITY;
    size_t i;
    for(i = 0 ; i < n; i++ ) {
        if(isnan(data[i]) || isnan(max)) {
            max = max + data[i];  /* propagate NA and NaN */
        }else{
            max = (max > data[i]) ? max : data[i];
            res[i] = max;
        }
    }
    return res;
}


double *cumsum(const double *data, size_t n, double *res) {
    double sum = 0;
    size_t i;
    for(i=0; i < n; i++) {
        sum += data[i];
        res[i] = sum;
    }
    return res;
}


double *pmin(const double *data, size_t n, double min, double *res) {
    size_t i;
    for(i = 0; i < n; i++) {
        res[i] = (min < data[i]) ? min : data[i];
    }
    return res;
}


double *pmax(const double *data, size_t n, double max, double *res) {
    size_t i;
    for(i = 0; i < n; i++) {
        res[i] = (max > data[i]) ? max : data[i];
    }
    return res;
}
