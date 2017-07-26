
#include "statistics.h"

/* **************************************************************
 *    		Public functions implementations 		*
 * *************************************************************/

/* **********************************************
 *    		Statistic functions  		*
 * *********************************************/

double stats_mean(const double* values, size_t length) {
    return gsl_stats_mean(values, 1, length);
}

double stats_variance(const double* values, size_t length) {
    return gsl_stats_variance(values, 1, length);
}

double stats_median(double* values, size_t length) {
    gsl_sort(values, 1, length);
    return stats_median_sorted_values(values, length);
}

double stats_median_sorted_values(const double* values, size_t length) {
    return gsl_stats_median_from_sorted_data(values, 1, length);
}

double stats_percentile(double* values, size_t length, double percentile) {
    gsl_sort(values, 1, length);    
    return stats_percentile_sorted_values(values, length, percentile);
}

double stats_percentile_sorted_values(const double* values, size_t length, double percentile) {
    assert ((percentile < 100) && (percentile > 0));
      
    return gsl_stats_quantile_from_sorted_data(values, 1, length, (percentile / 100));  
}

/* **********************************************
 *    		Vector functions  		*
 * *********************************************/

vector_t* vector_new(double *values, size_t length) {
    vector_t* vector_p = (vector_t*) calloc(1, sizeof(vector_t));
    
    vector_p->length = length;
    vector_p->values = (double*) calloc(length, sizeof(double));
    memcpy(vector_p->values, values, length * sizeof(double));
    
    return vector_p;
}

void vector_free(vector_t* v) {
    if (v == NULL) {
        return;
    }
    
    if (v->values != NULL) {
        free(v->values);
    }    
    free(v);
    
    return;
}

vector_t* vector_sort(vector_t* v) {
    gsl_sort(v->values, 1, v->length);
    return v;    
}
