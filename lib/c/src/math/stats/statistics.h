#ifndef STATS_H
#define STATS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>

//#include "log.h"

/* **************************************
 *    		Structures  		*
 * *************************************/

/**
* @brief Structure for storing t test results
*
* Structure for storing t test results
*/
typedef struct vector {
    size_t length;   	/**< Length of the vector. */
    double* values;   	/**< Values in the vector. */
} vector_t;

/* **************************************
 *    		Functions  		*
 * *************************************/

/**
*  @brief Calculates mean of double values
*  @param values pointer to the double values
*  @param length num of values
*  @return mean value
*  
*  Calculates mean of double values
*/
double stats_mean(const double* values, size_t length);

/**
*  @brief Calculates variance of double values
*  @param values pointer to the double values
*  @param length num of values
*  @return variance value
*  
*  Calculates variance of double values
*/
double stats_variance(const double* values, size_t length);

/**
*  @brief Calculates median of double values
*  @param values pointer to the double values
*  @param length num of values
*  @return median value
*  
*  Calculates median of double values
*/
double stats_median(double* values, size_t length);

/**
*  @brief Calculates median of SORTED double values
*  @param values pointer to the double values
*  @param length num of values
*  @return median value
*  
*  Calculates median of SORTED double values
*/
double stats_median_sorted_values(const double* values, size_t length);

/**
*  @brief Calculates percentile of double values
*  @param values pointer to the double values
*  @param length num of values
*  @param percentile value of the percentile to calculate
*  @return value meeting the percentile
*  
*  Calculates percentile of double values
*/
double stats_percentile(double* values, size_t length, double percentile);

/**
*  @brief Calculates percentile of SORTED double values
*  @param values pointer to the double values
*  @param length num of values
*  @param percentile value of the percentile to calculate
*  @return value meeting the percentile
*  
*  Calculates percentile of SORTED double values
*/
double stats_percentile_sorted_values(const double* values, size_t length, double percentile);

/**
*  @brief Creates a new vector from a double pointer
*  @param values pointer to the double values
*  @param length length of the vector
*  @return created vector 
*  
*  Creates a new vector from a double pointer and its length
*/
vector_t* vector_new(double *values, size_t length);

/**
*  @brief Frees a vector
*  @param v pointer to the vector structure
*  @return void
*  
*  Frees a vector
*/
void vector_free(vector_t* v);

/**
*  @brief Sorts double values of a vector
*  @param v pointer to the vector structure
*  @param length length of the vector
*  @return sorted vector 
*  
*  Sorts double values of a vector
*/
vector_t* vector_sort(vector_t* v);


#endif  /* STATS_H */
