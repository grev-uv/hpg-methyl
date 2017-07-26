
#ifndef CORRELATION_H
#define CORRELATION_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_cdf.h>

/* **************************************
 *    		Structures  		*
 * *************************************/

/* **************************************
 *    		Functions  		*
 * *************************************/

/**
*  @brief Computes correlation for two arrays of double 
*  @param values1 first array of data
*  @param values2 second array of data
*  @param size size of arrays
*  @return correlation value
*  
*  Computes correlation for two arrays of double
*/
double pearson_correlation(double* values1, double* values2, double size);

#endif /* CORRELATION_H */
