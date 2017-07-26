
#ifndef COX_H
#define COX_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "data/matrix.h"

#define COX_ERROR_LIMIT		0.00001

/* **************************************
 *    		Structures  		*
 * *************************************/

/* **************************************
 *    		Functions  		*
 * *************************************/

/**
*  @brief Computes Cox test for a given covariates
*  @param covariates pointer to matrix of covariates
*  @param num_features_in_covariate number of features (rows) in the covariates matrix
*  @param num_samples number of samples (columns) in the dataset
*  @param time pointer to the array with time variable values
*  @param censor pointer to the array with censor variable values
*  @param coefficients[in,out] pointer to coefficient array
*  @param variance[in,out] pointer variance matrix
*  @return void
*  
*  Computes Cox test for a given covariates (usually with 1 row, but not necessarily)
*/
void cox_test(double** covariates, size_t num_features_in_covariate, size_t num_samples, double* time, double* censor, double* coefficients, double** variance);

#endif /* COX_H */
