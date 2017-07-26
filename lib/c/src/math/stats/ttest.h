
#ifndef TTEST_H
#define TTEST_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#define P_VALUE_THRESHOLD	0.05

/* **************************************
 *    		Structures  		*
 * *************************************/

/**
* @brief Structure for storing t test results
*
* Structure for storing t test results
*/

typedef struct ttest_result {
    double statistic;   	/**< Statistic of the distribution. */
    double degrees_of_freedom;  /**< Degrees of freedom. */
    double p_value;		/**< p value. */
} ttest_result_t;

/* **************************************
 *    		Functions  		*
 * *************************************/

/**
*  @brief Creates a ttest_result array
*  @param num_rows number of rows (number of features)
*  @return ttest_result array
*
*  Creates and returns a ttest_result array
*/
ttest_result_t** ttest_results_new(int num_rows);

/**
*  @brief Frees a given ttest_result array
*  @param num_rows number of rows (number of features)
*  @param[in,out] ttest_results pointer to ttest_result array
*  @return void
*
*  Frees a given ttest_result array
*/
void ttest_results_free(int num_rows, ttest_result_t** ttest_results);

/**
*  @brief Gets a double statistics vector from a collection of ttest_result
*  @param num_rows number of rows (number of features)
*  @param anova_results pointer to ttest_result array
*  @param[in,out] statistics pointer to the statistics vector
*  @return void
*
*  Gets a double statistics vector from a collection of ttest_result
*/
double* ttest_results_get_statistics(int num_rows, ttest_result_t** ttest_results, double* statistics);

/**
*  @brief Gets a double p_values vector from a collection of ttest_result
*  @param num_rows number of rows (number of features)
*  @param ttest_results pointer to ttest_result array
*  @param[in,out] p_values pointer to the p_values vector
*  @return void
*
*  Gets a double p_values vector from a collection of ttest_result
*/
double* ttest_results_get_p_values(int num_rows, ttest_result_t** ttest_results, double* p_values);




/**
*  @brief Computes approximate degrees of freedom for 2-sample t-test
*  @param var1 first sample variance
*  @param var2 second sample variance
*  @param n1 first sample n
*  @param n2 second sample n
*  @return approximate degrees of freedom
*  
*  Computes approximate degrees of freedom for 2-sample t-test
*/
double degrees_of_freedom(double var1, double var2, double n1, double n2);

/**
*  @brief Computes t test statistic for 1-sample t-test
*  @param mean sample mean
*  @param mu constant to test again
*  @param var sample variance
*  @param n sample n
*  @return t test statistic
*  
*  Computes t test statistic for 1-sample t-test
*/
double ttest_one_sample(double mean, double mu, double var, double n);

/**
*  @brief Computes p-value for 2-sample t-test
*  @param mean1 first sample mean
*  @param mean2 second sample mean
*  @param var1 first sample variance
*  @param var2 second sample variance
*  @param n1 first sample n
*  @param n2 second sample n
*  @return t test statistic
*  
*  Computes p-value for 2-sided, 2-sample t-test.
*/
double ttest_two_samples(double mean1, double mean2, double var1, double var2, double n1, double n2);

/**
*  @brief Computes p-value for 2-sided, 1-sample t-test
*  @param mean sample mean
*  @param mu constant to test again
*  @param var sample variance
*  @param n sample n
*  @param ttest_results pointer to t test results structures
*  @return pointer to t test result structure (statistic, degrees of freedom and p-value)
*  
*  Computes p-value for 2-sided, 1-sample t-test
*/
ttest_result_t* ttest_two_sided_one_sample(double mean, double mu, double var, double n, ttest_result_t* ttest_result);

/**
*  @brief Computes p-value for 2-sided, 2-sample t-test
*  @param mean1 first sample mean
*  @param mean2 second sample mean
*  @param var1 first sample variance
*  @param var2 second sample variance
*  @param n1 first sample n
*  @param n2 second sample n
*  @param ttest_results pointer to t test results structures
*  @return pointer to t test result structure (statistic, degrees of freedom and p-value)
*  
*  Computes p-value for 2-sided, 2-sample t-test.
*/
ttest_result_t* ttest_two_sided_two_samples(double mean1, double mean2, double var1, double var2, double n1, double n2, ttest_result_t* ttest_result);

#endif /* TTEST_H */
