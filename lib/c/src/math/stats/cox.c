
#include "cox.h"

/* **************************************************************
 *    		Public functions implementations 		*
 * *************************************************************/

void cox_test(double** covariates, size_t num_features_in_covariate, size_t num_samples, double* time, double* censor, double* coefficients, double** variance) {
    // declare variables, init values and allocate memory
    // gsl matrices
    matrix_t* coefficients_matrix_p =  NULL;
    matrix_t* information_matrix_p =  NULL;
    matrix_t* information_matrix_inverse_p =  NULL;
    matrix_t* score_matrix_p =  NULL;
    matrix_t* error_matrix_p =  NULL;
    matrix_t* variance_matrix_p =  NULL;

    // other variables
    double denominator = 0, numerator = 0;
    double error1 = 1, error2 = 1;
    double* risk_factor = (double*) calloc(num_samples, sizeof(double));
    double* score = (double*) calloc(num_features_in_covariate, sizeof(double));
    double** expected_covariate = (double**) calloc(num_features_in_covariate, sizeof(double*));
    double** information = (double**) calloc(num_features_in_covariate, sizeof(double*));
    
    for (size_t i = 0; i < num_features_in_covariate; i++) {
        coefficients[i] = 0.1;
        expected_covariate[i] = (double*) calloc(num_samples, sizeof(double));
        information[i] = (double*) calloc(num_features_in_covariate, sizeof(double));
    }

    // create gsl matrices
    coefficients_matrix_p =  matrix_new(num_features_in_covariate, 1);
    matrix_init(0.1, coefficients_matrix_p);
    information_matrix_p = matrix_new(num_features_in_covariate, num_features_in_covariate);
    score_matrix_p = matrix_new(num_features_in_covariate, 1);
    error_matrix_p = matrix_new(num_features_in_covariate, 1);
    information_matrix_inverse_p = matrix_new(num_features_in_covariate, num_features_in_covariate);

    while((error1 > COX_ERROR_LIMIT) || (error2 > COX_ERROR_LIMIT)) {
        for (size_t i = 0; i < num_samples; ++i) {
            risk_factor[i] = 1.0;
            for (size_t s = 0; s < num_features_in_covariate; ++s) {
                risk_factor[i] *= exp(coefficients[s] * covariates[s][i]);
            }
        }

        for (size_t j = 0; j < num_features_in_covariate; j++) {
            score[j] = 0.0;
            for (size_t i = 0; i < num_samples; i++) {
                for (size_t k = 0; k < num_samples; k++) {
                    if (time[k] >= time[i]) {
                        denominator += risk_factor[k];
                        numerator += covariates[j][k] * risk_factor[k];
                    }
                }
                 
                expected_covariate[j][i] = numerator / denominator;
                score[j] += censor[i] * (covariates[j][i] - expected_covariate[j][i]);

                numerator = 0.0;
                denominator = 0.0;
            }
        }

        for (size_t r = 0; r < num_features_in_covariate; r++) {
            for (size_t s = 0; s < num_features_in_covariate; s++) {
                information[r][s] = 0.0;
                for (size_t i = 0; i < num_samples; i++) {
                    for (size_t k = 0; k < num_samples; k++) {
                        if (time[k] >= time[i]) {
                            denominator += risk_factor[k];
                            numerator += (covariates[r][k] * covariates[s][k] * risk_factor[k]);
                        }
                    }
                    information[r][s] +=  censor[i] * (expected_covariate[r][i] * expected_covariate[s][i] - (numerator / denominator));

                    numerator = 0.0;
                    denominator = 0.0;
                }
            }
        }

        // fill information_matrix
        matrix_fill(information, num_features_in_covariate, num_features_in_covariate, information_matrix_p);  // fill the matrix with data

        // fill score_matrix from score array
        for (size_t i = 0; i < num_features_in_covariate; i++) {
            matrix_set(i, 0, score[i], score_matrix_p);
        }

        // calculate error matrix: inv(information_matrix) * score_matrix
        matrix_inv(information_matrix_p, information_matrix_inverse_p);
        matrix_mul(information_matrix_inverse_p, score_matrix_p, error_matrix_p);

        // calculate coefficients matrix
        coefficients_matrix_p = matrix_sub(coefficients_matrix_p, error_matrix_p);

        // fill coefficientes
        for (size_t i = 0; i < num_features_in_covariate; i++) {
            coefficients[i] = matrix_get(i, 0, coefficients_matrix_p);
        }

        error1 = sqrt(matrix_Fnorm(error_matrix_p));
        error2 = sqrt(matrix_Fnorm(score_matrix_p));
    }  // end of while
    
    // calculate variance: (-1 * inv(information_matrix))
    variance_matrix_p = matrix_scale(information_matrix_inverse_p, -1.0);
    for (size_t i = 0; i < num_features_in_covariate; i++) {
        for (size_t j = 0; j < num_features_in_covariate; j++) {
            variance[i][j] = matrix_get(i, j, variance_matrix_p);
        }
    }  

    // free gsl matrices
    matrix_free(coefficients_matrix_p);
    matrix_free(information_matrix_p);
    matrix_free(information_matrix_inverse_p);
    matrix_free(score_matrix_p);
    matrix_free(error_matrix_p);
    variance_matrix_p = NULL;  // points to information_matrix_inverse_p previously freed
    matrix_free(variance_matrix_p); 

    // free other resources
    free(risk_factor);
    free(score);
    for (size_t i = 0; i < num_features_in_covariate; i++) {
        free(expected_covariate[i]);
        free(information[i]);
    }    
    free(expected_covariate);
    free(information); 
    
    return;
}
