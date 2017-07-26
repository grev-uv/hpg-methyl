
#ifndef CORRELATION_MATRIX_H
#define CORRELATION_MATRIX_H

/* **************************************
 *    		Structures  		*
 * *************************************/

/**
* @brief Matrix for storing correlation results
*
* Matrix for storing correlation results
*/
typedef struct correlation_matrix {
    float** values;   	/**< Values of correlation. */
} correlation_matrix_t;

/* **************************************
 *    		Functions  		*
 * *************************************/

/**
*  @brief Creates a correlation matrix
*  @param num_rows number of rows and cols of the square matrix
*  @return pointer to the correlation matrix
*  
*  Creates a correlation matrix
*/
correlation_matrix_t* correlation_matrix_new(size_t num_rows);

/**
*  @brief Frees a given correlation matrix
*  @param num_rows number of rows and cols of the square matrix
*  @return void
*  
*  Frees a given correlation matrix
*/
void correlation_matrix_free(size_t num_rows, correlation_matrix_t* correlation_matrix_p);

#endif /* CORRELATION_MATRIX_H */
