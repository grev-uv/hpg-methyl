
#ifndef MATRIX_H
#define MATRIX_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include "stats/statistics.h"

/* **************************************
 *    		Structures  		*
 * *************************************/

/**
* @brief Structure for the matrix
*
* Structure for the matrix
*/
typedef struct matrix {
    gsl_matrix *gsl_matrix;	/**< Pointer to GSL matrix structure. */
} matrix_t;

/* **************************************
 *    		Functions  		*
 * *************************************/

/**
*  @brief Creates a new matrix and allocates memory for the data
*  @param num_rows number of rows of the matrix
*  @param num_cols number of cols of the matrix
*  @return pointer to the created matrix
*  
*  Creates a new matrix and allocates memory for the data
*/
matrix_t* matrix_new(size_t num_rows, size_t num_cols);

/**
*  @brief Frees a given matrix
*  @param m pointer to the matrix
*  @return void
*  
*   Frees a given matrix
*/
void matrix_free(matrix_t *m);

/**
*  @brief Inits a matrix with a given constant value
*  @param value double constant value to init matrix elements
*  @param m[in,out] pointer to the matrix
*  @return pointer to the matrix
*  
*  Inits a matrix with a given constant value
*/
matrix_t* matrix_init(double value, matrix_t *m);

/**
*  @brief Fills a matrix with the content of a data double pointer array
*  @param data pointer to the double pointer array
*  @param num_rows number of rows of the matrix
*  @param num_cols number of cols of the matrix
*  @param m pointer to the matrix to fill
*  @return pointer to the matrix
*  
*  Fills a matrix with the content of a data double pointer array
*/
matrix_t* matrix_fill(double **data, size_t num_rows, size_t num_cols, matrix_t *m);

/**
*  @brief Gets a double value from a given coordinate in the matrix
*  @param row number of the row of the element
*  @param col number of the columns of the element
*  @param m pointer to the matrix
*  @return double value in the given coordinate
*  
*  Gets a double value from a given coordinate in the matrix
*/
double matrix_get(size_t row, size_t col, matrix_t *m);

/**
*  @brief Sets a double value in a given coordinate in the matrix
*  @param row number of the row of the element
*  @param col number of the columns of the element
*  @param value double value to insert in the matrix
*  @param m pointer to the matrix
*  @return double value in the given coordinate
*  
*  Sets a double value in a given coordinate in the matrix
*/
void matrix_set(size_t row, size_t col, double value, matrix_t *m);


int matrix_printf(const matrix_t *m, const char *format);

int matrix_fprintf(const matrix_t *m, const char *format, FILE *stream);


int matrix_memcpy(const matrix_t *src, matrix_t *dest);


int matrix_get_row(vector_t *v, const matrix_t *m, size_t i);

int matrix_get_col(vector_t *v, const matrix_t *m, size_t j);

int matrix_set_row(size_t i, const vector_t *v, matrix_t *m);

int matrix_set_col(size_t j, const vector_t *v, matrix_t *m);

/* ******************************************************
 *    		Operations with matrices  		*
 * *****************************************************/

/**
*  @brief Scales/multiplies by a factor all elements of a matrix
*  @param m pointer to the matrix to scale
*  @param factor double factor to scale the matrix
*  @return pointer to the scaled matrix
*  
*  Scales/multiplies by a factor all elements of a matrix
*/
matrix_t* matrix_scale(matrix_t *m, double factor);

/**
*  @brief Adds two matrices (m1 = m1 + m2)
*  @param m1 pointer to the matrix 1
*  @param m2 pointer to the matrix 2 
*  @return pointer to the result matrix 1
*  
*  Adds two matrices (m1 = m1 + m2)
*/
matrix_t* matrix_add(matrix_t *m1, matrix_t *m2);

/**
*  @brief Substracts two matrices (m1 = m1 - m2)
*  @param m1 pointer to the matrix 1
*  @param m2 pointer to the matrix 2  
*  @return pointer to the result matrix 1
*  
*  Substracts two matrices (m1 = m1 - m2)
*/
matrix_t* matrix_sub(matrix_t *m1, matrix_t *m2);

/**
*  @brief Multiply two matrices (m_result = m1 * m2)
*  @param m1 pointer to the matrix 1
*  @param m2 pointer to the matrix 2
*  @param m_result pointer to the result matrix 
*  @return pointer to the result matrix
*  
*  Multiply two matrices (m1 = m1 * m2)
*/
matrix_t* matrix_mul(matrix_t *m1, matrix_t *m2, matrix_t *m_result);

/**
*  @brief Calculate the inverse of a matrix (m_result = 1 / m)
*  @param m pointer to the matrix
*  @param m_result pointer to the inverse matrix
*  @return pointer to the inverse matrix
*  
*  Calculate the inverse of a matrix (m_result = 1 / m)
*/
matrix_t* matrix_inv(matrix_t *m, matrix_t* m_result);

/**
*  @brief Calculates the norm of a given matrix
*  @param m pointer to the matrix
*  @return matrix norm double value
*  
*  Calculates the norm of a given matrix
*/
double matrix_Fnorm(matrix_t *m);


#endif  /* MATRIX_H */
