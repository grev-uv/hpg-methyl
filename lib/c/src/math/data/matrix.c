
#include "matrix.h"

/* ******************************************************
 *    		Functions implementations 		*
 * *****************************************************/

matrix_t* matrix_new(size_t num_rows, size_t num_cols) {
    matrix_t* matrix_p = (matrix_t*) calloc(1, sizeof(matrix_t));
    matrix_p->gsl_matrix = gsl_matrix_alloc(num_rows, num_cols);
    
    return matrix_p;
}

void matrix_free(matrix_t *m) {
    if (m == NULL) {
        return;
    }
    
    if (m->gsl_matrix != NULL) {
        gsl_matrix_free(m->gsl_matrix);
        m->gsl_matrix = NULL;
    }

    free(m);
    m = NULL;
    
    return;
}

matrix_t* matrix_init(double value, matrix_t *m) {
    gsl_matrix_set_all(m->gsl_matrix, value);
    
    return m;
}

matrix_t* matrix_fill(double **data, size_t num_rows, size_t num_cols, matrix_t *m) {
    for (size_t i = 0; i < num_rows; i++) {
        for (size_t j = 0; j < num_cols; j++) {
            matrix_set(i, j, data[i][j], m);
        }
    }   
    
    return m;
}

double matrix_get(size_t row, size_t col, matrix_t *m) {
    return gsl_matrix_get(m->gsl_matrix, row, col);
}

void matrix_set(size_t row, size_t col, double value, matrix_t *m) {
    gsl_matrix_set(m->gsl_matrix, row, col, value);
    
    return;
}


int matrix_printf (const matrix_t * m, const char * format);

int matrix_fprintf (const matrix_t * m, const char * format, FILE * stream);


int matrix_memcpy (const matrix_t *src, matrix_t *dest);


int matrix_get_row (vector_t *v, const matrix_t * m, size_t i);

int matrix_get_col (vector_t *v, const matrix_t * m, size_t j);

int matrix_set_row (size_t i, const vector_t * v, matrix_t *m);

int matrix_set_col (size_t j, const vector_t * v, matrix_t *m);

/* ******************************************************
 *    		Operations with matrices  		*
 * *****************************************************/

matrix_t* matrix_scale(matrix_t *m, double factor) {
    gsl_matrix_scale(m->gsl_matrix, factor);

    return m;
}

matrix_t* matrix_add(matrix_t *m1, matrix_t *m2) {
    gsl_matrix_add(m1->gsl_matrix, m2->gsl_matrix);
    
    return m1;
}

matrix_t* matrix_sub(matrix_t *m1, matrix_t *m2) {
    gsl_matrix_sub(m1->gsl_matrix, m2->gsl_matrix);
    
    return m1;    
}

matrix_t* matrix_mul(matrix_t *m1, matrix_t *m2, matrix_t *m_result) {
    gsl_linalg_matmult(m1->gsl_matrix, m2->gsl_matrix, m_result->gsl_matrix);

    return m_result;    
}

matrix_t* matrix_inv(matrix_t *m, matrix_t* m_result) {
    int signum;
    gsl_permutation *p = gsl_permutation_alloc(m->gsl_matrix->size1);
  
    gsl_linalg_LU_decomp(m->gsl_matrix, p, &signum);
    gsl_linalg_LU_invert(m->gsl_matrix, p, m_result->gsl_matrix);
  
    // free permutation
    gsl_permutation_free(p);
    
    return m_result;
}

double matrix_Fnorm(matrix_t *m) {
    size_t num_rows = m->gsl_matrix->size1;
    size_t num_cols = m->gsl_matrix->size2;
    double norm = 0;
    
    for (size_t i = 0; i < num_rows; i++) {
        for (size_t j = 0; j < num_cols; j++) {
            norm = hypot(norm, matrix_get(i, j, m));
        }
    }

    return norm;
}

double hypot(double a, double b) {
    double r;

    if (fabs(a) > fabs(b)) {
        r = b / a;
        r = fabs(a) * sqrt(1 + r * r);
    } else if (b != 0) {
        r = a / b;
        r = fabs(b) * sqrt(1 + r * r);
    } else {
        r = 0.0;
    }
    
    return r;
}
