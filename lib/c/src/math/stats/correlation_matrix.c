
#include <stdlib.h>
#include <string.h>

#include "correlation_matrix.h"

/* **************************************************************
 *    		Public functions implementations 		*
 * *************************************************************/

correlation_matrix_t* correlation_matrix_new(size_t num_rows) {
    correlation_matrix_t* correlation_matrix_p = (correlation_matrix_t*) calloc(1, sizeof(correlation_matrix_t));
    
    correlation_matrix_p->values = (float**) calloc(num_rows, sizeof(float*));
    
    for (size_t i = 0; i < num_rows; i++) {
        correlation_matrix_p->values[i] = (float*) calloc(num_rows, sizeof(float));
    }
    
    return correlation_matrix_p;
}

void correlation_matrix_free(size_t num_rows, correlation_matrix_t* correlation_matrix_p) {
    if ((correlation_matrix_p == NULL) || (correlation_matrix_p->values == NULL)) {
        return;
    }
    
    for (size_t i = 0; i < num_rows; i++) {
        if (correlation_matrix_p->values[i] != NULL) {
            free(correlation_matrix_p->values[i]);
            correlation_matrix_p->values[i] = NULL;
        }
    }
    
    free(correlation_matrix_p->values);
    correlation_matrix_p->values = NULL;
    
    free(correlation_matrix_p);
    correlation_matrix_p = NULL;
 
    return;
}