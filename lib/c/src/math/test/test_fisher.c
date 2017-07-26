#include <stdio.h>
#include <stdlib.h>

#include <omp.h>

#include "../fisher.h"


int main(int argc, char *argv[]) {
    double start, end;
    int loop = 1e7;
    int mode = GREATER;
    
    int a = 3, b = 1, c = 1, d = 3;
    double *factorial_logarithms = init_logarithm_array((a+b+c+d)*2);
    double result = 0.0;
    double *result_v = NULL;
    
    printf("Factorial logarithms precalculated\n");
    
    // Benchmark sequential implementation
    start = omp_get_wtime();
    for (int i = 0; i < loop * 10; i++) {
        result = fisher_test(a, b, c, d, mode, factorial_logarithms);
    }
    end = omp_get_wtime();

    printf("Sequential implementation: %.3f s, result = %f\n", end - start, result);
    
    // Benchmark vectorized implementation
    int *a_v = (int*) malloc (loop * sizeof(int));
    int *b_v = (int*) malloc (loop * sizeof(int));
    int *c_v = (int*) malloc (loop * sizeof(int));
    int *d_v = (int*) malloc (loop * sizeof(int));
    for (int i = 0; i < loop; i++) {
        a_v[i] = a; b_v[i] = b;
        c_v[i] = c; d_v[i] = d;
    }
    
    start = omp_get_wtime();
    for (int i = 0; i < 10; i++) {
        result_v = fisher_test_vectorized(a_v, b_v, c_v, d_v, loop, mode, factorial_logarithms);
        if (i < 9) { free(result_v); }
    }
    end = omp_get_wtime();
    
    printf("Vectorized implementation: %.3f s, result_v[0] = %f\n", end - start, result_v[0]);
    for (int i = 0; i < loop; i++) {
        if (result != result_v[i]) {
            printf("Position %d is not correct (%e)\n", i, result_v[i]);
        }
    }
    free(result_v);
    
    // Benchmark vectorized implementation using OpenMP
    start = omp_get_wtime();
    for (int i = 0; i < 10; i++) {
        result_v = fisher_test_omp(a_v, b_v, c_v, d_v, loop, mode, factorial_logarithms);
        if (i < 9) { free(result_v); }
    }
    end = omp_get_wtime();
    
    printf("OpenMP implementation: %.3f s, result_v[0] = %f\n", end - start, result_v[0]);
    for (int i = 0; i < loop; i++) {
        if (result != result_v[i]) {
            printf("Position %d is not correct (%e)\n", i, result_v[i]);
        }
    }
    free(result_v);
    
    free(a_v); free(b_v); free(c_v); free(d_v);
    free(factorial_logarithms);
    
    return 0;
}
