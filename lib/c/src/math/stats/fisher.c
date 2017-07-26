#include "fisher.h"

double fisher_test(int a, int b, int c, int d, enum Fisher_mode mode, double *factorial_logarithms) {
    int init_a = 0, init_b = 0, init_c = 0, init_d = 0;
    double result;
    int steps = 0;
    switch (mode) {
        case 1:
            if(a > d) {
                init_a = a - d;
                init_b = b + d;
                init_c = c + d;
                init_d = 0;
                steps = d;
            }
            else {
                init_a = 0;
                init_b = b + a;
                init_c = c + a;
                init_d = d - a;
                steps = a;
            }
            break;
        case 2:
            init_a = a;
            init_b = b;
            init_c = c;
            init_d = d;
            if (b > c) {
                steps = c;
            } else {
                steps = b;
            }
            break;
        case 3:
            if(a > d) {
                init_a = a - d;
                init_b = b + d;
                init_c = c + d;
                init_d = 0;
                steps = d;
            } else {
                init_a = 0;
                init_b = b + a;
                init_c = c + a;
                init_d = d - a;
                steps = a;
            }
            if (b > c) {
                steps += c;
            } else {
                steps += b;
            }
            break;
    }
    
    double n = factorial_logarithms[a+b+c+d];
    double numerator = factorial_logarithms[init_a+init_b] + factorial_logarithms[init_b+init_d] + factorial_logarithms[init_a+init_c] + factorial_logarithms[init_c+init_d];
    double denominator = n + factorial_logarithms[init_a] + factorial_logarithms[init_b] + factorial_logarithms[init_c] + factorial_logarithms[init_d];
    double current_p = numerator - denominator;

    if(mode == 3) {
        numerator = factorial_logarithms[a+b] + factorial_logarithms[b+d] + factorial_logarithms[a+c] + factorial_logarithms[c+d];
        denominator = n + factorial_logarithms[a] + factorial_logarithms[b] + factorial_logarithms[c] + factorial_logarithms[d];
        double init_p = numerator - denominator;

        if (current_p <= init_p) {
            result = exp(current_p);
        } else {
            result = 0;
        }

        while(steps-- > 0) {
            init_a++;
            init_d++;
            init_c--;
            init_b--;
            numerator = factorial_logarithms[init_a+init_b] + factorial_logarithms[init_b+init_d] + factorial_logarithms[init_a+init_c] + factorial_logarithms[init_c+init_d];
            denominator = n + factorial_logarithms[init_a] + factorial_logarithms[init_b] + factorial_logarithms[init_c] + factorial_logarithms[init_d];
            current_p = numerator - denominator;
            if(current_p <= init_p) {
                result += exp(current_p);
            }
        }
    } else {
        result = exp(current_p);
        while(steps-- > 0) {
            init_a++;
            init_d++;
            init_c--;
            init_b--;
            numerator = factorial_logarithms[init_a+init_b] + factorial_logarithms[init_b+init_d] + factorial_logarithms[init_a+init_c] + factorial_logarithms[init_c+init_d];
            denominator = n + factorial_logarithms[init_a] + factorial_logarithms[init_b] + factorial_logarithms[init_c] + factorial_logarithms[init_d];
            current_p = numerator - denominator;
            result += exp(current_p);
        }
    }
    return result;
}


double *fisher_test_vectorized(int* a, int* b, int* c, int* d, int len, enum Fisher_mode mode, double* factorial_logarithms) {
    double *result = (double*) malloc (len * sizeof(double));
    
    int *init_a = (int*) malloc (len * sizeof(int));
    int *init_b = (int*) malloc (len * sizeof(int));
    int *init_c = (int*) malloc (len * sizeof(int));
    int *init_d = (int*) malloc (len * sizeof(int));
    
    int *steps = (int*) malloc (len * sizeof(int));
    double *n = (double*) malloc (len * sizeof(double));
    double *numerator = (double*) malloc (len * sizeof(double));
    double *denominator = (double*) malloc (len * sizeof(double));
    double *current_p = (double*) malloc (len * sizeof(double));
    
    double *init_p = (mode == 3) ? (double*) malloc (len * sizeof(double)) : NULL;
    
    for (int i = 0; i < len; i++) {
        switch (mode) {
            case 1:
                if(a[i] > d[i]) {
                    init_a[i] = a[i] - d[i];
                    init_b[i] = b[i] + d[i];
                    init_c[i] = c[i] + d[i];
                    init_d[i] = 0;
                    steps[i] = d[i];
                }
                else {
                    init_a[i] = 0;
                    init_b[i] = b[i] + a[i];
                    init_c[i] = c[i] + a[i];
                    init_d[i] = d[i] - a[i];
                    steps[i] = a[i];
                }
                break;
            case 2:
                init_a[i] = a[i];
                init_b[i] = b[i];
                init_c[i] = c[i];
                init_d[i] = d[i];
                if (b[i] > c[i]) {
                    steps[i] = c[i];
                } else {
                    steps[i] = b[i];
                }
                break;
            case 3:
                if(a[i] > d[i]) {
                    init_a[i] = a[i] - d[i];
                    init_b[i] = b[i] + d[i];
                    init_c[i] = c[i] + d[i];
                    init_d[i] = 0;
                    steps[i] = d[i];
                } else {
                    init_a[i] = 0;
                    init_b[i] = b[i] + a[i];
                    init_c[i] = c[i] + a[i];
                    init_d[i] = d[i] - a[i];
                    steps[i] = a[i];
                }
                if (b[i] > c[i]) {
                    steps[i] += c[i];
                } else {
                    steps[i] += b[i];
                }
                break;
        }
        
        n[i] = factorial_logarithms[a[i]+b[i]+c[i]+d[i]];
        numerator[i] = factorial_logarithms[init_a[i]+init_b[i]] + factorial_logarithms[init_b[i]+init_d[i]] + 
                       factorial_logarithms[init_a[i]+init_c[i]] + factorial_logarithms[init_c[i]+init_d[i]];
        denominator[i] = n[i] + factorial_logarithms[init_a[i]] + factorial_logarithms[init_b[i]] + factorial_logarithms[init_c[i]] + factorial_logarithms[init_d[i]];
        current_p[i] = numerator[i] - denominator[i];

        if(mode == 3) {
            numerator[i] = factorial_logarithms[a[i]+b[i]] + factorial_logarithms[b[i]+d[i]] + factorial_logarithms[a[i]+c[i]] + factorial_logarithms[c[i]+d[i]];
            denominator[i] = n[i] + factorial_logarithms[a[i]] + factorial_logarithms[b[i]] + factorial_logarithms[c[i]] + factorial_logarithms[d[i]];
            init_p[i] = numerator[i] - denominator[i];

            if (current_p[i] <= init_p[i]) {
                result[i] = exp(current_p[i]);
            } else {
                result[i] = 0;
            }

            while(steps[i]-- > 0) {
                init_a[i]++;
                init_d[i]++;
                init_c[i]--;
                init_b[i]--;
                numerator[i] = factorial_logarithms[init_a[i]+init_b[i]] + factorial_logarithms[init_b[i]+init_d[i]] + 
                               factorial_logarithms[init_a[i]+init_c[i]] + factorial_logarithms[init_c[i]+init_d[i]];
                denominator[i] = n[i] + factorial_logarithms[init_a[i]] + factorial_logarithms[init_b[i]] + factorial_logarithms[init_c[i]] + factorial_logarithms[init_d[i]];
                current_p[i] = numerator[i] - denominator[i];
                if(current_p[i] <= init_p[i]) {
                    result[i] += exp(current_p[i]);
                }
            }
        } else {
            result[i] = exp(current_p[i]);
            while(steps[i]-- > 0) {
                init_a[i]++;
                init_d[i]++;
                init_c[i]--;
                init_b[i]--;
                numerator[i] = factorial_logarithms[init_a[i]+init_b[i]] + factorial_logarithms[init_b[i]+init_d[i]] + 
                               factorial_logarithms[init_a[i]+init_c[i]] + factorial_logarithms[init_c[i]+init_d[i]];
                denominator[i] = n[i] + factorial_logarithms[init_a[i]] + factorial_logarithms[init_b[i]] + factorial_logarithms[init_c[i]] + factorial_logarithms[init_d[i]];
                current_p[i] = numerator[i] - denominator[i];
                result[i] += exp(current_p[i]);
            }
        }
    }
    
    free(init_a);
    free(init_b);
    free(init_c);
    free(init_d);
    
    free(steps);
    free(n);
    free(numerator);
    free(denominator);
    free(init_p);
    free(current_p);
    
    
    return result;
}


double *fisher_test_omp(int* a, int* b, int* c, int* d, int len, enum Fisher_mode mode, double* factorial_logarithms) {
    double *result = (double*) malloc (len * sizeof(double));
   
    int init_a, init_b, init_c, init_d;
    int steps;
    double n, numerator, denominator, current_p, init_p;
    
#pragma omp parallel for private(init_a, init_b, init_c, init_d, n, steps, numerator, denominator, current_p, init_p)
    for (int i = 0; i < len; i++) {
        switch (mode) {
            case 1:
                if(a[i] > d[i]) {
                    init_a = a[i] - d[i];
                    init_b = b[i] + d[i];
                    init_c = c[i] + d[i];
                    init_d = 0;
                    steps = d[i];
                }
                else {
                    init_a = 0;
                    init_b = b[i] + a[i];
                    init_c = c[i] + a[i];
                    init_d = d[i] - a[i];
                    steps = a[i];
                }
                break;
            case 2:
                init_a = a[i];
                init_b = b[i];
                init_c = c[i];
                init_d = d[i];
                if (b[i] > c[i]) {
                    steps = c[i];
                } else {
                    steps = b[i];
                }
                break;
            case 3:
                if(a[i] > d[i]) {
                    init_a = a[i] - d[i];
                    init_b = b[i] + d[i];
                    init_c = c[i] + d[i];
                    init_d = 0;
                    steps = d[i];
                } else {
                    init_a = 0;
                    init_b = b[i] + a[i];
                    init_c = c[i] + a[i];
                    init_d = d[i] - a[i];
                    steps = a[i];
                }
                if (b[i] > c[i]) {
                    steps += c[i];
                } else {
                    steps += b[i];
                }
                break;
        }
        
        n = factorial_logarithms[a[i]+b[i]+c[i]+d[i]];
        numerator = factorial_logarithms[init_a+init_b] + factorial_logarithms[init_b+init_d] + 
                       factorial_logarithms[init_a+init_c] + factorial_logarithms[init_c+init_d];
        denominator = n + factorial_logarithms[init_a] + factorial_logarithms[init_b] + factorial_logarithms[init_c] + factorial_logarithms[init_d];
        current_p = numerator - denominator;

        if(mode == 3) {
            numerator = factorial_logarithms[a[i]+b[i]] + factorial_logarithms[b[i]+d[i]] + factorial_logarithms[a[i]+c[i]] + factorial_logarithms[c[i]+d[i]];
            denominator = n + factorial_logarithms[a[i]] + factorial_logarithms[b[i]] + factorial_logarithms[c[i]] + factorial_logarithms[d[i]];
            init_p = numerator - denominator;

            if (current_p <= init_p) {
                result[i] = exp(current_p);
            } else {
                result[i] = 0;
            }

            while(steps-- > 0) {
                init_a++;
                init_d++;
                init_c--;
                init_b--;
                numerator = factorial_logarithms[init_a+init_b] + factorial_logarithms[init_b+init_d] + 
                               factorial_logarithms[init_a+init_c] + factorial_logarithms[init_c+init_d];
                denominator = n + factorial_logarithms[init_a] + factorial_logarithms[init_b] + factorial_logarithms[init_c] + factorial_logarithms[init_d];
                current_p = numerator - denominator;

                if(current_p <= init_p) {
                    result[i] += exp(current_p);
                }

            }
        } else {
            result[i] = exp(current_p);
            while(steps-- > 0) {
                init_a++;
                init_d++;
                init_c--;
                init_b--;
                numerator = factorial_logarithms[init_a+init_b] + factorial_logarithms[init_b+init_d] + 
                               factorial_logarithms[init_a+init_c] + factorial_logarithms[init_c+init_d];
                denominator = n + factorial_logarithms[init_a] + factorial_logarithms[init_b] + factorial_logarithms[init_c] + factorial_logarithms[init_d];
                current_p = numerator - denominator;
                result[i] += exp(current_p);
            }
        }
    }
    
    return result;
}



double *init_logarithm_array(int nmax) {
    double *factorial_logarithms = (double*) calloc (nmax, sizeof(double));
    factorial_logarithms[0] = 0;
    factorial_logarithms[1] = 0;
    
    for(int i = 2; i < nmax; i++) {
        factorial_logarithms[i] = factorial_logarithms[i-1] + log(i);
    }
    
    return factorial_logarithms;
}
