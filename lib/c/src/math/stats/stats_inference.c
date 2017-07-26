#include "stats_inference.h"

inline double chi_square(int observed, int expected, double* p_value) {
    double stat = ((observed - expected) * (observed - expected)) / expected;
    
    // TODO calculate p-value
    
    return stat;
}

double chi_square_test_1df(int a, int b, int c, int d) {
    double total = a + c + b + d;
    
    double total_category_1 = a + c;
    double total_category_2 = b + d;
    
    double total_observed = a + b;
    double total_expected = c + d;
    
    double expected_a = (total_category_1 * total_observed) / total;
    double expected_c = (total_category_1 * total_expected) / total;
    double expected_b = (total_category_2 * total_observed) / total;
    double expected_d = (total_category_2 * total_expected) / total;

    return ((a - expected_a) * (a - expected_a)) / expected_a + 
           ((c - expected_c) * (c - expected_c)) / expected_c +
           ((b - expected_b) * (b - expected_b)) / expected_b + 
           ((d - expected_d) * (d - expected_d)) / expected_d ;
}
