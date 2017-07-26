#include "math_utils.h"

double log_base(double value, double base) {
    return log(value) / log(base);
}

double round_digits(double a, double exp) {
    return round(a * pow(10, exp)) / pow(10, exp);
}
