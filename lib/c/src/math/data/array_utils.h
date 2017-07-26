#ifndef ARRAY_UTILS_H
#define ARRAY_UTILS_H

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "math/math_utils.h"

int array_scalar_multiply(double *data, size_t n, double value);

int array_scalar_sum(double *data, size_t n, double value);


int array_sum(double *data1, const double *data2, size_t n);

int array_substract(double *data1, const double *data2, size_t n);

int array_dotproduct(double *data1, const double *data2, size_t n, double *res);


int array_log(double *values, size_t n);

int array_log10(double *values, size_t n);

int array_log_base(double *values, size_t n, double base);


int array_accum(const double *values, size_t n, double *res);

int array_accum_range(const double *values, size_t begin, size_t end, double *res);


int array_order(double *values, size_t n, int asc, size_t *indices);

int array_ordered(const double *values, size_t n, const size_t *indices, double *ordered);

/**
 * @brief Shuffles an array of doubles
 * @details Shuffles an array of doubles using the Fisher–Yates shuffle algorithm, which guarantees that 
 * every permutation is equally likely (unbiased)
 * @see http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
 * 
 * @param[in,out] values The array to shuffle
 * @param n Number of elements in the array
 **/
void array_shuffle(double *values, size_t n);

/**
 * @brief Shuffles an array of integers
 * @details Shuffles an array of integers using the Fisher–Yates shuffle algorithm, which guarantees that 
 * every permutation is equally likely (unbiased)
 * @see http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
 * 
 * @param[in,out] values The array to shuffle
 * @param n Number of elements in the array
 **/
void array_shuffle_int(int *values, size_t n);


size_t array_printf(const double *values, size_t n, char *format);

size_t array_fprintf(const double *values, size_t n, char *format, FILE *file);

void array_fread(FILE *file, double *values, size_t n);


#endif
