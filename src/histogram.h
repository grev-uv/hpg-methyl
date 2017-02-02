#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <nmmintrin.h> // SSE4.2 support
#include "buffers.h"

//====================================================================================
//  Structures and prototypes
//====================================================================================

/**
 * @brief Structure used to store the parameters needed to compute the histogram
 * of the read batch and to return the result.
 * 
 */

typedef struct histogram_input {
  fastq_read_t* in_read; /**< FastQ read with the batch data */
  size_t out_nc; /**< Number of cytosines in the batch */
  size_t out_ng; /**< Number of guanines in the batch */
  float out_ncg; /**< Proportion of cytosines with respect to guanines */
  float out_ngc; /**< Proportion of guanines with respect to cytosines */
  unsigned char priv_use_simd; /**< (Private member) SIMD is supported in this runtime */
} histogram_input_t;

/**
 * @brief Initializes the histogram input structure.
 * @param h_input Pointer to the histogram input structure.
 * @param read Pointer to the FastQ read.
 *
 */
void histogram_input_init(histogram_input_t* h_input, const fastq_read_t* read);

/**
 * @brief Calculates the histogram for a given FastQ batch.
 * @param h_input Pointer to the histogram input structure.
 *
 */
void histogram_apply(histogram_input_t* h_input);

/**
 * @brief (Private function) Counts the number of occurrences
 * of cytosines and guanines.
 * @param h_input Pointer to the histogram input structure.
 *
 */
void __histogram_count_cg(histogram_input_t* h_input);

/**
 * @brief (Private function) Counts the number of occurrences
 * of cytosines and guanines using vector extensions (SSE4.2).
 * @param h_input Pointer to the histogram input structure.
 *
 */
void __histogram_count_cg_simd(histogram_input_t* h_input);

#endif // HISTOGRAM_H