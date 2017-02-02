#include "histogram.h"


//------------------------------------------------------------------------------------

void histogram_input_init(histogram_input_t* h_input, const fastq_read_t* read) {
  // Check if SSE4.2 and POPCNT CPU extensions are supported.
  //
  // If the application is compiled using GCC, it will
  // choose the best version at runtime. If other compiler
  // is used, the default scalar version is used.
  h_input->priv_use_simd = 0;

#ifdef __GNUC__
  if (__builtin_cpu_supports("sse4.2") && __builtin_cpu_supports("popcnt")) {
    h_input->priv_use_simd = 1;
  }
#endif // __GNUC__

  h_input->in_read = (fastq_read_t*)read;
}

//------------------------------------------------------------------------------------

void __histogram_count_cg(histogram_input_t* h_input) {
  const fastq_read_t* read = h_input->in_read;
  size_t nc = 0, ng = 0;

  // Count the C and G occurrences
  for (size_t i = 0; i < read->length; ++i) {
		if (read->sequence[i] == 'C') {
			++nc;
		} else if (read->sequence[i] == 'G') {
			++ng;
		}
	}

  // Store the occurrence count
  h_input->out_nc = nc;
  h_input->out_ng = ng;
}

//------------------------------------------------------------------------------------

void __histogram_count_cg_simd(histogram_input_t* h_input) {
  const fastq_read_t* read = h_input->in_read;
  size_t nc = 0, ng = 0;

  // Calculate how many bases will be counted using vector
  // instructions and how many using a scalar loop, since
  // in the case of using vector instructions these must
  // be counted in multiples of 16 characters
  size_t vec_count = read->length - (read->length % 16);

  // Store the reference strings for comparing with all C's
  // and all G's in the SIMD registers
  const char ref_c[16] __attribute__ ((aligned (16))) = {
    'C','C','C','C','C','C','C','C',
    'C','C','C','C','C','C','C','C'
  };

  const char ref_g[16] __attribute__ ((aligned (16))) = {
    'G','G','G','G','G','G','G','G',
    'G','G','G','G','G','G','G','G'
  };

  const __m128i* ref_c_simd = (__m128i*)ref_c;
  const __m128i* ref_g_simd = (__m128i*)ref_g;
  __m128i sequence, mask_nc, mask_ng;

  // Count the C and G in multiples of 16
  // Settings for the string comparation:
  // - Assume 8 bit unsigned integer characters.
  // - Do the comparation character by character.
  // - Return a mask with a bit per character.
  for (size_t i = 0; i < vec_count; i += 16) {
		// Load the next 16 characters
    sequence = _mm_loadu_si128((__m128i*)&read->sequence[i]);

    // Compare the sequence with G's
    mask_nc = _mm_cmpestrm(*ref_c_simd, 16, sequence, 16, _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | 
                          _SIDD_BIT_MASK);

    // Compare the sequence with C's
    mask_ng = _mm_cmpestrm(*ref_g_simd, 16, sequence, 16, _SIDD_UBYTE_OPS | _SIDD_CMP_EQUAL_EACH | 
                          _SIDD_BIT_MASK);

    // Perform the reduction on both results.
    // The reduction is performed on the lower 64 bits
    // of the SIMD register, which are unpacked with the
    // conversion intrinsic.
    nc += _mm_popcnt_u64(_mm_cvtsi128_si64(mask_nc));
    ng += _mm_popcnt_u64(_mm_cvtsi128_si64(mask_ng));
	}

  // Count the remaining Cs and Gs outside of the
  // 16 character boundary
  for (int i = vec_count; i < read->length; i++) {
    if (read->sequence[i] == 'C') {
			++nc;
		} else if (read->sequence[i] == 'G') {
			++ng;
		}
  }

  // Store the occurrence count
  h_input->out_nc = nc;
  h_input->out_ng = ng;
}

//------------------------------------------------------------------------------------

void histogram_apply(histogram_input_t* h_input) {
  // Initialize histogram values
  h_input->out_nc = 0;
  h_input->out_ng = 0;
  h_input->out_ncg = 0.0f;
  h_input->out_ngc = 0.0f;

  // Count ocurrences of C and G, using the scalar or
  // vector version if extensions are available.
  if (h_input->priv_use_simd) {
    __histogram_count_cg_simd(h_input);
  } else {
    __histogram_count_cg(h_input);
  }

  // Calculate the proportion of C and G
  h_input->out_ncg = (h_input->out_nc + h_input->out_ng == 0) ? (0.5) : 
                     (1.0 * h_input->out_nc / (h_input->out_nc + h_input->out_ng));
  h_input->out_ngc = 1.0 - h_input->out_ncg;

  LOG_DEBUG_F("========= END OF HISTOGRAM OF %s =========\n", h_input->in_read->sequence);
  LOG_DEBUG_F("========= VALUES Ncg %f =========\n", h_input->out_ncg);
  LOG_DEBUG_F("========= VALUES Ngc %f =========\n", h_input->out_ngc);
}
