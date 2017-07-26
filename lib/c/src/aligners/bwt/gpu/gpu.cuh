#ifndef _SEARCH_GPU_
#define _SEARCH_GPU_

#include <cuda.h>

#include "../search/csafm.h"
#include "../search/io.h"

#if defined FM_COMP_32 || FM_COMP_64

#define BWiterationGPU(k_in, l_in, k_out, l_out, b, C, C1, O)\
  do {\
    (k_out) = (C1)[(b)] + getOcompValueGPU((b), (k_in)  , (O));\
    (l_out) = (C)[(b)]  + getOcompValueGPU((b), (l_in)+1, (O));\
  } while(0);

#else

#define BWiterationGPU(k_in, l_in, k_out, l_out, b, C, C1, O)\
  do {\
    (k_out) = (C1)[(b)] + (O).desp[(b)][(k_in)  ];\
    (l_out) = (C)[(b)]  + (O).desp[(b)][(l_in)+1];\
  } while(0);

#endif

#define manageCudaError()\
  error = cudaGetLastError();\
  if (error != cudaSuccess) {\
    fprintf(stderr, "Error kernel: %s\n", cudaGetErrorString(error));\
    exit(1);\
  }

#ifdef __cplusplus
extern "C" {
#endif

void copy_vector_gpu(vector *device, vector *host);

void read_comp_matrix_gpu(comp_matrix *matrix, const char *directory, const char *name);
void copy_comp_matrix_gpu(comp_matrix *device, comp_matrix *host);
void free_comp_matrix_gpu_host(comp_matrix *reverse, comp_matrix *strand); 
void free_comp_matrix_gpu_device(comp_matrix *reverse, comp_matrix *strand);

void reverse_strand_gpu_O(comp_matrix *r_O, comp_matrix *s_O);

void BWExactSearchBackwardGPUWrapper(uintmax_t num_bloques, uintmax_t tam_bloques, uint8_t* W, uint64_t* nW, intmax_t* k, intmax_t* l, intmax_t k_ini, intmax_t l_ini, vector* C, vector* C1, comp_matrix* O);
void BWExactSearchForwardGPUWrapper(uintmax_t num_bloques, uintmax_t tam_bloques, uint8_t* W, uint64_t* nW, intmax_t* k, intmax_t* l, intmax_t k_ini, intmax_t l_ini, vector* C, vector* C1, comp_matrix* O);

void BWExactSearchBackwardVectorGPUWrapper(uintmax_t num_bloques, uintmax_t tam_bloques, uint8_t* W, uint64_t* nW, intmax_t* k, intmax_t* l, intmax_t k_ini, intmax_t l_ini, vector* C, vector* C1, comp_matrix* O);
void BWExactSearchForwardVectorGPUWrapper(uintmax_t num_bloques, uintmax_t tam_bloques, uint8_t* W, uint64_t* nW, intmax_t* k, intmax_t* l, intmax_t k_ini, intmax_t l_ini, vector* C, vector* C1, comp_matrix* O);

#ifdef __cplusplus
}
#endif

/*
typedef struct {
  result *results;
  uintmax_t *num_results;
} blocked_results_lists;

#define add_resultGPU(_orig, _r_list_results, _num_results, _max_results)\
  do {\
    if ((_num_results) < (_max_results)) {\
      result *dest;\
      dest = (_r_list_results) + (_num_results);\
      (_num_results)++;\
      copy_result(dest, (_orig));\
    }\
  } while(0);
*/
//TODO: Add error detection in CUDA kernel and notify which string (offset) has got problems

//void declare_blocked_results_list_cpu(blocked_results_lists *lists, uintmax_t max_results, uintmax_t num_lists); 
//void declare_blocked_results_list_gpu(blocked_results_lists *lists, uintmax_t max_results, uintmax_t num_lists);
//void copy_blocked_results_list_gpu(blocked_results_lists *lists_gpu, blocked_results_lists *lists_cpu, uintmax_t max_results, uintmax_t num_lists);
//void copy_blocked_results_list_cpu(blocked_results_lists *lists_cpu, blocked_results_lists *lists_gpu, uintmax_t max_results, uintmax_t num_lists);
//void write_blocked_results(blocked_results_lists *r_list, exome* ex, comp_vector *S, comp_vector *Si, vector *C, comp_matrix *O, comp_matrix *Oi, char *mappings, SA_TYPE nW, bool type, FILE *fp, uintmax_t max_results, uintmax_t num_lists, uintmax_t block_read_index);

//void BWSearchGPU(uintmax_t num_bloques, uintmax_t tam_bloques, uint8_t *W, char *h_W, SA_TYPE nW, vector *C, vector *h_C, vector *C1, vector *h_C1, comp_matrix *O, comp_matrix *h_O, comp_matrix *Oi, comp_matrix *h_Oi, comp_vector *S, comp_vector *R, comp_vector *Si, comp_vector *Ri, blocked_results_lists *rl_prev_cpu, blocked_results_lists *rl_next_cpu, blocked_results_lists *rl_prev_i_cpu, blocked_results_lists *rl_next_i_cpu, blocked_results_lists *rl_final_cpu, blocked_results_lists *rl_prev_gpu, blocked_results_lists *rl_next_gpu, blocked_results_lists *rl_prev_i_gpu, blocked_results_lists *rl_next_i_gpu, blocked_results_lists *rl_final_gpu, int16_t fragsize, uintmax_t max_results);

#endif
