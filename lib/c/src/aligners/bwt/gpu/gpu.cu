#include "gpu.cuh"

__device__ SA_TYPE getOcompValueGPU(SA_TYPE n, SA_TYPE m, comp_matrix O) {
	SA_TYPE pos, desp;
	pos  = m / FM_COMP_VALUE;
	desp = m % FM_COMP_VALUE;
#if   defined FM_COMP_32
	return O.desp[n][pos] + __popc  ( O.count[n][pos] << (FM_COMP_VALUE - (desp + 1)) );
#elif defined FM_COMP_64
	return O.desp[n][pos] + __popcll( O.count[n][pos] << (FM_COMP_VALUE - (desp + 1)) );
#endif
}

void reverse_strand_gpu_O(comp_matrix *r_O, comp_matrix *s_O) {

	cudaError_t error;

	r_O->siz = s_O->siz;

	r_O->n_desp = s_O->n_desp;
	r_O->m_desp = s_O->m_desp;

	SA_TYPE **r_desp = (SA_TYPE **) malloc(r_O->n_desp * sizeof(SA_TYPE *));
	check_malloc(r_desp, "reverse_strand_O");
	SA_TYPE **s_desp = (SA_TYPE **) malloc(s_O->n_desp * sizeof(SA_TYPE *));
	check_malloc(s_desp, "reverse_strand_O");

	cudaMemcpy(s_desp, s_O->desp, s_O->n_desp * sizeof(SA_TYPE *), cudaMemcpyDeviceToHost);
	manageCudaError();

	if (AA != ((uint8_t) -1) && TT !=((uint8_t) -1)) {
		r_desp[AA] = s_desp[TT];
		r_desp[TT] = s_desp[AA];
	} else if (AA != ((uint8_t) -1)) {
		r_desp[AA] = s_desp[AA];
	} else if (TT != ((uint8_t) -1)) {
		r_desp[TT] = s_desp[TT];
	}

	if (CC != ((uint8_t) -1) && GG !=((uint8_t) -1)) {
		r_desp[CC] = s_desp[GG];
		r_desp[GG] = s_desp[CC];
	} else if (CC != ((uint8_t) -1)) {
		r_desp[CC] = s_desp[CC];
	} else if (GG != ((uint8_t) -1)) {
		r_desp[GG] = s_desp[GG];
	}

	cudaMalloc((void**) &r_O->desp, r_O->n_desp * sizeof(SA_TYPE *));
	manageCudaError();
	cudaMemcpy(r_O->desp, r_desp, r_O->n_desp * sizeof(SA_TYPE *), cudaMemcpyHostToDevice);
	manageCudaError();

	free(r_desp);
	free(s_desp);

#if defined FM_COMP_32 || FM_COMP_64

	r_O->n_count = s_O->n_count;
	r_O->m_count = s_O->m_count;

	FM_COMP_TYPE **r_count = (FM_COMP_TYPE **) malloc(r_O->n_count * sizeof(FM_COMP_TYPE *));
	check_malloc(r_count, "reverse_strand_O");
	FM_COMP_TYPE **s_count = (FM_COMP_TYPE **) malloc(s_O->n_count * sizeof(FM_COMP_TYPE *));
	check_malloc(s_count, "reverse_strand_O");

	cudaMemcpy(s_count, s_O->count, s_O->n_count * sizeof(FM_COMP_TYPE *), cudaMemcpyDeviceToHost);
	manageCudaError();

	if (AA != ((uint8_t) -1) && TT !=((uint8_t) -1)) {
		r_count[AA] = s_count[TT];
		r_count[TT] = s_count[AA];
	} else if (AA != ((uint8_t) -1)) {
		r_count[AA] = s_count[AA];
	} else if (TT != ((uint8_t) -1)) {
		r_count[TT] = s_count[TT];
	}

	if (CC != ((uint8_t) -1) && GG !=((uint8_t) -1)) {
		r_count[CC] = s_count[GG];
		r_count[GG] = s_count[CC];
	} else if (CC != ((uint8_t) -1)) {
		r_count[CC] = s_count[CC];
	} else if (GG != ((uint8_t) -1)) {
		r_count[GG] = s_count[GG];
	}

	cudaMalloc((void**) &r_O->count, r_O->n_count * sizeof(FM_COMP_TYPE *));
	manageCudaError();
	cudaMemcpy(r_O->count, r_count, r_O->n_count * sizeof(FM_COMP_TYPE *), cudaMemcpyHostToDevice);
	manageCudaError();

	free(r_count);
	free(s_count);

#endif

}

void read_comp_matrix_gpu(comp_matrix *matrix, const char *directory, const char *name) {

	size_t err=0;
	cudaError_t error;
	FILE *fp;

	char path[500];

	path[0]='\0';
	strcat(path, directory);
	strcat(path, "/");
	strcat(path, name);
	strcat(path, ".desp");

	fp  = fopen(path,  "rb+");
	check_file_open(fp, path);

	err = fread(&matrix->siz,    sizeof(SA_TYPE),  1, fp);
	check_file_read(err, 1, path);

	err = fread(&matrix->n_desp, sizeof(SA_TYPE),  1, fp);
	check_file_read(err, 1, path);

	err = fread(&matrix->m_desp, sizeof(SA_TYPE),  1, fp);
	check_file_read(err, 1, path);

	cudaMallocHost((void**) &matrix->desp, matrix->n_desp * sizeof(SA_TYPE *));
	manageCudaError();

	for (SA_TYPE i=0; i<matrix->n_desp; i++) {
		cudaMallocHost((void**) &matrix->desp[i], matrix->m_desp * sizeof(SA_TYPE));
		manageCudaError();
		err = fread(matrix->desp[i], sizeof(SA_TYPE), matrix->m_desp, fp);
		check_file_read(err, matrix->m_desp, path);
	}

	fclose(fp);

#if defined FM_COMP_32 || FM_COMP_64
	path[0]='\0';
	strcat(path, directory);
	strcat(path, "/");
	strcat(path, name);
	strcat(path, ".count");

	fp  = fopen(path,  "rb+");
	check_file_open(fp, path);

	err = fread(&matrix->n_count,   sizeof(SA_TYPE),  1, fp);
	check_file_read(err, 1, path);

	err = fread(&matrix->m_count,   sizeof(SA_TYPE),  1, fp);
	check_file_read(err, 1, path);

	cudaMallocHost((void**) &matrix->count, matrix->n_count * sizeof(FM_COMP_TYPE *));
	manageCudaError();

	for (SA_TYPE i=0; i<matrix->n_count; i++){
		cudaMallocHost((void**) &matrix->count[i], matrix->m_count * sizeof(FM_COMP_TYPE));
		manageCudaError();
		err = fread(matrix->count[i], sizeof(FM_COMP_TYPE), matrix->m_count, fp);
		check_file_read(err, matrix->m_count, path);
	}

	fclose(fp);
#endif

}

void copy_vector_gpu(vector *device, vector *host) {

	cudaError_t error;

	device->n = host->n;

	cudaMalloc((void**) &device->vector,  device->n * sizeof(SA_TYPE));
	manageCudaError();
	cudaMemcpy(device->vector, host->vector, device->n * sizeof(SA_TYPE), cudaMemcpyHostToDevice);
	manageCudaError();

}

void copy_comp_matrix_gpu(comp_matrix *device, comp_matrix *host) {

	cudaError_t error;

	device->siz    = host->siz;
	device->n_desp = host->n_desp;
	device->m_desp = host->m_desp;

	SA_TYPE **desp = (SA_TYPE **) malloc(host->n_desp * sizeof(SA_TYPE *)); 
	check_malloc(desp, "reverse_strand_O");

	for (SA_TYPE i=0; i<device->n_desp; i++) {
		cudaMalloc((void**) &desp[i], host->m_desp * sizeof(SA_TYPE));
		manageCudaError();
		cudaMemcpy(desp[i], host->desp[i], host->m_desp * sizeof(SA_TYPE), cudaMemcpyHostToDevice);
		manageCudaError();
	}

	cudaMalloc((void**) &device->desp, host->n_desp * sizeof(SA_TYPE *));
	manageCudaError();
	cudaMemcpy(device->desp, desp, host->n_desp * sizeof(SA_TYPE *), cudaMemcpyHostToDevice);
	manageCudaError();
	free(desp);

#if defined FM_COMP_32 || FM_COMP_64

	device->n_count = host->n_count;
	device->m_count = host->m_count;

	FM_COMP_TYPE **count = (FM_COMP_TYPE **) malloc(host->n_count * sizeof(FM_COMP_TYPE *)); 
	check_malloc(count, "reverse_strand_O");

	for (SA_TYPE i=0; i<device->n_count; i++) {
		cudaMalloc((void**) &count[i], host->m_count * sizeof(FM_COMP_TYPE));
		manageCudaError();
		cudaMemcpy(count[i], host->count[i], host->m_count * sizeof(FM_COMP_TYPE), cudaMemcpyHostToDevice);
		manageCudaError();
	}

	cudaMalloc((void**) &device->count, host->n_count * sizeof(FM_COMP_TYPE *));
	manageCudaError();
	cudaMemcpy(device->count, count, host->n_count * sizeof(FM_COMP_TYPE *), cudaMemcpyHostToDevice);
	manageCudaError();
	free(count);

#endif

}

void free_comp_matrix_gpu_host(comp_matrix *reverse, comp_matrix *strand) {

	for (SA_TYPE i=0; i<strand->n_desp; i++) {
		cudaFreeHost(strand->desp[i]);
#if defined FM_COMP_32 || FM_COMP_64
		cudaFreeHost(strand->count[i]);
#endif
	}

	cudaFreeHost(strand->desp);
	if (reverse != NULL) cudaFreeHost(reverse->desp);
#if defined FM_COMP_32 || FM_COMP_64
	cudaFreeHost(strand->count);
	if (reverse != NULL) cudaFreeHost(reverse->count);
#endif

}

void free_comp_matrix_gpu_device(comp_matrix *reverse, comp_matrix *strand) {

	cudaError_t error;

	SA_TYPE **desp = (SA_TYPE **) malloc(strand->n_desp * sizeof(SA_TYPE *)); 
	check_malloc(desp, "reverse_strand_O");
	cudaMemcpy(desp, strand->desp, strand->n_desp * sizeof(SA_TYPE *), cudaMemcpyDeviceToHost);
	manageCudaError();

#if defined FM_COMP_32 || FM_COMP_64
	FM_COMP_TYPE **count = (FM_COMP_TYPE **) malloc(strand->n_count * sizeof(FM_COMP_TYPE *)); 
	check_malloc(count, "reverse_strand_O");
	cudaMemcpy(count, strand->count, strand->n_count * sizeof(FM_COMP_TYPE *), cudaMemcpyDeviceToHost);
	manageCudaError();
#endif

	for (SA_TYPE i=0; i<strand->n_desp; i++) {
		cudaFree(desp[i]);
#if defined FM_COMP_32 || FM_COMP_64
		cudaFree(count[i]);
#endif
	}

	cudaFree(strand->desp);
	if (reverse != NULL){
		cudaFree(reverse->desp);	
	}
	free(desp);

#if defined FM_COMP_32 || FM_COMP_64
	cudaFree(strand->count);
	if (reverse != NULL) {
		cudaFree(reverse->count);

	}
	free(count);
#endif

}

//void declare_blocked_results_list_cpu(blocked_results_lists *lists, uintmax_t max_results, uintmax_t num_lists) {
//	cudaMallocHost((void**) &(lists->results), max_results * num_lists * sizeof(result));
//	cudaMallocHost((void**) &(lists->num_results), num_lists * sizeof(uintmax_t));
//}
//
//void declare_blocked_results_list_gpu(blocked_results_lists *lists, uintmax_t max_results, uintmax_t num_lists) {
//	cudaMalloc((void**) &(lists->results), max_results * num_lists * sizeof(result));
//	cudaMalloc((void**) &(lists->num_results), num_lists * sizeof(uintmax_t));
//}
//
//void copy_blocked_results_list_gpu(blocked_results_lists *lists_gpu, blocked_results_lists *lists_cpu, uintmax_t max_results, uintmax_t num_lists) {
//	cudaMemcpy(lists_gpu->results, lists_cpu->results, max_results * num_lists * sizeof(result), cudaMemcpyHostToDevice);
//	cudaMemcpy(lists_gpu->num_results, lists_cpu->num_results, num_lists * sizeof(uintmax_t), cudaMemcpyHostToDevice);
//}
//
//void copy_blocked_results_list_cpu(blocked_results_lists *lists_cpu, blocked_results_lists *lists_gpu, uintmax_t max_results, uintmax_t num_lists) {
//	cudaMemcpy(lists_cpu->results, lists_gpu->results, max_results * num_lists * sizeof(result), cudaMemcpyDeviceToHost);
//	cudaMemcpy(lists_cpu->num_results, lists_gpu->num_results, num_lists * sizeof(uintmax_t), cudaMemcpyDeviceToHost);
//}
//
//void write_blocked_results(blocked_results_lists *r_list, exome* ex, comp_vector *S, comp_vector *Si, vector *C, comp_matrix *O, comp_matrix *Oi, char *mappings, uintmax_t nW, bool type, FILE *fp, uintmax_t max_results, uintmax_t num_lists, uintmax_t block_read_index) {
//
//	result *r;
//	bool found;
//
//	char search[MAXLINE+1];
//
//	for (uintmax_t k=0; k<num_lists; k++) {
//
//		found = false;
//
//		search[0] = '\0';
//		strncat(search, mappings + k*MAXLINE, nW);
//
//		for (uintmax_t i=0; i<r_list->num_results[k]; i++) {
//			r = &r_list->results[k*max_results + i];
//			manage_single_result(r, ex, S, Si, C, O, Oi, search, type, fp, block_read_index + k, &found);
//		}
//
//	}
//
//}
//

__global__ void BWExactSearchBackwardGPU(uint8_t *W, uint64_t *nW, intmax_t *k, intmax_t *l, intmax_t k_ini, intmax_t l_ini, SA_TYPE *C, SA_TYPE *C1, comp_matrix O) {

	intmax_t i;
	intmax_t k2, l2;
	uintmax_t offset  = blockIdx.x * blockDim.x + threadIdx.x;

	__shared__ SA_TYPE Cshared[4];
	__shared__ SA_TYPE C1shared[4];

	if (threadIdx.x<4) {
		Cshared[threadIdx.x]  = C[threadIdx.x];
		C1shared[threadIdx.x] = C1[threadIdx.x];
	}

	__syncthreads();

	k2 = k_ini; l2 = l_ini;

	for (i=nW[offset]-1; (k2<=l2) && (i>=0); i--)
		BWiterationGPU(k2, l2, k2, l2, W[offset*MAXLINE+i], Cshared, C1shared, O);

	k[offset] = k2;
	l[offset] = l2;

}

__global__ void BWExactSearchForwardGPU(uint8_t *W, uint64_t *nW, intmax_t *k, intmax_t *l, intmax_t k_ini, intmax_t l_ini, SA_TYPE *C, SA_TYPE *C1, comp_matrix O) {

	intmax_t i;
	intmax_t k2, l2;
	uintmax_t offset  = blockIdx.x * blockDim.x + threadIdx.x;

	__shared__ SA_TYPE Cshared[4];
	__shared__ SA_TYPE C1shared[4];

	if (threadIdx.x<4) {
		Cshared[threadIdx.x] = C[threadIdx.x];
		C1shared[threadIdx.x] = C1[threadIdx.x];
	}

	__syncthreads();

	k2 = k_ini;  l2 = l_ini;

	for (i=0; (k2<=l2) && (i<nW[offset]); i++)
		BWiterationGPU(k2, l2, k2, l2, W[offset*MAXLINE+i], Cshared, C1shared, O);

	k[offset] = k2;
	l[offset] = l2;

}

__global__ void BWExactSearchBackwardVectorGPU(uint8_t *W, uint64_t *nW, intmax_t *k, intmax_t *l, intmax_t k_ini, intmax_t l_ini, SA_TYPE *C, SA_TYPE *C1, comp_matrix O) {

	intmax_t i;
	intmax_t k2, l2;
	uintmax_t offset  = blockIdx.x * blockDim.x + threadIdx.x;

	__shared__ SA_TYPE Cshared[4];
	__shared__ SA_TYPE C1shared[4];

	if (threadIdx.x<4) {
		Cshared[threadIdx.x] = C[threadIdx.x];
		C1shared[threadIdx.x] = C1[threadIdx.x];
	}

	__syncthreads();

	k2 = k_ini;  l2 = l_ini;

	for (i=nW[offset]-1; i>=0; i--) {

		BWiterationGPU(k2, l2, k2, l2, W[offset*MAXLINE + i], Cshared, C1shared, O);

		k[offset*MAXLINE+i] = k2;
		l[offset*MAXLINE+i] = l2;

		if (k2 > l2) {
			i--;
			break;
		}

	}

	for(;i>=0; i--) {
		k[offset*MAXLINE+i] = k2;
		l[offset*MAXLINE+i] = l2;
	}

}

__global__ void BWExactSearchForwardVectorGPU(uint8_t *W, uint64_t *nW, intmax_t *k, intmax_t *l, intmax_t k_ini, intmax_t l_ini, SA_TYPE *C, SA_TYPE *C1, comp_matrix O) {

	intmax_t i;
	intmax_t k2, l2;
	uintmax_t offset  = blockIdx.x * blockDim.x + threadIdx.x;

	__shared__ SA_TYPE Cshared[4];
	__shared__ SA_TYPE C1shared[4];

	if (threadIdx.x<4) {
		Cshared[threadIdx.x] = C[threadIdx.x];
		C1shared[threadIdx.x] = C1[threadIdx.x];
	}

	__syncthreads();

	k2 = k_ini;  l2 = l_ini;

	for (i=0; i<nW[offset]; i++) {

		BWiterationGPU(k2, l2, k2, l2, W[offset*MAXLINE + i], Cshared, C1shared, O);

		k[offset*MAXLINE+i] = k2;
		l[offset*MAXLINE+i] = l2;

	}

}

//SOME TESTS:
/* __global__ void BWExactFinalResultsBackwardGPU(uint8_t *W, SA_TYPE *C, SA_TYPE *C1, comp_matrix O, results_list rl_prev, results_list rl_next, SA_TYPE chunk_size, SA_TYPE *stack_size) { */
/*   SA_TYPE k, l; */
/*   int16_t start, pos, pos_start, end; */
/*   unsigned read_index, read_offset; */
/*   SA_TYPE offset  = blockIdx.x * blockDim.x + threadIdx.x; */
/*   result *r_prev, *r_next; */

/*   __shared__ SA_TYPE Cshared[4]; */
/*   __shared__ SA_TYPE C1shared[4]; */

/*   if (offset==0) */
/*     *stack_size=0; */

/*   if (threadIdx.x<4) { */
/*     Cshared[threadIdx.x] = C[threadIdx.x]; */
/*     C1shared[threadIdx.x] = C1[threadIdx.x]; */
/*   } */

/*   __syncthreads(); */

/*   r_prev = &rl_prev.list[offset]; */

/*   start      = r_prev->start; */
/*   pos        = r_prev->pos; */
/*   end        = r_prev->end; */
/*   k          = r_prev->k; */
/*   l          = r_prev->l; */
/*   read_index = r_prev->read_index; */
/*   read_offset = /\*read_index*\/offset*MAXLINE; */

/*   pos_start = pos - chunk_size; */
/*   if (pos_start < start) pos_start = start; */

/*   for(; pos>=pos_start; pos--) { */
/*     BWiterationGPU(k, l, k, l, W[read_offset + pos], Cshared, C1shared, O); */
/*     if (k > l) { */
/*       pos=start-1; break; */
/*     } */
/*   } */

/*   r_next = &rl_next.list[/\*atomicAdd(stack_size,1)*\/offset]; */

/*   r_next->start = start; */
/*   r_next->pos = pos; */
/*   r_next->end = end; */
/*   r_next->k = k; */
/*   r_next->l = l; */
/*   r_next->read_index = read_index; */

/* } */

///////////////////////////////////////////////////MULTI-ERROR///////////////////////////////////////////////////////////
//__global__ void init_listsGPU(blocked_results_lists rl_prev, blocked_results_lists rl_next, blocked_results_lists rl_prev_i, blocked_results_lists rl_next_i, blocked_results_lists rl_final) {
//
//	SA_TYPE offset  = blockIdx.x * blockDim.x + threadIdx.x;
//
//	rl_prev.num_results[offset]   = 0;
//	rl_next.num_results[offset]   = 0;
//	rl_prev_i.num_results[offset] = 0;
//	rl_next_i.num_results[offset] = 0;
//	rl_final.num_results[offset]  = 0;
//
//}
//
//__global__ void BWExactSearchBackwardBlockedGPU(uint8_t *W, SA_TYPE nW, SA_TYPE *C, SA_TYPE *C1, comp_matrix O, result r, blocked_results_lists rl_prev, uintmax_t max_results) {
//
//	SA_TYPE offset  = blockIdx.x * blockDim.x + threadIdx.x;
//
//	SA_TYPE k2, l2;
//	result *r_iterator;
//	intmax_t i;
//
//	__shared__ SA_TYPE Cshared[4];
//	__shared__ SA_TYPE C1shared[4];
//
//	if (threadIdx.x<4) {
//		Cshared[threadIdx.x] = C[threadIdx.x];
//		C1shared[threadIdx.x] = C1[threadIdx.x];
//	}
//
//	__syncthreads();
//
//	k2 = r.k;
//	l2 = r.l;
//
//	for (i=r.pos; (k2<=l2) && (i>=r.start); i--)
//		BWiterationGPU(k2, l2, k2, l2, W[offset*MAXLINE+i], Cshared, C1shared, O);
//
//	if (k2 <= l2) {
//		r_iterator = rl_prev.results + offset * max_results;
//
//		init_result(r_iterator, 0);
//		bound_result(r_iterator, 0, r.end);
//		change_result(r_iterator, k2, l2, i);
//
//		rl_prev.num_results[offset] = 1; //Init the number of results
//	} else {
//		rl_prev.num_results[offset] = 0; //Init the number of results
//	}
//
//}
//
//__global__ void BWExactSearchForwardBlockedGPU(uint8_t *W, SA_TYPE nW, SA_TYPE *C, SA_TYPE *C1, comp_matrix O, result r, blocked_results_lists rl_prev, uintmax_t max_results) {
//
//	SA_TYPE offset  = blockIdx.x * blockDim.x + threadIdx.x;
//
//	SA_TYPE k2, l2;
//	result *r_iterator;
//	intmax_t i;
//
//	__shared__ SA_TYPE Cshared[4];
//	__shared__ SA_TYPE C1shared[4];
//
//	if (threadIdx.x<4) {
//		Cshared[threadIdx.x] = C[threadIdx.x];
//		C1shared[threadIdx.x] = C1[threadIdx.x];
//	}
//
//	__syncthreads();
//
//	k2 = r.k;
//	l2 = r.l;
//
//	for (i=r.pos; (k2<=l2) && (i<=r.end); i++)
//		BWiterationGPU(k2, l2, k2, l2, W[offset*MAXLINE+i], Cshared, C1shared, O);
//
//	if (k2 <= l2) {
//		r_iterator = rl_prev.results + offset * max_results;
//
//		init_result(r_iterator, 1);
//		bound_result(r_iterator, r.start, nW-1);
//		change_result(r_iterator, k2, l2, i);
//
//		rl_prev.num_results[offset] = 1; //Init the number of results
//	} else {
//		rl_prev.num_results[offset] = 0; //Init the number of results
//	}
//
//}
//
//__global__ void BWExactFinalResultsBackwardBlockedGPU(uint8_t *W, SA_TYPE *C, SA_TYPE *C1, comp_matrix O, blocked_results_lists rl_prev, blocked_results_lists rl_next_i, int16_t block_size, int16_t last_block, uintmax_t max_results) {
//
//	SA_TYPE offset  = blockIdx.x * blockDim.x + threadIdx.x;
//	result *rl_next_i_results = rl_next_i.results + offset * max_results;
//	SA_TYPE rl_next_i_num_results = rl_next_i.num_results[offset];
//
//	SA_TYPE k, l;
//	int16_t start, pos, current_block;
//
//	__shared__ SA_TYPE Cshared[4];
//	__shared__ SA_TYPE C1shared[4];
//
//	if (threadIdx.x<4) {
//		Cshared[threadIdx.x] = C[threadIdx.x];
//		C1shared[threadIdx.x] = C1[threadIdx.x];
//	}
//
//	result *r_iterator;
//
//	for (uintmax_t ii=0; ii < rl_prev.num_results[offset]; ii++) {
//
//		r_iterator = rl_prev.results + offset * max_results + ii;
//
//		start = r_iterator->start;
//		pos   = r_iterator->pos;
//
//		k = r_iterator->k;
//		l = r_iterator->l;
//
//		current_block = pos / block_size;
//
//		if ((current_block < last_block) || (pos == start-1)) { // Current block will be always >= start and previous results are propagated
//
//		} else {
//
//			if (current_block > last_block) { //Not in last previsited block
//
//				continue;
//
//			} else { //I am in the last previsited block
//
//				if ((pos + 1) % block_size) { //I am not in the first element of the block
//				} else { //I am in the first element in the block (all the block must be processed)
//					continue;
//				}
//
//			}
//
//		}
//
//		__syncthreads();
//
//		for(int16_t i=pos; i>=start; i--) {
//			BWiterationGPU(k, l, k, l, W[offset*MAXLINE+i], Cshared, C1shared, O);
//			if (k > l) break;
//		}
//
//		if (k <= l) {
//			change_result(r_iterator, k, l, start-1);
//			add_resultGPU(r_iterator, rl_next_i_results, rl_next_i_num_results, max_results);
//		}
//
//	} //r_prev
//
//	rl_next_i.num_results[offset] = rl_next_i_num_results;
//
//}
//
//__global__ void BWExactFinalResultsForwardBlockedGPU(uint8_t *W, SA_TYPE *C, SA_TYPE *C1, comp_matrix O, blocked_results_lists rl_prev, blocked_results_lists rl_next_i, int16_t block_size, int16_t last_block, uintmax_t max_results) {
//
//	SA_TYPE offset  = blockIdx.x * blockDim.x + threadIdx.x;
//	result *rl_next_i_results = rl_next_i.results + offset * max_results;
//	SA_TYPE rl_next_i_num_results = rl_next_i.num_results[offset];
//
//	SA_TYPE k, l;
//	int16_t pos, end;
//	int16_t current_block;
//
//	__shared__ SA_TYPE Cshared[4];
//	__shared__ SA_TYPE C1shared[4];
//
//	if (threadIdx.x<4) {
//		Cshared[threadIdx.x] = C[threadIdx.x];
//		C1shared[threadIdx.x] = C1[threadIdx.x];
//	}
//
//	result *r_iterator;
//
//	for (uintmax_t ii=0; ii < rl_prev.num_results[offset]; ii++) {
//
//		r_iterator = rl_prev.results + offset * max_results + ii;
//
//		pos   = r_iterator->pos;
//		end   = r_iterator->end;
//
//		k = r_iterator->k;
//		l = r_iterator->l;
//
//		current_block = pos / block_size;
//
//		if ( (current_block > last_block) || (pos == end+1) ) { // Current block will be always <= end and previous results are propagated
//
//		} else {
//
//			if (current_block < last_block) { //Not in last previsited block
//
//				continue;
//
//			} else { //I am in the last previsited block
//
//				if (pos % block_size) { //I am not in the first element of the block
//				} else { //I am in the first element in the block (all the block must be processed)
//					continue;
//				}
//
//			}
//
//		}
//
//		__syncthreads();
//
//		for(int16_t i=pos; i<=end; i++) {
//			BWiterationGPU(k, l, k, l, W[offset*MAXLINE+i], Cshared, C1shared, O);
//			if (k > l) break;
//		}
//
//		if (k <= l) {
//			change_result(r_iterator, k, l, end+1);
//			add_resultGPU(r_iterator, rl_next_i_results, rl_next_i_num_results, max_results);
//		}
//
//	} //r_prev
//
//	rl_next_i.num_results[offset] = rl_next_i_num_results;
//
//}
//
//__device__ uintmax_t BWExactFinalResultForwardGPU(uint8_t *W, SA_TYPE *C, SA_TYPE *C1, comp_matrix *O, result *r_iterator, result *rl_next_i_results, uintmax_t rl_next_i_num_results, int16_t block_size, int16_t last_block, uintmax_t max_results) {
//
//	SA_TYPE k, l;
//	int16_t pos, end;
//	int16_t current_block;
//
//	pos   = r_iterator->pos;
//	end   = r_iterator->end;
//
//	k = r_iterator->k;
//	l = r_iterator->l;
//
//	current_block = pos / block_size;
//
//	if ( (current_block > last_block) || (pos == end+1) ) { // Current block will be always <= end and previous results are propagated
//
//	} else {
//
//		if (current_block < last_block) { //Not in last previsited block
//
//			return rl_next_i_num_results;
//
//		} else { //I am in the last previsited block
//
//			if (pos % block_size) { //I am not in the first element of the block
//			} else { //I am in the first element in the block (all the block must be processed)
//				return rl_next_i_num_results;
//			}
//
//		}
//
//	}
//
//	for(int16_t i=pos; i<=end; i++) {
//		BWiterationGPU(k, l, k, l, W[i], C, C1, *O);
//		if (k > l) break;
//	}
//
//	if (k <= l) {
//		change_result(r_iterator, k, l, end+1);
//		add_resultGPU(r_iterator, rl_next_i_results, rl_next_i_num_results, max_results);
//	}
//
//	return rl_next_i_num_results;
//
//}
//
//__global__ void BWBranchFinalResultsForwardBlockedGPU(uint8_t *W, SA_TYPE *C, SA_TYPE *C1, comp_matrix O, blocked_results_lists rl_prev, blocked_results_lists rl_next, int16_t block_size, int16_t last_block, intmax_t max_results, uint8_t nA) {
//
//	SA_TYPE offset  = blockIdx.x * blockDim.x + threadIdx.x;
//	result *rl_next_results = rl_next.results + offset * max_results;
//	uintmax_t rl_next_num_results = rl_next.num_results[offset];
//
//	SA_TYPE k, l, k_aux, l_aux;
//	int16_t end, pos;
//	int16_t r_num_mismatches;
//	bool no_previous;
//	int16_t last_err_pos;
//	uint8_t last_err_kind;
//	uint8_t last_err_base;
//
//	result *r_iterator;
//
//	for (uintmax_t ii=0; ii < rl_prev.num_results[offset]; ii++) {
//
//		r_iterator = rl_prev.results + offset * max_results + ii;
//
//		if (r_iterator->dir != 1) continue;
//
//		end = r_iterator->end;
//		pos = r_iterator->pos;
//
//		if (pos > end) {
//			add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//			continue;
//		}
//
//		no_previous = true;
//		r_num_mismatches = r_iterator->num_mismatches-1;
//		if (r_num_mismatches>-1) {
//			last_err_pos  = r_iterator->err_pos[r_num_mismatches];
//			last_err_kind = r_iterator->err_kind[r_num_mismatches];
//			last_err_base = r_iterator->err_base[r_num_mismatches];
//		} else {
//			last_err_pos  = -10;
//			last_err_kind = 0;
//			last_err_base = (uint8_t) -1;
//		}
//
//		k = r_iterator->k;
//		l = r_iterator->l;
//
//		add_mismatch(r_iterator, DELETION, (uint8_t) -1, pos);
//
//		if (last_err_pos == pos - 1) { //Previous MISMATCH or DELETION
//
//			if (last_err_kind == MISMATCH) { //Previous MISMATCH
//
//				//Deletion
//				if (W[offset*MAXLINE+pos]!=last_err_base) {
//					change_result(r_iterator, k, l, pos+1);
//					rl_next_num_results = BWExactFinalResultForwardGPU(W, C, C1, &O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//				}
//
//				for (uint8_t b=0;b<nA;b++) {
//
//					BWiterationGPU(k, l, k_aux, l_aux, b, C, C1, O);
//
//					if (k_aux > l_aux) continue;
//
//					//Insertion
//					if (b!=W[offset*MAXLINE+last_err_pos]) {
//						change_result(r_iterator, k_aux, l_aux, pos);
//						modify_last_mismatch2(r_iterator, INSERTION, b);
//						rl_next_num_results = BWExactFinalResultForwardGPU(W, C, C1, &O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//					}
//
//					//Mismatch
//					if (b!=W[offset*MAXLINE+pos]) {
//						change_result(r_iterator, k_aux, l_aux, pos+1);
//						modify_last_mismatch2(r_iterator, MISMATCH, b);
//						rl_next_num_results = BWExactFinalResultForwardGPU(W, C, C1, &O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//					}
//
//				}
//
//				no_previous = false;
//
//			} else if (last_err_kind == DELETION) { //Previous DELETION
//
//				//Deletion
//				change_result(r_iterator, k, l, pos+1);
//				rl_next_num_results = BWExactFinalResultForwardGPU(W, C, C1, &O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//
//				for (uint8_t b=0;b<nA;b++) {
//
//					BWiterationGPU(k, l, k_aux, l_aux, b, C, C1, O);
//
//					if (k_aux > l_aux) continue;
//
//					// NO INSERTION
//
//					//Mismatch
//					if (b!=W[offset*MAXLINE+pos]) {
//
//						if (b!=W[offset*MAXLINE+last_err_pos]) {
//							change_result(r_iterator, k_aux, l_aux, pos+1);
//							modify_last_mismatch2(r_iterator, MISMATCH, b);
//							rl_next_num_results = BWExactFinalResultForwardGPU(W, C, C1, &O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//						}
//
//					}
//
//				}
//
//				no_previous = false;
//
//			}
//
//		} else if (last_err_pos == pos) { //Previous INSERTION
//
//			//NO DELETION
//
//			for (uint8_t b=0;b<nA;b++) {
//
//				BWiterationGPU(k, l, k_aux, l_aux, b, C, C1, O);
//
//				if (k_aux > l_aux) continue;
//
//				//Insertion
//				change_result(r_iterator, k_aux, l_aux, pos);
//				modify_last_mismatch2(r_iterator, INSERTION, b);
//				rl_next_num_results = BWExactFinalResultForwardGPU(W, C, C1, &O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//
//				//Mismatch
//				if (b!=W[offset*MAXLINE+pos]) {
//
//					if (W[offset*MAXLINE+pos]!=last_err_base) {
//						r_iterator->pos = pos+1;
//						modify_last_mismatch1(r_iterator, MISMATCH);
//						rl_next_num_results = BWExactFinalResultForwardGPU(W, C, C1, &O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//					}
//
//				}
//
//			}
//
//			no_previous = false;
//
//		}
//
//		if (no_previous) { //Previous MATCH
//
//			//Deletion
//			change_result(r_iterator, k, l, pos+1);
//			rl_next_num_results = BWExactFinalResultForwardGPU(W, C, C1, &O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//
//			for (uint8_t b=0;b<nA;b++) {
//
//				BWiterationGPU(k, l, k_aux, l_aux, b, C, C1, O);
//
//				if (k_aux > l_aux) continue;
//
//				//Insertion
//				change_result(r_iterator, k_aux, l_aux, pos);
//				modify_last_mismatch2(r_iterator, INSERTION, b);
//				rl_next_num_results = BWExactFinalResultForwardGPU(W, C, C1, &O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//
//				if (b!=W[offset*MAXLINE+pos]) { //Mismatch
//					r_iterator->pos = pos+1;
//					modify_last_mismatch1(r_iterator, MISMATCH);
//					rl_next_num_results = BWExactFinalResultForwardGPU(W, C, C1, &O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//				}
//
//			}
//
//		}
//
//	}
//
//	rl_prev.num_results[offset] = 0;
//	rl_next.num_results[offset] = rl_next_num_results;
//
//}
//
//uintmax_t BWExactFinalResultForwardCPU(uint8_t *W, vector *C, vector *C1, comp_matrix *O, result *r_iterator, result *rl_next_i_results, uintmax_t rl_next_i_num_results, int16_t block_size, int16_t last_block, uintmax_t max_results) {
//
//	SA_TYPE k, l;
//	int16_t pos, end;
//	int16_t current_block;
//
//	pos   = r_iterator->pos;
//	end   = r_iterator->end;
//
//	k = r_iterator->k;
//	l = r_iterator->l;
//
//	current_block = pos / block_size;
//
//	if ( (current_block > last_block) || (pos == end+1) ) { // Current block will be always <= end and previous results are propagated
//
//	} else {
//
//		if (current_block < last_block) { //Not in last previsited block
//
//			return rl_next_i_num_results;
//
//		} else { //I am in the last previsited block
//
//			if (pos % block_size) { //I am not in the first element of the block
//			} else { //I am in the first element in the block (all the block must be processed)
//				return rl_next_i_num_results;
//			}
//
//		}
//
//	}
//
//	for(int16_t i=pos; i<=end; i++) {
//		BWiteration(k, l, k, l, W[i], C, C1, O);
//		if (k > l) break;
//	}
//
//	if (k <= l) {
//		change_result(r_iterator, k, l, end+1);
//		add_resultGPU(r_iterator, rl_next_i_results, rl_next_i_num_results, max_results);
//	}
//
//	return rl_next_i_num_results;
//
//}
//
//void BWBranchFinalResultsForwardBlockedCPU(uintmax_t num_bloques, uintmax_t tam_bloques, uint8_t *W, vector *C, vector *C1, comp_matrix *O, blocked_results_lists rl_prev, blocked_results_lists rl_next, int16_t block_size, int16_t last_block, uintmax_t max_results, uint8_t nA) {
//
//	//#pragma omp parallel for
//	for(SA_TYPE offset=0; offset < num_bloques*tam_bloques; offset++) {
//
//		result *rl_next_results = rl_next.results + offset * max_results;
//		uintmax_t rl_next_num_results = rl_next.num_results[offset];
//
//		SA_TYPE k, l, k_aux, l_aux;
//		int16_t end, pos;
//		int16_t r_num_mismatches;
//		bool no_previous;
//		int16_t last_err_pos;
//		uint8_t last_err_kind;
//		uint8_t last_err_base;
//
//		result *r_iterator;
//
//		for (uintmax_t ii=0; ii < rl_prev.num_results[offset]; ii++) {
//
//			r_iterator = rl_prev.results + offset * max_results + ii;
//
//			if (r_iterator->dir != 1) continue;
//
//			end = r_iterator->end;
//			pos = r_iterator->pos;
//
//			if (pos > end) {
//				add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//				continue;
//			}
//
//			no_previous = true;
//			r_num_mismatches = r_iterator->num_mismatches-1;
//			if (r_num_mismatches>-1) {
//				last_err_pos  = r_iterator->err_pos[r_num_mismatches];
//				last_err_kind = r_iterator->err_kind[r_num_mismatches];
//				last_err_base = r_iterator->err_base[r_num_mismatches];
//			} else {
//				last_err_pos  = -10;
//				last_err_kind = 0;
//				last_err_base = (uint8_t) -1;
//			}
//
//			k = r_iterator->k;
//			l = r_iterator->l;
//
//			add_mismatch(r_iterator, DELETION, (uint8_t) -1, pos);
//
//			if (last_err_pos == pos - 1) { //Previous MISMATCH or DELETION
//
//				if (last_err_kind == MISMATCH) { //Previous MISMATCH
//
//					//Deletion
//					if (W[offset*MAXLINE+pos]!=last_err_base) {
//						change_result(r_iterator, k, l, pos+1);
//						rl_next_num_results = BWExactFinalResultForwardCPU(W, C, C1, O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//					}
//
//					for (uint8_t b=0;b<nA;b++) {
//
//						BWiteration(k, l, k_aux, l_aux, b, C, C1, O);
//
//						if (k_aux > l_aux) continue;
//
//						//Insertion
//						if (b!=W[offset*MAXLINE+last_err_pos]) {
//							change_result(r_iterator, k_aux, l_aux, pos);
//							modify_last_mismatch2(r_iterator, INSERTION, b);
//							rl_next_num_results = BWExactFinalResultForwardCPU(W, C, C1, O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//						}
//
//						//Mismatch
//						if (b!=W[offset*MAXLINE+pos]) {
//							change_result(r_iterator, k_aux, l_aux, pos+1);
//							modify_last_mismatch2(r_iterator, MISMATCH, b);
//							rl_next_num_results = BWExactFinalResultForwardCPU(W, C, C1, O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//						}
//
//					}
//
//					no_previous = false;
//
//				} else if (last_err_kind == DELETION) { //Previous DELETION
//
//					//Deletion
//					change_result(r_iterator, k, l, pos+1);
//					rl_next_num_results = BWExactFinalResultForwardCPU(W, C, C1, O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//
//					for (uint8_t b=0;b<nA;b++) {
//
//						BWiteration(k, l, k_aux, l_aux, b, C, C1, O);
//
//						if (k_aux > l_aux) continue;
//
//						// NO INSERTION
//
//						//Mismatch
//						if (b!=W[offset*MAXLINE+pos]) {
//
//							if (b!=W[offset*MAXLINE+last_err_pos]) {
//								change_result(r_iterator, k_aux, l_aux, pos+1);
//								modify_last_mismatch2(r_iterator, MISMATCH, b);
//								rl_next_num_results = BWExactFinalResultForwardCPU(W, C, C1, O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//							}
//
//						}
//
//					}
//
//					no_previous = false;
//
//				}
//
//			} else if (last_err_pos == pos) { //Previous INSERTION
//
//				//NO DELETION
//
//				for (uint8_t b=0;b<nA;b++) {
//
//					BWiteration(k, l, k_aux, l_aux, b, C, C1, O);
//
//					if (k_aux > l_aux) continue;
//
//					//Insertion
//					change_result(r_iterator, k_aux, l_aux, pos);
//					modify_last_mismatch2(r_iterator, INSERTION, b);
//					rl_next_num_results = BWExactFinalResultForwardCPU(W, C, C1, O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//
//					//Mismatch
//					if (b!=W[offset*MAXLINE+pos]) {
//
//						if (W[offset*MAXLINE+pos]!=last_err_base) {
//							r_iterator->pos = pos+1;
//							modify_last_mismatch1(r_iterator, MISMATCH);
//							rl_next_num_results = BWExactFinalResultForwardCPU(W, C, C1, O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//						}
//
//					}
//
//				}
//
//				no_previous = false;
//
//			}
//
//			if (no_previous) { //Previous MATCH
//
//				//Deletion
//				change_result(r_iterator, k, l, pos+1);
//				rl_next_num_results = BWExactFinalResultForwardCPU(W, C, C1, O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//
//				for (uint8_t b=0;b<nA;b++) {
//
//					BWiteration(k, l, k_aux, l_aux, b, C, C1, O);
//
//					if (k_aux > l_aux) continue;
//
//					//Insertion
//					change_result(r_iterator, k_aux, l_aux, pos);
//					modify_last_mismatch2(r_iterator, INSERTION, b);
//					rl_next_num_results = BWExactFinalResultForwardCPU(W, C, C1, O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//
//					if (b!=W[offset*MAXLINE+pos]) { //Mismatch
//						r_iterator->pos = pos+1;
//						modify_last_mismatch1(r_iterator, MISMATCH);
//						rl_next_num_results = BWExactFinalResultForwardCPU(W, C, C1, O, r_iterator, rl_next_results, rl_next_num_results, block_size, last_block, max_results);
//					}
//
//				}
//
//			}
//
//		}
//
//		rl_prev.num_results[offset] = 0;
//		rl_next.num_results[offset] = rl_next_num_results;
//
//	}
//
//}
//
//__global__ void BWExactPartialResultsBackwardBlockedGPU(uint8_t *W, SA_TYPE *C, SA_TYPE *C1, comp_matrix O, blocked_results_lists rl_prev, blocked_results_lists rl_next, blocked_results_lists rl_next_i, int16_t block_size, int16_t last_block, uintmax_t max_results) {
//
//	SA_TYPE offset  = blockIdx.x * blockDim.x + threadIdx.x;
//	result *rl_next_results = rl_next.results + offset * max_results;
//	uintmax_t rl_next_num_results = 0;
//	result *rl_next_i_results = rl_next_i.results + offset * max_results;
//	uintmax_t rl_next_i_num_results = rl_next_i.num_results[offset];
//
//	SA_TYPE k, l, k_next, l_next;
//	int16_t start, pos;
//	int16_t current_block, last_block_pos;
//	bool complete_search;
//
//	__shared__ SA_TYPE Cshared[4];
//	__shared__ SA_TYPE C1shared[4];
//
//	if (threadIdx.x<4) {
//		Cshared[threadIdx.x] = C[threadIdx.x];
//		C1shared[threadIdx.x] = C1[threadIdx.x];
//	}
//
//	result *r_iterator;
//	SA_TYPE results, results_next;
//
//	for (uintmax_t ii=0; ii < rl_prev.num_results[offset]; ii++) {
//
//		r_iterator = rl_prev.results + offset * max_results + ii;
//
//		start  = r_iterator->start;
//		pos    = r_iterator->pos;
//
//		k_next = r_iterator->k;
//		l_next = r_iterator->l;
//		results_next = l_next - k_next;
//
//		current_block = pos / block_size;
//
//		if ((current_block < last_block) || (pos == start-1)) { // Current block will be always >= start and previous results are propagated
//
//			last_block_pos = start;
//			complete_search = true;
//
//		} else {
//
//			if (current_block > last_block) { //Not in last previsited block
//
//				if ((pos + 1) % block_size) { //Not in first element of the block
//					last_block_pos = (current_block-1) * block_size;
//				} else { //I am in the first element in the block (all the block must be processed)
//					last_block_pos = current_block * block_size;
//				}
//
//				complete_search = false;
//
//			} else { //I am in the last previsited block
//
//				if ((pos + 1) % block_size) { //I am not in the first element of the block
//					last_block_pos = start;
//					complete_search = true;
//				} else { //I am in the first element of the block (all the block must be processed)
//					last_block_pos = current_block * block_size;
//					complete_search = false;
//				}
//
//			}
//
//		}
//
//		__syncthreads();
//
//		for(int16_t i=pos; i>=last_block_pos; i--) {
//
//			k = k_next;
//			l = l_next;
//
//			if (k > l) break;
//
//			BWiterationGPU(k, l, k_next, l_next, W[offset*MAXLINE+i], Cshared, C1shared, O);
//			results      = results_next;
//			results_next = l_next - k_next;
//			if (results == results_next) continue;
//
//			change_result(r_iterator, k, l, i);
//			add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//
//		}
//
//		if (complete_search && k_next <= l_next) {
//			change_result(r_iterator, k_next, l_next, start-1);
//			add_resultGPU(r_iterator, rl_next_i_results, rl_next_i_num_results, max_results);
//		}
//
//	} //r_prev
//
//	rl_next.num_results[offset] = rl_next_num_results;
//	rl_next_i.num_results[offset] = rl_next_i_num_results;
//
//}
//
//__global__ void BWExactPartialResultsForwardBlockedGPU(uint8_t *W, SA_TYPE *C, SA_TYPE *C1, comp_matrix O, blocked_results_lists rl_prev, blocked_results_lists rl_next, blocked_results_lists rl_next_i, int16_t block_size, int16_t last_block, uintmax_t max_results) {
//
//	SA_TYPE offset  = blockIdx.x * blockDim.x + threadIdx.x;
//	result *rl_next_results = rl_next.results + offset * max_results;
//	uintmax_t rl_next_num_results = 0;
//	result *rl_next_i_results = rl_next_i.results + offset * max_results;
//	uintmax_t rl_next_i_num_results = rl_next_i.num_results[offset];
//
//	SA_TYPE k, l, k_next, l_next;
//	int16_t pos, end;
//	int16_t current_block, last_block_pos;
//	SA_TYPE results, results_next;
//	bool complete_search;
//
//	__shared__ SA_TYPE Cshared[4];
//	__shared__ SA_TYPE C1shared[4];
//
//	if (threadIdx.x<4) {
//		Cshared[threadIdx.x] = C[threadIdx.x];
//		C1shared[threadIdx.x] = C1[threadIdx.x];
//	}
//
//	result *r_iterator;
//
//
//	for (uintmax_t ii=0; ii < rl_prev.num_results[offset]; ii++) {
//
//		r_iterator = rl_prev.results + offset * max_results + ii;
//
//		pos   = r_iterator->pos;
//		end   = r_iterator->end;
//
//		k_next = r_iterator->k;
//		l_next = r_iterator->l;
//		results_next = l_next - k_next;
//
//		current_block = pos / block_size;
//
//		if ( (current_block > last_block) || (pos == end+1) ) { // Current block will be always <= end and previous results are propagated
//
//			last_block_pos = end;
//			complete_search = true;
//
//		} else {
//
//			if (current_block < last_block) { //Not in last previsited block
//
//				if (pos % block_size) { //Not in first element of the block
//					last_block_pos = (current_block+2) * block_size - 1;
//				} else { //I am in the first element in the block (all the block must be processed)
//					last_block_pos = (current_block+1) * block_size - 1;
//				}
//
//				complete_search = false;
//
//			} else { //I am in the last previsited block
//
//				if (pos % block_size) { //I am not in the first element of the block
//					last_block_pos = end;
//					complete_search = true;
//				} else { //I am in the first element in the block (all the block must be processed)
//					last_block_pos = (current_block+1) * block_size - 1;
//					complete_search = false;
//				}
//
//			}
//
//		}
//
//		__syncthreads();
//
//		for(int16_t i=pos; i<=last_block_pos; i++) {
//
//			k = k_next;
//			l = l_next;
//
//			if (k > l) break;
//
//			BWiterationGPU(k, l, k_next, l_next, W[offset*MAXLINE+i], Cshared, C1shared, O);
//			results      = results_next;
//			results_next = l_next - k_next;
//			if (results == results_next) continue;
//
//			change_result(r_iterator, k, l, i);
//			add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//		}
//
//		if (complete_search && k_next <= l_next) {
//			change_result(r_iterator, k_next, l_next, end+1);
//			add_resultGPU(r_iterator, rl_next_i_results, rl_next_i_num_results, max_results);
//		}
//
//	} //r_prev
//
//	rl_next.num_results[offset] = rl_next_num_results;
//	rl_next_i.num_results[offset] = rl_next_i_num_results;
//
//}
//
//__global__ void BWBranchPartialResultsBackwardBlockedGPU(uint8_t *W, SA_TYPE *C, SA_TYPE *C1, comp_matrix O, blocked_results_lists rl_prev, blocked_results_lists rl_next, uintmax_t max_results, uint8_t nA) {
//
//	SA_TYPE offset  = blockIdx.x * blockDim.x + threadIdx.x;
//	result *rl_next_results = rl_next.results + offset * max_results;
//	uintmax_t rl_next_num_results = 0;
//
//	SA_TYPE k, l, k_aux, l_aux;
//	int16_t start, pos;
//	int16_t r_num_mismatches;
//	bool no_previous;
//	int16_t last_err_pos;
//	uint8_t last_err_kind;
//	uint8_t last_err_base;
//
//	__shared__ SA_TYPE Cshared[4];
//	__shared__ SA_TYPE C1shared[4];
//
//	if (threadIdx.x<4) {
//		Cshared[threadIdx.x] = C[threadIdx.x];
//		C1shared[threadIdx.x] = C1[threadIdx.x];
//	}
//
//	result *r_iterator;
//
//	for (uintmax_t ii=0; ii < rl_prev.num_results[offset]; ii++) {
//
//		r_iterator = rl_prev.results + offset * max_results + ii;
//
//		if (r_iterator->dir != 0) continue;
//
//		start = r_iterator->start;
//		pos   = r_iterator->pos;
//
//		if (pos < start) {
//			add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//			continue;
//		}
//
//		no_previous = true;
//		r_num_mismatches = r_iterator->num_mismatches-1;
//		if (r_num_mismatches>-1) {
//			last_err_pos  = r_iterator->err_pos[r_num_mismatches];
//			last_err_kind = r_iterator->err_kind[r_num_mismatches];
//			last_err_base = r_iterator->err_base[r_num_mismatches];
//		} else {
//			last_err_pos  = -10;
//			last_err_kind = 0;
//			last_err_base = (uint8_t) -1;
//		}
//
//		k = r_iterator->k;
//		l = r_iterator->l;
//
//		add_mismatch(r_iterator, DELETION, (uint8_t) -1, pos);
//
//		__syncthreads();
//
//		if (last_err_pos == pos + 1) { //Previous MISMATCH or DELETION
//
//			if (last_err_kind == MISMATCH) { //Previous MISMATCH
//
//				//Deletion
//				if (W[offset*MAXLINE+pos]!=last_err_base) {
//					change_result(r_iterator, k, l, pos-1);
//					add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//				}
//
//				for (uint8_t b=0;b<nA;b++) {
//
//					BWiterationGPU(k, l, k_aux, l_aux, b, Cshared, C1shared, O);
//
//					if (k_aux > l_aux) continue;
//
//					//Insertion
//					if (b!=W[offset*MAXLINE+last_err_pos]) {
//						change_result(r_iterator, k_aux, l_aux, pos);
//						modify_last_mismatch2(r_iterator, INSERTION, b);
//						add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//					}
//
//					//Mismatch
//					if (b!=W[offset*MAXLINE+pos]) {
//						change_result(r_iterator, k_aux, l_aux, pos-1);
//						modify_last_mismatch2(r_iterator, MISMATCH, b);
//						add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//					}
//
//				}
//
//				no_previous = false;
//
//			} else if (last_err_kind == DELETION) { //Previous DELETION
//
//				//Deletion
//				change_result(r_iterator, k, l, pos-1);
//				add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//
//				for (uint8_t b=0;b<nA;b++) {
//
//					BWiterationGPU(k, l, k_aux, l_aux, b, Cshared, C1shared, O);
//
//					if (k_aux > l_aux) continue;
//
//					// NO INSERTION
//
//					if (b!=W[offset*MAXLINE+pos]) { //Mismatch
//
//						if (b!=W[offset*MAXLINE+last_err_pos]) {
//							change_result(r_iterator, k_aux, l_aux, pos-1);
//							modify_last_mismatch2(r_iterator, MISMATCH, b);
//							add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//						}
//
//					}
//
//				}
//
//				no_previous = false;
//
//			}
//
//		} else if (last_err_pos == pos) { //Previous INSERTION
//
//			//NO DELETION
//
//			for (uint8_t b=0;b<nA;b++) {
//
//				BWiterationGPU(k, l, k_aux, l_aux, b, Cshared, C1shared, O);
//
//				if (k_aux > l_aux) continue;
//
//				//Insertion
//				change_result(r_iterator, k_aux, l_aux, pos);
//				modify_last_mismatch2(r_iterator, INSERTION, b);
//				add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//
//				//Mismatch
//				if (b!=W[offset*MAXLINE+pos]) {
//
//					if (W[offset*MAXLINE+pos]!=last_err_base) {
//						r_iterator->pos = pos-1;
//						modify_last_mismatch1(r_iterator, MISMATCH);
//						add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//					}
//
//				}
//
//			}
//
//			no_previous = false;
//
//		}
//
//		if (no_previous) { //Previous MATCH
//
//			//Deletion
//			change_result(r_iterator, k, l, pos-1);
//			add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//
//			for (uint8_t b=0;b<nA;b++) {
//
//				BWiterationGPU(k, l, k_aux, l_aux, b, Cshared, C1shared, O);
//
//				if (k_aux > l_aux) continue;
//
//				//Insertion
//				change_result(r_iterator, k_aux, l_aux, pos);
//				modify_last_mismatch2(r_iterator, INSERTION, b);
//				add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//
//				//Mismatch
//				if (b!=W[offset*MAXLINE+pos]) {
//					r_iterator->pos = pos-1;
//					modify_last_mismatch1(r_iterator, MISMATCH);
//					add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//				}
//
//			}
//
//		}
//
//	}
//
//	rl_next.num_results[offset] = rl_next_num_results;
//
//}
//
//__global__ void BWBranchPartialResultsForwardBlockedGPU(uint8_t *W, SA_TYPE *C, SA_TYPE *C1, comp_matrix O, blocked_results_lists rl_prev, blocked_results_lists rl_next, uintmax_t max_results, uint8_t nA) {
//
//	SA_TYPE offset  = blockIdx.x * blockDim.x + threadIdx.x;
//	result *rl_next_results = rl_next.results + offset * max_results;
//	uintmax_t rl_next_num_results = 0;
//
//	SA_TYPE k, l, k_aux, l_aux;
//	int16_t end, pos;
//	int16_t r_num_mismatches;
//	bool no_previous;
//	int16_t last_err_pos;
//	uint8_t last_err_kind;
//	uint8_t last_err_base;
//
//	__shared__ SA_TYPE Cshared[4];
//	__shared__ SA_TYPE C1shared[4];
//
//	if (threadIdx.x<4) {
//		Cshared[threadIdx.x] = C[threadIdx.x];
//		C1shared[threadIdx.x] = C1[threadIdx.x];
//	}
//
//	result *r_iterator;
//
//	for (uintmax_t ii=0; ii < rl_prev.num_results[offset]; ii++) {
//
//		r_iterator = rl_prev.results + offset * max_results + ii;
//
//		__syncthreads();
//
//		if (r_iterator->dir != 1) continue;
//
//		end = r_iterator->end;
//		pos = r_iterator->pos;
//
//		if (pos > end) {
//			add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//			continue;
//		}
//
//		no_previous = true;
//		r_num_mismatches = r_iterator->num_mismatches-1;
//		if (r_num_mismatches>-1) {
//			last_err_pos  = r_iterator->err_pos[r_num_mismatches];
//			last_err_kind = r_iterator->err_kind[r_num_mismatches];
//			last_err_base = r_iterator->err_base[r_num_mismatches];
//		} else {
//			last_err_pos  = -10;
//			last_err_kind = 0;
//			last_err_base = (uint8_t) -1;
//		}
//
//		k = r_iterator->k;
//		l = r_iterator->l;
//
//		add_mismatch(r_iterator, DELETION, (uint8_t) -1, pos);
//
//		__syncthreads();
//
//		if (last_err_pos == pos - 1) { //Previous MISMATCH or DELETION
//
//			if (last_err_kind == MISMATCH) { //Previous MISMATCH
//
//				//Deletion
//				if (W[offset*MAXLINE+pos]!=last_err_base) {
//					change_result(r_iterator, k, l, pos+1);
//					add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//				}
//
//				for (uint8_t b=0;b<nA;b++) {
//
//					BWiterationGPU(k, l, k_aux, l_aux, b, Cshared, C1shared, O);
//
//					if (k_aux > l_aux) continue;
//
//					//Insertion
//					if (b!=W[offset*MAXLINE+last_err_pos]) {
//						change_result(r_iterator, k_aux, l_aux, pos);
//						modify_last_mismatch2(r_iterator, INSERTION, b);
//						add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//					}
//
//					//Mismatch
//					if (b!=W[offset*MAXLINE+pos]) {
//						change_result(r_iterator, k_aux, l_aux, pos+1);
//						modify_last_mismatch2(r_iterator, MISMATCH, b);
//						add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//					}
//
//				}
//
//				no_previous = false;
//
//			} else if (last_err_kind == DELETION) { //Previous DELETION
//
//				//Deletion
//				change_result(r_iterator, k, l, pos+1);
//				add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//
//				for (uint8_t b=0;b<nA;b++) {
//
//					BWiterationGPU(k, l, k_aux, l_aux, b, Cshared, C1shared, O);
//
//					if (k_aux > l_aux) continue;
//
//					// NO INSERTION
//
//					//Mismatch
//					if (b!=W[offset*MAXLINE+pos]) {
//
//						if (b!=W[offset*MAXLINE+last_err_pos]) {
//							change_result(r_iterator, k_aux, l_aux, pos+1);
//							modify_last_mismatch2(r_iterator, MISMATCH, b);
//							add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//						}
//
//					}
//
//				}
//
//				no_previous = false;
//
//			}
//
//		} else if (last_err_pos == pos) { //Previous INSERTION
//
//			//NO DELETION
//
//			for (uint8_t b=0;b<nA;b++) {
//
//				BWiterationGPU(k, l, k_aux, l_aux, b, Cshared, C1shared, O);
//
//				if (k_aux > l_aux) continue;
//
//				//Insertion
//				change_result(r_iterator, k_aux, l_aux, pos);
//				modify_last_mismatch2(r_iterator, INSERTION, b);
//				add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//
//				//Mismatch
//				if (b!=W[offset*MAXLINE+pos]) {
//
//					if (W[offset*MAXLINE+pos]!=last_err_base) {
//						r_iterator->pos = pos+1;
//						modify_last_mismatch1(r_iterator, MISMATCH);
//						add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//					}
//
//				}
//
//			}
//
//			no_previous = false;
//
//		}
//
//		if (no_previous) { //Previous MATCH
//
//			//Deletion
//			change_result(r_iterator, k, l, pos+1);
//			add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//
//			for (uint8_t b=0;b<nA;b++) {
//
//				BWiterationGPU(k, l, k_aux, l_aux, b, Cshared, C1shared, O);
//
//				if (k_aux > l_aux) continue;
//
//				//Insertion
//				change_result(r_iterator, k_aux, l_aux, pos);
//				modify_last_mismatch2(r_iterator, INSERTION, b);
//				add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//
//				if (b!=W[offset*MAXLINE+pos]) { //Mismatch
//					r_iterator->pos = pos+1;
//					modify_last_mismatch1(r_iterator, MISMATCH);
//					add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//
//				}
//
//			}
//
//		}
//
//	}
//
//	rl_next.num_results[offset] = rl_next_num_results;
//
//}
//
//void BWBranchPartialResultsForwardBlockedCPU(uintmax_t num_bloques, uintmax_t tam_bloques, uint8_t *W, vector *C, vector *C1, comp_matrix *O, blocked_results_lists rl_prev, blocked_results_lists rl_next, uintmax_t max_results, uint8_t nA) {
//
//	//#pragma omp parallel for
//	for(SA_TYPE offset=0; offset < num_bloques*tam_bloques; offset++) {
//
//		SA_TYPE k, l, k_aux, l_aux;
//		int16_t end, pos;
//		int16_t r_num_mismatches;
//		bool no_previous;
//		int16_t last_err_pos;
//		uint8_t last_err_kind;
//		uint8_t last_err_base;
//
//		result *r_iterator;
//
//		result *rl_next_results = rl_next.results + offset * max_results;
//		uintmax_t rl_next_num_results = 0;
//
//		for (uintmax_t ii=0; ii < rl_prev.num_results[offset]; ii++) {
//
//			r_iterator = rl_prev.results + offset * max_results + ii;
//
//			if (r_iterator->dir != 1) continue;
//
//			end = r_iterator->end;
//			pos = r_iterator->pos;
//
//			if (pos > end) {
//				add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//				continue;
//			}
//
//			no_previous = true;
//			r_num_mismatches = r_iterator->num_mismatches-1;
//			if (r_num_mismatches>-1) {
//				last_err_pos  = r_iterator->err_pos[r_num_mismatches];
//				last_err_kind = r_iterator->err_kind[r_num_mismatches];
//				last_err_base = r_iterator->err_base[r_num_mismatches];
//			} else {
//				last_err_pos  = -10;
//				last_err_kind = 0;
//				last_err_base = (uint8_t) -1;
//			}
//
//			k = r_iterator->k;
//			l = r_iterator->l;
//
//			add_mismatch(r_iterator, DELETION, (uint8_t) -1, pos);
//
//			if (last_err_pos == pos - 1) { //Previous MISMATCH or DELETION
//
//				if (last_err_kind == MISMATCH) { //Previous MISMATCH
//
//					//Deletion
//					if (W[offset*MAXLINE+pos]!=last_err_base) {
//						change_result(r_iterator, k, l, pos+1);
//						add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//					}
//
//					for (uint8_t b=0;b<nA;b++) {
//
//						BWiteration(k, l, k_aux, l_aux, b, C, C1, O);
//
//						if (k_aux > l_aux) continue;
//
//						//Insertion
//						if (b!=W[offset*MAXLINE+last_err_pos]) {
//							change_result(r_iterator, k_aux, l_aux, pos);
//							modify_last_mismatch2(r_iterator, INSERTION, b);
//							add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//						}
//
//						//Mismatch
//						if (b!=W[offset*MAXLINE+pos]) {
//							change_result(r_iterator, k_aux, l_aux, pos+1);
//							modify_last_mismatch2(r_iterator, MISMATCH, b);
//							add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//						}
//
//					}
//
//					no_previous = false;
//
//				} else if (last_err_kind == DELETION) { //Previous DELETION
//
//					//Deletion
//					change_result(r_iterator, k, l, pos+1);
//					add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//
//					for (uint8_t b=0;b<nA;b++) {
//
//						BWiteration(k, l, k_aux, l_aux, b, C, C1, O);
//
//						if (k_aux > l_aux) continue;
//
//						// NO INSERTION
//
//						//Mismatch
//						if (b!=W[offset*MAXLINE+pos]) {
//
//							if (b!=W[offset*MAXLINE+last_err_pos]) {
//								change_result(r_iterator, k_aux, l_aux, pos+1);
//								modify_last_mismatch2(r_iterator, MISMATCH, b);
//								add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//							}
//
//						}
//
//					}
//
//					no_previous = false;
//
//				}
//
//			} else if (last_err_pos == pos) { //Previous INSERTION
//
//				//NO DELETION
//
//				for (uint8_t b=0;b<nA;b++) {
//
//					BWiteration(k, l, k_aux, l_aux, b, C, C1, O);
//
//					if (k_aux > l_aux) continue;
//
//					//Insertion
//					change_result(r_iterator, k_aux, l_aux, pos);
//					modify_last_mismatch2(r_iterator, INSERTION, b);
//					add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//
//					//Mismatch
//					if (b!=W[offset*MAXLINE+pos]) {
//
//						if (W[offset*MAXLINE+pos]!=last_err_base) {
//							r_iterator->pos = pos+1;
//							modify_last_mismatch1(r_iterator, MISMATCH);
//							add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//						}
//
//					}
//
//				}
//
//				no_previous = false;
//
//			}
//
//			if (no_previous) { //Previous MATCH
//
//				//Deletion
//				change_result(r_iterator, k, l, pos+1);
//				add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//
//				for (uint8_t b=0;b<nA;b++) {
//
//					BWiteration(k, l, k_aux, l_aux, b, C, C1, O);
//
//					if (k_aux > l_aux) continue;
//
//					//Insertion
//					change_result(r_iterator, k_aux, l_aux, pos);
//					modify_last_mismatch2(r_iterator, INSERTION, b);
//					add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//
//					if (b!=W[offset*MAXLINE+pos]) { //Mismatch
//						r_iterator->pos = pos+1;
//						modify_last_mismatch1(r_iterator, MISMATCH);
//						add_resultGPU(r_iterator, rl_next_results, rl_next_num_results, max_results);
//					}
//
//				}
//
//			}
//
//		}
//
//		rl_next.num_results[offset] = rl_next_num_results;
//	}
//
//}
//
//void BWSearchGPU(uintmax_t num_bloques, uintmax_t tam_bloques, uint8_t *W, char *h_W, SA_TYPE nW, vector *C, vector *h_C, vector *C1, vector *h_C1, comp_matrix *O, comp_matrix *h_O, comp_matrix *Oi, comp_matrix *h_Oi, comp_vector *S, comp_vector *R, comp_vector *Si, comp_vector *Ri, blocked_results_lists *rl_prev_cpu, blocked_results_lists *rl_next_cpu, blocked_results_lists *rl_prev_i_cpu, blocked_results_lists *rl_next_i_cpu, blocked_results_lists *rl_final_cpu, blocked_results_lists *rl_prev_gpu, blocked_results_lists *rl_next_gpu, blocked_results_lists *rl_prev_i_gpu, blocked_results_lists *rl_next_i_gpu, blocked_results_lists *rl_final_gpu, int16_t fragsize, uintmax_t max_results) {
//
//	result r;
//
//	int16_t fragments = nW / fragsize;
//	int16_t half = fragments / 2;
//	if (fragments % 2) half++;
//	int err_count;
//
//	timevars();
//
//	//printf("\n----> Tamaño: %ju, Fragmentos: %ju Errores: %ju\n", nW, fragments, fragments-1);
//
//	/* //////////////////////////////FORWARD/////////////////////////////////////////// */
//
//	/* for (int16_t i = half-1; i>0; i--) { */
//
//	/* err_count = fragments-1; */
//	/* init_result(&r, 1); */
//	/* change_result(&r, 0, O->siz-2, fragsize*i); */
//	/* bound_result(&r, fragsize*i, fragsize*(i+1) - 1); */
//	/* BWExactSearchForwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, nW, C->vector, C1->vector, *Oi, r, *rl_prev_gpu, max_results); */
//
//	/*   while (err_count > 0) { */
//
//	/*     BWExactPartialResultsForwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, C->vector, C1->vector, Oi, rl_prev, rl_next, rl_prev_i, fragsize, half-1); */
//	/*     BWChangeDirectionForwardBlockedGPU<<<num_bloques,tam_bloques>>>(Si, R, C, Oi, O, rl_prev_i, 0); */
//	/*     BWExactPartialResultsBackwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, C->vector, C1->vector, O, rl_prev_i, rl_next_i, rl_final, fragsize, half-1); */
//	/*     BWBranchPartialResultsForwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, C->vector, C1->vector, Oi, rl_next, rl_prev); */
//	/*     BWBranchPartialResultsBackwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, C->vector, C1->vector, O, rl_next_i, rl_prev_i); */
//
//	/*     err_count--; */
//
//	/*   } */
//
//	/*   BWExactFinalResultsForwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, C->vector, C1, Oi, rl_prev, rl_prev_i, fragsize, half-1); */
//	/*   BWChangeDirectionForwardBlockedGPU<<<num_bloques,tam_bloques>>>(Si, R, C, Oi, O, rl_prev_i, 0); */
//	/*   BWExactFinalResultsBackwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, C->vector, C1, O, rl_prev_i, rl_final, fragsize, half-1); */
//
//	/* } */
//
//	///////BLOCK 0/////////////////////////////////////
//	err_count = fragments-1;
//
//	init_listsGPU<<<num_bloques,tam_bloques>>>(*rl_prev_gpu, *rl_next_gpu, *rl_prev_i_gpu, *rl_next_i_gpu, *rl_final_gpu);
//
//	init_result(&r, 1);
//	change_result(&r, 0, O->siz-2, 0);
//	bound_result(&r, 0, fragsize - 1);
//
//	// cudaThreadSynchronize();
//	// tic("Initial");
//	BWExactSearchForwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, nW, C->vector, C1->vector, *Oi, r, *rl_prev_gpu, max_results);
//	// cudaThreadSynchronize();
//	// toc();
//
//	while (err_count > 0) {
//
//		// cudaThreadSynchronize();
//		// tic("Forward");
//		BWExactPartialResultsForwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, C->vector, C1->vector, *Oi, *rl_prev_gpu, *rl_next_gpu, *rl_final_gpu, fragsize, /*half-1*/0, max_results);
//		// cudaThreadSynchronize();
//		// toc();
//
//		// cudaThreadSynchronize();
//		// tic("To CPU");
//		// copy_blocked_results_list_cpu(rl_next_cpu, rl_next_gpu, max_results, num_bloques * tam_bloques);
//		// cudaThreadSynchronize();
//		// toc();
//		// cudaThreadSynchronize();
//		// tic("Branch");
//		// BWBranchPartialResultsForwardBlockedCPU(num_bloques, tam_bloques, h_W, h_C, h_C1, h_Oi, *rl_next_cpu, *rl_prev_cpu, max_results, nA);
//		// cudaThreadSynchronize();
//		// toc();
//		// cudaThreadSynchronize();
//		// tic("To GPU");
//		// copy_blocked_results_list_gpu(rl_prev_gpu, rl_prev_cpu, max_results, num_bloques * tam_bloques);
//		// cudaThreadSynchronize();
//		// toc();
//		if (err_count == 1) {
//
//			//  copy_blocked_results_list_cpu(rl_next_cpu, rl_next_gpu, max_results, num_bloques * tam_bloques);
//
//			//  cudaThreadSynchronize();
//			//  toc();
//			break;
//		}
//
//		// cudaThreadSynchronize();
//		// tic("Branch");
//		BWBranchPartialResultsForwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, C->vector, C1->vector, *Oi, *rl_next_gpu, *rl_prev_gpu, max_results, nA);
//		// cudaThreadSynchronize();
//		// toc();
//
//		err_count--;
//
//	}
//
//	cudaThreadSynchronize();
//	tic("FinalForward");
//	//BWBranchFinalResultsForwardBlockedCPU(num_bloques, tam_bloques, h_W, h_C, h_C1, h_Oi, *rl_next_cpu, *rl_final_cpu, fragsize, /*half-1*/0, max_results, nA);
//	BWBranchFinalResultsForwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, C->vector, C1->vector, *Oi, *rl_next_gpu, *rl_final_gpu, fragsize, /*half-1*/0, max_results, nA);
//	//BWExactFinalResultsForwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, C->vector, C1->vector, *Oi, *rl_prev_gpu, *rl_final_gpu, fragsize, /*half-1*/0, max_results);
//	cudaThreadSynchronize();
//	toc();
//
//	//copy_blocked_results_list_gpu(rl_final_gpu, rl_final_cpu, max_results, num_bloques * tam_bloques);
//
//
//	//////////////////////////////BACKWARD///////////////////////////////////////////
//
//	/* for (int16_t i = half; i<fragments-1; i++) { */
//
//	/* /\* printf("\n****BLOCK %d****\n", i); *\/ */
//	/* err_count = fragments-1; */
//
//	/* init_result(&r, 0); */
//	/* change_result(&r, 0, O->siz-2, fragsize*(i+1) - 1); */
//	/* bound_result(&r, fragsize*i, fragsize*(i+1) - 1); */
//	/* BWExactSearchBackwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, nW, C->vector, C1->vector, *O, r, *rl_prev_gpu, max_results); */
//
//	/*   while (err_count > 0) { */
//	/*     BWExactPartialResultsBackwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, C, C1, O, rl_prev, rl_next, rl_prev_i, fragsize, 0); */
//	/*     BWChangeDirectionBackwardBlockedGPU<<<num_bloques,tam_bloques>>>(S, Ri, C, O, Oi, rl_prev_i, nW-1); */
//	/*     BWExactPartialResultsForwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, C, C1, Oi, rl_prev_i, rl_next_i, rl_final, fragsize, 0); */
//	/*     BWBranchPartialResultsBackwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, C, C1, O, rl_next, rl_prev); */
//	/*     BWBranchPartialResultsForwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, C, C1, Oi, rl_next_i, rl_prev_i); */
//	/*     err_count--; */
//	/*   } */
//
//	/*   BWExactFinalResultsBackwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, C, C1, O, rl_prev, rl_prev_i, fragsize, 0); */
//	/*   BWChangeDirectionBackwardBlockedGPU<<<num_bloques,tam_bloques>>>(S, Ri, C, O, Oi, rl_prev_i, nW-1); */
//	/*   BWExactFinalResultsForwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, C, C1, Oi, rl_prev_i, rl_final, fragsize, 0); */
//
//	/* } */
//
//	/* ///////BLOCK FRAGMENTS-1///////////////////////////////////// */
//	/* printf("\n****BLOCK %d****\n", fragments-1); */
//
//	/* err_count = fragments-1; */
//
//	/* init_result(&r, 0); */
//	/* change_result(&r, 0, O->siz-2, /\*fragsize*fragments - 1 Last block is larger*\/nW-1); */
//	/* bound_result(&r, fragsize*(fragments-1), /\*fragsize*fragments - 1 Last block is larger*\/nW-1); */
//	/* BWExactSearchBackwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, nW, C->vector, C1->vector, *O, r, *rl_prev_gpu, max_results); */
//
//	/* while (err_count > 0) { */
//	/*   BWExactPartialResultsBackwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, C->vector, C1->vector, *O, *rl_prev_gpu, *rl_next_gpu, *rl_final_gpu, fragsize, 0, max_results); */
//	/*   BWBranchPartialResultsBackwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, C->vector, C1->vector, *O, *rl_next_gpu, *rl_prev_gpu, max_results, nA); */
//	/*   err_count--; */
//	/* } */
//
//	/* BWExactFinalResultsBackwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, C->vector, C1->vector, *O, *rl_prev_gpu, *rl_final_gpu, fragsize, 0, max_results); */
//
//}
//
///* void BWExactSearchBackwardBlockedGPUWrapper(uintmax_t num_bloques, uintmax_t tam_bloques, uint8_t *W, SA_TYPE nW, vector *C, vector *C1, comp_matrix *O, result *r, blocked_results_lists *rl_prev_gpu, uintmax_t max_results) { */
///*   BWExactSearchBackwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, nW, C->vector, C1->vector, *O, *r, *rl_prev_gpu, max_results); */
///* } */
//
///* void BWExactSearchForwardBlockedGPUWrapper(uintmax_t num_bloques, uintmax_t tam_bloques, uint8_t *W, SA_TYPE nW, vector *C, vector *C1, comp_matrix *O, result *r, blocked_results_lists *rl_prev_gpu, uintmax_t max_results) { */
///*   BWExactSearchForwardBlockedGPU<<<num_bloques,tam_bloques>>>(W, nW, C->vector, C1->vector, *O, *r, *rl_prev_gpu, max_results); */
///* } */
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BWExactSearchBackwardGPUWrapper(uintmax_t num_bloques, uintmax_t tam_bloques, uint8_t* W, uint64_t* nW, intmax_t* k, intmax_t* l, intmax_t k_ini, intmax_t l_ini, vector* C, vector* C1, comp_matrix* O) {
	BWExactSearchBackwardGPU<<<num_bloques,tam_bloques>>>(W, nW, k, l, k_ini, l_ini, C->vector, C1->vector, *O);
}

void BWExactSearchForwardGPUWrapper(uintmax_t num_bloques, uintmax_t tam_bloques, uint8_t* W, uint64_t* nW, intmax_t* k, intmax_t* l, intmax_t k_ini, intmax_t l_ini, vector* C, vector* C1, comp_matrix* O) {
	BWExactSearchForwardGPU<<<num_bloques,tam_bloques>>>(W, nW, k, l, k_ini, l_ini, C->vector, C1->vector, *O);
}

void BWExactSearchBackwardVectorGPUWrapper(uintmax_t num_bloques, uintmax_t tam_bloques, uint8_t* W, uint64_t* nW, intmax_t* k, intmax_t* l, intmax_t k_ini, intmax_t l_ini, vector* C, vector* C1, comp_matrix* O) {
	BWExactSearchBackwardVectorGPU<<<num_bloques,tam_bloques>>>(W, nW, k, l, k_ini, l_ini, C->vector, C1->vector, *O);
}

void BWExactSearchForwardVectorGPUWrapper(uintmax_t num_bloques, uintmax_t tam_bloques, uint8_t* W, uint64_t* nW, intmax_t* k, intmax_t* l, intmax_t k_ini, intmax_t l_ini, vector* C, vector* C1, comp_matrix* O) {
	BWExactSearchForwardVectorGPU<<<num_bloques,tam_bloques>>>(W, nW, k, l, k_ini, l_ini, C->vector, C1->vector, *O);
}
