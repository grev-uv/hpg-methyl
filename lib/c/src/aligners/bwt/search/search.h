#ifndef _SEARCH_SEARCH_
#define _SEARCH_SEARCH_

#include<inttypes.h>

#include "results.h"
#include "runtime.h"

static inline void BWExactSearchBackward(uint8_t *W, bwt_index *index, result *r, result *r_anchor) {

        intmax_t k = 0, l = 0, k2, l2;
	int16_t i;
	
	k2 = r->k;
	l2 = r->l;
	
	//printf("B1º -> %lu - %lu\n", k2, l2);

	for(i=r->pos; i>=r->start; i--) {
	  BWiteration(k2, l2, k2, l2, W[i], index);
	  //printf("B -> %d -> %lu - %lu\n", i, k2, l2);
	  if (k2 > l2) break;

	  k = k2; 
	  l = l2;
	}

	if (r_anchor != NULL) {
	  init_result(r_anchor, 0);
	  bound_result(r_anchor, i + 1, r->pos);
	  change_result(r_anchor, k, l, i);
	}

	change_result(r, k2, l2, i);

}

static inline void BWExactSearchForward(uint8_t *W, bwt_index *index, result *r, result *r_anchor) {

        intmax_t k = 0, l = 0, k2, l2;
	int16_t i;

	k2 = r->k;
	l2 = r->l;

	//printf("F1º-> %d -> %lu - %lu\n", i, k2, l2);
	for(i=r->pos; i<=r->end; i++) {

	  BWiteration(k2, l2, k2, l2, W[i], index);
	  //printf("F-> %d -> %lu - %lu\n", i, k2, l2);
	  if (k2 > l2) break;

	  k = k2; 
	  l = l2;
	}

	if (r_anchor != NULL) {
	  init_result(r_anchor, 1);
	  bound_result(r_anchor, r->pos, i - 1);
	  change_result(r_anchor, k, l, i);
	}

	change_result(r, k2, l2, i);

}

static inline bool BWExactFinalResultBackward(uint8_t *W, bwt_index *index, result *r_iterator, results_list *rl_final, int16_t block_size, int16_t last_block) {

	intmax_t k, l;
	int16_t start, pos;
	int16_t current_block;

	start = r_iterator->start;
	pos   = r_iterator->pos;

	k = r_iterator->k;
	l = r_iterator->l;

	current_block = pos / block_size;

	if ((current_block < last_block) || (pos == start-1)) { // Current block will be always >= start and previous results are propagated

	} else {

		if (current_block > last_block) { //Not in last previsited block

			return false;

		} else { //I am in the last previsited block

			if ((pos + 1) % block_size) { //I am not in the first element of the block
			} else { //I am in the first element in the block (all the block must be processed)
				return false;
			}

		}

	}

	for(int16_t i=pos; i>=start; i--) {
		BWiteration(k, l, k, l, W[i], index);
		if (k > l) break;
	}

	if (k <= l) {
		change_result(r_iterator, k, l, start-1);
		add_result(r_iterator, rl_final);
	}

	return false;

}

static inline bool BWExactFinalResultForward(uint8_t *W, bwt_index *index, result *r_iterator, results_list *rl_final, int16_t block_size, int16_t last_block) {

	intmax_t k, l;
	int16_t pos, end;
	int16_t current_block;

	pos   = r_iterator->pos;
	end   = r_iterator->end;

	k = r_iterator->k;
	l = r_iterator->l;

	current_block = pos / block_size;

	if ( (current_block > last_block) || (pos == end+1) ) { // Current block will be always <= end and previous results are propagated

	} else {

		if (current_block < last_block) { //Not in last previsited block

			return 0;

		} else { //I am in the last previsited block

			if (pos % block_size) { //I am not in the first element of the block
			} else { //I am in the first element in the block (all the block must be processed)
				return 0;
			}

		}

	}

	for(int16_t i=pos; i<=end; i++) {
		BWiteration(k, l, k, l, W[i], index);
		if (k > l) break;
	}

	if (k <= l) {
		change_result(r_iterator, k, l, end+1);
		add_result(r_iterator, rl_final);
	}

	return 0;

}

static inline void change_direction(bwt_index *backward, bwt_index *forward, result *res) {

	intmax_t k, l, ki, li, aux, aux2;
	int16_t start, end, err_offset;

	k  = res->k;
	l  = res->l;
	ki = size_SA(backward)-1;
	li = 0;

	start = res->start;
	end   = res->end;

	err_offset=0;

	for (uint8_t rr=0; rr<res->num_mismatches; rr++) {

		if      (res->err_kind[rr]==DELETION)
			err_offset--;
		else if (res->err_kind[rr]==INSERTION)
			err_offset++;

	}

	for (intmax_t i = k; i <= l; i++) {

		aux  = size_SA(backward) - get_SA(i, backward) - (end - start + 2) - err_offset;
		aux2 = get_ISA(aux, forward);

		if (aux2 < ki) ki = aux2;
		if (aux2 > li) li = aux2;

	}

	res->k = ki;
	res->l = li;

}

bool BWSearchCPU(uint8_t *W, uint64_t nW, bwt_index *backward, bwt_index *forward, results_list *rl_prev, results_list *rl_next, results_list *rl_prev_i, results_list *rl_next_i, results_list *rl_final, results_list *rl_anchor, int16_t fragsize, bool type, uint8_t nA);

int16_t BWExactSearchVectorBackward(uint8_t *W, int16_t start, int16_t end, intmax_t k, intmax_t l, intmax_t *vec_k, intmax_t *vec_l, bwt_index *index);
int16_t BWExactSearchVectorForward(uint8_t *W, int16_t start, int16_t end, intmax_t k, intmax_t l, intmax_t *vec_k, intmax_t *vec_l, bwt_index *index);

bool BWSearch1VectorHelper(uint8_t *W, int16_t start, int16_t end, intmax_t *vec_k, intmax_t *vec_l, intmax_t *vec_ki, intmax_t *vec_li, bwt_index *backward, bwt_index *forward, results_list *r_list, uint8_t nA);

bool BWSearch1CPU(uint8_t *W, bwt_index *backward, bwt_index *forward, result *res, results_list *r_list, results_list *r_list_anchor, uint8_t nA);

#endif
