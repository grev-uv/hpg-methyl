#include "search.h"

bool BWExactFinalResultsBackward(uint8_t *W, bwt_index *index, results_list *rl_prev, results_list *rl_next_i, int16_t block_size, int16_t last_block) {

	intmax_t k, l;
	int16_t start, pos;
	int16_t current_block;
	result *r_iterator;

	for (uintmax_t ii=0; ii < rl_prev->num_results; ii++) {

		r_iterator = &rl_prev->list[ii];

		start = r_iterator->start;
		pos   = r_iterator->pos;

		k = r_iterator->k;
		l = r_iterator->l;

		current_block = pos / block_size;

		if ((current_block < last_block) || (pos == start-1)) { // Current block will be always >= start and previous results are propagated
	
		} else {

			if (current_block > last_block) { //Not in last previsited block

				continue;

			} else { //I am in the last previsited block

				if ((pos + 1) % block_size) { //I am not in the first element of the block
				} else { //I am in the first element in the block (all the block must be processed)
					continue;
				}

			}

		}

		for(int16_t i=pos; i>=start; i--) {
			BWiteration(k, l, k, l, W[i], index);
			if (k > l) break;
		}

		if (k <= l) {
			change_result(r_iterator, k, l, start-1);
			add_result(r_iterator, rl_next_i);
		}

	} //r_prev

	rl_prev->num_results = 0;

	return false;

}

bool BWExactFinalResultsForward(uint8_t *W, bwt_index *index, results_list *rl_prev, results_list *rl_next_i, int16_t block_size, int16_t last_block) {

	intmax_t k, l;
	int16_t pos, end;
	int16_t current_block;

	result *r_iterator;

	for (uintmax_t ii=0; ii < rl_prev->num_results; ii++) {

		r_iterator = &rl_prev->list[ii];

		pos   = r_iterator->pos;
		end   = r_iterator->end;

		k = r_iterator->k;
		l = r_iterator->l;

		current_block = pos / block_size;

		if ( (current_block > last_block) || (pos == end+1) ) { // Current block will be always <= end and previous results are propagated

		} else {

			if (current_block < last_block) { //Not in last previsited block

				continue;

			} else { //I am in the last previsited block

				if (pos % block_size) { //I am not in the first element of the block
				} else { //I am in the first element in the block (all the block must be processed)
					continue;
				}

			}

		}

		for(int16_t i=pos; i<=end; i++) {
			BWiteration(k, l, k, l, W[i], index);
			if (k > l) break;
		}

		if (k <= l) {
			change_result(r_iterator, k, l, end+1);
			add_result(r_iterator, rl_next_i);
		}

	} //r_prev

	rl_prev->num_results = 0;

	return false;

}

bool BWBranchFinalResultsBackward(uint8_t *W, bwt_index *index, results_list *rl_prev, results_list *rl_final, int16_t block_size, int16_t last_block, uint8_t nA) {

	intmax_t k, l, k_aux, l_aux;
	int16_t start, pos, last_err_pos, r_num_mismatches;
	uint8_t last_err_kind;
	uint8_t last_err_base;
	bool no_previous;

	result *r_iterator;

	for (uintmax_t ii=0; ii < rl_prev->num_results; ii++) {

		r_iterator = &rl_prev->list[ii];

		if (r_iterator->dir != 0) continue;

		start = r_iterator->start;
		pos   = r_iterator->pos;

		if (pos < start) {
			if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
			continue;
		}

		no_previous = true;
		r_num_mismatches = r_iterator->num_mismatches-1;
		if (r_num_mismatches>-1) {
			last_err_pos  = r_iterator->err_pos[r_num_mismatches];
			last_err_kind = r_iterator->err_kind[r_num_mismatches];
			last_err_base = r_iterator->err_base[r_num_mismatches];
		} else {
			last_err_pos  = -10;
			last_err_kind = 0;
			last_err_base = (uint8_t) -1;
		}

		k = r_iterator->k;
		l = r_iterator->l;

		add_mismatch(r_iterator, DELETION, -1, pos);

		if (last_err_pos == pos + 1) { //Previous MISMATCH or DELETION

			if (last_err_kind == MISMATCH) { //Previous MISMATCH

				//Deletion
				if (W[pos]!=last_err_base) {
					change_result(r_iterator, k, l, pos-1);
					if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
				}

				for (uint8_t b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, index);

					if (k_aux > l_aux) continue;

					if (b!=W[pos]) {

						//Insertion
						if (b!=W[last_err_pos]) {
							change_result(r_iterator, k_aux, l_aux, pos);
							modify_last_mismatch2(r_iterator, INSERTION, b);
							if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
						}

						//Mismatch
						change_result(r_iterator, k_aux, l_aux, pos-1);
						modify_last_mismatch2(r_iterator, MISMATCH, b);
						if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

					}

				}

				no_previous = false;

			} else if (last_err_kind == DELETION) { //Previous DELETION

				//Deletion
				change_result(r_iterator, k, l, pos-1);
				if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

				for (uint8_t b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, index);

					if (k_aux > l_aux) continue;

					// NO INSERTION

					//Mismatch
					if (b!=W[pos]) {

						if (b!=W[last_err_pos]) {
							change_result(r_iterator, k_aux, l_aux, pos-1);
							modify_last_mismatch2(r_iterator, MISMATCH, b);
							if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
						}

					}

				}

				no_previous = false;

			}

		} else if (last_err_pos == pos) { //Previous INSERTION

			//NO DELETION

			for (uint8_t b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, index);

				if (k_aux > l_aux) continue;

				if (b!=W[pos]) {

					//Insertion
					change_result(r_iterator, k_aux, l_aux, pos);
					modify_last_mismatch2(r_iterator, INSERTION, b);
					if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

					//Mismatch
					if (W[pos]!=last_err_base) {
						r_iterator->pos = pos-1;
						modify_last_mismatch1(r_iterator, MISMATCH);
						if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
					}

				}

			}

			no_previous = false;

		}

		if (no_previous) { //Previous MATCH

			//Deletion
			change_result(r_iterator, k, l, pos-1);
			if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

			for (uint8_t b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, index);

				if (k_aux > l_aux) continue;

				if (b!=W[pos]) {

					//Insertion
					change_result(r_iterator, k_aux, l_aux, pos);
					modify_last_mismatch2(r_iterator, INSERTION, b);
					if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

					//Mismatch
					r_iterator->pos = pos-1;
					modify_last_mismatch1(r_iterator, MISMATCH);
					if (BWExactFinalResultBackward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

				}

			}

		}

	}

	rl_prev->num_results = 0;

	return false;

}

bool BWBranchFinalResultsForward(uint8_t *W, bwt_index *index, results_list *rl_prev, results_list *rl_final, int16_t block_size, int16_t last_block, uint8_t nA) {

	intmax_t k, l, k_aux, l_aux;
	int16_t end, pos, last_err_pos, r_num_mismatches;
	uint8_t last_err_kind;
	uint8_t last_err_base;
	bool no_previous;

	result *r_iterator;

	for (uintmax_t ii=0; ii < rl_prev->num_results; ii++) {

		r_iterator = &rl_prev->list[ii];

		if (r_iterator->dir != 1) continue;

		end = r_iterator->end;
		pos = r_iterator->pos;

		if (pos > end) {
			add_result(r_iterator, rl_final);
			continue;
		}

		no_previous = true;
		r_num_mismatches = r_iterator->num_mismatches-1;
		if (r_num_mismatches>-1) {
			last_err_pos  = r_iterator->err_pos[r_num_mismatches];
			last_err_kind = r_iterator->err_kind[r_num_mismatches];
			last_err_base = r_iterator->err_base[r_num_mismatches];
		} else {
			last_err_pos  = -10;
			last_err_kind = 0;
			last_err_base = (uint8_t) -1;
		}

		k = r_iterator->k;
		l = r_iterator->l;

		add_mismatch(r_iterator, DELETION, -1, pos);

		if (last_err_pos == pos - 1) { //Previous MISMATCH or DELETION

			if (last_err_kind == MISMATCH) { //Previous MISMATCH

				//Deletion
				if (W[pos]!=last_err_base) {
					change_result(r_iterator, k, l, pos+1);
					if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
				}

				for (uint8_t b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, index);

					if (k_aux > l_aux) continue;

					if (b!=W[pos]) {

						//Insertion
						if (b!=W[last_err_pos]) {
							change_result(r_iterator, k_aux, l_aux, pos);
							modify_last_mismatch2(r_iterator, INSERTION, b);
							if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
						}

						//Mismatch
						change_result(r_iterator, k_aux, l_aux, pos+1);
						modify_last_mismatch2(r_iterator, MISMATCH, b);
						if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

					}

				}

				no_previous = false;

			} else if (last_err_kind == DELETION) { //Previous DELETION

				//Deletion
				change_result(r_iterator, k, l, pos+1);
				if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

				for (uint8_t b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, index);

					if (k_aux > l_aux) continue;

					// NO INSERTION

					if (b!=W[pos]) { //Mismatch

						if (b!=W[last_err_pos]) {
							change_result(r_iterator, k_aux, l_aux, pos+1);
							modify_last_mismatch2(r_iterator, MISMATCH, b);
							if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
						}

					}

				}

				no_previous = false;

			}

		} else if (last_err_pos == pos) { //Previous INSERTION

			//NO DELETION

			for (uint8_t b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, index);

				if (k_aux > l_aux) continue;

				if (b!=W[pos]) {

					//Insertion
					change_result(r_iterator, k_aux, l_aux, pos);
					modify_last_mismatch2(r_iterator, INSERTION, b);
					if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

					//Mismatch
					if (W[pos]!=last_err_base) {
						r_iterator->pos = pos+1;
						modify_last_mismatch1(r_iterator, MISMATCH);
						if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;
					}

				}

			}

			no_previous = false;

		}

		if (no_previous) { //Previous MATCH

			//Deletion
			change_result(r_iterator, k, l, pos+1);
			if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

			for (uint8_t b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, index);

				if (k_aux > l_aux) continue;

				if (b!=W[pos]) {

					//Insertion
					change_result(r_iterator, k_aux, l_aux, pos);
					modify_last_mismatch2(r_iterator, INSERTION, b);
					if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

					//Mismatch
					r_iterator->pos = pos+1;
					modify_last_mismatch1(r_iterator, MISMATCH);
					if (BWExactFinalResultForward(W, index, r_iterator, rl_final, block_size, last_block)) return true;

				}

			}

		}

	}

	rl_prev->num_results = 0;

	return false;

}

bool BWExactPartialResultsBackward(uint8_t *W, uint8_t *D, uint8_t max_errors, bwt_index *index, results_list *rl_prev, results_list *rl_next, results_list *rl_next_i, int16_t block_size, int16_t last_block) {

	intmax_t k, l, k_next, l_next;
	int16_t start, pos, current_block, last_block_pos;
	uint8_t num_mismatches, error_bounding;
	bool complete_search;

	result *r_iterator;
	intmax_t results, results_next;

	for (uintmax_t ii=0; ii < rl_prev->num_results; ii++) {

		r_iterator = &rl_prev->list[ii];

		start  = r_iterator->start;
		pos    = r_iterator->pos;

		num_mismatches = r_iterator->num_mismatches;

		k_next = r_iterator->k;
		l_next = r_iterator->l;
		results_next = l_next - k_next;

		current_block = pos / block_size;

		error_bounding = max_errors - num_mismatches;

		if ((current_block < last_block) || (pos == start-1)) { // Current block will be always >= start and previous results are propagated

			last_block_pos = start;
			complete_search = true;

		} else {

			if (current_block > last_block) { //Not in last previsited block

				if ((pos + 1) % block_size) { //Not in first element of the block
					last_block_pos = (current_block-1) * block_size;
				} else { //I am in the first element in the block (all the block must be processed)
					last_block_pos = current_block * block_size;
				}

				complete_search = false;

			} else { //I am in the last previsited block

				if ((pos + 1) % block_size) { //I am not in the first element of the block
					last_block_pos = start;
					complete_search = true;
				} else { //I am in the first element of the block (all the block must be processed)
					last_block_pos = current_block * block_size;
					complete_search = false;
				}

			}

		}

		for(int16_t i=pos; i>=last_block_pos; i--) {

			k = k_next;
			l = l_next;

			if (k > l) break;

			BWiteration(k, l, k_next, l_next, W[i], index);
			results      = results_next;
			results_next = l_next - k_next;

			if (results == results_next || error_bounding < D[i]) continue;

			change_result(r_iterator, k, l, i);
			add_result(r_iterator, rl_next);

		}

		if (complete_search && k_next <= l_next) {
			change_result(r_iterator, k_next, l_next, start-1);
			add_result(r_iterator, rl_next_i);
		}

	} //r_prev

	rl_prev->num_results = 0;

	return false;

}

bool BWExactPartialResultsForward(uint8_t *W, uint8_t *D, uint8_t max_errors, bwt_index *index, results_list *rl_prev, results_list *rl_next, results_list *rl_next_i, int16_t block_size, int16_t last_block) {

	intmax_t k, l, k_next, l_next;
	int16_t pos, end, current_block, last_block_pos;
	uint8_t num_mismatches, error_bounding;
	bool complete_search;

	result *r_iterator;
	intmax_t results, results_next;

	for (uintmax_t ii=0; ii < rl_prev->num_results; ii++) {

		r_iterator = &rl_prev->list[ii];

		pos   = r_iterator->pos;
		end   = r_iterator->end;

		num_mismatches = r_iterator->num_mismatches;

		k_next = r_iterator->k;
		l_next = r_iterator->l;
		results_next = l_next - k_next;

		current_block = pos / block_size;

		error_bounding = max_errors - num_mismatches;

		if ( (current_block > last_block) || (pos == end+1) ) { // Current block will be always <= end and previous results are propagated
			last_block_pos = end;
			complete_search = true;

		} else {

			if (current_block < last_block) { //Not in last previsited block

				if (pos % block_size) { //Not in first element of the block
					last_block_pos = (current_block+2) * block_size - 1;
				} else { //I am in the first element in the block (all the block must be processed)
					last_block_pos = (current_block+1) * block_size - 1;
				}

				complete_search = false;

			} else { //I am in the last previsited block

				if (pos % block_size) { //I am not in the first element of the block
					last_block_pos = end;
					complete_search = true;
				} else { //I am in the first element in the block (all the block must be processed)
					last_block_pos = (current_block+1) * block_size - 1;
					complete_search = false;
				}

			}

		}

		for(int16_t i=pos; i<=last_block_pos; i++) {

			k = k_next;
			l = l_next;

			if (k > l) break;

			BWiteration(k, l, k_next, l_next, W[i], index);
			results      = results_next;
			results_next = l_next - k_next;

			if (results == results_next || error_bounding < D[i]) continue;

			change_result(r_iterator, k, l, i);
			add_result(r_iterator, rl_next);

		}

		if (complete_search && k_next <= l_next) {
			change_result(r_iterator, k_next, l_next, end+1);
			add_result(r_iterator, rl_next_i);
		}

	} //r_prev

	rl_prev->num_results = 0;

	return false;

}

bool BWBranchPartialResultsBackward(uint8_t *W, bwt_index *index, results_list *rl_prev, results_list *rl_next, uint8_t nA) {

	intmax_t k, l, k_aux, l_aux;
	int16_t start, pos, last_err_pos, r_num_mismatches;
	uint8_t last_err_kind;
	uint8_t last_err_base;
	bool no_previous;

	result *r_iterator;

	for (uintmax_t ii=0; ii < rl_prev->num_results; ii++) {

		r_iterator = &rl_prev->list[ii];

		if (r_iterator->dir != 0) continue;

		start = r_iterator->start;
		pos   = r_iterator->pos;

		if (pos < start) {
			add_result(r_iterator, rl_next);
			continue;
		}

		no_previous = true;
		r_num_mismatches = r_iterator->num_mismatches-1;
		if (r_num_mismatches>-1) {
			last_err_pos  = r_iterator->err_pos[r_num_mismatches];
			last_err_kind = r_iterator->err_kind[r_num_mismatches];
			last_err_base = r_iterator->err_base[r_num_mismatches];
		} else {
			last_err_pos  = -10;
			last_err_kind = 0;
			last_err_base = (uint8_t) -1;
		}

		k = r_iterator->k;
		l = r_iterator->l;

		add_mismatch(r_iterator, DELETION, -1, pos);

		if (last_err_pos == pos + 1) { //Previous MISMATCH or DELETION

			if (last_err_kind == MISMATCH) { //Previous MISMATCH

				//Deletion
				if (W[pos]!=last_err_base) {
					change_result(r_iterator, k, l, pos-1);
					add_result(r_iterator, rl_next);
				}

				for (uint8_t b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, index);

					if (k_aux > l_aux) continue;

					if (b!=W[pos]) {

						//Insertion
						if (b!=W[last_err_pos]) {
							change_result(r_iterator, k_aux, l_aux, pos);
							modify_last_mismatch2(r_iterator, INSERTION, b);
							add_result(r_iterator, rl_next);
						}

						//Mismatch
						change_result(r_iterator, k_aux, l_aux, pos-1);
						modify_last_mismatch2(r_iterator, MISMATCH, b);
						add_result(r_iterator, rl_next);

					}

				}

				no_previous = false;

			} else if (last_err_kind == DELETION) { //Previous DELETION

				//Deletion
				change_result(r_iterator, k, l, pos-1);
				add_result(r_iterator, rl_next);

				for (uint8_t b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, index);

					if (k_aux > l_aux) continue;

					// NO INSERTION

					if (b!=W[pos]) { //Mismatch

						if (b!=W[last_err_pos]) {
							change_result(r_iterator, k_aux, l_aux, pos-1);
							modify_last_mismatch2(r_iterator, MISMATCH, b);
							add_result(r_iterator, rl_next);
						}

					}

				}

				no_previous = false;

			}

		} else if (last_err_pos == pos) { //Previous INSERTION

			//NO DELETION

			for (uint8_t b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, index);

				if (k_aux > l_aux) continue;

				if (b!=W[pos]) {
					//Insertion
					change_result(r_iterator, k_aux, l_aux, pos);
					modify_last_mismatch2(r_iterator, INSERTION, b);
					add_result(r_iterator, rl_next);

					//Mismatch
					if (W[pos]!=last_err_base) {
						r_iterator->pos = pos-1;
						modify_last_mismatch1(r_iterator, MISMATCH);
						add_result(r_iterator, rl_next);
					}

				}

			}

			no_previous = false;

		}

		if (no_previous) { //Previous MATCH

			//Deletion
			change_result(r_iterator, k, l, pos-1);
			add_result(r_iterator, rl_next);

			for (uint8_t b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, index);

				if (k_aux > l_aux) continue;

				if (b!=W[pos]) {

					//Insertion
					change_result(r_iterator, k_aux, l_aux, pos);
					modify_last_mismatch2(r_iterator, INSERTION, b);
					add_result(r_iterator, rl_next);

					//Mismatch
					r_iterator->pos = pos-1;
					modify_last_mismatch1(r_iterator, MISMATCH);
					add_result(r_iterator, rl_next);

				}

			}

		}

	}

	rl_prev->num_results = 0;

	return false;

}

bool BWBranchPartialResultsForward(uint8_t *W, bwt_index *index, results_list *rl_prev, results_list *rl_next, uint8_t nA) {

	intmax_t k, l, k_aux, l_aux;
	int16_t end, pos=0, last_err_pos, r_num_mismatches;
	uint8_t last_err_kind;
	uint8_t last_err_base;
	bool no_previous;

	result *r_iterator;

	for (uintmax_t ii=0; ii < rl_prev->num_results; ii++) {

		r_iterator = &rl_prev->list[ii];

		if (r_iterator->dir != 1) continue;

		end = r_iterator->end;
		pos = r_iterator->pos;

		if (pos > end) {
			add_result(r_iterator, rl_next);
			continue;
		}

		no_previous = true;
		r_num_mismatches = r_iterator->num_mismatches-1;
		if (r_num_mismatches>-1) {
			last_err_pos  = r_iterator->err_pos[r_num_mismatches];
			last_err_kind = r_iterator->err_kind[r_num_mismatches];
			last_err_base = r_iterator->err_base[r_num_mismatches];
		} else {
			last_err_pos  = -10;
			last_err_kind = 0;
			last_err_base = (uint8_t) -1;
		}

		k = r_iterator->k;
		l = r_iterator->l;

		add_mismatch(r_iterator, DELETION, -1, pos);

		if (last_err_pos == pos - 1) { //Previous MISMATCH or DELETION


			if (last_err_kind == MISMATCH) { //Previous MISMATCH

				//Deletion
				if (W[pos]!=last_err_base) {
					change_result(r_iterator, k, l, pos+1);
					add_result(r_iterator, rl_next);
				}

				for (uint8_t b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, index);

					if (k_aux > l_aux) continue;

					if (b!=W[pos]) {

						//Insertion
						if (b!=W[last_err_pos]) {
							change_result(r_iterator, k_aux, l_aux, pos);
							modify_last_mismatch2(r_iterator, INSERTION, b);
							add_result(r_iterator, rl_next);
						}

						//Mismatch
						change_result(r_iterator, k_aux, l_aux, pos+1);
						modify_last_mismatch2(r_iterator, MISMATCH, b);
						add_result(r_iterator, rl_next);

					}

				}

				no_previous = false;

			} else if (last_err_kind == DELETION) { //Previous DELETION

				//Deletion
				change_result(r_iterator, k, l, pos+1);
				add_result(r_iterator, rl_next);

				for (uint8_t b=0;b<nA;b++) {

					BWiteration(k, l, k_aux, l_aux, b, index);

					if (k_aux > l_aux) continue;

					// NO INSERTION

					//Mismatch
					if (b!=W[pos]) {

						if (b!=W[last_err_pos]) {
							change_result(r_iterator, k_aux, l_aux, pos+1);
							modify_last_mismatch2(r_iterator, MISMATCH, b);
							add_result(r_iterator, rl_next);
						}

					}

				}

				no_previous = false;

			}

		} else if (last_err_pos == pos) { //Previous INSERTION

			//NO DELETION

			for (uint8_t b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, index);

				if (k_aux > l_aux) continue;

				if (b!=W[pos]) {

					//Insertion
					change_result(r_iterator, k_aux, l_aux, pos);
					modify_last_mismatch2(r_iterator, INSERTION, b);
					add_result(r_iterator, rl_next);

					//Mismatch
					if (W[pos]!=last_err_base) {
						r_iterator->pos = pos+1;
						modify_last_mismatch1(r_iterator, MISMATCH);
						add_result(r_iterator, rl_next);
					}

				}

			}

			no_previous = false;

		}

		if (no_previous) { //Previous MATCH

			//Deletion
			change_result(r_iterator, k, l, pos+1);
			add_result(r_iterator, rl_next);

			for (uint8_t b=0;b<nA;b++) {

				BWiteration(k, l, k_aux, l_aux, b, index);

				if (k_aux > l_aux) continue;

				if (b!=W[pos]) { //Mismatch

					//Insertion
					change_result(r_iterator, k_aux, l_aux, pos);
					modify_last_mismatch2(r_iterator, INSERTION, b);
					add_result(r_iterator, rl_next);

					r_iterator->pos = pos+1;
					modify_last_mismatch1(r_iterator, MISMATCH);
					add_result(r_iterator, rl_next);

				}

			}

		}

	}

	rl_prev->num_results = 0;

	return false;

}

void BWChangeDirectionBackward(bwt_index *backward, bwt_index *forward, results_list *r_list, int16_t end) {

	result *r_iterator;

	for (uintmax_t ii=r_list->num_results; ii > 0; ii--) {
		r_iterator = &r_list->list[ii-1];

		if (r_iterator->dir == 1) break;

		if (r_iterator->pos == (r_iterator->start-1)) {
			change_direction(backward, forward, r_iterator);
			r_iterator->pos = r_iterator->end+1;
			r_iterator->end = end;
			r_iterator->dir = 1;
		}

	}

}

void BWChangeDirectionForward(bwt_index *backward, bwt_index *forward, results_list *r_list, int16_t start) {

	result *r_iterator;

	for (uintmax_t ii=r_list->num_results; ii > 0; ii--) {
		r_iterator = &r_list->list[ii-1];

		if (r_iterator->dir == 0) break;

		if (r_iterator->pos == (r_iterator->end+1)) {

			change_direction(backward, forward, r_iterator);

			r_iterator->pos = r_iterator->start-1;
			r_iterator->start = start;
			r_iterator->dir = 0;

		}

	}

}

void calculateDBackward(uint8_t *D, uint8_t *W, uint64_t nW, bwt_index *backward, bwt_index *forward, result r) {

	intmax_t k, l;

	uint8_t z;
	int8_t last;

	z = 0;

	k = r.k;
	l = r.l;

	for (int16_t i=r.pos; i>= 0; i--) {

		BWiteration(k, l, k, l, W[i], backward);

		if (k>l) {
      k=0;
      l=size_SA(backward)-1;
      z++;
    }

		D[i] = z;

  }

	if (r.end == (int16_t) nW-1) {

		last = D[0];

	} else {

    k=0;
    l=size_SA(backward)-1;
 
		for (uint64_t i=r.end+1; i < nW; i++) {

			BWiteration(k, l, k, l, W[i], forward);

			if (k>l) {
				k=0;
				l=size_SA(forward)-1;
				z++;
			}

			D[i] = z;

		}

		last = D[nW-1];

		for (uint64_t i=r.end+1; i < nW; i++) {
			D[i] = last - D[i];
		}

	}

	for (int16_t i=r.pos; i>= 0; i--) {
		D[i] = last - D[i];
	}

}

void calculateDForward(uint8_t *D, uint8_t *W, uint64_t nW, bwt_index *backward, bwt_index *forward, result r) {

	intmax_t k,l;

	uint8_t z;
	int8_t last;

	z = 0;

	k = r.k;
	l = r.l;

	for (uint64_t i=r.pos; i < nW; i++) {

		BWiteration(k, l, k, l, W[i], forward);

		if (k>l) {
      k=0;
      l=size_SA(forward)-1;
      z++;
    }

		D[i] = z;

  }

	if (r.start == 0) {

		last = D[nW-1];

	} else {

    k=0;
    l=size_SA(backward)-1;
 
		for (int16_t i=r.start-1; i >= 0; i--) {

			BWiteration(k, l, k, l, W[i], backward);

			if (k>l) {
				k=0;
				l=size_SA(backward)-1;
				z++;
			}

			D[i] = z;

		}

		last = D[0];

		for (int16_t i=r.start-1; i >= 0; i--) {
			D[i] = last - D[i];
		}

	}

	for (uint64_t i=r.pos; i< nW; i++) {
		D[i] = last - D[i];
	}

}

bool BWSearchCPU(uint8_t *W, uint64_t nW, bwt_index *backward, bwt_index *forward, results_list *rl_prev, results_list *rl_next, results_list *rl_prev_i, results_list *rl_next_i, results_list *rl_final, results_list *rl_anchor, int16_t fragsize, bool type, uint8_t nA) {

        result r, r_anchor;

	int16_t fragments = nW / fragsize;
	int16_t half = fragments / 2;
	if (fragments % 2) half++;
	int err_count;

	bool flow;

	uint8_t *D = (uint8_t *) malloc (nW * sizeof(uint8_t));
	
	if (rl_anchor != NULL) rl_anchor->num_results = 0;

	//printf("\n***** Size: %d, Fragments: %d Errors: %d\n", nW, fragments, fragments-1);

	if (fragments <= 1) {

		if (type == 0) {
			init_result(&r, 1);
			change_result(&r, 0, size_SA(forward)-1, 0);
			bound_result(&r, 0, nW-1);
			if (rl_anchor != NULL) {
			  BWExactSearchForward(W, forward, &r, &r_anchor);
			  add_result(&r_anchor, rl_anchor);
			} else {
			  BWExactSearchForward(W, forward, &r, NULL);
			}
			if (r.k<=r.l)
				add_result(&r, rl_final);
		} else {
			init_result(&r, 0);
			change_result(&r, 0, size_SA(backward)-1, nW-1);
			bound_result(&r, 0, nW-1);
			if (rl_anchor != NULL) {
			  BWExactSearchBackward(W, backward, &r, &r_anchor);
			  add_result(&r_anchor, rl_anchor);
			} else {
			  BWExactSearchBackward(W, backward, &r, NULL);
			}
			if (r.k<=r.l)
			  add_result(&r, rl_final);
		}

		return false;

	}

	//////////////////////////////FORWARD///////////////////////////////////////////

	for (int16_t i = half-1; i > 0; i--) {
		//printf("\n****BLOCK %d****\n", i);

		flow = true;

		err_count = fragments-1;

		rl_prev->num_results = 0; rl_prev_i->num_results = 0;
		rl_next->num_results = 0; rl_next_i->num_results = 0;

		init_result(&r, 1);
		change_result(&r, 0, size_SA(forward)-1, fragsize*i);
		bound_result(&r, fragsize*i, fragsize*(i+1) - 1);
		if (rl_anchor != NULL) {
		  BWExactSearchForward(W, forward, &r, &r_anchor);
		  add_result(&r_anchor, rl_anchor);
		} else {
		  BWExactSearchForward(W, forward, &r, NULL);
		}
		r.end = nW-1;

		if (r.k <= r.l) {

			calculateDForward(D, W, nW, backward, forward, r);

			add_result(&r, rl_prev);

			while (err_count > 0) {

				if (BWExactPartialResultsForward(W, D, fragments-1, forward, rl_prev, rl_next, rl_prev_i, fragsize, half-1)) {flow = false; break;}
				BWChangeDirectionForward(forward, backward, rl_prev_i, 0);
				if (BWExactPartialResultsBackward(W, D, fragments-1, backward, rl_prev_i, rl_next_i, rl_final, fragsize, half-1)) {flow = false; break;}
				if (err_count==1) break;
				if (BWBranchPartialResultsForward(W, forward, rl_next, rl_prev, nA)) {flow = false; break;}
				if (BWBranchPartialResultsBackward(W, backward, rl_next_i, rl_prev_i, nA)) {flow = false; break;}
				err_count--;

			}

			if (flow) {
				BWBranchFinalResultsForward(W, forward, rl_next, rl_prev_i, fragsize, half-1, nA);
				BWChangeDirectionForward(forward, backward, rl_prev_i, 0);
				BWExactFinalResultsBackward(W, backward, rl_prev_i, rl_final, fragsize, half-1);

				BWBranchFinalResultsBackward(W, backward, rl_next_i, rl_final, fragsize, half-1, nA);
			}

		}

	}

	///////BLOCK 0/////////////////////////////////////
	//printf("\n****BLOCK %d****\n", 0);

	flow = true;

	err_count = fragments-1;

	rl_prev->num_results = 0; rl_prev_i->num_results = 0;
	rl_next->num_results = 0; rl_next_i->num_results = 0;

	init_result(&r, 1);
	change_result(&r, 0, size_SA(forward)-1, 0);
	bound_result(&r, 0, fragsize - 1);
	if (rl_anchor!=NULL) {
	  BWExactSearchForward(W, forward, &r, &r_anchor);
	add_result(&r_anchor, rl_anchor);
	} else {
	  BWExactSearchForward(W, forward, &r, NULL);
	}
	r.end = nW-1;

	if (r.k <= r.l) {

		calculateDForward(D, W, nW, backward, forward, r);

		add_result(&r, rl_prev);

		while (err_count > 0) {
			if (BWExactPartialResultsForward(W, D, fragments-1, forward, rl_prev, rl_next, rl_final, fragsize, half-1)) {flow = false; break;}
			if (err_count==1) break;
			if (BWBranchPartialResultsForward(W, forward, rl_next, rl_prev, nA)) {flow = false; break;}
			err_count--;
		}

		if (flow) {
			BWBranchFinalResultsForward(W, forward, rl_next, rl_final, fragsize, half-1, nA);
		}

	}

	//////////////////////////////BACKWARD///////////////////////////////////////////

	for (int16_t i = half; i<fragments-1; i++) {
		//printf("\n****BLOCK %d****\n", i);
		flow = true;

		err_count = fragments-1;

		rl_prev->num_results = 0; rl_prev_i->num_results = 0;
		rl_next->num_results = 0; rl_next_i->num_results = 0;

		init_result(&r, 0);
		change_result(&r, 0, size_SA(backward)-1, fragsize*(i+1) - 1);
		bound_result(&r, fragsize*i, fragsize*(i+1) - 1);
		if (rl_anchor != NULL) {
		  BWExactSearchBackward(W, backward, &r, &r_anchor);
		  add_result(&r_anchor, rl_anchor);
		} else {
		  BWExactSearchBackward(W, backward, &r, NULL);
		}
		r.start = 0;

		if (r.k <= r.l) {

			calculateDBackward(D, W, nW, backward, forward, r);

			add_result(&r, rl_prev);

			while (err_count > 0) {
				if (BWExactPartialResultsBackward(W, D, fragments-1, backward, rl_prev, rl_next, rl_prev_i, fragsize, 0)) {flow = false; break;}
				BWChangeDirectionBackward(backward, forward, rl_prev_i, nW-1);
				if (BWExactPartialResultsForward(W, D, fragments-1, forward, rl_prev_i, rl_next_i, rl_final, fragsize, 0)) {flow = false; break;}
				if (err_count==1) break;
				if (BWBranchPartialResultsBackward(W, backward, rl_next, rl_prev, nA)) {flow = false; break;}
				if (BWBranchPartialResultsForward(W, forward, rl_next_i, rl_prev_i, nA)) {flow = false; break;}
				err_count--;
			}

			if (flow) {
				BWBranchFinalResultsBackward(W, backward, rl_next, rl_prev_i, fragsize, 0, nA);
				BWChangeDirectionBackward(backward, forward, rl_prev_i, nW-1);
				BWExactFinalResultsForward(W, forward, rl_prev_i, rl_final, fragsize, 0);

				BWBranchFinalResultsForward(W, forward, rl_next_i, rl_final, fragsize, 0, nA);
			}

		}

	}

	///////BLOCK FRAGMENTS-1/////////////////////////////////////
	//printf("\n****BLOCK %d****\n", fragments-1);

	flow = true;

	err_count = fragments-1;

	rl_prev->num_results = 0; rl_prev_i->num_results = 0;
	rl_next->num_results = 0; rl_next_i->num_results = 0;

	init_result(&r, 0);
	change_result(&r, 0, size_SA(backward)-1, /*fragsize*fragments - 1 Last block is larger*/nW-1);
	bound_result(&r, fragsize*(fragments-1), /*fragsize*fragments - 1 Last block is larger*/nW-1);
	if (rl_anchor != NULL) {
	  BWExactSearchBackward(W, backward, &r, &r_anchor);
	  add_result(&r_anchor, rl_anchor);
	} else {
	  BWExactSearchBackward(W, backward, &r, NULL);
	}
	r.start = 0;

	if (r.k <= r.l) {

		calculateDBackward(D, W, nW, backward, forward, r);

		add_result(&r, rl_prev);

		while (err_count > 0) {
			if (BWExactPartialResultsBackward(W, D, fragments-1, backward, rl_prev, rl_next, rl_final, fragsize, 0)) {flow = false; break;}
			if(err_count==1) break;
			if (BWBranchPartialResultsBackward(W, backward, rl_next, rl_prev, nA)) {flow = false; break;}
			err_count--;
		}

		if (flow) {
			BWBranchFinalResultsBackward(W, backward, rl_next, rl_final, fragsize, 0, nA);
		}

	}

	return false;

}

int16_t BWExactSearchVectorBackward(uint8_t *W, int16_t start, int16_t end, intmax_t k, intmax_t l, intmax_t *vec_k, intmax_t *vec_l, bwt_index *index) {

  int16_t last_good = 0;

  if (k > l)       return 0;
  if (start > end) return 0;

  intmax_t k2, l2;
  int16_t last, i, j;

  last = end-start;

  k2 = k;
  l2 = l;

  //printf("B -> %d -> %lu - %lu\n", i, k2, l2);

  for(i=end, j=last; i>=start; i--, j--) {

    BWiteration(k2, l2, k2, l2, W[i], index);
    //printf("B -> %d -> %lu - %lu\n", i, k2, l2);

    vec_k[j] = k2;
    vec_l[j] = l2;

    if (k2 > l2) {
      i--; j--;
      break;
    }

    last_good++;

  }

  for(;i>=start; i--, j--) {
    vec_k[j] = k2;
    vec_l[j] = l2;
  }

  return last_good;

}

int16_t BWExactSearchVectorForward(uint8_t *W, int16_t start, int16_t end, intmax_t k, intmax_t l, intmax_t *vec_k, intmax_t *vec_l, bwt_index *index) {

  int16_t last_good = 0;

  if (k > l)       return 0;
  if (start > end) return 0;

  intmax_t k2, l2;
  int16_t i, j;

  k2 = k;
  l2 = l;

  //printf("F -> %d -> %lu - %lu\n", i, k2, l2);

  for(i=start, j=0; i<=end; i++, j++) {

    BWiteration(k2, l2, k2, l2, W[i], index);
    //printf("F -> %d -> %lu - %lu\n", i, k2, l2);

    vec_k[j] = k2;
    vec_l[j] = l2;

    if (k2 > l2) {
      i++; j++;
      break;
    }

    last_good++;

  }

  for(; i<=end; i++, j++) {
    vec_k[j] = k2;
    vec_l[j] = l2;
  }

  return last_good;

}

bool BWSearch1VectorHelper(uint8_t *W, int16_t start, int16_t end, intmax_t *vec_k, intmax_t *vec_l, intmax_t *vec_ki, intmax_t *vec_li, bwt_index *backward, bwt_index *forward, results_list *r_list, uint8_t nA) {

	intmax_t _k, _l, _ki, _li, _k_aux, _l_aux, _ki_aux, _li_aux;
	intmax_t results, results_last;

	int16_t i, j, half, n;

	result r;

	n = end - start;
	half = n / 2;

	init_result(&r, 0);

	bound_result(&r, start, end);

	if (vec_k[0] <= vec_l[0]) {
		change_result(&r, vec_k[0], vec_l[0], -1);
		add_result(&r, r_list);
	}

	add_mismatch(&r, MATCH, -1, start);

	results = vec_l[0] - vec_k[0];

	results_last = results;
	_k  = vec_k[1];
	_l  = vec_l[1];
	results = _l  - _k;

	//printf("B -> %d: %d -> %d, %u, %u\n", 0, results, results_last, _k, _l);

	if (results != results_last) {

		//printf("*B -> %d: %d -> %d, %u, %u\n", 0, results, results_last, _k, _l);

		//printf("%d: %d -> %d\n", 0, results, results_last);

		//Deletion
		change_result(&r, _k, _l, -1);
		modify_last_mismatch3(&r, DELETION, -1, start);
		add_result(&r, r_list);

		for (uint8_t b=0;b<nA;b++) {

			BWiteration(_k, _l, _k_aux, _l_aux, b, backward);
			//printf("W -> %d, %d - %d\n", b, _k_aux, _l_aux);

			if (_k_aux > _l_aux) continue;
			//printf("*W -> %d, %d - %d\n", b, _k_aux, _l_aux);

			uint8_t b_w = W[start];

			//Missmatch
			if (b!=b_w) {
				change_result(&r, _k_aux, _l_aux, -1);
				modify_last_mismatch2(&r, MISMATCH, b);
				add_result(&r, r_list);
			}

			//Insertion
			BWiteration(_k_aux, _l_aux, _k_aux, _l_aux, b_w, backward);

			if (_k_aux <= _l_aux) {
				change_result(&r, _k_aux, _l_aux, -1);
				modify_last_mismatch3(&r, 3, INSERTION, b);
				add_result(&r, r_list);
			}

		}

	}

	for (i=start+2, j=2; j<=half; i++, j++) {

		results_last = results;
		_k = vec_k[j];
		_l = vec_l[j];
		results = _l  - _k;

		//printf("B -> %d: %d -> %d, %u, %u\n", j-1, results, results_last, _k, _l);

		if (results == results_last) continue;

		//printf("*B -> %d: %d -> %d, %u, %u\n", j-1, results, results_last, _k, _l);

		//Deletion
		change_result(&r, _k, _l, i-2);
		modify_last_mismatch3(&r, DELETION, -1, i-1);
		BWExactSearchBackward(W, backward, &r, NULL);
		if (r.k<=r.l) add_result(&r, r_list);

		for (uint8_t b=0;b<nA;b++) {

			BWiteration(_k, _l, _k_aux, _l_aux, b, backward);

			if (_k_aux > _l_aux) continue;

			//Insertion
			change_result(&r, _k_aux, _l_aux, i-1);
			modify_last_mismatch2(&r, INSERTION, b);
			BWExactSearchBackward(W, backward, &r, NULL);
			if (r.k<=r.l) add_result(&r, r_list);

			//Mismatch
			if (b!=W[i-1]) {
				change_result(&r, _k_aux, _l_aux, i-2);
				modify_last_mismatch1(&r, MISMATCH);
				BWExactSearchBackward(W, backward, &r, NULL);
				if (r.k<=r.l) add_result(&r, r_list);
			}

		}

	}

	//printf("\n");

	half--;
	results = vec_li[n] - vec_ki[n];

	r.dir=1; //Change direction

	results_last = results;
	_ki  = vec_ki[n-1];
	_li  = vec_li[n-1];
	results = _li - _ki;

	//printf("F-> %d: %d -> %d, %u, %u\n", n, results, results_last, _ki, _li);

	if (results != results_last) {

		//printf("*F -> %d: %d -> %d, %u - %u\n", i+1, results, results_last, _ki, _li);

		//Deletion
		change_result(&r, _ki, _li, -1);
		modify_last_mismatch3(&r, DELETION, -1, end);
		add_result(&r, r_list);

		for (uint8_t b=0;b<nA;b++) {

			BWiteration(_ki, _li, _ki_aux, _li_aux, b, forward);

			if (_ki_aux > _li_aux) continue;

			uint8_t b_w = W[end];

			//Mismatch
			if (b!=b_w) {
				change_result(&r, _ki_aux, _li_aux, -1);
				modify_last_mismatch2(&r, MISMATCH, b);
				add_result(&r, r_list);
			}

			//Insertion
			BWiteration(_ki_aux, _li_aux, _ki_aux, _li_aux, b_w, forward);

			//printf("\tI -> %d - %d\n", _ki_aux, _li_aux);

			if (_ki_aux <= _li_aux){
				change_result(&r, _ki_aux, _li_aux, -1);
				modify_last_mismatch2(&r, INSERTION, b);
				add_result(&r, r_list);
			}

		}

	}

	for(i=end-2,j=n-2; j>=half; i--, j--) {

		results_last = results;
		_ki  = vec_ki[j];
		_li  = vec_li[j];
		results = _li - _ki;

		//printf("F -> %d: %d -> %d, %u - %u\n", i+1, results, results_last, _ki, _li);

		if (results == results_last) continue;

		//printf("*F -> %d: %d -> %d, %u - %u\n", i+1, results, results_last, _ki, _li);

		//Deletion
		change_result(&r, _ki, _li, i+2);
		modify_last_mismatch3(&r, DELETION, -1, i+1);
		BWExactSearchForward(W, forward, &r, NULL);
		if (r.k<=r.l) add_result(&r, r_list);

		for (uint8_t b=0;b<nA;b++) {

			BWiteration(_ki, _li, _ki_aux, _li_aux, b, forward);

			//printf("W -> %d, %d - %d\n", b, _ki_aux, _li_aux);

			if (_ki_aux > _li_aux) continue;

			//Insertion
			change_result(&r, _ki_aux, _li_aux, i+1);
			modify_last_mismatch2(&r, INSERTION, b);
			BWExactSearchForward(W, forward, &r, NULL);
			if (r.k<=r.l) add_result(&r, r_list);

			//Mismatch
			if (b!= W[i+1]) {
				change_result(&r, _ki_aux, _li_aux, i+2);
				modify_last_mismatch1(&r, MISMATCH);
				BWExactSearchForward(W, forward, &r, NULL);
				if (r.k<=r.l) add_result(&r, r_list);
			}

		}

	}

	//printf("\n");

	return false;

}

bool BWSimpleSearch1Backward(uint8_t *W, bwt_index *index, result *res, results_list *r_list, uint8_t nA) {

	intmax_t _k,_l, _k_next, _l_next, _k_aux, _l_aux;
	intmax_t results, results_next;
	int16_t start, end, i;

	result r;

	start   = res->start;
	end     = res->pos;

	_k_next = res->k;
	_l_next = res->l;
	results_next = _l_next - _k_next;

	init_result(&r, 0);
	bound_result(&r, start, end);
	add_mismatch(&r, MATCH, -1, start);

	for(i=end; i>=start; i--) {

		_k = _k_next;
		_l = _l_next;

		//printf("%d:\n", i);

		if (_k > _l) {
			change_result(res, _k, _l, -1);
			return false;
		}

		BWiteration(_k, _l, _k_next, _l_next, W[i], index);
		results      = results_next;
		results_next = _l_next - _k_next;
		//printf("(%lu, %lu, %lu)\t", results, _k, _l);

		if (results == results_next) continue;

		//Deletion
		change_result(&r, _k, _l, i-1);
		BWExactSearchBackward(W, index, &r, NULL);
		if (r.k<=r.l) {
			modify_last_mismatch3(&r, DELETION, -1, i);
			add_result(&r, r_list);
		}

		for (uint8_t b=0;b<nA;b++) {

			BWiteration(_k, _l, _k_aux, _l_aux, b, index);

			//printf("W -> %d, %d - %d\n", b, _k_aux, _l_aux);

			if (_k_aux > _l_aux) continue;

			//Insertion
			change_result(&r, _k_aux, _l_aux, i);
			BWExactSearchBackward(W, index, &r, NULL);
			if (r.k<=r.l) {
				modify_last_mismatch3(&r, INSERTION, b, i);
				add_result(&r, r_list);
			}

			//Mismatch
			if (b!=(int)W[i]) {
				change_result(&r, _k_aux, _l_aux, i-1);
				BWExactSearchBackward(W, index, &r, NULL);
				if (r.k<=r.l) {
					modify_last_mismatch3(&r, MISMATCH, b, i);
					add_result(&r, r_list);
				}

			}

		}

	}

	//Match at exit in res
	change_result(res, _k_next, _l_next, -1);

	return false;

}

bool BWSimpleSearch1Forward(uint8_t *W, bwt_index *index, result *res, results_list *r_list, uint8_t nA) {

	intmax_t _k, _l, _k_next, _l_next, _k_aux, _l_aux;
	intmax_t results, results_next;
	int16_t start, end, i;

	result r;

	start   = res->pos;
	end     = res->end;

	_k_next = res->k;
	_l_next = res->l;
	results_next = _l_next - _k_next;

	init_result(&r, 1);
	bound_result(&r, start, end);
	add_mismatch(&r, MATCH, -1, start);

	for(i=start; i<=end; i++) {

		_k = _k_next;
		_l = _l_next;

		//printf("%d:\n", i);

		if (_k > _l) {
			change_result(res, _k, _l, -1);
			return false;
		}

		BWiteration(_k, _l, _k_next, _l_next, W[i], index);
		results      = results_next;
		results_next = _l_next - _k_next;
		if (results == results_next) continue;

		//Deletion
		change_result(&r, _k, _l, i+1);
		BWExactSearchForward(W, index, &r, NULL);
		if (r.k<=r.l) {
			modify_last_mismatch3(&r, DELETION, -1, i);
			add_result(&r, r_list);
		}

		for (uint8_t b=0;b<nA;b++) {

			BWiteration(_k, _l, _k_aux, _l_aux, b, index);

			//printf("W -> %d, %d - %d\n", b, _k_aux, _l_aux);

			if (_k_aux > _l_aux) continue;

			//Insertion
			change_result(&r, _k_aux, _l_aux, i);
			BWExactSearchForward(W, index, &r, NULL);
			if (r.k<=r.l) {
				modify_last_mismatch3(&r, INSERTION, b, i);
				add_result(&r, r_list);
			}

			//Mismatch
			if (b!=(int)W[i]) {
				change_result(&r, _k_aux, _l_aux, i+1);
				BWExactSearchForward(W, index, &r, NULL);
				if (r.k<=r.l) {
					modify_last_mismatch3(&r, MISMATCH, b, i);
					add_result(&r, r_list);
				}
			}

		}

	}

	//Match at exit in res
	change_result(res, _k_next, _l_next, -1);

	return false;

}

bool BWSearch1CPU(uint8_t *W, bwt_index *backward, bwt_index *forward, result *res, results_list *r_list, results_list *r_list_anchor, uint8_t nA) {

	int16_t start, end, half, n;;
	intmax_t _k, _l;

	_k = res->k;
	_l = res->l;

	start = res->start;
	end   = res->end;

	n = end - start + 1;
	half = n / 2;

	result r, r_anchor;

	init_result(&r, 0);
	bound_result(&r, half, end);
	change_result(&r, _k, _l, end);

	if (r_list_anchor == NULL) {
	  BWExactSearchBackward(W, backward, &r, NULL);
	} else {
	  r_list_anchor->num_results = 0; 
	  BWExactSearchBackward(W, backward, &r, &r_anchor);
	  add_result(&r_anchor, r_list_anchor);
	}

	if (r.k <= r.l) {
		r.start = start;
		r.pos = half-1;
		BWSimpleSearch1Backward(W, backward, &r, r_list, nA);

		if (r.k <= r.l) add_result(&r, r_list); //Match
	}

	half--;

	init_result(&r, 1);
	bound_result(&r, start, half);
	change_result(&r, _k, _l, start);

	if (r_list_anchor == NULL) {
	  BWExactSearchForward(W, forward, &r, NULL);
	} else {
	  BWExactSearchForward(W, forward, &r, &r_anchor);
	  add_result(&r_anchor, r_list_anchor);
	}

	if (r.k <= r.l) {
		r.pos = half+1;
		r.end = end;
		BWSimpleSearch1Forward(W, forward, &r, r_list, nA);
	}

	return false;

}
