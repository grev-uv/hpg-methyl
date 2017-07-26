#ifndef _SEARCH_RESULTS_
#define _SEARCH_RESULTS_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "runtime.h"
#include "io.h"

#define MAX_MISMATCHES 20

typedef struct {
	bool dir;                      //0 - Backward, 1 - Forward
	intmax_t k, l;
	int16_t start, pos, end;
	uint8_t num_mismatches;
	uint8_t err_kind[MAX_MISMATCHES];
	uint8_t err_base[MAX_MISMATCHES];
	int16_t err_pos[MAX_MISMATCHES];
} result;

typedef struct {
	result *list;
	uintmax_t num_results;
	uintmax_t max_results;
	uintmax_t read_index;
} results_list;

#define new_results_list(_r_list, _max_results)\
	do {\
		(_r_list)->list = (result *) malloc((_max_results) * sizeof(result));\
		check_malloc((_r_list)->list, "new_result_list");\
		(_r_list)->num_results = 0;\
		(_r_list)->max_results = (_max_results);\
	} while(0);

#define init_result(_r, _dir)\
	do {\
		(_r)->dir = (_dir);\
		(_r)->num_mismatches = 0;\
	} while(0);

#define bound_result(_r, _start, _end)\
	do {\
		(_r)->start = (_start);\
		(_r)->end = (_end);\
	} while(0);

#define change_result(_r, _k, _l, _pos)\
	do {\
		(_r)->k = (_k);\
		(_r)->l = (_l);\
		(_r)->pos = (_pos);\
	} while(0);

#define add_mismatch(_r, _err_kind, _base, _position)\
	do {\
		(_r)->err_kind[(_r)->num_mismatches] = (_err_kind);\
		(_r)->err_base[(_r)->num_mismatches] = (_base);\
		(_r)->err_pos[(_r)->num_mismatches] = (_position);\
		(_r)->num_mismatches++;\
	} while(0);

#define modify_last_mismatch3(_r, _err_kind, _base, _position)\
	do {\
		(_r)->err_kind[(_r)->num_mismatches-1] = (_err_kind);\
		(_r)->err_base[(_r)->num_mismatches-1] = (_base);\
		(_r)->err_pos[(_r)->num_mismatches-1] = (_position);\
	} while(0);

#define modify_last_mismatch2(_r, _err_kind, _base)\
	do {\
		(_r)->err_kind[(_r)->num_mismatches-1] = (_err_kind);\
		(_r)->err_base[(_r)->num_mismatches-1] = (_base);\
	} while(0);

#define modify_last_mismatch1(_r, _err_kind)\
	do {\
		(_r)->err_kind[(_r)->num_mismatches-1] = (_err_kind);\
	} while(0);

#define copy_result(_dest, _orig)\
	do {\
		init_result((_dest), (_orig)->dir);\
		bound_result((_dest), (_orig)->start, (_orig)->end);\
		change_result((_dest), (_orig)->k, (_orig)->l, (_orig)->pos);\
		for(uint8_t rr=0; rr<(_orig)->num_mismatches; rr++) {\
			add_mismatch((_dest), (_orig)->err_kind[rr], (_orig)->err_base[rr], (_orig)->err_pos[rr]);\
		}\
	} while(0);

#define add_result(orig, _r_list)\
	do {\
		if ((_r_list)->num_results < (_r_list)->max_results) {\
			result *dest;\
			dest = &(_r_list)->list[(_r_list)->num_results];\
			(_r_list)->num_results++;\
			copy_result(dest, (orig));\
		} else {\
			return true;\
		}\
	} while(0);

static inline void concat_error_string(char *mask, char *mask_aux, result *r, uint8_t rr, uintmax_t *enW, bwt_config_t *bwt_config) {

	if      (r->err_kind[rr]==DELETION)
		(*enW)--;
	else if (r->err_kind[rr]==INSERTION)
		(*enW)++;

	if      (r->err_kind[rr]==DELETION) {
		sprintf(mask_aux, "_%d%c",   r->err_pos[rr], 'd');
		strcat(mask, mask_aux);
	} else if (r->err_kind[rr]==MISMATCH) {
		sprintf(mask_aux, "_%d%c%c", r->err_pos[rr], 'm', bwt_config->nucleotides[r->err_base[rr]]);
		strcat(mask, mask_aux);
	} else {
		sprintf(mask_aux, "_%d%c%c", r->err_pos[rr], 'i', bwt_config->nucleotides[r->err_base[rr]]);
		strcat(mask, mask_aux);
	}

}

static inline void manage_single_result(result *r, exome* ex, bwt_index *backward, bwt_index *forward, char *search, unsigned int type, FILE *fp, uintmax_t read_index, bool *found, bwt_config_t *bwt_config) {

	bool direction;
	uintmax_t enW;
	char plusminus[] = "-+";
	uintmax_t index, index2, key, key2;
	int type2=0;

	char mask[6*(MAXLINE+1)];
	char mask_aux[6];

	if        (type==0) {
		direction = !r->dir;
		type2 = 0;
	} else if (type==1) {
		direction = r->dir;
		type2 = 1;
	} else {
		direction = r->dir;
	}

	enW = r->end - r->start + 1;

	mask[0] = '\0';

	for (int rr=0; rr<r->num_mismatches; rr++) {
	  concat_error_string(mask, mask_aux, r, rr, &enW, bwt_config);
	}

	//printf("%d %d %d %u %u\n", r->start, r->pos, r->end, r->k, r->l);

	for (intmax_t j=r->k; j<=r->l; j++) {

		if (direction)
			key = size_SA(forward) - get_SA(j, forward) - enW -1;
		else
			key = get_SA(j, backward);

		index = binsearch(ex->offset, ex->size, key);

		if(key + enW <= ex->offset[index]) {

			if (type>1) {

				if (key >= (uintmax_t) (size_SA(backward) / 2)) {
					index2 = ex->size - index;
					key2 = ex->start[index2] + (size_SA(forward) - key - enW - 1) - ex->offset[index2];
					type2 = 0; 
				} else {
		  		key2 = ex->start[index-1] + key - ex->offset[index-1];
					type2 = 1;
				}

			} else {

				key2 = ex->start[index-1] + key - ex->offset[index-1];

			}

			fprintf(fp, "read_%ju %c %s %ju %d%s %s\n", (uintmax_t) read_index, plusminus[type2], ex->chromosome + (index-1)*IDMAX, key2, r->num_mismatches, mask, search);
			*found=1;

		}

	}

}

bool write_results(results_list *r_list, intmax_t *k, intmax_t *l, exome* ex, bwt_index *backward, bwt_index *forward, char *mapping, uintmax_t nW, int type, FILE *fp, bwt_config_t *bwt_config);

#endif
