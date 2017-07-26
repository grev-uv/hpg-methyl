/*
    bwa_gpu a set of tools which allow short sequence alignment using the Burrows-Wheeler
    transform usign both CPU and GPU approaches.
    Copyright (C) 2011  Jose Salavert Torres, Ignacio Blanquer Espert,
                        Andres Tomas Dominguez, Vicente Hernandez Garcia,
	 		Ignacio Medina Castello, Joaquin Tarraga Gimenez,
			Joaquin Dopazo Blazquez

    Contact e-mail: josator@fiv.upv.es, iblanque@dsic.upv.es, atomas@dsic.upv.es,
                    vhernand@dsic.upv.es, imedina@cipf.es, jtarraga@cipf.es,
                    jdopazo@cipf.es

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _BW_SEARCH_
#define _BW_SEARCH_

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "BW_io.h"

#define BWiterationVariables() size_t aux_b;

#if defined VECTOR_O_32BIT_COMPRESSION || VECTOR_O_64BIT_COMPRESSION

#define BWiteration(k_in, l_in, k_out, l_out, b, C, C1, O)		 \
  {									 \
    aux_b = (b);							 \
    (k_out) = (C1)->vector[aux_b] + getOcompValue(aux_b, (k_in)  , (O)); \
    (l_out) = (C)->vector[aux_b]  + getOcompValue(aux_b, (l_in)+1, (O)); \
  }
//printf("k-> %lu, l-> %lu, O -> %u\n", k_out, l_out, getOcompValue(aux_b, (l_in)+1, (O)));

#else

#define BWiteration(k_in, l_in, k_out, l_out, b, C, C1, O)		\
  {									\
    aux_b = (b);							\
    (k_out) = (C1)->vector[aux_b] + (O)->desp[aux_b][(k_in)];		\
    (l_out) = (C)->vector[aux_b]  + (O)->desp[aux_b][(l_in)+1];		\
  }
//printf("k-> %lu, l-> %lu, O -> %u\n", k_out, l_out, (O)->desp[aux_b][(l_in)+1]);

#endif

inline void change_direction(comp_vector *S, comp_vector *Ri, vector *C, comp_matrix *O, comp_matrix *Oi, result *res) {

  size_t k, l, ki, li, aux, aux2;
  int start, end;

  k  = res->k;
  l  = res->l;
  ki = O->siz-2;
  li = 0;

  start = res->start;
  end   = res->end;

  if (S->ratio == 1) {

    for (size_t i = k; i <= l; i++) {

      aux  = S->siz - S->vector[i]  - (end - start + 2);
      aux2 = Ri->vector[aux];

      if (aux2 < ki) ki = aux2;
      if (aux2 > li) li = aux2;

    }

  } else {

    for (size_t i = k; i <= l; i++) {

      aux  = S->siz - getScompValue(i, S, C, O) - (end - start + 2);
      aux2 = getRcompValue(aux, Ri, C, Oi);

      if (aux2 < ki) ki = aux2;
      if (aux2 > li) li = aux2;

    }

  }

  res->k = ki;
  res->l = li;

}

//--------------------------------------------------------------------------------------
// In case the build is being compiled with all optimizations disabled (ie: debug build)
// inline functions must be forward declared to enable the compiler to find the appropriate
// simbols to link into them on all the translation units.
// - Date: 14 / 11 / 2016
// - Who: Cesar
#ifdef __GNUC__
#ifdef __NO_INLINE__
static inline void BWExactSearchBackward(char *W, vector *C, vector *C1, comp_matrix *O, result *r) {
#else
inline void BWExactSearchBackward(char *W, vector *C, vector *C1, comp_matrix *O, result *r) {
#endif /* __NO_INLINE__ */
#else
inline void BWExactSearchBackward(char *W, vector *C, vector *C1, comp_matrix *O, result *r) {
#endif /* __GNUC__ */
  BWiterationVariables();
  size_t k2, l2;
  int i;

  k2 = r->k;
  l2 = r->l;

  //printf("B1º -> %lu - %lu\n", k2, l2);

  for(i=r->pos; i>=r->start; i--) {

    BWiteration(k2, l2, k2, l2, W[i], C, C1, O);
    //printf("B -> %d -> %lu - %lu\n", i, k2, l2);
    if (k2 > l2) break;

  }

  r->k = k2;
  r->l = l2;
  r->pos = i;

}

//--------------------------------------------------------------------------------------
// In case the build is being compiled with all optimizations disabled (ie: debug build)
// inline functions must be forward declared to enable the compiler to find the appropriate
// simbols to link into them on all the translation units.
// - Date: 14 / 11 / 2016
// - Who: Cesar
#ifdef __GNUC__
#ifdef __NO_INLINE__
static inline void BWExactSearchForward(char *W, vector *C, vector *C1, comp_matrix *Oi, result *r) {
#else
inline void BWExactSearchForward(char *W, vector *C, vector *C1, comp_matrix *Oi, result *r) {
#endif /* __NO_INLINE__ */
#else
inline void BWExactSearchForward(char *W, vector *C, vector *C1, comp_matrix *Oi, result *r) {
#endif /* __GNUC__ */
  BWiterationVariables();
  size_t k2, l2;
  int i;

  k2 = r->k;
  l2 = r->l;

  //printf("F1º -> %lu - %lu\n", k2, l2);

  for(i=r->pos; i<=r->end; i++) {
    BWiteration(k2, l2, k2, l2, W[i], C, C1, Oi);
    //printf("F-> %d -> %lu - %lu\n", i, k2, l2);
    if (k2 > l2) break;

  }

  r->k = k2;
  r->l = l2;
  r->pos = i;

}

//--------------------------------------------------------------------------------------
// In case the build is being compiled with all optimizations disabled (ie: debug build)
// inline functions must be forward declared to enable the compiler to find the appropriate
// simbols to link into them on all the translation units.
// - Date: 14 / 11 / 2016
// - Who: Cesar
#ifdef __GNUC__
#ifdef __NO_INLINE__
static void change_direction(comp_vector*, comp_vector*, vector*, comp_matrix*, comp_matrix*, result*);
#endif
#endif

void BWExactFinalResultsBackward(char *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next);
void BWExactFinalResultsForward(char *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next);
void BWExactPartialResultsBackward(char *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next);
void BWExactPartialResultsForward(char *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next);
void BWBranchPartialResultsBackward(char *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next);
void BWBranchPartialResultsForward(char *W, vector *C, vector *C1, comp_matrix *O, results_list *rl_prev, results_list *rl_next);

size_t BWExactSearchVectorBackward(char *W, int start, int end, size_t k, size_t l,
				 size_t *vec_k, size_t *vec_l, vector *C, vector *C1, 
				 comp_matrix *O, size_t *last_k, size_t *last_l, int *nt);
size_t BWExactSearchVectorForward(char *W, int start, int end, size_t k, size_t l,
				size_t *vec_k, size_t *vec_l, vector *C, vector *C1,
				comp_matrix *O, size_t *last_k, size_t *last_l, int *nt);

void BWSearch1(char *W, int start, int end, size_t *vec_k, size_t *vec_l, size_t *vec_ki, size_t *vec_li, vector *C, vector *C1, comp_matrix *O, comp_matrix *Oi, results_list *r_list);

void BWSearch1CPU(char *W, vector *C, vector *C1, comp_matrix *O, comp_matrix *Oi, result * res, results_list *r_list);

void BWSearchCPUBackward(char *W, vector *C, vector *C1, comp_matrix *O, result *res, results_list *rl_prev, results_list *rl_next, int num_errors);
void BWSearchCPUForward(char *W, vector *C, vector *C1, comp_matrix *O, result *res, results_list *rl_prev, results_list *rl_next, int num_errors);

void BWSimpleSearch1Backward(char *W, vector *C, vector *C1, comp_matrix *O, result *res, results_list *r_list);
void BWSimpleSearch1Forward(char *W, vector *C, vector *C1, comp_matrix *O, result *res, results_list *r_list);

#endif
