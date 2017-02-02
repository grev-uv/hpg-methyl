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
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _BW_PREPROCESS_
#define _BW_PREPROCESS_

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<limits.h>

#include"BW_io.h"

#define MAX_SEARCH 200

//--------------------------------------------------------------------------------------
// In case the build is being compiled with all optimizations disabled (ie: debug build)
// inline functions must be forward declared to enable the compiler to find the appropriate
// simbols to link into them on all the translation units.
// - Date: 14 / 11 / 2016
// - Who: Cesar
#ifdef __GNUC__
#ifdef __NO_INLINE__
unsigned int ternary_quicksort_start(unsigned int*, char*, unsigned int*, unsigned int, unsigned int, size_t);
void ternary_quicksort(unsigned int*, unsigned int*, int*, unsigned int, unsigned int, size_t);
#endif
#endif

void calculateBWT(byte_vector *B, comp_vector *S, byte_vector *X, int reverse, exome *ex, char *reference);
void calculateBWTdebug(byte_vector *B, comp_vector *S, byte_vector *X, int reverse);

void calculateR(comp_vector *S, comp_vector *R);
void calculateSRcomp(comp_vector *SR, comp_vector *SRcomp, unsigned int ratio);

void calculateC(vector *C, vector *C1, byte_vector *B, size_t offset);
void calculateO(comp_matrix *O, byte_vector *B);

int my_sort(const void *a, const void *b);

#endif
