#ifndef _SEARCH_PREPROCESS_
#define _SEARCH_PREPROCESS_

#include <omp.h>

#include "../dbwt/dbwt.h"
#include "csafm.h"
#include "io.h"

void calculate_and_save_B(ref_vector *X, const char *directory, const char *name);
void calculate_C(vector *C, vector *C1, ref_vector *B, uint8_t nA);
void calculate_O(comp_matrix *O, ref_vector *B, uint8_t nA);
void calculate_S_and_R(comp_vector *S, comp_vector *R, ref_vector *B, vector *C, comp_matrix *O, SA_TYPE ratio);

void compress_S_or_R(comp_vector *SRcomp, comp_vector *SR, SA_TYPE ratio);

#endif

/*
This file is part of GNU BWT Aligner.

		Copyright José Salavert Torres, Ignacio Blanquer Espert
							Universitat Politècnica of València
							 - Institut de Instrumentació per a la Imatge Molecular (GRyCAP)
							with partial support of Bull (S.A) Informatique
							2010 - 2013

    GNU BWT Aligner is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GNU BWT Aligner is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with GNU BWT Aligner.  If not, see <http://www.gnu.org/licenses/>.
*/
