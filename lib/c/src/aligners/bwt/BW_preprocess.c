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
    GNU General Public License for more details.eo

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "BW_preprocess.h"
#include "BW_io.h"

char *BW_aux;
size_t nB_aux;
size_t chunk;
size_t max_nX_aux;
size_t reads_aux;

int my_sort(const void *a, const void *b) {

  // Casting pointer types

  const unsigned int *ia = (const unsigned int *)a;
  const unsigned int *ib = (const unsigned int *)b;

  char *pa = BW_aux + *ia;
  char *pb = BW_aux + *ib;

  for (size_t i=0; /*i<MAX_SEARCH &&*/ i<nB_aux; i++) {
    if (pa[i] == pb[i]) continue;
    if (pa[i]  < pb[i]) return -1;
    else                return  1;
  }

  return 0;

}

int my_sort_reverse(const void *a, const void *b) {

  // Casting pointer types

  const unsigned int *ia = (const unsigned int *)a;
  const unsigned int *ib = (const unsigned int *)b;

  char *pa = BW_aux + *ia;
  char *pb = BW_aux + *ib;

  for (size_t i=nB_aux - 1; /*i<MAX_SEARCH &&*/ i>=0; i--) {
    if (pa[i] == pb[i]) continue;
    if (pa[i] < pb[i])  return -1;
    else                return  1;
  }

  return 0;

}

void calculateBWTdebug(byte_vector *B, comp_vector *S, byte_vector *X, int reverse) {

  size_t i;

  B->n     = X->n+1;  /*We add the $ character*/
  nB_aux   = B->n;

  S->siz   = B->n;
  S->n     = S->siz;
  S->ratio = 1;

  B->vector  = ( char *) malloc(B->n * sizeof(char));
  checkMalloc(B->vector, "calculateBWTdebug");
  S->vector  = ( unsigned int *) malloc(S->n * sizeof(unsigned int));
  checkMalloc(S->vector, "calculateBWTdebug");

  BW_aux  = X->vector;

  /*Generation of B and S*/

  for(i=0; i<S->n; i++)
    S->vector[i] = i;

  if (reverse) {

    qsort(S->vector, S->n, sizeof(unsigned int), my_sort_reverse);

    for(i=0; i<B->n; i++) {
      B->vector[i] = BW_aux[S->vector[i]];
      S->vector[i] = S->n - (S->vector[i]+1);
    }

  } else {

    qsort(S->vector, S->n, sizeof(unsigned int), my_sort);

    size_t nBlast = X->n; //B.n-1

    for(i=0; i<B->n; i++)
      B->vector[i] = BW_aux[S->vector[i] + nBlast];

  }

#ifdef VERBOSE_DBG
  size_t j;

  for(i=0; i<B->n; i++) {
    for(j=0; j<B->n; j++) {

      if (reverse) {

        if((int)BW_aux[S->n - (S->vector[i] + 1) + j] == -1)
          printf("-");
        else
          printf("%d", BW_aux[S->n - (S->vector[i] + 1) + j]);

      } else {

        if((int)BW_aux[S->vector[i]+j]==-1)
          printf("-");
        else
          printf("%d", BW_aux[S->vector[i]+j]);

      }

    }
    printf("\n");
  }
#endif

}

inline unsigned int ternary_quicksort_start(unsigned int *S, char *X, unsigned int *ranges, unsigned int start, unsigned int end, size_t h) {

  if (h>1) {
    fprintf(stderr, "Ternary quicksort start only works in the two first iterations\n");
    exit(1);
  }

  if (end<=start) {
    fprintf(stderr, "Start and end range is wrong\n");
    exit(1);
  }

  unsigned int start_pivot_pos, start_pivot, end_pivot_pos, end_pivot, l, p, r, value;

  // Put $ suffix at the beginning and ?$ in the end and keep it until sorting ends, just to make second phase faster
  if (h==0) {
    end_pivot_pos = S[end];
    S[end] = S[end - 1];
    S[end - 1] = S[start];
    S[start] = end_pivot_pos;
    start++;
    end--;
  }

  //Select an efficient pivot and put it at the end
  for(value=start; value < end; value++) {
    if (X[S[value]+h] == 1) {end_pivot_pos = S[value]; break;}
  }

  if (value < end) {
    S[value] = S[end];
    S[end] = end_pivot_pos;
  } else {
    printf("ERROR: Pivot not found\n");
    exit(1);
  }

  //Init pivots and perform first pass
  start_pivot_pos = S[start];
  start_pivot = X[start_pivot_pos + h];

  end_pivot_pos = S[end];
  end_pivot = X[end_pivot_pos + h];

  l = start;
  p = start+1;
  r = end;

  while (p != r) {

    value = X[S[p] + h];

    if (value == end_pivot) { //Ternary quicksort (if = does nothing)
      p++;
    } else if (value > end_pivot) {
      S[r] = S[p];
      r--;
      S[p] = S[r];
    } else {
      S[l] = S[p];
      l++;
      S[p] = S[l];
      p++;
    }

  }

  if (start_pivot < end_pivot) {
    S[l] = start_pivot_pos;
    S[r] = end_pivot_pos;
    if (r < end  ) r++;
  } else if (start_pivot > end_pivot) {
    S[r] = start_pivot_pos;
    S[l] = end_pivot_pos;
    if (l > start) l--;
  } else {
    S[r] = start_pivot_pos;
    S[l] = end_pivot_pos;
    if (l > start) l--;
    if (r < end  ) r++;
  }

  ranges[0] = l - start + 1;
  ranges[1] = r - l - 1;

  if (nA==4) {
    //Init pivots for second pass
    start = r;
    
    start_pivot_pos = S[start];
    start_pivot = X[start_pivot_pos + h];

    end_pivot_pos = S[end];
    end_pivot = X[end_pivot_pos + h];
    
    l = start;
    p = start+1;
    r = end;
    
    while (p != r) {
      
      value = X[S[p] + h];
      
      if (value == end_pivot) { //Ternary quicksort (if = does nothing)
	p++;
      } else if (value > end_pivot)  {
	S[r] = S[p];
	r--;
	S[p] = S[r];
      } else {
	S[l] = S[p];
	l++;
	S[p] = S[l];
	p++;
      }
      
    }

    if (start_pivot < end_pivot) {
      S[l] = start_pivot_pos;
      S[r] = end_pivot_pos;
      if (r < end  ) r++;
    } else if (start_pivot > end_pivot) {
      S[r] = start_pivot_pos;
      S[l] = end_pivot_pos;
      if (l > start) l--;
    } else {
      S[r] = start_pivot_pos;
      S[l] = end_pivot_pos;
      if (l > start) l--;
      if (r < end  ) r++;
    }

    if (end_pivot==2) {
      ranges[2] = r - l;
      ranges[3] = end - r + 1;
    } else {
      ranges[2] = l - start + 1;
      ranges[3] = r - l;
    }

  } else if (nA==3) {
    ranges[2] = end - r + 1;
  } else if (nA==2) {
  }

  //Put ?$ on place
  if (h==0) {
    end++;
    value=1;
    int i;
    for (i=0; i < X[S[end] + h]; i++) {
      value += ranges[i];
    }

    end_pivot_pos=S[value];
    end_pivot=X[S[value] + h];
    S[value]=S[end];

    for (; i <nA; i++) {
      value += ranges[i];
      start_pivot_pos = S[value];
      S[value]=end_pivot_pos;
      end_pivot_pos=start_pivot_pos;
    }

    ranges[end_pivot]++;

  } else {

    end_pivot=-1;

  }

  return end_pivot;

}

inline void ternary_quicksort(unsigned int *S, unsigned int *V, int *L, unsigned int start, unsigned int end, size_t h) {

  unsigned int stack[1024*1024];
  int i=0;

  unsigned int start_pivot_pos, start_pivot, end_pivot_pos, end_pivot, left, right, l, p, r, value;

  stack[i++] = start;
  stack[i++] = end;

  if (end<=start) {
    fprintf(stderr, "Start and end range is wrong\n");
    exit(1);
  }

  while (i > 0) {

    right = stack[--i];
    left  = stack[--i];

    if (h>256) {
      value = left + (random() % (right-left));
      end_pivot_pos = S[value];
      S[value] = S[right];
      S[right] = end_pivot_pos;
    }

    //Init pivots
    start_pivot_pos = S[left];
    start_pivot = V[start_pivot_pos];

    end_pivot_pos = S[right];
    end_pivot = V[end_pivot_pos];

    l = left;
    p = left+1;
    r = right;

    while (p != r) {

      value = V[S[p]];

      if (value == end_pivot) { //Ternary quicksort (if = does nothing)
        p++;
      } else if (value > end_pivot) {
        S[r] = S[p];
        r--;
        S[p] = S[r];
      } else {
        S[l] = S[p];
        l++;
        S[p] = S[l];
        p++;
      }

    }

    if (start_pivot < end_pivot) {
      S[l] = start_pivot_pos;
      S[r] = end_pivot_pos;
      if (r < right ) r++;
    } else if (start_pivot > end_pivot) {
      S[l] = end_pivot_pos;
      S[r] = start_pivot_pos;
      if (l > left) l--;
    } else {
      S[r] = start_pivot_pos;
      S[l] = end_pivot_pos;
      if (l > left ) l--;
      if (r < right) r++;
    }

    if (l > left) {
      stack[i++] = left;
      stack[i++] = l;
    }

    if (r < right) {
      stack[i++] = r;
      stack[i++] = right;
    }

  }

}

void calculateBWT(byte_vector *B, comp_vector *S, byte_vector *X, int reverse, exome *ex, char *reference) {

  vector V;
  int **L;
  unsigned int ranges[nA];
  unsigned int offsets[nA];
  unsigned int ranges2[nA*nA];
  unsigned int offsets2[nA*nA];

  load_reference(X, 1, NULL, reference);

  if (reverse) revstring(X->vector, X->n*2 + 1);
  X->n++;

  printUIntVector(X->vector, X->n);

  S->siz    = X->n;
  S->n      = S->siz;
  S->ratio  = 1;
  S->vector = (unsigned int *) malloc(S->n * sizeof(unsigned int));

  for (unsigned int i=0; i < X->n; i++) { //TODO: Paralelizar bucle con OpenMP
    S->vector[i] = i;
  }

  printUIntVector(S->vector, S->n);

  size_t h = 0;
  int dollar;

  printf("**** Suborden: %lu ****\n", h);
  fflush(stdout);

  dollar = ternary_quicksort_start(S->vector, X->vector, ranges, 0, S->n-1, h);

  for (int i=0, group=1; i < nA; i++) {
    offsets[i] = group;
    group += ranges[i];
  }

  printUIntVector(S->vector, S->n);
  printUIntVector(ranges, nA);
  printUIntVector(offsets, nA);

  h = 1;

  printf("**** Suborden: %lu ****\n", h);
  fflush(stdout);

  L = (int**) malloc(nA * nA * sizeof(int*));

#pragma omp parallel for
  for (int r=0; r < nA; r++) {

    if (r==dollar) // Manage ?$
      ternary_quicksort_start(S->vector + offsets[r] + 1, X->vector, ranges2 + nA*r, 0, ranges[r]-2, h);
    else
      ternary_quicksort_start(S->vector + offsets[r]    , X->vector, ranges2 + nA*r, 0, ranges[r]-1, h);

    //Initialize L vectors
    
    unsigned int aux_desp;
    for (int j=0; j < nA; j++) {

      aux_desp=0;
      if(r==0 && j==0)      aux_desp++;
      if(r==dollar && j==0) aux_desp++;

      L[r*nA + j] = (int *) malloc((ranges2[r*nA + j] + aux_desp) * sizeof(int));
      if (aux_desp) L[r*nA + j][0] = -aux_desp;
      L[r*nA + j][aux_desp] = (ranges2[r*nA + j]==1)?-1:ranges2[r*nA + j];

      ranges2[r*nA + j] += aux_desp;

      //printIntVector(L[i*nA + j], ranges2[i*nA + j]);

    }

  }

  for (int i=0, group=0; i < nA * nA; i++) {
    offsets2[i] = group;
    group += ranges2[i];
  }

  free(X->vector);

  V.n = X->n;
  V.vector  = (unsigned int *) malloc(V.n * sizeof(unsigned int));

  printUIntVector(S->vector, S->n);

  printUIntVector(ranges2, nA*nA);
  printUIntVector(offsets2, nA*nA);

  int end;

  do {

    end = 0;
    h = h * 2;

    printf("**** Suborden: %lu ****\n", h);
    fflush(stdout);

#pragma omp parallel for
    for(int r = 0; r < nA*nA; r++) {

      unsigned int offset=offsets2[r], offset_last;
      for(unsigned int i = 0; i<ranges2[r]; i += abs(L[r][i])) {

        offset_last=offset;
        offset += abs(L[r][i]);
 
        if (L[r][i]<0) {
          for(size_t j = offset_last; j < offset; j++)
            V.vector[(S->vector[j] + X->n - h) % X->n] = j;
        } else {
          for(size_t j = offset_last; j < offset; j++)
            V.vector[(S->vector[j] + X->n - h) % X->n] = offset - 1;
        }

      }

    }

    printUIntVector(V.vector, V.n);

#pragma omp parallel for
    for(int r = 0; r < nA*nA; r++) {

      int pos_neg=-1, pos_pos=-1;
      unsigned int pos_size=0;
      unsigned int offset=offsets2[r], offset_last, increment;
      for(unsigned int i = 0; i<ranges2[r]; i += increment) {

        increment = abs(L[r][i]);
        offset_last=offset;
        offset += abs(L[r][i]);

        if (L[r][i]<0) {
          
          if (pos_neg==-1)
            pos_neg=i;
          else
            L[r][pos_neg] += L[r][i];
          
          continue;
          
        }

        //printf("++ Sorting segment %lu - %lu\n", offset_last, offset - 1);
        ternary_quicksort(S->vector + offset_last, V.vector, L[r] + i, 0, increment - 1, h);
        
        //printIntVector(L[r], ranges2[r]);

        pos_pos=-1, pos_size=0;

        unsigned int j;
        for(j = i; j < i + increment - 1; j++) {

          if (V.vector[S->vector[offsets2[r] + j]]==V.vector[S->vector[offsets2[r] + j+1]]) {

            pos_neg = -1;

            if (pos_pos==-1) {
              pos_pos = j;
              pos_size = 2;
            } else {
              pos_size++;
            }

          } else {

            if (pos_pos != -1) {
              L[r][pos_pos] = pos_size;
              pos_pos = -1;
              continue;
            }

            if(pos_neg == -1) {
              L[r][j] = -1;
              pos_neg = j;
            } else {
              L[r][pos_neg]--;
            }

            //printf("Condicion restar %lu,%lu %d\n", j, i, pos_neg != -1);

          }

        }

        //printf(":%u %u\n", V.vector[d_S[j-1]], V.vector[d_S[j]]);

        if (V.vector[S->vector[offsets2[r] + j-1]]!=V.vector[S->vector[offsets2[r] + j]]) {

          if(pos_neg == -1) {
            L[r][j] = -1;
            pos_neg = j;
          } else {
            L[r][pos_neg]--;
          }

        } else {

          if (pos_pos != -1) {
            L[r][pos_pos] = pos_size;
            pos_pos = -1;
          }

        }

      }

      //printIntVector(L[r], ranges2[r]);
      if(L[r][0] != (int) -ranges2[r]) {
        //printf("Ranges => %d != %d\n", L[r][0], (int) -ranges2[r]);
        end=1;
      }

    }

    printUIntVector(S->vector, S->n);

  } while(end);

  free(V.vector);
  for(int i=0; i<nA*nA; i++)
    free(L[i]);

  printf("**** Terminado ****\n");
  fflush(stdout);

  load_reference(X, 1, ex, reference);
  if (reverse) revstring(X->vector, X->n*2 + 1);
  B->n = X->n + 1;
  B->vector = (char *) malloc(B->n * sizeof(char));

  size_t nBlast = B->n - 1;

#pragma omp parallel for //TODO: Comprobar si va bien aqui OpenMP
  for(size_t i=0; i<B->n; i++)
    B->vector[i] = X->vector[S->vector[i] + nBlast];

#ifdef VERBOSE_DBG
  size_t i, j;

  for(i=0; i<B->n; i++) {
    for(j=0; j<B->n; j++) {

      if((int)X->vector[S->vector[i]+j]==-1)
        printf("-");
      else
        printf("%d", X->vector[S->vector[i]+j]);

    }
    printf("\n");
  }
#endif

}

void calculateR(comp_vector *S, comp_vector *R) {

  if (S->ratio != 1) {
    fprintf(stderr, "calculateR: The S vector must be uncompressed\n");
    exit(1);
  }

  R->siz   = S->siz;
  R->n     = S->n;
  R->ratio = S->ratio;

  R->vector = (unsigned int *) malloc(S->n * sizeof(unsigned int));

#pragma omp parallel for
  for(size_t i=0; i<S->n; i++) {
    R->vector[S->vector[i]] = i;
  }

}

void calculateC(vector *C, vector *C1, byte_vector *B, size_t offset) {

  size_t i;

  C->n = nA;
  C1->n = nA;

  C->vector  = (unsigned int *) malloc(C->n * sizeof(unsigned int));
  checkMalloc(C->vector, "calculateC");
  C1->vector = (unsigned int *) malloc(C1->n * sizeof(unsigned int));
  checkMalloc(C1->vector, "calculateC");

  for (i = 0; i < C->n; i++)
    C->vector[i] = 0;

  for (i = 0; i<B->n; i++) {
    if (B->vector[i] == -1) continue;
    C->vector[B->vector[i] + 1]++;
  }

  for (i = 1; i < C->n; i++)
    C->vector[i] += C->vector[i-1];

  for (i = 0; i < C->n; i++)
    C1->vector[i] = C->vector[i] + 1;

}

void calculateO(comp_matrix *O, byte_vector *B) {

#if defined VECTOR_O_32BIT_COMPRESSION

  if (sizeof(int)*CHAR_BIT != 32) {
    fprintf(stderr, "Integer size is not 32: %lu\n", (long unsigned int) sizeof(int)*CHAR_BIT);
    exit(1);
  }

  O->siz = B->n+1;          // Position 0 is -1, so I add one element

  O->n_desp = nA;
  O->m_desp = (O->siz / 32);
  if (O->siz % 32) O->m_desp++;

  O->n_count = nA;
  O->m_count = O->m_desp;

#pragma omp parallel for
  for (size_t i=0; i<O->n_count; i++) {

    O->desp[i]  = (unsigned int *) malloc(O->m_desp * sizeof(unsigned int));
    checkMalloc(O->desp[i], "calculateO");
    O->count[i] = (unsigned int *) malloc(O->m_count * sizeof(unsigned int));
    checkMalloc(O->count[i], "calculateO");

    size_t k=0;
    size_t pos=0;
    unsigned int bit;

    O->desp[i][pos]  = 0;
    O->count[i][pos] = 0;

    /* printf("%d", 0); */

    for (size_t j=0; j<B->n; j++) {

      bit = (j+1) % 32;                //First column is -1 index

      if (!bit) {
        pos++;
        O->desp[i][pos]  = k;
        O->count[i][pos] = 0;          //Initialize to 0 bit-vector
      }

      if (((size_t)B->vector[j]) == i) {
        k++;
        O->count[i][pos] |= 1 << bit;
      }

    }

  }

#elif defined VECTOR_O_64BIT_COMPRESSION

  if (sizeof(unsigned long long) * CHAR_BIT != 64) {
    fprintf(stderr, "Long long size is not 64: %lu\n", (long unsigned int) sizeof(unsigned long long)*CHAR_BIT);
    exit(1);
  }

  O->siz = B->n+1;          // Position 0 is -1, so I add one element

  O->n_desp = nA;
  O->m_desp = (O->siz / 64);
  if (O->siz % 64) O->m_desp++;

  O->n_count = nA;
  O->m_count = O->m_desp;

#pragma omp parallel for
  for (size_t i=0; i<O->n_count; i++) {

    O->desp[i]  = (unsigned int *) malloc(O->m_desp * sizeof(unsigned int));
    checkMalloc(O->desp[i], "calculateO");
    O->count[i] = (unsigned long long *) malloc(O->m_count * sizeof(unsigned long long));
    checkMalloc(O->count[i], "calculateO");

    size_t k=0;
    size_t pos=0;
    unsigned int bit;

    O->desp[i][pos]  = 0;
    O->count[i][pos] = 0;

    /* printf("%d", 0); */
    for (size_t j=0; j<B->n; j++) {

      bit = (j+1) % 64;                //First column is -1 index

      if (!bit) {
        pos++;
        O->desp[i][pos]  = k;
        O->count[i][pos] = 0;          //Initialize to 0 bit-vector
      }

      if (((size_t)B->vector[j]) == i) {
        k++;
        O->count[i][pos] |= 1L << bit;
      }

    }

  }

#else

  O->siz    = B->n+1;
  O->n_desp = nA;
  O->m_desp = O->siz;       // Position 0 is -1, so I add one element

#pragma omp parallel for
  for (size_t i=0; i<O->n_desp; i++) {

    O->desp[i] = (unsigned int *) malloc(O->m_desp * sizeof(unsigned int));
    checkMalloc(O->desp[i], "calculateO");

    size_t k=0;
    O->desp[i][0] = 0;     // First column is -1 index
    for (size_t j=0; j<B->n; j++) {
      if (((size_t)B->vector[j]) == i)
        k++;
      O->desp[i][j+1] = k; // First column is -1 index
    }

  }

#endif

}

void calculateSRcomp(comp_vector *SR, comp_vector *SRcomp, unsigned int ratio) {

  SRcomp->n = (SR->siz / ratio);
  if (SR->siz % ratio) SRcomp->n++;

  SRcomp->siz = SR->siz;
  SRcomp->ratio = ratio;

  SRcomp->vector = (unsigned int *) malloc(SRcomp->n * sizeof(unsigned int));

#pragma omp parallel for
  for (size_t i=0; i < SRcomp->siz; i+=ratio)
    SRcomp->vector[i/ratio] = SR->vector[i];

}
