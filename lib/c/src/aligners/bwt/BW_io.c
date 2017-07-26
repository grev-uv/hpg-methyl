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

#include "BW_io.h"

char init_mask[MAXLINE+1];
char plusminus[] = "-+";

void freeCompMatrix(comp_matrix *matrix) {

for (size_t i=0; i<matrix->n_desp; i++) {
    free(matrix->desp[i]);
#if defined VECTOR_O_32BIT_COMPRESSION || VECTOR_O_64BIT_COMPRESSION
    free(matrix->count[i]);
#endif
  }

}

//TODO: Arreglar para mapear el complementario reverso cuando trabajo con
//el análisis con menos de 4 bases.
void reverseStrandC(vector *r_C, vector *s_C, vector *r_C1, vector *s_C1) {

  r_C->n  = s_C->n; r_C1->n = s_C1->n;

  //printf("%d %d %d %d\n", AA, CC, GG, TT);

  r_C->vector   = (unsigned int*)malloc(r_C->n * sizeof(unsigned int));
  checkMalloc(r_C->vector,  "reverseStrandC");
  r_C1->vector  = (unsigned int*)malloc(r_C1->n * sizeof(unsigned int));
  checkMalloc(r_C1->vector,  "reverseStrandC");

  if (AA != -1 && TT !=-1) {
    r_C->vector[AA] = s_C->vector[TT]; r_C1->vector[AA] = s_C1->vector[TT];
    r_C->vector[TT] = s_C->vector[AA]; r_C1->vector[TT] = s_C1->vector[AA];
  } else if (AA != -1) {
    r_C->vector[AA] = s_C->vector[AA]; r_C1->vector[AA] = s_C1->vector[AA];
  } else if (TT != -1) {
    r_C->vector[TT] = s_C->vector[TT]; r_C1->vector[TT] = s_C1->vector[TT];
  }

  if (CC != -1 && GG !=-1) {
    r_C->vector[CC] = s_C->vector[GG]; r_C1->vector[CC] = s_C1->vector[GG];
    r_C->vector[GG] = s_C->vector[CC]; r_C1->vector[GG] = s_C1->vector[CC];
  } else if (CC != -1) {
    r_C->vector[CC] = s_C->vector[CC]; r_C1->vector[CC] = s_C1->vector[CC];
  } else if (GG != -1) {
    r_C->vector[GG] = s_C->vector[GG]; r_C1->vector[GG] = s_C1->vector[GG];
  }

}

void reverseStrandO(comp_matrix *r_O, comp_matrix *s_O) {

  r_O->siz = s_O->siz;

  r_O->n_desp = s_O->n_desp;
  r_O->m_desp = s_O->m_desp;

  if (AA != -1 && TT !=-1) {
    r_O->desp[AA] = s_O->desp[TT];
    r_O->desp[TT] = s_O->desp[AA];
  } else if (AA != -1) {
    r_O->desp[AA] = s_O->desp[AA];
  } else if (TT != -1) {
    r_O->desp[TT] = s_O->desp[TT];
  }

  if (CC != -1 && GG !=-1) {
    r_O->desp[CC] = s_O->desp[GG];
    r_O->desp[GG] = s_O->desp[CC];
  } else if (CC != -1) {
    r_O->desp[CC] = s_O->desp[CC];
  } else if (GG != -1) {
    r_O->desp[GG] = s_O->desp[GG];
  }

#if defined VECTOR_O_32BIT_COMPRESSION || VECTOR_O_64BIT_COMPRESSION

  r_O->n_count = s_O->n_count;
  r_O->m_count = s_O->m_count;

  if (AA != -1 && TT !=-1) {
    r_O->count[AA] = s_O->count[TT];
    r_O->count[TT] = s_O->count[AA];
  } else if (AA != -1) {
    r_O->count[AA] = s_O->count[AA];
  } else if (TT != -1) {
    r_O->count[TT] = s_O->count[TT];
  }

  if (CC != -1 && GG !=-1) {
    r_O->count[CC] = s_O->count[GG];
    r_O->count[GG] = s_O->count[CC];
  } else if (CC != -1) {
    r_O->count[CC] = s_O->count[CC];
  } else if (GG != -1) {
    r_O->count[GG] = s_O->count[GG];
  }

#endif

}

void readUIntVector(vector *vector, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];
  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".vec");

  fp  = fopen(path,  "rb+");
  checkFileOpen(fp, path);

  err = fread(&vector->n, sizeof(size_t),  1, fp);
  checkFileRead(err, 1, path);

  vector->vector = (unsigned int *) malloc(vector->n * sizeof(unsigned int));
  checkMalloc(vector->vector, path);

  err = fread(vector->vector, sizeof(unsigned int), vector->n, fp);
  checkFileRead(err, vector->n, path);

  fclose(fp);

}

void readUIntCompVector(comp_vector *vector, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];
  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".vec");

  fp  = fopen(path,  "rb+");
  checkFileOpen(fp, path);

  err = fread(&vector->n, sizeof(size_t),  1, fp);
  checkFileRead(err, 1, path);

  err = fread(&vector->siz, sizeof(size_t),  1, fp);
  checkFileRead(err, 1, path);

  err = fread(&vector->ratio, sizeof(size_t),  1, fp);
  checkFileRead(err, 1, path);

  vector->vector = (unsigned int *) malloc(vector->n * sizeof(unsigned int));
  //int page_size = sysconf(_SC_PAGESIZE);
  //vector->vector = (int *) memalign(page_size, vector->n * sizeof(unsigned int));

  checkMalloc(vector->vector, path);

  err = fread(vector->vector, sizeof(unsigned int), vector->n, fp);
  checkFileRead(err, vector->n, path);

  fclose(fp);

}

void readCharVector(byte_vector *vector, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".vec");

  fp  = fopen(path,  "rb+");
  checkFileOpen(fp, path);

  err = fread(&vector->n, sizeof(size_t),  1, fp);
  checkFileRead(err, 1, path);

  vector->vector = (char *) malloc(vector->n * sizeof(char));
  //int page_size = sysconf(_SC_PAGESIZE);
  //vector->vector = (int *) memalign(page_size, vector->n * sizeof(unsigned int));

  checkMalloc(vector->vector, path);

  err = fread(vector->vector, sizeof(char), vector->n, fp);
  checkFileRead(err, vector->n, path);

  fclose(fp);

}

void readCompMatrix(comp_matrix *matrix, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".desp");

  fp  = fopen(path,  "rb+");
  checkFileOpen(fp, path);

  err = fread(&matrix->siz,      sizeof(size_t),  1, fp);
  checkFileRead(err, 1, path);

  err = fread(&matrix->n_desp,   sizeof(size_t),  1, fp);
  checkFileRead(err, 1, path);

  err = fread(&matrix->m_desp,   sizeof(size_t),  1, fp);
  checkFileRead(err, 1, path);

  for (size_t i=0; i<matrix->n_desp; i++) {
    matrix->desp[i] = (unsigned int *) malloc(matrix->m_desp * sizeof(unsigned int));
    checkMalloc(matrix->desp[i], path);
    err = fread(matrix->desp[i], sizeof(unsigned int), matrix->m_desp, fp);
    checkFileRead(err, matrix->m_desp, path);
  }

  fclose(fp);

#if defined VECTOR_O_32BIT_COMPRESSION

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".count");

  fp  = fopen(path,  "rb+");
  checkFileOpen(fp, path);

  err = fread(&matrix->n_count,   sizeof(size_t),  1, fp);
  checkFileRead(err, 1, path);

  err = fread(&matrix->m_count,   sizeof(size_t),  1, fp);
  checkFileRead(err, 1, path);

  for (size_t i=0; i<matrix->n_count; i++){
    matrix->count[i] = (unsigned int *) malloc(matrix->m_count * sizeof(unsigned int));
    checkMalloc(matrix->count[i], path);
    err = fread(matrix->count[i], sizeof(unsigned int), matrix->m_count, fp);
    checkFileRead(err, matrix->m_count, path);
  }

  fclose(fp);

#elif defined VECTOR_O_64BIT_COMPRESSION

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".count");

  fp  = fopen(path,  "rb+");
  checkFileOpen(fp, path);

  err = fread(&matrix->n_count,   sizeof(size_t),  1, fp);
  checkFileRead(err, 1, path);

  err = fread(&matrix->m_count,   sizeof(size_t),  1, fp);
  checkFileRead(err, 1, path);

  for (size_t i=0; i<matrix->n_count; i++){
    matrix->count[i] = (unsigned long long *) malloc(matrix->m_count * sizeof(unsigned long long));
    checkMalloc(matrix->count[i], path);
    err = fread(matrix->count[i], sizeof(unsigned long long), matrix->m_count, fp);
    checkFileRead(err, matrix->m_count, path);
  }

  fclose(fp);

#endif

}

void saveUIntVector(vector *vector, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".vec");

  fp  = fopen(path,  "wb+");
  checkFileOpen(fp, path);

  err = fwrite(&vector->n,     sizeof(size_t), 1, fp);
  checkFileWrite(err, 1, path);

  err = fwrite(vector->vector, sizeof(unsigned int), vector->n, fp);
  checkFileWrite(err, vector->n, path);

  fclose(fp);

}

void saveUIntCompVector(comp_vector *vector, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".vec");

  fp  = fopen(path,  "wb+");
  checkFileOpen(fp, path);

  err = fwrite(&vector->n,     sizeof(size_t), 1, fp);
  checkFileWrite(err, 1, path);

  err = fwrite(&vector->siz,   sizeof(size_t), 1, fp);
  checkFileWrite(err, 1, path);

  err = fwrite(&vector->ratio, sizeof(size_t), 1, fp);
  checkFileWrite(err, 1, path);

  err = fwrite(vector->vector, sizeof(unsigned int), vector->n, fp);
  checkFileWrite(err, vector->n, path);

  fclose(fp);

}

void saveCharVector(byte_vector *vector, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".vec");

  fp  = fopen(path,  "wb+");
  checkFileOpen(fp, path);

  err = fwrite(&vector->n,     sizeof(size_t), 1, fp);
  checkFileWrite(err, 1, path);

  err = fwrite(vector->vector, sizeof(char), vector->n, fp);
  checkFileWrite(err, vector->n, path);

  fclose(fp);

}

void saveCompMatrix(comp_matrix *matrix, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".desp");

  fp  = fopen(path,  "wb+");
  checkFileOpen(fp, path);

  err = fwrite(&matrix->siz,    sizeof(size_t), 1, fp);
  checkFileWrite(err, 1, path);

  err = fwrite(&matrix->n_desp, sizeof(size_t), 1, fp);
  checkFileWrite(err, 1, path);

  err = fwrite(&matrix->m_desp, sizeof(size_t), 1, fp);
  checkFileWrite(err, 1, path);

  for (size_t i=0; i<matrix->n_desp; i++) {
    err = fwrite(matrix->desp[i], sizeof(unsigned int), matrix->m_desp, fp);
    checkFileWrite(err, matrix->m_desp, path);
  }

  fclose(fp);

#if defined VECTOR_O_32BIT_COMPRESSION

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".count");

  fp  = fopen(path,  "wb+");
  checkFileOpen(fp, path);

  err = fwrite(&matrix->n_count, sizeof(size_t), 1, fp);
  checkFileWrite(err, 1, path);

  err = fwrite(&matrix->m_count, sizeof(size_t), 1, fp);
  checkFileWrite(err, 1, path);

  for (size_t i=0; i<matrix->n_count; i++) {
    err = fwrite(matrix->count[i], sizeof(unsigned int), matrix->m_count, fp);
    checkFileWrite(err, matrix->m_count, path);
  }

  fclose(fp);

#elif defined VECTOR_O_64BIT_COMPRESSION

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".count");

  fp  = fopen(path,  "wb+");
  checkFileOpen(fp, path);

  err = fwrite(&matrix->n_count, sizeof(size_t), 1, fp);
  checkFileWrite(err, 1, path);

  err = fwrite(&matrix->m_count, sizeof(size_t), 1, fp);
  checkFileWrite(err, 1, path);

  for (size_t i=0; i<matrix->n_count; i++) {
    err = fwrite(matrix->count[i], sizeof(unsigned long long), matrix->m_count, fp);
    checkFileWrite(err, matrix->m_count, path);
  }

  fclose(fp);

#endif

}

int nextFASTAToken(FILE *queries_file, char *uncoded, char *coded, unsigned int *nquery, char *compressed, unsigned int*ncompress) {

  char line[MAXLINE];
  unsigned int length;

  *nquery=0;
  uncoded[0]='\0';

  while ( fgets(line, MAXLINE, queries_file) ) {

    if (line[0] == '>') {
      if (*nquery) break;
      else continue;
    }

    length=strlen(line);
    if (line[length-1]=='\n')
      length--;

    uncoded[*nquery] = '\0';
    strncpy(uncoded + *nquery, line , length);

    *nquery += length;

  }

  if (*nquery) {

    //replaceBases(uncoded, coded, *nquery);
    encodeBases(coded, uncoded, *nquery);

    if (compressed != NULL)
      *ncompress = comp4basesInChar(coded, *nquery, compressed);

    return 1;

  } else {

    return 0;

  }

}

size_t comp4basesInChar(char *X, size_t nX, char *Y) {

  size_t desp=2, iaux=0;
  char aux=0;
  size_t i;

  if (nX>0) aux = X[0];

  for (i=1; i<nX; i++) {

    if (desp==8) {
      desp=0;
      Y[iaux] = aux;
      iaux++ ;
      aux = aux & 0;
    }

    aux = aux | (X[i] << desp);
    desp = desp +2;

  }

  Y[iaux] = aux;

  unsigned int nY = nX / 4;
  if (nX % 4) nY++;

  return nY;
}

void revstring(char *X, size_t nX) {

  char tmp;
  size_t i, j;

  if (nX <= 1) return;

  for (i=0, j=nX-1; i<=j; i++, j--) {
    tmp = X[i];
    X[i] = X[j];
    X[j] = tmp;
  }

}

unsigned int binsearch(unsigned *array, unsigned int size, size_t key) {

  if( !array ) return 0;

  unsigned int *p = array;
  unsigned int w;

  while( size > 0 ) {

    w=size/2;
    
    if ( p[w] <= key ) {
      p+=w+1;
      size-=w+1;
    } else {
      size=w;
    }

  }

  return p - array;

}

void duplicate_reference(byte_vector *X) {

  size_t naux = (X->n)*2 + 1;

  for (size_t i=(X->n)+1, j=0; i<naux; i++,j++)
    X->vector[i] = X->vector[j];

}

void load_reference(byte_vector *X, int duplicate, exome *ex, const char *path) {

  FILE *ref_file;
  ref_file = fopen(path,  "r");
  checkFileOpen(ref_file, path);

  size_t read=0, size;

  fseek(ref_file, 0, SEEK_END);
  read = ftell(ref_file);
  fseek(ref_file, 0, SEEK_SET);

  if (duplicate) size = 2*read + 1000;
  else           size = read + 1000;

  X->vector = (char *) malloc( size * sizeof(char) );
  checkMalloc(X->vector, path);

  if (ex !=NULL) ex->size=0;

  //char line[MAXLINE];
  unsigned int length, partial_length, total_length;

  total_length=0;
  partial_length=0;

  printf("START FUNC\n");

  while ( fgets(X->vector + total_length, MAXLINE, ref_file) ) {

    if ( (X->vector + total_length)[0] == '>') {

      if (ex!=NULL) {

	if (total_length == 0) {

	  sscanf(X->vector + total_length, ">%s ", ex->chromosome + ex->size * IDMAX);
	  ex->start[ex->size] = 0;

	} else {

	  ex->end[ex->size] = partial_length - 1;
	  partial_length=0;

	  if (ex->size==0)
	    ex->offset[0] = 0;
	  else
	    ex->offset[ex->size] = ex->offset[ex->size-1] + (ex->end[ex->size-1] - ex->start[ex->size-1] + 1);
	  ex->size++;

	  sscanf(X->vector + total_length, ">%s ", ex->chromosome + ex->size * IDMAX);
	  ex->start[ex->size] = 0;

	}

      }

      continue;

    }

    length = strlen(X->vector + total_length);
    if ((X->vector + total_length)[length-1]=='\n')
      length--;

    partial_length += length;
    total_length += length;

  }

  if (ex != NULL) {
    ex->end[ex->size] = partial_length - 1;
    partial_length=0;
    
    if (ex->size==0)
      ex->offset[0] = 0;
    else
      ex->offset[ex->size] = ex->offset[ex->size-1] + (ex->end[ex->size-1] - ex->start[ex->size-1] + 1);
    ex->size++;
  }

  //replaceBases(X->vector, X->vector, total_length);
  encodeBases(X->vector, X->vector, total_length);

  X->vector[total_length] = -1; // dollar symbol
  X->n = total_length;

  fclose(ref_file);

  if (duplicate) duplicate_reference(X);
  
  printf("FIN FUNC\n");

}

void load_exome_file(exome *ex, const char *directory) {

  FILE *fp;

  char path[500];
  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/index");

  fp  = fopen(path,  "r");
  checkFileOpen(fp, path);

  int j;
  char c;
  char line[MAXLINE];

  ex->offset[0]=0;
  ex->size = 0;

  while (fgets(line, MAXLINE, fp) ) {

    if (line[0]=='>') {

      c = 0;
      for(j = 0; j < IDMAX-1; j++) {
	c = line[j+1];
	if (c==' ') break;
	ex->chromosome[ex->size*IDMAX+j] = c;
      }

      ex->chromosome[ex->size*IDMAX+j] = '\0';

      sscanf(line + j + 2,"%u %u %*s", &ex->start[ex->size], &ex->end[ex->size]);
      //printf(">%u %s %u %u\n", ex->offset[ex->size], ex->chromosome + ex->size*IDMAX, ex->start[ex->size], ex->end[ex->size]);
      ex->size++;
      ex->offset[ex->size] = ex->offset[ex->size-1] + (ex->end[ex->size-1] - ex->start[ex->size-1]+1);
    }
  }

  fclose(fp);
}

void save_exome_file(exome *ex, const char *directory) {

  FILE *fp;

  char path[500];
  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/index");

  fp  = fopen(path, "w");
  checkFileOpen(fp, path);

  for(size_t i=0; i<ex->size; i++) {
    fprintf(fp, ">%s %u %u\n", ex->chromosome + i*IDMAX, ex->start[i], ex->end[i]);
  }

}

void initialize_init_mask() {

  size_t i;
  for (i=0; i<MAXLINE; i++)
    init_mask[i] = '-';
  init_mask[i] = '\0';

}

void concat_error_string(char *mask, char *mask_aux, result *r, int rr, int *enW) {

  if      (r->err_kind[rr]==DELETION)
    (*enW)--;
  else if (r->err_kind[rr]==INSERTION)
    (*enW)++;

  if      (r->err_kind[rr]==DELETION) {
    sprintf(mask_aux, "_%d%c",   r->err_pos[rr], 'd');
    strcat(mask, mask_aux);
  } else if (r->err_kind[rr]==MISMATCH) {
    sprintf(mask_aux, "_%d%c%c", r->err_pos[rr], 'm', rev_table[(unsigned char)r->err_base[rr]]);
    strcat(mask, mask_aux);
  } else {
    sprintf(mask_aux, "_%d%c%c", r->err_pos[rr], 'i', rev_table[(unsigned char)r->err_base[rr]]);
    strcat(mask, mask_aux);
  }

}

int write_results(results_list *r_list, exome* ex, comp_vector *S, comp_vector *Si, vector *C, comp_matrix *O, comp_matrix *Oi, char *mapping, int nW, int type, FILE *fp) {

  result *r;

  unsigned int index, key;
  int enW;
  int found=0;

  char search[MAXLINE+1];
  char mask[6*(MAXLINE+1)];
  char mask_aux[6];

  int direction;

  search[0] = '\0';
  strncat(search, mapping, nW);

  for (size_t i=0;i<r_list->num_results; i++) {

    r = &r_list->list[i];

    if (type) {
      direction = r->dir;
    } else {
      direction = !r->dir;
    }

    //enW = nW;

    enW = r->end - r->start + 1; //Partial results

    mask[0] = '\0';
    
    if (direction) {

      for (int rr=0; rr<r->num_mismatches; rr++) {
	concat_error_string(mask, mask_aux, r, rr, &enW);
      }

    } else {
 
      for (int rr=r->num_mismatches-1; rr>=0; rr--) {
	concat_error_string(mask, mask_aux, r, rr, &enW);
      }

    }

    //printf("%lu %lu -> %d\n", r->k, r->l, r->err_kind[0]);
    //printf("%d %d %d\n", r->start, r->pos, r->end);
    //printf("\n");

    for (size_t j=r->k; j<=r->l; j++) {

      if (S->ratio==1) {

      	if (direction)
      	  key = Si->siz - Si->vector[j] - enW - 1;
      	else
      	  key = S->vector[j];
	
      } else {

      	if (direction)
      	  key = Si->siz - getScompValue(j, Si, C, Oi) - enW -1;
      	else
      	  key = getScompValue(j, S, C, O);
 
      }

      index = binsearch(ex->offset, ex->size, key);

      if(key + enW <= ex->offset[index]) {
    	found = 1;
    	//printf("%lu\n", r_list->read_index);
	fprintf(fp, "read_%u %c %s %u %d%s\n", r_list->read_index, plusminus[type], ex->chromosome + (index-1)*IDMAX, ex->start[index-1] + (key - ex->offset[index-1]), r->num_mismatches, mask);
	//printf("read_%u %c %s %u %d%s\n", r_list->read_index, plusminus[type], ex->chromosome + (index-1)*IDMAX, ex->start[index-1] + (key - ex->offset[index-1]), r->num_mismatches, mask);
      }

    }

  }

  return found;

}
