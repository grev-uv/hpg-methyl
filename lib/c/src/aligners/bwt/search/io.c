#include "io.h"

exome *exome_new() {
  exome *ex = (exome *)malloc(sizeof(exome));
  
  ex->max_chromosomes = INDEX_EXOME;
  ex->size = 0;

  ex->chromosome = (char *)malloc(sizeof(char)*INDEX_EXOME*IDMAX);
  ex->start      = (uintmax_t *)malloc(sizeof(uintmax_t)*INDEX_EXOME);
  ex->end        = (uintmax_t *)malloc(sizeof(uintmax_t)*INDEX_EXOME);
  ex->offset     = (uintmax_t *)malloc(sizeof(uintmax_t)*INDEX_EXOME);

  return ex;

}

void exome_free(exome *ex) {
  free(ex->chromosome);
  free(ex->start);
  free(ex->end);
  free(ex->offset);
  free(ex);
}

void load_exome_file(exome *ex, const char *directory) {

  FILE *fp;

  char path[500];
  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/index");

  fp  = fopen(path, "r");
  check_file_open(fp, path);

  char line[MAXLINE];

  ex->offset[0]=0;
  ex->size = 0;

  while (fgets(line, MAXLINE, fp) ) {

    if (line[0]=='>') {

      int j;
      char c = 0;

      for(j=0; j<IDMAX-1; j++) {
	c = line[j+1];
	if (c==' ') break;
	ex->chromosome[ex->size*IDMAX+j] = c;
      }

      ex->chromosome[ex->size*IDMAX+j] = '\0';

      sscanf(line + j + 2, "%ju %ju %*s", &ex->start[ex->size], &ex->end[ex->size]);
      ex->size++;

      if (ex->size >= ex->max_chromosomes) {
	ex->max_chromosomes *= 2;
	ex->chromosome = (char *)realloc(ex->chromosome, 
					 sizeof(char)*ex->max_chromosomes*IDMAX);
	ex->start      = (uintmax_t *)realloc(ex->start, 
					      sizeof(uintmax_t)*ex->max_chromosomes);
	ex->end        = (uintmax_t *)realloc(ex->end, 
					      sizeof(uintmax_t)*ex->max_chromosomes);
	ex->offset     = (uintmax_t *)realloc(ex->offset, 
					      sizeof(uintmax_t)*ex->max_chromosomes); 
      }

      ex->offset[ex->size] = ex->offset[ex->size-1] + (ex->end[ex->size-1] - ex->start[ex->size-1]+1);

    }

  }

  fclose(fp);

}

void save_exome_file(exome *ex, bool reverse, const char *directory) {

  FILE *fp;

  char path[500];
  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/index");

  fp  = fopen(path, "w");
  check_file_open(fp, path);

  for(uintmax_t i=0; i<ex->size; i++) {
    fprintf(fp, ">%s %ju %ju\n", ex->chromosome + i*IDMAX, (uintmax_t) ex->start[i], (uintmax_t) ex->end[i]);
  }

  if(reverse) {
    for(int i=ex->size-1; i>=0; i--) {
      fprintf(fp, ">%s %ju %ju\n", ex->chromosome + i*IDMAX, (uintmax_t) ex->start[i], (uintmax_t) ex->end[i]);
    }
  }

}

void encode_reference(ref_vector *X, exome *ex, const char *ref_path, bwt_config_t *bwt_config) {

  FILE *ref_file;
  ref_file = fopen(ref_path, "r");
  check_file_open(ref_file, ref_path);

  size_t read, size;

  fseek(ref_file, 0, SEEK_END);
  read = ftell(ref_file); //Valgrind errors on dbwt
  fseek(ref_file, 0, SEEK_SET);

  if (bwt_config->duplicate_strand) size = read*2 + 1;
  else         size = read   + 1;

  X->vector = (uint8_t *) malloc( size * sizeof(uint8_t) );
  check_malloc(X->vector, ref_path);

  char *reference = (char *) X->vector;

  if (ex !=NULL) ex->size=0;

  uintmax_t partial_length=0, total_length=0;

  while ( fgets(reference + total_length, MAXLINE, ref_file) ) {

    if ( (reference + total_length)[0] == '>') {

      if (ex!=NULL) {

	if (total_length == 0) {

	  sscanf(reference + total_length, ">%s ", ex->chromosome + ex->size * IDMAX);
	  ex->start[ex->size] = 0;

	} else {

	  ex->end[ex->size] = partial_length - 1;
	  partial_length=0;

	  if (ex->size == 0)
	    ex->offset[0] = 0;
	  else
	    ex->offset[ex->size] = ex->offset[ex->size-1] + (ex->end[ex->size-1] - ex->start[ex->size-1] + 1);
	  ex->size++;
	  if (ex->size >= ex->max_chromosomes) {
	    ex->max_chromosomes *= 2;
	    ex->chromosome = (char *)realloc(ex->chromosome, 
					     sizeof(char)*ex->max_chromosomes*IDMAX);
	    ex->start      = (uintmax_t *)realloc(ex->start, 
						  sizeof(uintmax_t)*ex->max_chromosomes);
	    ex->end        = (uintmax_t *)realloc(ex->end, 
						  sizeof(uintmax_t)*ex->max_chromosomes);
	    ex->offset     = (uintmax_t *)realloc(ex->offset, 
						  sizeof(uintmax_t)*ex->max_chromosomes); 
	  }

	  sscanf(reference + total_length, ">%s ", ex->chromosome + ex->size * IDMAX);
	  ex->start[ex->size] = 0;

	}

      }

      continue;

    }

    size_t length = strlen(reference + total_length);
    if ((reference + total_length)[length-1]=='\n')
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

  encode_bases(X->vector, reference, total_length, bwt_config->table);

  if (bwt_config->duplicate_strand) {
    duplicate_reverse(X->vector, total_length, bwt_config->reverse);
    X->n = total_length * 2;
  } else {
    X->n = total_length;
  }

  X->dollar = 0;
  X->vector[X->n] = 0; //Valgrind errors on dbwt

  fclose(ref_file);

}

bool nextFASTAToken(FILE *queries_file, char *uncoded, uint8_t *coded, uintmax_t *nquery, bwt_config_t *bwt_config) {

	char line[MAXLINE];
  size_t length=0;

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

    encode_bases(coded, uncoded, *nquery, bwt_config->table);

    return true;

  } else {

    return false;

  }

}

void read_ref_vector(ref_vector *vector, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".vec");

  fp  = fopen(path,  "rb+");
  check_file_open(fp, path);

  err = fread(&vector->n, sizeof(uint64_t),  1, fp);
  check_file_read(err, 1, path);

  err = fread(&vector->dollar, sizeof(uint64_t),  1, fp);
  check_file_read(err, 1, path);

  check_file_read(err, 1, path);
  vector->vector = (uint8_t *) malloc((vector->n + 1) * sizeof(uint8_t)); //Valgrind errors on dbwt
  check_malloc(vector->vector, path);

  err = fread(vector->vector, sizeof(uint8_t), vector->n, fp);
  check_file_read(err, vector->n, path);

  vector->vector[vector->n] = 0; //Valgrind errors on dbwt

  fclose(fp);

}

void save_ref_vector(ref_vector *vector, const char *directory, const char *name) {

  size_t err=0;
  FILE *fp;

  char path[500];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/");
  strcat(path, name);
  strcat(path, ".vec");

  fp  = fopen(path,  "wb+");
  check_file_open(fp, path);

  err = fwrite(&vector->n,      sizeof(uint64_t), 1, fp);
  check_file_write(err, 1, path);

  err = fwrite(&vector->dollar, sizeof(uint64_t), 1, fp);
  check_file_write(err, 1, path);

  err = fwrite(vector->vector, sizeof(uint8_t), vector->n, fp);
  check_file_write(err, vector->n, path);

  fclose(fp);

}

//==============================================================================

void save_config(char *nucleotides, bool duplicate_strand, const char *directory) {
  FILE *fp;

  char path[500];
  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/config.txt");

  fp  = fopen(path,  "w");
  check_file_open(fp, path);

  fputs(nucleotides, fp);
  fputc('\n', fp);
  fputs((duplicate_strand)?"1":"0", fp);
  fputc('\n', fp);

  fclose(fp);
}

void read_config(char *nucleotides, bool *duplicate_strand, const char *directory) {

  if (nucleotides == NULL) {
    fprintf(stderr, "%s -> Nucleotides NULL\n", __func__);
    exit(EXIT_FAILURE);
  }

  FILE *fp;
  char path[500];
  char ds_str[32];

  path[0]='\0';
  strcat(path, directory);
  strcat(path, "/config.txt");

  fp  = fopen(path,  "r");
  check_file_open(fp, path);

  char *res = fgets(nucleotides, 128, fp);
  nucleotides[strlen(nucleotides) - 1] = '\0';

  res = fgets(ds_str, 1, fp);
  *duplicate_strand = atoi(ds_str);

  fclose(fp);

}

//==============================================================================
