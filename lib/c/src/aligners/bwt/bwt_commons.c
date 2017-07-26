#include "bwt_commons.h"

size_t cur_alloc=0, max_alloc=0;

void *mymalloc(size_t n)
{
  void *p;

  p = malloc(n);
  if (p == NULL) {
    printf("malloc failed.\n");
    exit(1);
  }
  cur_alloc += n;
  if (cur_alloc > max_alloc) {
    max_alloc = cur_alloc;
    //printf("allocated %ld\n",max_alloc);
  }

  if (n == 0) {
    printf("warning: 0 bytes allocated p=%p\n",p);
  }

  return p;
}

void *myrealloc(void *ptr, size_t next, size_t last)
{
  void *p;

  p = realloc(ptr, next);
  if (next > 0 && p == NULL) {
    printf("realloc failed. ptr=%p new=%zd old=%zd\n",ptr,next,last);
    exit(1);
  }
  cur_alloc += next - last;
  if (cur_alloc > max_alloc) {
    max_alloc = cur_alloc;
    //printf("allocated %ld\n",max_alloc);
  }

  return p;
}

void report_mem(const char *s)
{
  puts(s);
  printf("allocated total %zu   max %zu\n",cur_alloc,max_alloc);
}


void myfree(void *p, size_t s)
{
  free(p);
  cur_alloc -= s;
}

//-----------------------------------------------------------------------------

void bwt_init_replace_table(bwt_config_t *config, char *nucleotides) {
  //printf("****Bases = %s\n", config->nucleotides);
  
  if (nucleotides == NULL) {
    fprintf(stderr, "%s -> Nucleotides NULL\n", __func__);
    exit(EXIT_FAILURE);
  }

  config->nA          = strlen(nucleotides);
  config->table       = (uint8_t *) malloc(128 * sizeof(uint8_t));
  config->nucleotides = (char *) malloc(config->nA * sizeof(char));
  config->reverse     = (char *) malloc(config->nA * sizeof(char));

  for (int i = 0; i < config->nA; i++) {

    config->nucleotides[i] = toupper(nucleotides[i]);

    config->table[toupper(nucleotides[i])] = i;
    config->table[tolower(nucleotides[i])] = i;

    if      (toupper(nucleotides[i]) == 'A') config->AA = i;
    else if (toupper(nucleotides[i]) == 'C') config->CC = i;
    else if (toupper(nucleotides[i]) == 'G') config->GG = i;
    else if (toupper(nucleotides[i]) == 'T') config->TT = i;

  }

  config->reverse[config->AA] = config->TT;
  config->reverse[config->CC] = config->GG;
  config->reverse[config->GG] = config->CC;
  config->reverse[config->TT] = config->AA;

}
void bwt_free_replace_table(bwt_config_t *config) {
  free(config->table);
  free(config->nucleotides);
  free(config->reverse);
}

void encode_bases(uint8_t* dest, char* src, uintmax_t length, uint8_t *table) {
  for (uintmax_t i=0; i<length; i++)
    dest[i] = table[(uintmax_t)src[i]];
}

void decode_bases(char* dest, uint8_t* src, uintmax_t length, char *rev_table) {
  for (uintmax_t i=0; i<length; i++)
    dest[i] = rev_table[(uintmax_t)src[i]];
  dest[length] = '\0';
}

void revstring(uint8_t *X, uintmax_t nX) {

  char tmp;
  uintmax_t i, j;

  if (nX <= 1) return;

  for (i=0, j=nX-1; i<=j; i++, j--) {
    tmp = X[i];
    X[i] = X[j];
    X[j] = tmp;
  }

}

void revstrand(uint8_t *X, uintmax_t nX, char *reverse) {

  char tmp;
  uintmax_t i, j;

  if (nX <= 1) return;

  for (i=0, j=nX-1; i<=j; i++, j--) {
    tmp = reverse[X[i]];
    X[i] = reverse[X[j]];
    X[j] = tmp;
  }

}

void duplicate_reverse(uint8_t *X, uintmax_t nX, char *reverse) {

  uintmax_t i, j;

  if (nX == 0) return;

  for (i=0, j=nX*2-1; i<nX; i++, j--) {
    X[j] = reverse[X[i]];
  }

}
