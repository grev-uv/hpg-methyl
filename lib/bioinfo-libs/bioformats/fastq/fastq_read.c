#include <string.h>
#include <stdlib.h>

#include "fastq_read.h"


fastq_read_t *fastq_read_new(char *id, char *sequence, char *quality) {
	fastq_read_t *fq_read = (fastq_read_t*) calloc(1, sizeof(fastq_read_t));

	//size_t id_length = strlen(id); // + 1;
	size_t seq_length = strlen(sequence); // + 1;
	
	fq_read->id = strdup(id);
	//printf("%s\n", sequence);
	fq_read->sequence = strdup(sequence);
	//fq_read->sequence = strndup(sequence, seq_length + 1);
	//printf("%s\n", fq_read->sequence);
	fq_read->quality = strdup(quality);
	fq_read->length = seq_length;

	return fq_read;
}

fastq_read_t *fastq_read_dup(fastq_read_t *fq) {
  fastq_read_t *fq_read = (fastq_read_t*) calloc(1, sizeof(fastq_read_t));
  fq_read->id = strdup(fq->id);
  fq_read->sequence = strdup(fq->sequence);
  fq_read->quality = strdup(fq->quality);
  fq_read->length = fq->length;

  return fq_read;
}

fastq_read_pe_t *fastq_read_pe_new(char *id1, char *id2, char *sequence1, char *quality1, char *sequence2, char *quality2, int mode) {
	fastq_read_pe_t *fq_read_pe = (fastq_read_pe_t *) calloc(1, sizeof(fastq_read_pe_t));

	size_t id_length1 = strlen(id1); // + 1;
	size_t id_length2 = strlen(id2); // + 1;
	size_t seq_length1 = strlen(sequence1); // + 1;
	size_t seq_length2 = strlen(sequence2); // + 1;

	fq_read_pe->id1 = strndup(id1, id_length1);
	fq_read_pe->id2 = strndup(id2, id_length2);
	fq_read_pe->sequence1 = strndup(sequence1, seq_length1);
	fq_read_pe->sequence2 = strndup(sequence2, seq_length2);
	fq_read_pe->quality1 = strndup(quality1, seq_length1);
	fq_read_pe->quality2 = strndup(quality2, seq_length2);
	fq_read_pe->length1 = seq_length1;
	fq_read_pe->length2 = seq_length2;
	fq_read_pe->mode = mode;

	return fq_read_pe;
}



void fastq_read_free(fastq_read_t *fq_read) {
	if (fq_read == NULL) return;

	if (fq_read->id != NULL) free(fq_read->id);
	if (fq_read->sequence != NULL) free(fq_read->sequence);
	if (fq_read->quality != NULL) free(fq_read->quality);

	free(fq_read);
}

void fastq_read_pe_free(fastq_read_pe_t *fq_read_pe) {
	if (fq_read_pe == NULL) return;

	if (fq_read_pe->id1 != NULL) free(fq_read_pe->id1);
	if (fq_read_pe->id2 != NULL) free(fq_read_pe->id2);
	if (fq_read_pe->sequence1 != NULL) free(fq_read_pe->sequence1);
	if (fq_read_pe->sequence2 != NULL) free(fq_read_pe->sequence2);
	if (fq_read_pe->quality1 != NULL) free(fq_read_pe->quality1);
	if (fq_read_pe->quality2 != NULL) free(fq_read_pe->quality2);

	free(fq_read_pe);
}



void fastq_read_print(fastq_read_t *read) {
	printf("%s\n", read->id);
	printf("%s\n", read->sequence);
	printf("+\n");
	printf("%s\n", read->quality);
}

void fastq_read_pe_print(fastq_read_pe_t *fq_read_pe) {
	printf("%s\n", fq_read_pe->id1);
	printf("%s\n", fq_read_pe->sequence1);
	printf("+\n");
	printf("%s\n", fq_read_pe->quality1);
	printf("%s\n", fq_read_pe->id2);
	printf("%s\n", fq_read_pe->sequence2);
	printf("+\n");
	printf("%s\n", fq_read_pe->quality2);
}

