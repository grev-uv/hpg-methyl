#ifndef FASTQ_READ_H
#define FASTQ_READ_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "containers/array_list.h"

#define ID_MAX_BYTES		(1 * 1024 * 1024)
#define INDEX_MAX		100000
#define SEQUENCE_MAX_BYTES	(2 * 1024 * 1024)

/* **************************************
 *  		Structures		*
 * *************************************/

/**
* @brief Fastq read structure 
* 
* Structure for storing a fastq read
*/
typedef struct fastq_read {
    char *id;			/**< Id of the read. */
    char *sequence;		/**< Sequence of nts. */
    char *quality;		/**< Qualities. */

    int length;
} fastq_read_t;

typedef struct fastq_read_pe {
    char *id1;			/**< Id of the read. */
    char *id2;
    char *sequence1;		/**< Sequence of nts. */
    char *sequence2;
    char *quality1;		/**< Qualities. */
    char *quality2;

    int length1;
    int length2;
    int mode;
} fastq_read_pe_t;


/* **************************************
 *  		Functions		*
 * *************************************/

fastq_read_t *fastq_read_dup(fastq_read_t *fq);

fastq_read_t *fastq_read_new(char *id, char *sequence, char *quality);

fastq_read_pe_t *fastq_read_pe_new(char *id1, char *id2, char *sequence1, char *quality1, char *sequence2, char *quality2, int mode);


void fastq_read_free(fastq_read_t *fq_read);

void fastq_read_pe_free(fastq_read_pe_t *fq_read);


void fastq_read_print(fastq_read_t *read);

void fastq_read_pe_print(fastq_read_pe_t *read);


//float fastq_quality_average(fastq_read_t *fq_read_t);

void fastq_nt_quality_average(fastq_read_t fq_read_t);
void fastq_nt_count(fastq_read_t fq_read_t);
void fastq_length(fastq_read_t fq_read_t);
void fastq_remove_barcode(fastq_read_t fq_read_t, char *barcode);

#endif	/*  FASTQ_READ_H  */
