#include "fastq_batch.h"


#define MAX_HEADER_LENGTH 128
#define MIN_READ_SIZE (2*50)   // 70 = 2 * 35
                               // size is splitted in sequence and quality, so divided 2
                               // minimum expected sequence length is 35

/* ******************************************************
 *    		Function implementations  		*
 * ******************************************************/

fastq_batch_t* fastq_batch_new(unsigned long size) {
    fastq_batch_t* fastq_batch_p = (fastq_batch_t*) calloc(1, sizeof(fastq_batch_t));
    fastq_batch_p->num_reads = 0;
    fastq_batch_p->data_size = size;
    fastq_batch_p->data_indices_size = (sizeof(int) * fastq_batch_p->data_size) / MIN_READ_SIZE;

    fastq_batch_p->header_indices = (int*) calloc(fastq_batch_p->data_indices_size, 1);
    fastq_batch_p->header = (char*) calloc(MAX_HEADER_LENGTH * fastq_batch_p->data_size / MIN_READ_SIZE, sizeof(char));

    fastq_batch_p->data_indices = (int*) calloc(fastq_batch_p->data_indices_size, 1);

    fastq_batch_p->seq = (char*) calloc(fastq_batch_p->data_size / 2, sizeof(char));
    fastq_batch_p->quality = (char*) calloc(fastq_batch_p->data_size / 2, sizeof(char));

    return fastq_batch_p;
}
/*
fastq_batch_t* fastq_batch_init(fastq_batch_t* fastq_batch_p, unsigned long size) {
    fastq_batch_p->num_reads = 0;
    fastq_batch_p->data_size = size;
    fastq_batch_p->data_indices_size = sizeof(int) * fastq_batch_p->data_size / MIN_READ_SIZE;

    fastq_batch_p->header_indices = (int*) calloc(fastq_batch_p->data_indices_size, 1);
    fastq_batch_p->header = (char*) calloc(MAX_HEADER_LENGTH * fastq_batch_p->data_size / MIN_READ_SIZE, sizeof(char));

    fastq_batch_p->data_indices = (int*) calloc(fastq_batch_p->data_indices_size, 1);

    fastq_batch_p->seq = (char*) calloc(fastq_batch_p->data_size / 2, sizeof(char));
    fastq_batch_p->quality = (char*) calloc(fastq_batch_p->data_size / 2, sizeof(char));

    return fastq_batch_p;
}*/

void fastq_batch_free(fastq_batch_t* fastq_batch_p) {
    if (fastq_batch_p == NULL) {
        return;
    }

    if (fastq_batch_p->header_indices != NULL) {
        free(fastq_batch_p->header_indices);
        fastq_batch_p->header_indices = NULL;
    }
    if (fastq_batch_p->header != NULL) {
        free(fastq_batch_p->header);
        fastq_batch_p->header = NULL;
    }
    if (fastq_batch_p->data_indices != NULL) {
        free(fastq_batch_p->data_indices);
        fastq_batch_p->data_indices = NULL;
    }

    if (fastq_batch_p->seq != NULL) {
        free(fastq_batch_p->seq);
        fastq_batch_p->seq = NULL;
    }
    if (fastq_batch_p->quality != NULL) {
        free(fastq_batch_p->quality);
        fastq_batch_p->quality = NULL;
    }

    free(fastq_batch_p);
    fastq_batch_p = NULL;
}

void fastq_batch_print(fastq_batch_t* fastq_batch_p, FILE* fd) {
    for (int i = 0; i < fastq_batch_p->num_reads; i++) {
        fprintf(fd, "%s\n", &(fastq_batch_p->header[fastq_batch_p->header_indices[i]]));
        fprintf(fd, "%s\n", &(fastq_batch_p->seq[fastq_batch_p->data_indices[i]]));
        fprintf(fd, "+\n");
        fprintf(fd, "%s\n", &(fastq_batch_p->quality[fastq_batch_p->data_indices[i]]));
    }
}

void fprintf_read(FILE* fd, fastq_batch_t* batch_p, int index) {
    fprintf(fd, "%s\n", &(batch_p->header[batch_p->header_indices[index]]));
    fprintf(fd, "%s\n", &(batch_p->seq[batch_p->data_indices[index]]));
    fprintf(fd, "+\n");
    fprintf(fd, "%s\n", &(batch_p->quality[batch_p->data_indices[index]]));
}

void fprintf_rtrim_read(FILE* fd, fastq_batch_t* batch_p, int index, int rtrim_length) {
    fprintf(fd, "%s\n", &(batch_p->header[batch_p->header_indices[index]]));
    fprintf(fd, "%s\n", rtrim(&(batch_p->seq[batch_p->data_indices[index]]), rtrim_length));
    fprintf(fd, "+\n");
    fprintf(fd, "%s\n", rtrim(&(batch_p->quality[batch_p->data_indices[index]]), rtrim_length));
}

void fprintf_ltrim_read(FILE* fd, fastq_batch_t* batch_p, int index, int ltrim_length) {
    fprintf(fd, "%s\n", &(batch_p->header[batch_p->header_indices[index]]));
    fprintf(fd, "%s\n", ltrim(&(batch_p->seq[batch_p->data_indices[index]]), ltrim_length));
    fprintf(fd, "+\n");
    fprintf(fd, "%s\n", ltrim(&(batch_p->quality[batch_p->data_indices[index]]), ltrim_length));
}

void fprintf_trim_read(FILE* fd, fastq_batch_t* batch_p, int index, int rtrim_length, int ltrim_length) {
    char* str;
    fprintf(fd, "%s\n", &(batch_p->header[batch_p->header_indices[index]]));

    str = &(batch_p->seq[batch_p->data_indices[index]]);

    if (rtrim_length > 0) {
        str = (char*) rtrim(str, rtrim_length);
    }
    if (ltrim_length > 0) {
        str = (char*) ltrim(str, ltrim_length);
    }

    fprintf(fd, "%s\n", str);
    fprintf(fd, "+\n");
    str = &(batch_p->quality[batch_p->data_indices[index]]);

    if (rtrim_length > 0) {
        str = (char*) rtrim(str, rtrim_length);
    }

    if (ltrim_length > 0) {
        str = (char*) ltrim(str, ltrim_length);
    }

    fprintf(fd, "%s\n", str);
}
