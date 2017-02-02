
#include "fastq_gzfile.h"


fastq_gzfile_t *fastq_gzopen(char *filename) {
	FILE *fd = fopen(filename, (char*)"r");
	char log_message[50];

	if (fd == NULL) {
		LOG_FATAL_F("Error opening file: %s\n", filename);
		return NULL;
	}

	fastq_gzfile_t* fq_gzfile = (fastq_gzfile_t*) malloc(sizeof(fastq_gzfile_t));

	fq_gzfile->filename = filename;
	fq_gzfile->fd = fd;

	fq_gzfile->strm.zalloc = Z_NULL;
	fq_gzfile->strm.zfree = Z_NULL;
	fq_gzfile->strm.opaque = Z_NULL;
	fq_gzfile->strm.avail_in = 0;
	fq_gzfile->strm.next_in = Z_NULL;
	fq_gzfile->ret = inflateInit2 (&fq_gzfile->strm, 15 + 32);    // Using inflateInit2 for GZIP support

	fq_gzfile->data = NULL;
	fq_gzfile->data_size = 0;

	return fq_gzfile;
}

void fastq_gzclose(fastq_gzfile_t *fq_gzfile) {
	fclose(fq_gzfile->fd);
	if(fq_gzfile->ret == Z_STREAM_END) {
		(void)inflateEnd(&fq_gzfile->strm);
	}
	free(fq_gzfile->data);
	free(fq_gzfile);
}



size_t fastq_gzread_se(array_list_t *reads, size_t num_reads, fastq_gzfile_t *fq_gzfile) {
	size_t count = 0;
	char header1[MAX_READ_ID_LENGTH];
	char sequence[MAX_READ_SEQUENCE_LENGTH];
	char header2[MAX_READ_ID_LENGTH];
	char qualities[MAX_READ_SEQUENCE_LENGTH];
	int header_length, sequence_length, quality_length;
	fastq_read_t *read;

	size_t num_lines_to_read = 4 * num_reads;	/* Each read consists of 4 lines */

	int max_data_len = CHUNK;
	int max_read_len = MAX_READ_SEQUENCE_LENGTH;	/* Each read is supposed to be shorter than MAX_READ_SEQUENCE_LENGTH */
	int eof_found = 0;
	int c = 0;
	int i = 0;
	//	fq_gzfile->i = 0;
	size_t lines = 0;
	char *aux;
	//	fq_gzfile->data = (char*) calloc (CHUNK, sizeof(char));
	char *data; // = (char*) calloc (CHUNK, sizeof(char));
	char *id = (char*) calloc (max_read_len, sizeof(char));
	char *seq = (char*) calloc (max_read_len, sizeof(char));
	char *qual = (char*) calloc (max_read_len, sizeof(char));

	// ZLIB variables
	unsigned have;
	unsigned char in[CHUNK];
	unsigned char out[CHUNK];


	// If there is some data from before calls
	if(fq_gzfile->data != NULL) {
		if(fq_gzfile->data_size > max_data_len) {
			data = (char*) calloc (fq_gzfile->data_size+max_data_len, sizeof(char));
			max_data_len = fq_gzfile->data_size+max_data_len;
		}else{
			data = (char*) calloc (max_data_len, sizeof(char));
		}
		strncpy(data, fq_gzfile->data, fq_gzfile->data_size);
		i = fq_gzfile->data_size;
	}else {
		// first time, no data has been saved before
		data = (char*) calloc (max_data_len, sizeof(char));
	}


	do {
		fq_gzfile->strm.avail_in = fread(in, 1, CHUNK, fq_gzfile->fd);
		//		printf("fq_gzfile->strm.avail_in: %i, CHUNK: %i\nnext_in: %s\n\n", fq_gzfile->strm.avail_in, CHUNK, fq_gzfile->strm.next_in);
		if (ferror(fq_gzfile->fd)) {
			(void)inflateEnd(&fq_gzfile->strm);
			return Z_ERRNO;
		}
		if (fq_gzfile->strm.avail_in == 0)
			break;
		fq_gzfile->strm.next_in = in;

		/* run inflate() on input until output buffer not full */
		do {
			fq_gzfile->strm.avail_out = CHUNK;
			fq_gzfile->strm.next_out = out;
			fq_gzfile->ret = inflate(&fq_gzfile->strm, Z_NO_FLUSH);
			assert(fq_gzfile->ret != Z_STREAM_ERROR);  /* state not clobbered */
			switch (fq_gzfile->ret) {
			case Z_NEED_DICT:
				fq_gzfile->ret = Z_DATA_ERROR;     /* and fall through */
			case Z_DATA_ERROR:
			case Z_MEM_ERROR:
				(void)inflateEnd(&fq_gzfile->strm);
				return fq_gzfile->ret;
			}
			have = CHUNK - fq_gzfile->strm.avail_out;
			for (int j = 0; j < have && !eof_found; j++) {
				c = out[j];

				if (c != EOF) {
					max_data_len = consume_input(c, &data, max_data_len, i);
					if (c == '\n') {
						lines++;
					}
					i++;
				} else {
					eof_found = 1;
				}
			}
		} while (fq_gzfile->strm.avail_out == 0);

		/* done when inflate() says it's done */
	} while (lines < num_lines_to_read && fq_gzfile->ret != Z_STREAM_END);

	//	printf("data: %s\n", data);
	//	LOG_DEBUG_F("lines: %i, num_lines_to_read: %i\n", lines, num_lines_to_read);

	// check if have read the expected number of lines
	size_t parsed_chars;
	size_t parsed_lines = 0;
	size_t data_size;
	//	if(lines > 0) { //= num_lines_to_read
	aux = data;
	for(parsed_chars = 0; parsed_chars < i && parsed_lines < num_lines_to_read; parsed_chars++) {
		if(data[parsed_chars] == '\n') {
//		printf(">>i: %i, parsed_chars: %i, %i, aux: %s\n", i, parsed_chars, data[i-1], aux);
			data[parsed_chars] = '\0';
			if(count % 4 == 0) {
				strcpy(id, aux);  //printf("%s\n", id);
			}
			if(count % 4 == 1) {
				strcpy(seq, aux);  //printf("%s\n", seq);
			}
			if(count % 4 == 2) {
			}
			if(count % 4 == 3) {
				strcpy(qual, aux);  //printf("%s\n", qual);
				read = fastq_read_new(id, seq, qual);
				array_list_insert(read, reads);
			}
			count++;
			aux = data + parsed_chars + 1;
			parsed_lines++;
		}
	}
	//		LOG_DEBUG_F("i: %lu, parsed_lines: %lu\n", i, parsed_lines);
	//		LOG_DEBUG_F("parsed_chars: %lu, parsed_lines: %lu\n", parsed_chars, parsed_lines);
	//		lines = 0;
	//		LOG_DEBUG_F("BEFORE memcpy: fq_gzfile->data_size: %lu, new size: %lu\n", fq_gzfile->data_size, data_size);
	data_size = i - parsed_chars;
	if(fq_gzfile->data == NULL) {
		fq_gzfile->data = (char*)malloc(data_size*sizeof(char));
	}
	if(fq_gzfile->data_size != 0 && fq_gzfile->data_size < data_size) {
		fq_gzfile->data = realloc(fq_gzfile->data, data_size);
	}
	if(data_size > 0) {
		memcpy(fq_gzfile->data, data+parsed_chars, data_size);
	}
	fq_gzfile->data_size = data_size;
	//	}

	free(data);
	free(id);
	free(seq);
	free(qual);

	//	if(fq_gzfile->ret == Z_STREAM_END) {
	//		(void)inflateEnd(&fq_gzfile->strm);
	//	}
	//		return fq_gzfile->ret == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
	//	printf(">>>>reads->size: %lu, num_reads: %lu\n", reads->size, num_reads);
	return reads->size;
}


size_t fastq_gzread_bytes_se(array_list_t *reads, size_t bytes_to_read, fastq_gzfile_t *fq_gzfile) {
	size_t count = 0;
	char header1[MAX_READ_ID_LENGTH];
	char sequence[MAX_READ_SEQUENCE_LENGTH];
	char header2[MAX_READ_ID_LENGTH];
	char qualities[MAX_READ_SEQUENCE_LENGTH];
	int header_length, sequence_length, quality_length;
	fastq_read_t *read;

//	size_t num_lines_to_read = bytes;	/* Each read consists of 4 lines */

	int max_data_len = CHUNK;
	int max_read_len = MAX_READ_SEQUENCE_LENGTH;	/* Each read is supposed to be shorter than MAX_READ_SEQUENCE_LENGTH */
	int eof_found = 0;
	int c = 0;
	int i = 0;
	size_t bytes_processed = 0;
	char *aux;

	char *data;
	char *id = (char*) calloc (max_read_len, sizeof(char));
	char *seq = (char*) calloc (max_read_len, sizeof(char));
	char *qual = (char*) calloc (max_read_len, sizeof(char));

	// ZLIB variables
	unsigned have;
	unsigned char in[CHUNK];
	unsigned char out[CHUNK];


	// If there is some data from before calls
	if(fq_gzfile->data != NULL) {
		if(fq_gzfile->data_size > max_data_len) {
			data = (char*) calloc (fq_gzfile->data_size+max_data_len, sizeof(char));
			max_data_len = fq_gzfile->data_size+max_data_len;
		}else{
			data = (char*) calloc (max_data_len, sizeof(char));
		}
		strncpy(data, fq_gzfile->data, fq_gzfile->data_size);
		i = fq_gzfile->data_size;
	}else {
		// first time, no data has been saved before
		data = (char*) calloc (max_data_len, sizeof(char));
	}


	do {
		fq_gzfile->strm.avail_in = fread(in, 1, CHUNK, fq_gzfile->fd);
		//		printf("fq_gzfile->strm.avail_in: %i, CHUNK: %i\nnext_in: %s\n\n", fq_gzfile->strm.avail_in, CHUNK, fq_gzfile->strm.next_in);
		if (ferror(fq_gzfile->fd)) {
			(void)inflateEnd(&fq_gzfile->strm);
			return Z_ERRNO;
		}
		if (fq_gzfile->strm.avail_in == 0)
			break;
		fq_gzfile->strm.next_in = in;

		/* run inflate() on input until output buffer not full */
		do {
			fq_gzfile->strm.avail_out = CHUNK;
			fq_gzfile->strm.next_out = out;
			fq_gzfile->ret = inflate(&fq_gzfile->strm, Z_NO_FLUSH);
			assert(fq_gzfile->ret != Z_STREAM_ERROR);  /* state not clobbered */
			switch (fq_gzfile->ret) {
			case Z_NEED_DICT:
				fq_gzfile->ret = Z_DATA_ERROR;     /* and fall through */
			case Z_DATA_ERROR:
			case Z_MEM_ERROR:
				(void)inflateEnd(&fq_gzfile->strm);
				return fq_gzfile->ret;
			}
			have = CHUNK - fq_gzfile->strm.avail_out;
			for (int j = 0; j < have && !eof_found; j++) {
				c = out[j];

				if (c != EOF) {
					max_data_len = consume_input(c, &data, max_data_len, i);
//					if (c == '\n') {
//						bytes_processed++;
//					}
					i++;
					bytes_processed++;
				} else {
					eof_found = 1;
				}
			}
		} while (fq_gzfile->strm.avail_out == 0);

		/* done when inflate() says it's done */
	} while (i < bytes_to_read && fq_gzfile->ret != Z_STREAM_END);

	// check if have read the expected number of lines
	size_t parsed_chars;
	size_t parsed_lines = 0;
	size_t data_size;
	aux = data;
	for(parsed_chars = 0; parsed_chars < i; parsed_chars++) {	//parsed_chars < bytes_to_read || parsed_lines % 4 == 0
		if(data[parsed_chars] == '\n') {
			data[parsed_chars] = '\0';
			if(count % 4 == 0) {
				strcpy(id, aux);  //printf("%s\n", id);
			}
			if(count % 4 == 1) {
				strcpy(seq, aux);  //printf("%s\n", seq);
			}
			if(count % 4 == 2) {
			}
			if(count % 4 == 3) {
				strcpy(qual, aux);  //printf("%s\n", qual);
				read = fastq_read_new(id, seq, qual);
				array_list_insert(read, reads);
				if(parsed_chars+1 > bytes_to_read) {
					parsed_chars++;
					break;
				}
			}
			count++;
			aux = data + parsed_chars + 1;
//			parsed_lines++;
		}
	}
	data_size = i - parsed_chars;
	if(fq_gzfile->data == NULL) {
		fq_gzfile->data = (char*)malloc(data_size*sizeof(char));
	}
	if(fq_gzfile->data_size != 0 && fq_gzfile->data_size < data_size) {
		fq_gzfile->data = realloc(fq_gzfile->data, data_size);
	}
	if(data_size > 0) {
		memcpy(fq_gzfile->data, data+parsed_chars, data_size);
	}
	fq_gzfile->data_size = data_size;

	free(data);
	free(id);
	free(seq);
	free(qual);

	return parsed_chars;
}

static size_t consume_input(int c, char **data, size_t max_len, int position_in_data) {
	assert(data);
	(*data)[position_in_data] = c;
	// Text too long to be stored in 'data', realloc
	if (position_in_data == max_len - 1) {
//		printf("reallocating data:  position_in_data: %i, max_len: %lu, data: %c\n", position_in_data, max_len, c);
		char *aux = realloc(*data, max_len * 2);
		if (aux) {
			*data = aux;
			return max_len * 2;
		} else {
			LOG_FATAL("Could not allocate enough memory for reading input VCF file\n");
		}
	}
	return max_len;
}
