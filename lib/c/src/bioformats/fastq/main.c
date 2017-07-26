/*
 * main.c
 *
 *  Created on: Sep 28, 2012
 *      Author: imedina
 */

#include "fastq_file.h"
#include "fastq_gzfile.h"
#include "containers/array_list.h"
#include "fastq_filter.h"

int main (int argc, char *argv[]) {

	if(!strcmp("count-lines", argv[1])) {
		fastq_file_t *file = fastq_fopen(argv[2]);
		array_list_t *reads = array_list_new(2000000, 1.8, COLLECTION_MODE_SYNCHRONIZED);
		size_t nread = 1;
		int count = 0;
		while((nread = fastq_fread_se(reads, 100000, file)) != 0) {
			count += nread;
			for(int i=0; i<reads->size; i++) {
				fastq_read_print(array_list_get(i, reads));
			}
			//			printf("Size: %i, Capacity: %i\n", reads->size, reads->capacity);
			array_list_clear(reads, (void *)fastq_read_free);
		}
		//		printf("Total num reads: %i\n", reads->size);
		//		fastq_read_print(array_list_get(0, reads));
		//		fastq_read_print(array_list_get(reads->size-1, reads));
		array_list_free(reads, (void *)fastq_read_free);
		fastq_fclose(file);
	}

	if(!strcmp("count-lines-pe", argv[1])) {
		fastq_file_t *file1 = fastq_fopen(argv[2]);
		fastq_file_t *file2 = fastq_fopen(argv[3]);
		array_list_t *reads = array_list_new(1000000, 1.8, COLLECTION_MODE_SYNCHRONIZED);
		size_t nread = 1;
		int count = 0;
		while((nread = fastq_fread_pe(reads, 100000, file1, file2, FASTQ_FILE_PAIRED_END_MODE)) != 0) {
			count += nread;
			for(int i=0; i<reads->size; i++) {
				fastq_read_pe_print(array_list_get(i, reads));
			}
			//			printf("Size: %i, Capacity: %i\n", reads->size, reads->capacity);
			array_list_clear(reads, (void *)fastq_read_pe_free);
		}
		//		printf("Total num reads: %i\n", reads->size);
		//		fastq_read_print(array_list_get(0, reads));
		//		fastq_read_print(array_list_get(reads->size-1, reads));
		array_list_free(reads, (void *)fastq_read_pe_free);
		fastq_fclose(file1);
		fastq_fclose(file2);
	}

	if(!strcmp("count-bytes", argv[1])) {
		fastq_file_t *file = fastq_fopen(argv[2]);
		//			printf("=>%i\n", file->ret);
		array_list_t *reads = array_list_new(100000, 1.8, COLLECTION_MODE_SYNCHRONIZED);
		size_t nread = 1;
		int count = 0;
		while((nread = fastq_fread_bytes_se(reads, 100000, file)) != 0) {
			//				nread = fastq_gzread_bytes_se(reads, 100000, file);
			count += reads->size;
			//				printf("Size: %i, Capacity: %i, count = %i, nread: %i\n", reads->size, reads->capacity, count, nread);
			for(int i=0; i<reads->size; i++) {
				fastq_read_print(array_list_get(i, reads));
			}
			//				fastq_read_print(array_list_get(reads->size-1, reads));
			array_list_clear(reads, (void *)fastq_read_free);
		}
		//			printf("Total num reads: %i\n", count);
		//		fastq_read_print(array_list_get(0, reads));
		array_list_free(reads, (void *)fastq_read_free);
		fastq_fclose(file);
	}

	if(!strcmp("count-bytes-pe", argv[1])) {
		fastq_file_t *file1 = fastq_fopen(argv[2]);
		fastq_file_t *file2 = fastq_fopen(argv[3]);
		//			printf("=>%i\n", file->ret);
		array_list_t *reads = array_list_new(100000, 1.8, COLLECTION_MODE_SYNCHRONIZED);
		size_t nread = 1;
		int count = 0;
		while((nread = fastq_fread_bytes_pe(reads, 100000, file1, file2, 1)) != 0) {
			//				nread = fastq_gzread_bytes_se(reads, 100000, file);
			count += reads->size;
			//				printf("Size: %i, Capacity: %i, count = %i, nread: %i\n", reads->size, reads->capacity, count, nread);
			for(int i=0; i<reads->size; i++) {
				fastq_read_pe_print(array_list_get(i, reads));
			}
			//				fastq_read_print(array_list_get(reads->size-1, reads));
			array_list_clear(reads, (void *)fastq_read_pe_free);
		}
		//			printf("Total num reads: %i\n", count);
		//		fastq_read_print(array_list_get(0, reads));
		array_list_free(reads, (void *)fastq_read_pe_free);
		fastq_fclose(file1);
		fastq_fclose(file2);
	}

	if(!strcmp("count-lines-gz", argv[1])) {
		fastq_gzfile_t *file = fastq_gzopen(argv[2]);
		//		printf("=>%i\n", file->ret);
		array_list_t *reads = array_list_new(1000000, 1.8, COLLECTION_MODE_SYNCHRONIZED);
		size_t nread = 1;
		int count = 0;
		while((nread = fastq_gzread_se(reads, 100000, file)) != 0) {
			//			nread = fastq_gzread_se(reads, 1000000, file);
			count += nread;
			//			printf("Size: %i, Capacity: %i, count = %i, nread: %i\n", reads->size, reads->capacity, count, nread);
			for(int i=0; i<reads->size; i++) {
				fastq_read_print(array_list_get(i, reads));
			}
			//			fastq_read_print((fastq_read_t*)array_list_get(reads->size-1, reads));
			array_list_clear(reads, (void *)fastq_read_free);
		}
		//		printf("Total num reads: %i\n", count);
		//		fastq_read_print(array_list_get(0, reads));
		array_list_free(reads, (void *)fastq_read_free);
		fastq_gzclose(file);
	}

	if(!strcmp("count-bytes-gz", argv[1])) {
		fastq_gzfile_t *file = fastq_gzopen(argv[2]);
		//			printf("=>%i\n", file->ret);
		array_list_t *reads = array_list_new(1000000, 1.8, COLLECTION_MODE_SYNCHRONIZED);
		size_t nread = 1;
		int count = 0;
		while((nread = fastq_gzread_bytes_se(reads, 10000000, file)) != 0) {
			//				nread = fastq_gzread_bytes_se(reads, 100000, file);
			count += reads->size;
			//				printf("Size: %i, Capacity: %i, count = %i, nread: %i\n", reads->size, reads->capacity, count, nread);
			for(int i=0; i<reads->size; i++) {
				fastq_read_print(array_list_get(i, reads));
			}
			//				fastq_read_print(array_list_get(reads->size-1, reads));
			array_list_clear(reads, (void *)fastq_read_free);
		}
		//			printf("Total num reads: %i\n", count);
		//		fastq_read_print(array_list_get(0, reads));
		array_list_free(reads, (void *)fastq_read_free);
		fastq_gzclose(file);
	}

	if(!strcmp("filter", argv[1])) {
		fastq_file_t *file = fastq_fopen(argv[2]);
		fastq_filter_options_t *fastq_filter_options = fastq_filter_options_new(50, 600, 10, 200, 5, 0, 0, 0, 0, 0, 0, 100);
		array_list_t *reads = array_list_new(200000, 1.8, COLLECTION_MODE_SYNCHRONIZED);
		array_list_t *passed_reads = array_list_new(200000, 1.8, COLLECTION_MODE_SYNCHRONIZED);
		array_list_t *failed_reads = array_list_new(200000, 1.8, COLLECTION_MODE_SYNCHRONIZED);
		size_t nread = 1;
		int count = 0;
		while((nread = fastq_fread_se(reads, 1000000, file)) != 0) {
			count += reads->size;
			//			passed_reads = array_list_new(200000, 1.8, COLLECTION_MODE_SYNCHRONIZED);
			//			failed_reads = array_list_new(200000, 1.8, COLLECTION_MODE_SYNCHRONIZED);
			//			for(int i=0; i<reads->size; i++) {
			//				fastq_read_print(array_list_get(i, reads));
			//			}
			fastq_filter(reads, passed_reads, failed_reads, fastq_filter_options);
			//fastq_read_print(array_list_get(0, passed_reads));
			for(int i=0; i<passed_reads->size; i++) {
				fastq_read_print(array_list_get(i, passed_reads));
			}
			//fastq_read_print(array_list_get(0, failed_reads));
			//printf("Total Reads: %lu, Passed Reads: %lu, Reads failed: %lu\n", reads->size, passed_reads->size, failed_reads->size);
			array_list_clear(reads, (void *)fastq_read_free);
			array_list_clear(passed_reads, (void *)NULL);
			array_list_clear(failed_reads, (void *)NULL);
			//			fastq_read_print(array_list_get(0, passed_reads));
			//			fastq_read_print(array_list_get(0, failed_reads));
			//			printf("Total Reads: %lu, Passed Reads: %lu, Reads filter: %lu\n", reads->size, passed_reads->size, failed_reads->size);
		}
		//		fastq_read_print(array_list_get(0, passed_reads));
		//		fastq_read_print(array_list_get(0, failed_reads));
		//		printf("Total Reads: %lu, Passed Reads: %lu, Reads filter: %lu\n", reads->size, passed_reads->size, failed_reads->size);

		fastq_filter_options_free(fastq_filter_options);
		array_list_free(reads, (void *)fastq_read_free);
		array_list_free(passed_reads, (void *)fastq_read_free);
		array_list_free(failed_reads, (void *)fastq_read_free);
		fastq_fclose(file);
	}

	return 0;
}

