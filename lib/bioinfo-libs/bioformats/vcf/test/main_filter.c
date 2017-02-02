#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include <containers/list.h>
#include <commons/log.h>

#include "vcf_file_structure.h"
#include "vcf_file.h"
#include "vcf_filters.h"
#include "vcf_read.h"
#include "vcf_reader.h"


int main(int argc, char *argv[]) {
    size_t max_batches = 20;
    size_t batch_size = 2000;
    list_t *read_list = (list_t*) malloc(sizeof(list_t));
    list_init("batches", 1, max_batches, read_list);

    int ret_code;
    double start, stop, total;
    vcf_file_t* file;

    char input[] = "1:1000000-2000000,22";
	filter_t *region_f = region_filter_new(input, 0, "http://localhost:8080", "hsa", "v1");
	filter_t *snp_f = snp_filter_new("exclude");
	filter_chain *chain = add_to_filter_chain(region_f, NULL);
	chain = add_to_filter_chain(snp_f, chain);
	
	int num_filters;
	filter_t **filter_list = sort_filter_chain(chain, &num_filters);

    init_log_custom(2, 1, NULL);

#pragma omp parallel sections private(start, stop, total) lastprivate(file)
    {
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d reads the VCF file\n", omp_get_thread_num());
            // Reading
            start = omp_get_wtime();

            file = vcf_open(argv[1]);
            ret_code = vcf_parse_batches(read_list, batch_size, file, 1);

            stop = omp_get_wtime();
            total = (stop - start);

            if (ret_code) { LOG_FATAL_F("[%dR] Error code = %d\n", omp_get_thread_num(), ret_code); }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            // Writing to a new file
            if (argc == 3) {
                start = omp_get_wtime();

                ret_code = vcf_write(file, argv[2]);

                stop = omp_get_wtime();
                total = (stop - start);

                if (ret_code) { LOG_ERROR_F("[%dW] Error code = %d\n", omp_get_thread_num(), ret_code); }

                LOG_INFO_F("[%dW] Time elapsed = %f s\n", omp_get_thread_num(), total);
                LOG_INFO_F("[%dW] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
            }

            list_decr_writers(read_list);

            vcf_close(file);

        }
#pragma omp section
        {
            int debug = 1;
            LOG_DEBUG_F("OMP num threads = %d\n", omp_get_num_threads());
            LOG_DEBUG_F("Thread %d prints info\n", omp_get_thread_num());

            start = omp_get_wtime();

            int i = 0;
            list_item_t* item = NULL;

            while ((item = list_remove_item(read_list)) != NULL) {
                vcf_batch_t *batch = (vcf_batch_t*) item->data_p;
                list_t *input_records = batch;
                list_t *passed_records;
                list_t *failed_records = (list_t*) malloc(sizeof(list_t));
                list_init("failed_records", 1, INT_MAX, failed_records);

                if (i % 200 == 0) {
                    LOG_DEBUG_F("Batch %d reached by thread %d - %zu/%zu records \n", 
                            i, omp_get_thread_num(),
                            batch->size, batch->capacity);
                }

                passed_records = run_filter_chain(input_records, failed_records, 
                                                  filter_list, num_filters);

                LOG_DEBUG_F("passed records = %zu, failed records = %zu\n", 
                        passed_records->length, failed_records->length);

//    vcf_batch_print(stdout, item->data_p);
                vcf_batch_free(item->data_p);
                list_item_free(item);
                i++;

                free(passed_records);
                free(failed_records);
            }

            stop = omp_get_wtime();

            total = (stop - start);

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

        }
    }

    free(filter_list);
    free(read_list);

    return 0;
}


