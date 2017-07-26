#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include <commons/log.h>
#include <containers/list.h>

#include "ped_batch.h"
#include "ped_file_structure.h"
#include "ped_file.h"
#include "ped_read.h"
#include "ped_reader.h"
#include "ped_write.h"


int main (int argc, char *argv[])
{
    size_t max_batches = 20;
    size_t batch_size = 2000;
    list_t *read_list = (list_t*) malloc (sizeof(list_t));
    list_init("batches", 1, max_batches, read_list);

    int ret_code;
    double start, stop, total;
    char *filename = (char*) malloc ((strlen(argv[1])+1) * sizeof(char));
    strncat(filename, argv[1], strlen(argv[1]));
    ped_file_t* file = ped_open(filename);
    
    init_log_custom(LOG_DEBUG_LEVEL, 1, NULL, "w");
    
#pragma omp parallel sections private(start, stop, total)
{
    #pragma omp section
    {
        LOG_DEBUG_F("Thread %d reads the PED file\n", omp_get_thread_num());
        // Reading
        start = omp_get_wtime();
        
        ret_code = ped_read_batches(read_list, batch_size, file);
        
        stop = omp_get_wtime();
        total = (stop - start);
        
        if (ret_code) {
            LOG_FATAL_F("[%dR] Error code = %d\n", omp_get_thread_num(), ret_code);
        }
        LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
        LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
        
        // Writing to a new file
        if (argc == 3) 
        {
            start = omp_get_wtime();
        
            ret_code = ped_write(file, argv[2]);
            
            stop = omp_get_wtime();
            total = (stop - start);
            
            if (ret_code) {
                LOG_ERROR_F("[%dW] Error code = %d\n", omp_get_thread_num(), ret_code);
            }
            LOG_INFO_F("[%dW] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dW] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
        }
        
        list_decr_writers(read_list);
        
    }
    #pragma omp section
    {
        LOG_DEBUG_F("OMP num threads = %d\n", omp_get_num_threads());
        LOG_DEBUG_F("Thread %d prints info\n", omp_get_thread_num());
        
        start = omp_get_wtime();
        
        int i = 0;
        int add_ret_code = 0;
        list_item_t* item = NULL;
        ped_batch_t *batch = NULL;
        list_item_t *batch_item = NULL;
        FILE *out = fopen("result.ped", "w");
        while ( (item = list_remove_item(read_list)) != NULL ) {
            batch = (ped_batch_t*) item->data_p;
            if (i % 200 == 0) 
            {
                int debug = 1;
                LOG_DEBUG_F("Batch %d reached by thread %d - %zu/%zu records \n", i, omp_get_thread_num(), 
                    ((ped_batch_t*) item->data_p)->length, ((ped_batch_t*) item->data_p)->max_length);
            }
            
            while ( (batch_item = list_remove_item_async(batch)) != NULL) {
                add_ret_code = add_ped_record(batch_item->data_p, file);
                if (add_ret_code > 0) {
                    LOG_ERROR_F("%s - %s\n", ((ped_record_t*) batch_item->data_p)->family_id, get_ped_semantic_error_msg(add_ret_code));
                }
            }
            
            ped_batch_free(item->data_p);
            list_item_free(item);
            i++;
        }
        
        ped_write_to_file(file, out);
        fclose(out);
        
        stop = omp_get_wtime();
        total = (stop - start);
        
        LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
        LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
    }
}

    ped_close(file, 1, 1);
    free(read_list);

    return 0;
}
