#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <limits.h>

#include <containers/list.h>
#include <commons/log.h>

#include "vcf_file_structure.h"
#include "vcf_file.h"
#include "vcf_util.h"


extern int mmap_vcf;

int main (int argc, char *argv[])
{
    if (argc < 2) {
        printf("Usage: main_batch <VCF file> [output-file] \n");
    }
    
    size_t max_batches = 20;
    size_t batch_size = 500;
//     list_t *read_list = (list_t*) malloc (sizeof(list_t));
//     list_init("batches", 1, max_batches, read_list);

    int ret_code;
    double start, stop, total;
    vcf_file_t* file = vcf_open(argv[1], INT_MAX);
    FILE *out_file;
    if (argc == 3) {
        out_file = fopen(argv[2], "w");
    }

    init_log_custom(LOG_DEBUG_LEVEL, 2, 1, NULL);
    
//     if (argc > 2 && strcmp(argv[2], "mmap-vcf") == 0) {
//         mmap_vcf = 1;
//     }
    
#pragma omp parallel sections private(start, stop, total)
{
    #pragma omp section
    {
        LOG_DEBUG_F("Thread %d reads the VCF file\n", omp_get_thread_num());
        // Reading
        start = omp_get_wtime();
        
        ret_code = vcf_parse_batches(batch_size, file);
        notify_end_reading(file);
        
        stop = omp_get_wtime();
        total = (stop - start);
        
        if (ret_code) { LOG_FATAL_F("[%dR] Error code = %d\n", omp_get_thread_num(), ret_code); }
        LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
        LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
    }
    
    #pragma omp section
    {
        LOG_DEBUG_F("OMP num threads = %d\n", omp_get_num_threads());
        LOG_DEBUG_F("Thread %d prints info\n", omp_get_thread_num());
        
        start = omp_get_wtime();
        
        int header_written = 0;
        
        int i = 0;
        vcf_batch_t *batch;
        while ( (batch = fetch_vcf_batch(file)) != NULL ) {
//             if (i % 20 == 0) 
//             {
                printf("Batch %d reached by thread %d - %zu/%zu records \n", i, omp_get_thread_num(), 
                    batch->records->size, batch->records->capacity);
//             }
            
            if (out_file) {
                // Writing to a new file
                if (!header_written) {
                    write_vcf_header(file, out_file);
                    header_written = 1;
                }
                    
//                 vcf_batch_print(stdout, item->data_p);
                write_vcf_batch(batch, out_file);
            }
            vcf_batch_free(batch);
            
            i++;
        }
        
        stop = omp_get_wtime();
        total = (stop - start);
        
        LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
        LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
    }
}

//     free(read_list);
    vcf_close(file);

    return 0;
}
