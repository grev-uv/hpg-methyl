#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>

#include "commons/commons.h"
#include "commons/log.h"
#include "commons/system_utils.h"
#include "containers/list.h"

#include "fastq_batch_reader_omp.h"
#include "fastq_batch.h"
#include "fastq_file.h"

/* ******************************************************
 *      	Function implementations    		*
 * ******************************************************/

fastq_batch_reader_t* fastq_batch_reader_new(char* filename, int source_id, list_t* list_p, size_t batch_size, int batch_list_max_length) {
    fastq_batch_reader_t* fastq_batch_reader_p = (fastq_batch_reader_t*) malloc(sizeof(fastq_batch_reader_t));

    // open the input file
    fastq_batch_reader_p->fastq_file_p = fastq_fopen_mode(filename, "r");
    fastq_batch_reader_p->source_id = source_id;

    fastq_batch_reader_p->batch_size = batch_size;
    fastq_batch_reader_p->batch_list_max_length = batch_list_max_length;

    fastq_batch_reader_p->eof = 0;
    pthread_mutex_init(&(fastq_batch_reader_p->eof_lock), NULL);

    fastq_batch_reader_p->alive = 1;
    pthread_mutex_init(&(fastq_batch_reader_p->alive_lock), NULL);

    fastq_batch_reader_p->batch_list_p = list_p;
    list_incr_writers(list_p);

    return fastq_batch_reader_p;
}

void fastq_batch_reader_free(fastq_batch_reader_t* fastq_batch_reader_p) {
    // close the input file and exit
    fastq_fclose(fastq_batch_reader_p->fastq_file_p);
}

void* fastq_batch_reader_thread_function(void* param_p) {
    unsigned int total_reads = 0;
    char log_message[100];

    // Cast param_p to 'fastq_batch_reader_t'
    fastq_batch_reader_t* fastq_batch_reader_p = (fastq_batch_reader_t*) param_p;

    sprintf(log_message, "Thread-READ: START, for file %s\n", fastq_batch_reader_p->fastq_file_p->filename);
    LOG_DEBUG(log_message);

    int num_reads = 0, num_batchs = 0;
    fastq_batch_t* fastq_batch_p;
    list_item_t *item_p = NULL;

    while (1) {
        // allocationg memory for the current fastq batch
        fastq_batch_p = (fastq_batch_t*) calloc(1, sizeof(fastq_batch_t));

        fastq_batch_init(fastq_batch_p, fastq_batch_reader_p->batch_size);
        fastq_batch_p->source_id = fastq_batch_reader_p->source_id;

        while ((list_get_length(fastq_batch_reader_p->batch_list_p) >= fastq_batch_reader_p->batch_list_max_length) ||
                (fastq_batch_reader_p->batch_list_p->max_length  > (fastq_batch_reader_p->batch_list_max_length / 2)))  {
            LOG_DEBUG("Thread-READ: go to sleep for a while...\n");

            // Delay for a bit
            sched_yield();
            sleep(1);
            break;
        }

        // read reads from file
        num_reads = fastq_read_batch_max_size(fastq_batch_p, fastq_batch_reader_p->batch_size, fastq_batch_reader_p->fastq_file_p);
        total_reads += num_reads;

        // if there is no read exit the loop
        if (num_reads == 0) {
            fastq_batch_free(fastq_batch_p);
            free(fastq_batch_p);
            break;
        }

        // otherwise, create a new batch object..
        item_p = list_item_new(num_batchs, 0, fastq_batch_p);

        // insert this batch to the list (synchronization is built in the list, see list.c)
        list_insert_item(item_p, fastq_batch_reader_p->batch_list_p);

        sprintf(log_message, "Thread-READ: ....reading batch %i (%i reads), done !!!!\n", num_batchs, num_reads);
        LOG_DEBUG(log_message);

        num_batchs++;
    } // end of batch loop

    list_print((list_t*)fastq_batch_reader_p->batch_list_p);

    list_decr_writers(fastq_batch_reader_p->batch_list_p);
    pthread_exit((void*)total_reads);
}
