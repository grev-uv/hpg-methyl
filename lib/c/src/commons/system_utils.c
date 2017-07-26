
#include <limits.h>

#include "system_utils.h"

/* **************************************************************
 *      Functions implementations     *
 * **************************************************************/

unsigned long int get_free_memory() {
    FILE* fd_memory = NULL;

    long int free_memory;
    char str_memory[15];
    size_t res;
    res = system("free | grep Mem | awk '{ print $4 }' > /tmp/free_memory");

    fd_memory = fopen("/tmp/free_memory", "r");

    if (fd_memory == NULL) {
        LOG_FATAL("Free memory cannot be obtained.\n");
    }

    char *res_f = fgets(str_memory, 15, fd_memory);
    if (!res_f) {
      return 0;
    }

    sscanf(str_memory, "%lu", &free_memory);

    LOG_DEBUG_F("free memory (KB): %li\n", free_memory);

    fclose(fd_memory);
    return free_memory;
}

unsigned long int get_estimated_memory_needed(int process, int batch_size, int max_list_length) {
    unsigned long int list_memory = max_list_length * batch_size;

    if (process == FASTQ_QC) {
        LOG_DEBUG_F("required memory (KB): %li\n", QC_MEMORY_USAGE_FACTOR * list_memory / 1000);
        return QC_MEMORY_USAGE_FACTOR * list_memory / 1000;
    } else if (process == FASTQ_PREPRO) {
        LOG_DEBUG_F("required memory (KB): %li\n", PREPRO_MEMORY_USAGE_FACTOR * list_memory / 1000);
        return PREPRO_MEMORY_USAGE_FACTOR * list_memory / 1000;
    } else { return 0; }
}

int get_max_estimated_alignments_by_chromosome(char* input_filename) {
    int max_estimated_alignments;
    struct stat st;

    stat(input_filename, &st);
    max_estimated_alignments = (int)(1.4 * (st.st_size / MEAN_COMPRESSED_ALIGNMENT_SIZE) / 12);  //security margin: 1.4, nts of chromosome 1: 1/12 *total nts

    return max_estimated_alignments;
}
/*
int get_optimal_cpu_num_threads() {
    FILE* fd_cpu_num_cores = NULL;

    int optimal_cpu_num_threads;
    char cpu_num_cores[2];
    char log_message[50];

    fd_cpu_num_cores = popen("cat /proc/cpuinfo | grep processor | wc -l", "r");

    if (fd_cpu_num_cores == NULL) {
        LOG_DEBUG("Num of CPU cores cannot be obtained.\n");
    }

    fgets(cpu_num_cores, 2, fd_cpu_num_cores);
    sscanf(cpu_num_cores, "%i", &optimal_cpu_num_threads);

    optimal_cpu_num_threads = max((optimal_cpu_num_threads - 2), 1);
    sprintf(log_message, "optimal_cpu_num_threads: %i\n", optimal_cpu_num_threads);
    LOG_DEBUG(log_message);

    pclose(fd_cpu_num_cores);

    return optimal_cpu_num_threads;

}
*/

size_t get_optimal_cpu_num_threads() {
  FILE* fd_cpu_num_cores = NULL;
  const unsigned int MAX_LINE = 1024;
  size_t optimal_cpu_num_threads = 0;

  char line[MAX_LINE];
  
  fd_cpu_num_cores = fopen("/proc/cpuinfo", "r");

  if (fd_cpu_num_cores == NULL) {
    printf("Num of CPU cores cannot be obtained.\n");
    return 0;
  }
  
  while (fgets(line, MAX_LINE, fd_cpu_num_cores)) {
    if (!strncmp("processor", line, 9)) {
      optimal_cpu_num_threads++;
    }
  }
  
  LOG_DEBUG_F("optimal_cpu_num_threads: %zu", optimal_cpu_num_threads);

  fclose(fd_cpu_num_cores);
  
  return optimal_cpu_num_threads;

}


int get_optimal_gpu_num_threads() {
    int optimal_gpu_num_threads = 0;

    #ifdef CUDA_VERSION
    optimal_gpu_num_threads = 16 * get_cuda_device_warp_size();
    #endif

    return optimal_gpu_num_threads;
}

int get_optimal_batch_size(int process, int max_list_length) {
    unsigned long int optimal_batch_size = 0;
    unsigned long int gpu_global_memory = ULLONG_MAX;

    #ifdef CUDA_VERSION
    gpu_global_memory = get_cuda_device_global_memory();
    #endif

    unsigned long int free_memory = get_free_memory();

    if (process == FASTQ_QC) {
        optimal_batch_size = (1000 * free_memory) / (QC_MEMORY_USAGE_FACTOR * max_list_length);
    } else if (process == FASTQ_PREPRO) {
        optimal_batch_size = (1000 * free_memory) / (PREPRO_MEMORY_USAGE_FACTOR * max_list_length);
    } else if (process == BAM_QC) {
        optimal_batch_size = 1000 * free_memory / BAM_BATCH_SIZE_TO_FREE_MEMORY_RATIO;
    }

    unsigned long int size = min(optimal_batch_size, (gpu_global_memory / 2));

    return (int) size;
}


