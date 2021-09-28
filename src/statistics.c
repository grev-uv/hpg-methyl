#include "statistics.h"


//------------------------------------------------------------------------------------

void basic_statistics_display(basic_statistics_t *statistics, int rna_mode, float alig_time, float load_time) {
  size_t total_reads = statistics->total_reads;
  size_t num_mapped_reads = statistics->num_mapped_reads;
  size_t total_mappings = statistics->total_mappings;

  size_t total_sp = statistics->total_sp;
  size_t uniq_sp = statistics->uniq_sp;

  printf("+--------------------------------------------------------------------------------------+\n");
  printf("|                                     GLOBAL STATISTICS                                |\n");
  printf("+--------------------------------------------------------------------------------------+\n");
  printf("| Loading Time (s)  : %-65.2f", load_time);
  printf("|\n");
  printf("| Alignment Time (s): %-65.2f", alig_time);
  printf("|\n");
  printf("| Total Time (s)    : %-65.2f", load_time + alig_time);
  printf("|\n");
  printf("========================================================================================\n");
  printf("| Total Reads Processed: %-62llu", total_reads);
  printf("|\n");
  printf("+-------------------------------------------+------------------------------------------+\n");
  printf("| Reads Mapped: %-18llu  %6.2f", num_mapped_reads, num_mapped_reads * 100.0 / total_reads);
  printf("% | ");
  printf(" Reads Unmapped: %-14llu  %6.2f", total_reads - num_mapped_reads, (total_reads - num_mapped_reads) * 100.0 / total_reads);
  printf("%  |\n");
  if (rna_mode) {
    printf("+-------------------------------------------+------------------------------------------+\n");
  } else {
    printf("+-------------------------------------------+------------------------------------------+\n");
  }

}

//------------------------------------------------------------------------------------------

basic_statistics_t *basic_statistics_new() {
  basic_statistics_t *basic = (basic_statistics_t *)malloc(sizeof(basic_statistics_t));
  pthread_mutex_init(&(basic->mutex), NULL);
  basic->total_reads = 0;
  basic->num_mapped_reads = 0;
  basic->total_mappings = 0;
  basic->total_sp = 0;
  basic->uniq_sp = 0;
  return basic;
}

//-------------------------------------------------------------------------------------------

void basic_statistics_add(size_t total_reads, size_t num_mapped_reads, size_t total_mappings, basic_statistics_t *basic) {
  pthread_mutex_lock(&basic->mutex);
  basic->total_reads += total_reads;
  basic->num_mapped_reads += num_mapped_reads;
  basic->total_mappings += total_mappings;
  pthread_mutex_unlock(&basic->mutex);
}
