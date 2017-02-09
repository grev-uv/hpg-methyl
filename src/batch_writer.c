#include "batch_writer.h"


//------------------------------------------------------------------------------------

bam_header_t *create_bam_header_by_genome(genome_t *genome) {

  bam_header_t *bam_header = (bam_header_t *) calloc(1, sizeof(bam_header_t));
  int num_targets = genome->num_chromosomes;

  bam_header->n_targets = num_targets;
  bam_header->target_name = (char **) calloc(num_targets, sizeof(char *));
  bam_header->target_len = (uint32_t*) calloc(num_targets, sizeof(uint32_t));

  // Allocate free space to store the global methylation statistics
  const char *header_pg = "@PG\tID:HPG-Methyl\tVN:3.1\n";
  const char *header_co = "@CO\tZM:00\tC:00000000000\n";

  size_t pg_len = strlen(header_pg);
  size_t co_len = strlen(header_co);
  size_t header_tx_len = 1 + strlen(header_pg) + num_targets * (1 + strlen(header_co));

  bam_header->text = calloc(header_tx_len, sizeof(char));
  memcpy(bam_header->text, header_pg, pg_len);

  char tg_str[4];

  for (size_t i = 0, k = pg_len; i < num_targets; ++i, k += co_len) {
    memcpy(bam_header->text + k, header_co, co_len);
    snprintf(tg_str, 4, "%02lu", i);

    bam_header->text[k + 7] = tg_str[0];
    bam_header->text[k + 8] = tg_str[1];
  }

  bam_header->text[header_tx_len - 1] = '\n';

  for (int i = 0; i < num_targets; i++) {
    bam_header->target_name[i] = strdup(genome->chr_name[i]);
    bam_header->target_len[i] = genome->chr_size[i] + 1;
  }

  bam_header->l_text = header_tx_len;

  return bam_header;
}

//------------------------------------------------------------------------------------

void batch_writer_input_init(char* match_filename, char* splice_exact_filename, 
			     char* splice_extend_filename, 
			     linked_list_t* list_p, genome_t* genome, 
			     batch_writer_input_t* input_p) {

  input_p->match_filename = match_filename;
  input_p->splice_exact_filename = splice_exact_filename;
  input_p->splice_extend_filename = splice_extend_filename;
  input_p->list_p = list_p;
  input_p->genome = genome;

  // Internal
  input_p->bam_file = NULL;
  input_p->total_batches = 0;
  input_p->total_reads = 0;
  input_p->total_mappings = 0;
  input_p->num_mapped_reads = 0;
  input_p->limit_print = 10000;
}
