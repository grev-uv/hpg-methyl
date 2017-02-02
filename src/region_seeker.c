#include "region_seeker.h"


void region_seeker_input_init(list_t *unmapped_read_list, cal_optarg_t *cal_optarg, 
			      bwt_optarg_t *bwt_optarg, bwt_index_t *bwt_index, 
			      list_t* region_list, unsigned int region_threads, 
			      unsigned int gpu_enable, int padding_left, int padding_right,
			      genome_t *genome, metaexons_t *metaexons,
			      region_seeker_input_t *input_p) {

  input_p->unmapped_read_list_p = unmapped_read_list;
  input_p->cal_optarg_p = cal_optarg;
  input_p->bwt_optarg_p = bwt_optarg;
  input_p->bwt_index_p = bwt_index;
  input_p->region_list_p = region_list;
  input_p->region_threads = region_threads;
  input_p->gpu_enable = gpu_enable;
  input_p->padding_left = padding_left;
  input_p->padding_right = padding_right;
  input_p->genome = genome;
  input_p->metaexons = metaexons;
}
