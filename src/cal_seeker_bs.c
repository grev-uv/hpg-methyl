#include "cal_seeker.h"

//------------------------------------------------------------------------------------
// functions to handle gaps between mapped seeds
//    - fill_gaps
//    - merge_seed_regions
//------------------------------------------------------------------------------------

void display_sr_lists_bs(char *msg, mapping_batch_t *mapping_batch, int bs_id) {
  fastq_read_t *read1;
  array_list_t *fq_batch1;

  size_t read_index;

  cal_t *cal;
  array_list_t *cal_list = NULL;
  array_list_t **mapping_lists;
  size_t num_cals, num_targets;
  size_t *targets;

  if (bs_id == 0) {
    mapping_lists = mapping_batch->mapping_lists;
    targets = mapping_batch->targets;
    num_targets = mapping_batch->num_targets;
    fq_batch1 = mapping_batch->GA_fq_batch;
  } else {
    mapping_lists = mapping_batch->mapping_lists2;
    targets = mapping_batch->targets2;
    num_targets = mapping_batch->num_targets2;
    fq_batch1 = mapping_batch->CT_fq_batch;
  }

  seed_region_t *s;
  linked_list_iterator_t *itr;

  LOG_DEBUG_F("%s\n", msg);

  // debugging....
  for (size_t i = 0; i < num_targets; i++) {
    read_index = targets[i];
    read1 = (fastq_read_t *)array_list_get(read_index, fq_batch1);

    LOG_DEBUG_F("Read %s\n", read1->id);

    cal_list = mapping_lists[read_index];
    num_cals = array_list_size(cal_list);

    if (num_cals <= 0) {
      continue;
    }

    // processing each CAL from this read
    for (size_t j = 0; j < num_cals; j++) {
      // get cal and read index
      cal = array_list_get(j, cal_list);
      LOG_DEBUG_F("\tCAL #%i of %i (strand %i, chr %i:%lu-%lu), sr_list size = %i\n",
		  						j, num_cals, cal->strand, cal->chromosome_id - 1, cal->start, cal->end, 
									cal->sr_list->size);
			itr = linked_list_iterator_new(cal->sr_list);

      s = (seed_region_t *)linked_list_iterator_curr(itr);

      while (s != NULL) {
				LOG_DEBUG_F("\t\t%s %x (dist. %i)\t[%i|%i - %i|%i]\n",
		    			(s->info ? new_cigar_code_string((cigar_code_t *)s->info) : ">>>>>> gap"),
		    			 s->info, (s->info ? ((cigar_code_t *)s->info)->distance : -1),
		    			 s->genome_start, s->read_start, s->read_end, s->genome_end);

				linked_list_iterator_next(itr);
				s = linked_list_iterator_curr(itr);
      }

      linked_list_iterator_free(itr);
    }
  }
}

//------------------------------------------------------------------------------------

void fill_gaps_bs(mapping_batch_t *mapping_batch, sw_optarg_t *sw_optarg,
		  genome_t *genome1, genome_t *genome2, int min_gap, int min_distance,
		  int bs_id, sw_optarg_t *sw_optarg1, sw_optarg_t *sw_optarg2,
			apply_sw_bs_stage_workspace_t *workspace) {
	array_list_t **mapping_lists;
	size_t num_targets;
	size_t *targets;
	array_list_t *fq_batch = mapping_batch->fq_batch;
	array_list_t *fq_batch1;
	array_list_t *fq_batch2;

	if (bs_id == 0) {
		mapping_lists = mapping_batch->mapping_lists;
		num_targets = mapping_batch->num_targets;
		targets = mapping_batch->targets;
		fq_batch1 = mapping_batch->GA_fq_batch;
		fq_batch2 = mapping_batch->GA_rev_fq_batch;
	} else {
		mapping_lists = mapping_batch->mapping_lists2;
		num_targets = mapping_batch->num_targets2;
		targets = mapping_batch->targets2;
		fq_batch1 = mapping_batch->CT_fq_batch;
		fq_batch2 = mapping_batch->CT_rev_fq_batch;
	}

	int sw_count = 0;
	int sw_count1 = 0;
	int sw_count2 = 0;

	fastq_read_t *read1;
	fastq_read_t *read2;

	fastq_read_t *read;
	genome_t *genome;

	size_t read_index, read_len;

	cal_t *cal;
	array_list_t *cal_list = NULL;
	size_t num_cals;

	char *revcomp_seq = NULL;

	seed_region_t *s, *prev_s, *new_s;
	linked_list_iterator_t *itr;

	cigar_code_t *cigar_code;

	size_t start, end;
	size_t gap_read_start, gap_read_end, gap_read_len;
	size_t gap_genome_start, gap_genome_end, gap_genome_len;

	int left_flank, right_flank;
	sw_prepare_t *sw_prepare;
	array_list_t *sw_prepare_list;

	if (workspace->fill_gaps_sw_prepare_list1 == NULL) {
		workspace->fill_gaps_sw_prepare_list1 = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
	} else {
		array_list_clear(workspace->fill_gaps_sw_prepare_list1, NULL);
	}

	if (workspace->fill_gaps_sw_prepare_list2 == NULL) {
		workspace->fill_gaps_sw_prepare_list2 = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
	} else {
		array_list_clear(workspace->fill_gaps_sw_prepare_list2, NULL);
	}

	array_list_t *sw_prepare_list1 = workspace->fill_gaps_sw_prepare_list1;
	array_list_t *sw_prepare_list2 = workspace->fill_gaps_sw_prepare_list2;

	char *query, *ref;
	char c1, c2;
	int distance, first, last;

	// initialize query and reference sequences to Smith-Waterman
	for (size_t i = 0; i < num_targets; i++) {
		read_index = targets[i];

		read = (fastq_read_t *)array_list_get(read_index, fq_batch);
		read1 = (fastq_read_t *)array_list_get(read_index, fq_batch1);
		read2 = (fastq_read_t *)array_list_get(read_index, fq_batch2);




		cal_list = mapping_lists[read_index];
		num_cals = array_list_size(cal_list);

		if (num_cals <= 0) {
			continue;
		}

		read_len = read->length;
		min_distance = read_len * 0.2;

		// processing each CAL from this read
		for (size_t j = 0; j < num_cals; j++) {
			// get cal and read index
			cal = array_list_get(j, cal_list);

			read = (cal->strand == 0) ? read1 : read2;
			genome = (cal->strand == 0) ? genome1 : genome2;

			if (cal->strand == 0 && bs_id == 0 || cal->strand == 1 && bs_id == 1) {
				c1 = 'C';
				c2 = 'T';
			} else {
				c1 = 'G';
				c2 = 'A';
			}

			prev_s = NULL;
			itr = linked_list_iterator_new(cal->sr_list);
			s = (seed_region_t *)linked_list_iterator_curr(itr);

			while (s != NULL) {
				// set the cigar for the current region
				gap_read_len = s->read_end - s->read_start + 1;
				cigar_code = cigar_code_new();
				cigar_code_add_op(cigar_op_new(gap_read_len, 'M'), cigar_code);
				
				s->info = (void *)cigar_code;
				cigar_code = NULL;
				sw_prepare = NULL;

				if ((prev_s == NULL && s->read_start != 0) || (prev_s != NULL)) {
					distance = 0;
					mapping_batch->num_gaps++;

					if (prev_s == NULL) {
						// gap at the first position
						gap_read_start = 0;
						gap_read_end = s->read_start - 1;

						gap_genome_start = s->genome_start - s->read_start;
						gap_genome_end = s->genome_start - 1;

						gap_read_len = gap_read_end - gap_read_start + 1;
						gap_genome_len = gap_genome_end - gap_genome_start + 1;

						if ((int)gap_genome_start >0)	//cal->start must be higher than 1
								cal->start = gap_genome_start;
						else
							;//RICARDO

						assert(gap_read_len != 0);
						assert(gap_genome_len != 0);

						if (gap_read_len > min_gap) {
							// the gap is too big, may be there's another CAL to cover it
							cigar_code = cigar_code_new();
							cigar_code_add_op(cigar_op_new(gap_read_len, 'H'), cigar_code);
						} else {
							left_flank = 0;
							right_flank = DOUBLE_FLANK;
						}
					} else {
						assert(prev_s->read_end < s->read_start);

						// gap at a middle position
						gap_read_start = prev_s->read_end + 1;
						gap_read_end = s->read_start - 1;

						gap_genome_start = prev_s->genome_end + 1;
						gap_genome_end = s->genome_start - 1;

						gap_read_len = gap_read_end - gap_read_start + 1;
						gap_genome_len = gap_genome_end - gap_genome_start + 1;

						assert(gap_genome_len != 0);

						if (gap_read_len == 0) {
							cigar_code = (cigar_code_t *)prev_s->info;

							cigar_code_append_op(cigar_op_new(gap_genome_len, 'D'), cigar_code);
							cigar_code->distance += gap_genome_len;

							cigar_code_append_op(cigar_op_new(s->read_end - s->read_start + 1, 'M'), cigar_code);
							cigar_code->distance += ((cigar_code_t *)s->info)->distance;

							prev_s->read_end = s->read_end;
							prev_s->genome_end = s->genome_end;

							// continue loop...
							linked_list_iterator_remove(itr);
							s = linked_list_iterator_curr(itr);
							continue;
						}

						left_flank = SINGLE_FLANK;
						right_flank = SINGLE_FLANK;
					}

					if (!cigar_code) {
						// we have to try to fill this gap and get a cigar
						if (gap_read_len == gap_genome_len) {
							// 1) first, for from  begin -> end, and begin <- end
							start = gap_genome_start; // + 1;
							end = gap_genome_end;     // + 1;
							first = -1;
							last = -1;
							ref = (char *)malloc((gap_genome_len + 5) * sizeof(char));

							genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1,
																								&start, &end, genome);

							query = &read->sequence[gap_read_start];

							for (int k = 0; k < gap_read_len; k++) {
								if (query[k] != ref[k]) {
									if (!(query[k] == c2 && ref[k] == c1)) {
										distance++;
									}

									if (first == -1) {
										first = k;
									}

									last = k;
								}
							}

							// free memory
							free(ref);

							if (distance < min_distance) {
								cigar_code = cigar_code_new();
								cigar_code_add_op(cigar_op_new(gap_read_len, 'M'), cigar_code);
								cigar_code_inc_distance(distance, cigar_code);
							}
						}

						if (!cigar_code) {
							// 2) second, prepare SW to run
							// get query sequence, revcomp if necessary
							size_t read_start = gap_read_start - left_flank;
							size_t read_end = gap_read_end + right_flank;
							int gap_read_len_ex = read_end - read_start + 1;

							query = (char *)malloc((gap_read_len_ex + 1) * sizeof(char));
							memcpy(query, &read->sequence[read_start], gap_read_len_ex);
							query[gap_read_len_ex] = '\0';

							// get ref. sequence
							size_t genome_start = gap_genome_start - left_flank; // + 1;
							size_t genome_end = gap_genome_end + right_flank;    // + 1;
							int gap_genome_len_ex = genome_end - genome_start + 1;

							ref = (char *)malloc((gap_genome_len_ex + 1) * sizeof(char));

							genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1,
							&genome_start, &genome_end, genome);

							ref[gap_genome_len_ex] = '\0';

							if (prev_s == NULL) {
								sw_prepare = sw_prepare_new(query, ref, left_flank, right_flank, FIRST_SW);
							} else {
								sw_prepare = sw_prepare_new(query, ref, left_flank, right_flank, MIDDLE_SW);
							}

							if (cal->strand == 0) {
								array_list_insert(sw_prepare, sw_prepare_list1);

								// increase counter
								sw_count1++;
							} else {
								array_list_insert(sw_prepare, sw_prepare_list2);

								// increase counter
								sw_count2++;
							}
						}
					}

					// insert gap in the list
					new_s = seed_region_new(gap_read_start, gap_read_end, gap_genome_start, gap_genome_end, 0);
					new_s->info = (void *)cigar_code;
					linked_list_iterator_insert(new_s, itr);

					if (sw_prepare) {
						sw_prepare->seed_region = new_s;
						sw_prepare->cal = cal;
						sw_prepare->read = read;
					}
				}

				// continue loop...
				prev_s = s;
				linked_list_iterator_next(itr);
				s = linked_list_iterator_curr(itr);
			}

			// check for a gap at the last position
			sw_prepare = NULL;

			if (prev_s != NULL && prev_s->read_end < read_len - 1) {
				cigar_code = NULL;
				mapping_batch->num_gaps++;

				// gap at the last position
				gap_read_start = prev_s->read_end + 1;
				gap_read_end = read_len - 1;
				gap_read_len = gap_read_end - gap_read_start + 1;

				assert(gap_read_len != 0);

				gap_genome_len = gap_read_len;
				gap_genome_start = prev_s->genome_end + 1;
				gap_genome_end = gap_genome_start + gap_genome_len - 1;

				cal->end = gap_genome_end;

				assert(gap_genome_len != 0);

				if (gap_read_len > min_gap) {
					// the gap is too big, may be there's another CAL to cover it
					cigar_code = cigar_code_new();
					cigar_code_add_op(cigar_op_new(gap_read_len, 'H'), cigar_code);
				} else {
					// we have to try to fill this gap and get a cigar
					//    1) first, for from  begin -> end, and begin <- end
					start = gap_genome_start; // + 1;
					end = gap_genome_end;     // + 1;
					first = -1;
					last = -1;
					ref = (char *)malloc((gap_genome_len + 1) * sizeof(char));

					genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1,
					&start, &end, genome);
					query = &read->sequence[gap_read_start];
					distance = 0;

					for (int k = 0; k < gap_read_len; k++) {
						if (query[k] != ref[k]) {
							if (!(query[k] == c2 && ref[k] == c1)) {
								distance++;
							}

							if (first == -1) {
								first = k;
							}

							last = k;
						}
					}

					// free memory
					free(ref);

					if (distance < min_distance) {
						cigar_code = cigar_code_new();
						cigar_code_add_op(cigar_op_new(gap_read_len, 'M'), cigar_code);
						cigar_code_inc_distance(distance, cigar_code);
					} else {
						//    2) second, prepare SW to run
						left_flank = DOUBLE_FLANK;
						right_flank = 0;

						// get query sequence, revcomp if necessary
						size_t read_start = gap_read_start - left_flank;
						size_t read_end = gap_read_end + right_flank;
						int gap_read_len_ex = read_end - read_start + 1;
						query = (char *)malloc((gap_read_len_ex + 1) * sizeof(char));

						// get ref. sequence
						size_t genome_start = gap_genome_start - left_flank; // + 1;
						size_t genome_end = gap_genome_end + right_flank;    // + 1;
						int gap_genome_len_ex = genome_end - genome_start + 1;
						ref = (char *)malloc((gap_genome_len_ex + 1) * sizeof(char));

						genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1,
																							&genome_start, &genome_end, genome);

						memcpy(query, &read->sequence[read_start], gap_read_len_ex);
						ref[gap_genome_len_ex] = '\0';
						query[gap_read_len_ex] = '\0';

						sw_prepare = sw_prepare_new(query, ref, left_flank, right_flank, LAST_SW);

						if (cal->strand == 0) {
							array_list_insert(sw_prepare, sw_prepare_list1);
							sw_count1++;
						} else {
							array_list_insert(sw_prepare, sw_prepare_list2);
							sw_count2++;
						}
					}
				}

				// insert gap in the list
				new_s = seed_region_new(gap_read_start, gap_read_end, gap_genome_start, gap_genome_end, 0);
				new_s->info = (void *)cigar_code;
				linked_list_insert_last(new_s, cal->sr_list);

				if (sw_prepare) {
					sw_prepare->seed_region = new_s;
					sw_prepare->cal = cal;
					sw_prepare->read = read;
				}
			}

			linked_list_iterator_free(itr);
		}

		// free memory
		if (revcomp_seq) {
			free(revcomp_seq);
			revcomp_seq = NULL;
		}
	}

	assert(sw_count1 == array_list_size(sw_prepare_list1));
	assert(sw_count2 == array_list_size(sw_prepare_list2));

	char *q[sw_count1], *r[sw_count1];

	for (int i = 0; i < sw_count1; i++) {
		sw_prepare = array_list_get(i, sw_prepare_list1);
		q[i] = sw_prepare->query;
		r[i] = sw_prepare->ref;
	}

	sw_multi_output_t *output1 = sw_multi_output_new(sw_count1);
	char *q2[sw_count2], *r2[sw_count2];

	for (int i = 0; i < sw_count2; i++) {
		sw_prepare = array_list_get(i, sw_prepare_list2);
		q2[i] = sw_prepare->query;
		r2[i] = sw_prepare->ref;
	}

	sw_multi_output_t *output2 = sw_multi_output_new(sw_count2);
	sw_multi_output_t *output;

	// run Smith-Waterman
	smith_waterman_mqmr(q, r, sw_count1, sw_optarg1, 1, output1);
	smith_waterman_mqmr(q2, r2, sw_count2, sw_optarg2, 1, output2);

	cigar_op_t *cigar_op;
	cigar_code_t *cigar_c;

	for (int type = 0; type < 2; type++){
		if (type == 0) {
			output = output1;
			sw_prepare_list = sw_prepare_list1;
			sw_count = sw_count1;
		} else {
			output = output2;
			sw_prepare_list = sw_prepare_list2;
			sw_count = sw_count2;
		}

		for (int i = 0; i < sw_count; i++) {
			sw_prepare = array_list_get(i, sw_prepare_list);
			s = sw_prepare->seed_region;

			int read_gap_len = s->read_end - s->read_start + 1;
			int genome_gap_len = s->genome_end - s->genome_start + 1;

			int read_gap_len_ex = read_gap_len_ex + sw_prepare->left_flank + sw_prepare->right_flank;
			int genome_gap_len_ex = genome_gap_len_ex + sw_prepare->left_flank + sw_prepare->right_flank;

			cigar_code_t *cigar_c = generate_cigar_code(output->query_map_p[i], output->ref_map_p[i],
																			strlen(output->query_map_p[i]), output->query_start_p[i],
																			output->ref_start_p[i], read_gap_len, genome_gap_len,
																			&distance, sw_prepare->ref_type);

			cigar_op = cigar_code_get_op(0, cigar_c);

			if (cigar_op) {
				if (cigar_op->name == 'H') {
					if (output->ref_start_p[i] == 0) { 
						cigar_op->name = 'I';
					} else {
						cigar_op->name = 'M';
					}
				} else if (cigar_op->name == '=') {
					cigar_op->name = 'M';
				}
			}

			cigar_op = cigar_code_get_last_op(cigar_c);

			if (cigar_op && cigar_op->name == 'H') {
				cigar_op->name = 'I';
			}

			// and now set the cigar for this gap
			s->info = (void *) cigar_c;

			// free
			sw_prepare_free(sw_prepare);
		}
	}

	// free memory
	sw_multi_output_free(output1);
	sw_multi_output_free(output2);
}

//------------------------------------------------------------------------------------

void fill_end_gaps_bs(mapping_batch_t *mapping_batch, sw_optarg_t *sw_optarg,
		      genome_t *genome1, genome_t *genome2, int min_H, int min_distance,
		      int bs_id, apply_sw_bs_stage_workspace_t *workspace) {
	array_list_t **mapping_lists;
	size_t num_targets;
	size_t *targets;
	array_list_t *fq_batch;
	array_list_t *fq_batch1;
	array_list_t *fq_batch2;

	fastq_read_t *read1, *read2;
	fastq_read_t *read;
	genome_t *genome;

	if (bs_id == 0) {
		mapping_lists = mapping_batch->mapping_lists;
		num_targets = mapping_batch->num_targets;
		targets = mapping_batch->targets;
		fq_batch1 = mapping_batch->GA_fq_batch;
		fq_batch2 = mapping_batch->GA_rev_fq_batch;
	} else {
		mapping_lists = mapping_batch->mapping_lists2;
		num_targets = mapping_batch->num_targets2;
		targets = mapping_batch->targets2;
		fq_batch1 = mapping_batch->CT_fq_batch;
		fq_batch2 = mapping_batch->CT_rev_fq_batch;
	}

	int sw_count, sw_count1 = 0, sw_count2 = 0;
	size_t read_index, read_len;

	cal_t *cal;
	array_list_t *cal_list = NULL;
	size_t num_cals;

	char *seq, *revcomp_seq = NULL;

	seed_region_t *s;

	cigar_op_t *cigar_op;
	cigar_code_t *cigar_code;

	size_t start, end;
	size_t gap_read_start, gap_read_end, gap_read_len;
	size_t gap_genome_start, gap_genome_end, gap_genome_len;

	int first, last, mode, distance, flank = 5;
	sw_prepare_t *sw_prepare;
	array_list_t *sw_prepare_list;

	if (workspace->fill_end_gaps_sw_prepare_list1 == NULL) {
		workspace->fill_end_gaps_sw_prepare_list1 = array_list_new(200, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
	} else {
		array_list_clear(workspace->fill_end_gaps_sw_prepare_list1, NULL);
	}

	if (workspace->fill_end_gaps_sw_prepare_list2 == NULL) {
		workspace->fill_end_gaps_sw_prepare_list2 = array_list_new(200, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
	} else {
		array_list_clear(workspace->fill_end_gaps_sw_prepare_list2, NULL);
	}

	array_list_t *sw_prepare_list1 = workspace->fill_end_gaps_sw_prepare_list1;
	array_list_t *sw_prepare_list2 = workspace->fill_end_gaps_sw_prepare_list2;

	char *ref, *query;

	float Ncg;
	float margen = mapping_batch->margin;
	char c1, c2;

	// initialize query and reference sequences to Smith-Waterman
	for (size_t i = 0; i < num_targets; i++) {
		read_index = targets[i];

		read1 = (fastq_read_t *)array_list_get(read_index, fq_batch1);
		read2 = (fastq_read_t *)array_list_get(read_index, fq_batch2);

		cal_list = mapping_lists[read_index];
		num_cals = array_list_size(cal_list);

		if (num_cals <= 0) {
			continue;
		}

		read_len = read1->length;
		revcomp_seq = NULL;

		// processing each CAL from this read
		for (size_t j = 0; j < num_cals; j++) {
			// get cal and read index
			cal = array_list_get(j, cal_list);

			if (cal->sr_list->size == 0) {
				continue;
			}

			read = (cal->strand == 0) ? read1 : read2;
			genome = (cal->strand == 0) ? genome1 : genome2;

			if (cal->strand == 0 && bs_id == 0 || cal->strand == 1 && bs_id == 1) {
				c1 = 'C';
				c2 = 'T';
			} else {
				c1 = 'G';
				c2 = 'A';
			}

			sw_prepare = NULL;
			s = (seed_region_t *)linked_list_get_first(cal->sr_list);
			cigar_code = (cigar_code_t *)s->info;

			LOG_DEBUG_F("CAL #%i of %i (strand %i), sr_list size = %i, cigar = %s (distance = %i)\n",
									j, num_cals, cal->strand, cal->sr_list->size, new_cigar_code_string(cigar_code), 
									cigar_code->distance);


			for (int k = 0; k < 2; k++) {
				mode = NONE_POS;

				if (k == 0) {
					if ((cigar_op = cigar_code_get_op(0, cigar_code)) &&
							 cigar_op->name == 'H' && cigar_op->number > min_H) {
						LOG_DEBUG_F("%i%c\n", cigar_op->number, cigar_op->name);

						mode = BEGIN_POS;
						gap_read_start = 0;
						gap_read_end = cigar_op->number - 1;
						gap_genome_start = s->genome_start;
						gap_genome_end = gap_genome_start + cigar_op->number - 1;
					}
				} else {
					if ((cigar_op = cigar_code_get_last_op(cigar_code)) &&
							 cigar_op->name == 'H' && cigar_op->number > min_H) {
						LOG_DEBUG_F("%i%c\n", cigar_op->number, cigar_op->name);

						mode = END_POS;
						gap_read_start = read_len - cigar_op->number;
						gap_read_end = read_len - 1;
						gap_genome_end = s->genome_end;
						gap_genome_start = gap_genome_end - cigar_op->number + 1;
					}
				}

				if (mode == NONE_POS) {
					continue;
				}

				// get query sequence, revcomp if necessary
				// get ref. sequence
				start = gap_genome_start; 
				end = gap_genome_end; 
				gap_genome_len = end - start + 1;

				ref = (char *)malloc((gap_genome_len + 1) * sizeof(char));

				seq = read->sequence;
				genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1,
																					&start, &end, genome);

				gap_read_len = gap_read_end - gap_read_start + 1;
				ref[gap_genome_len] = '\0';

				first = -1;
				last = -1;
				distance = 0;
				for (int k = 0, k1 = gap_read_start; k < gap_read_len; k++, k1++) {
					if (seq[k1] != ref[k]) {
						if (!(seq[k1] == c2 && ref[k] == c1)) {
							distance++;
						}

						if (first == -1) {
							first = k;
						}

						last = k;
					}
				}

				if (distance < min_distance) {
					cigar_op->name = 'M';
					cigar_code->distance += distance;
					free(ref);
					continue;
				}

				// we must run the SW algorithm, prepare this query-ref pair
				query = (char *)malloc((gap_read_len + 1) * sizeof(char));
				memcpy(query, &seq[gap_read_start], gap_read_len);
				query[gap_read_len] = 0;

				sw_prepare = sw_prepare_new(query, ref, 0, 0, (mode == BEGIN_POS ? FIRST_SW : LAST_SW));
				sw_prepare->seed_region = s;
				sw_prepare->cal = cal;
				sw_prepare->read = read;

				// increase counter
				if (cal->strand == 0) {
					array_list_insert(sw_prepare, sw_prepare_list1);
					sw_count1++;
				} else {
					array_list_insert(sw_prepare, sw_prepare_list2);
					sw_count2++;
				}
			}
		}
	}

	int num_ops;
	cigar_op_t *op;
	cigar_code_t *cigar_c;

	for (int type = 0; type < 2; type++) {
		if (type == 0) {
			sw_prepare_list = sw_prepare_list1;
			sw_count = sw_count1;
		} else {
			sw_prepare_list = sw_prepare_list2;
			sw_count = sw_count2;
		}

		// SW pre-processing
		char *q[sw_count], *r[sw_count];

		for (int i = 0; i < sw_count; i++) {
			sw_prepare = array_list_get(i, sw_prepare_list);
			q[i] = sw_prepare->query;
			r[i] = sw_prepare->ref;
		}

		sw_multi_output_t *output = sw_multi_output_new(sw_count);

		// SW execution
		smith_waterman_mqmr(q, r, sw_count, sw_optarg, 1, output);

		// SW post-processing
		for (int i = 0; i < sw_count; i++) {
			sw_prepare = array_list_get(i, sw_prepare_list);
			s = sw_prepare->seed_region;
			cal = sw_prepare->cal;

			if (sw_prepare->ref_type == FIRST_SW) {
				cigar_op = cigar_code_get_first_op(s->info);
			} else {
				cigar_op = cigar_code_get_last_op(s->info);
			}

			gap_read_len = cigar_op->number;
			gap_genome_len = cigar_op->number;

			cigar_code = generate_cigar_code(output->query_map_p[i], output->ref_map_p[i],
													strlen(output->query_map_p[i]), output->query_start_p[i],
													output->ref_start_p[i], gap_read_len, gap_genome_len,
													&distance, sw_prepare->ref_type);


			if (sw_prepare->ref_type == FIRST_SW) {
				cigar_c = (cigar_code_t *)s->info;
				num_ops = cigar_code_get_num_ops(cigar_c);
				cigar_op = cigar_code_get_last_op(cigar_c);

				if (cigar_op->name == 'H' && cigar_op->number > min_H) {
					num_ops--;
				}

				cigar_code->distance += cigar_c->distance;

				for (int i = 1; i < num_ops; i++) {
					op = cigar_code_get_op(i, cigar_c);
					cigar_code_append_new_op(op->number, op->name, cigar_code);
				}

				cal->info = cigar_code;
			} else {
				cigar_c = cal->info;

				if (cigar_c == NULL) {
					cigar_code_t *aux_cigar = (cigar_code_t *)s->info;
					cigar_c = cigar_code_new();
					num_ops = cigar_code_get_num_ops(aux_cigar) - 1;

					for (int i = 0; i < num_ops; i++) {
						op = cigar_code_get_op(i, aux_cigar);
						cigar_code_append_new_op(op->number, op->name, cigar_c);
					}

					cigar_c->distance += aux_cigar->distance;
					cal->info = cigar_c;
				}

				num_ops = cigar_code_get_num_ops(cigar_code);

				for (int i = 0; i < num_ops; i++) {
					op = cigar_code_get_op(i, cigar_code);
					cigar_code_append_new_op(op->number, op->name, cigar_c);
				}

				cigar_c->distance += cigar_code->distance;
				cigar_code_free(cigar_code);
			}

			// free
			sw_prepare_free(sw_prepare);
		}

		// free memory
		sw_multi_output_free(output);
	}
}

//------------------------------------------------------------------------------------

void merge_seed_regions_bs(mapping_batch_t *mapping_batch, int bs_id) {
	array_list_t **mapping_lists;
	register size_t num_targets;
	size_t *targets;

	if (bs_id == 0) {
		mapping_lists = mapping_batch->mapping_lists;
		num_targets = mapping_batch->num_targets;
		targets = mapping_batch->targets;
	} else {
		mapping_lists = mapping_batch->mapping_lists2;
		num_targets = mapping_batch->num_targets2;
		targets = mapping_batch->targets2;
	}

	cal_t *cal;
	seed_region_t *s, *s_first;
	cigar_code_t *cigar_code, *cigar_code_prev;
	cigar_op_t *cigar_op;
	int num_ops;
	int op;
	array_list_t *cals_list;
	linked_list_iterator_t itr;

	register size_t num_cals;
	register int i;
	register size_t t;

	for (t = 0; t < num_targets; t++) {
		cals_list = mapping_lists[targets[t]];
		num_cals = array_list_size(cals_list);

		for (i = 0; i < num_cals; i++) {
			cal = array_list_get(i, cals_list);
			cal->info = NULL;
			linked_list_iterator_init(cal->sr_list, &itr);

			s_first = linked_list_iterator_curr(&itr);

			if (s_first) {
				cigar_code_prev = (cigar_code_t *)s_first->info;
				s = linked_list_iterator_next(&itr);

				while (s) {
					cigar_code = (cigar_code_t *)s->info;

					if (cigar_code) {
						num_ops = array_list_size(cigar_code->ops);
						
						for (op = 0, cigar_op = array_list_get(op, cigar_code->ops); 
								 op < num_ops; 
								 op++, cigar_op = array_list_get(op, cigar_code->ops)) {
							cigar_code_append_new_op(cigar_op->number, cigar_op->name, cigar_code_prev);
						}

						cigar_code_prev->distance += cigar_code->distance;
						cigar_code_free(cigar_code);
					}

					s_first->read_end = s->read_end;
					s_first->genome_end = s->genome_end;

					seed_region_free(s);

					linked_list_iterator_remove(&itr);
					s = linked_list_iterator_curr(&itr);
				}
			}
		}
	}
}
