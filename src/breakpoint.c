#include "breakpoint.h"


//------------------------------------------------------------------------------------

cigar_op_t *cigar_op_new(int number, char name) {
  cigar_op_t *p = (cigar_op_t *) malloc(sizeof(cigar_op_t));

  p->number = number;
  p->name = name;

  return p;
}

//------------------------------------------------------------------------------------

void cigar_op_free(cigar_op_t *p) {
  if (p) {
    free(p);
  }
}

//--------------------------------------------------------------------------------------

cigar_code_t *cigar_code_new() {
  cigar_code_t *p = (cigar_code_t *) calloc(1, sizeof(cigar_code_t));

  p->distance = 0;
  p->cigar_str = NULL;
  p->ops = array_list_new(20, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  return p;
}

cigar_code_t *cigar_code_new_by_string(char *cigar_str) {
  cigar_code_t *p = cigar_code_new();
  int cigar_len = strlen(cigar_str);
  int c = 0;
  char op;
  char op_value[1024];

  for (int j = 0; j < cigar_len; j++) {
    op = cigar_str[j];

    if (op < 58) {
      op_value[c++] = op;
    } else {
      op_value[c] = '\0';
      cigar_code_append_new_op(atoi(op_value), op, p);
      c = 0;     
    }
  }
  
  return p;
}

//--------------------------------------------------------------------------------------

void cigar_code_free(cigar_code_t* p) {
  if (p) {
    if (p->ops) {
      array_list_free(p->ops, (void *) cigar_op_free);
    }

    if (p->cigar_str) {
      free(p->cigar_str);
    }

    free(p);
  }
}

//--------------------------------------------------------------------------------------

void cigar_code_update(cigar_code_t *p) {
  cigar_op_t *prev_op, *curr_op;

  size_t j = 0, num_ops = array_list_size(p->ops);
  prev_op = array_list_get(0, p->ops);

  for (size_t i = 1; i < num_ops; i++) {
    curr_op = array_list_get(i, p->ops);

    if (prev_op->name == curr_op->name) {
      prev_op->number += curr_op->number;
    } else {
      j++;
      prev_op = array_list_get(j, p->ops);
      prev_op->name = curr_op->name;
      prev_op->number = curr_op->number;
    }
  }

  for (size_t i = num_ops - 1; i > j; i--) {
    array_list_remove_at(i, p->ops);
  }
}

//--------------------------------------------------------------------------------------

int cigar_code_get_num_ops(cigar_code_t *p) {
  int num = 0;

  if (p && p->ops) { 
    return array_list_size(p->ops);
  }

  return num;
}

//--------------------------------------------------------------------------------------

cigar_op_t *cigar_code_get_first_op(cigar_code_t *p) {
  if (!p) { 
    return NULL; 
  }  

  return array_list_get(0, p->ops);
}

cigar_op_t *cigar_code_get_op(int index, cigar_code_t *p) {
  int num_ops = cigar_code_get_num_ops(p);

  if (num_ops > 0 && index < num_ops) {
    return array_list_get(index, p->ops);
  }

  return NULL;
}

cigar_op_t *cigar_code_get_last_op(cigar_code_t *p) {
  int num_ops = cigar_code_get_num_ops(p);
  return array_list_get(num_ops - 1, p->ops);
}

//--------------------------------------------------------------------------------------

void cigar_code_add_op(cigar_op_t *op, cigar_code_t *p) {
  if (p && p->ops && op) {
    array_list_insert(op, p->ops);
  }
}

//--------------------------------------------------------------------------------------

void cigar_code_append_op(cigar_op_t *op, cigar_code_t *p) {
  if (p && p->ops && op) {
    cigar_op_t *last = cigar_code_get_last_op(p);

    if (last && last->name == op->name) {
      last->number += op->number;
      cigar_op_free(op);
    } else {
      int num_ops = cigar_code_get_num_ops(p);

      if (last && num_ops > 1 && last->name == 'S') {
        printf("00 change S -> M : %s (insert %i%c)\n", new_cigar_code_string(p), 
               op->number, op->name);
        
        // A little bit tricky
        cigar_op_t *last1 = cigar_code_get_op(num_ops - 2, p);

        if (last1->name == 'M') {
          array_list_remove_at(num_ops - 1, p->ops);
          last1->number += last->number;
        } else {
          last->name = 'M';
        }

        printf("11 change S -> M : %s (insert %i%c)\n", new_cigar_code_string(p), 
               op->number, op->name);
      }
      
      array_list_insert(op, p->ops);
    }
  }
}

//--------------------------------------------------------------------------------------

void cigar_code_append_new_op(int value, char name, cigar_code_t *p) {
  cigar_code_append_op(cigar_op_new(value, name), p);
}

//--------------------------------------------------------------------------------------

void cigar_code_inc_distance(int distance, cigar_code_t *p) {
  if (p && p->ops) {
    p->distance += distance;
  }
}

//--------------------------------------------------------------------------------------

char *new_cigar_code_string(cigar_code_t *p) {
  if (!p) { 
    return "\0"; 
  }

  if (p->cigar_str) {
    free(p->cigar_str);
  }

  int num_ops = cigar_code_get_num_ops(p);

  if (num_ops == 0) { 
    return NULL; 
  }

  // Using sprintf to append data to the same buffer produces a
  // memory overlapping memcpy! This is undefined behavior under
  // the C standard
  char *str = (char *)calloc(num_ops * 16 + 1, sizeof(char));
  char temp[16];

  cigar_op_t *op;
  
  // Use safe functions to prevent buffer overflows if the CIGAR
  // op length overflowed and the resulting string was too long
  // to fit
  for (int i = 0; i < num_ops; i++) {
    op = array_list_get(i, p->ops);
    
    snprintf(temp, 16, "%i%c", op->number, op->name);
    strncat(str, temp, 16);
  }

  p->cigar_str = str;
  return p->cigar_str;
}

//--------------------------------------------------------------------------------------

int cigar_code_nt_length(cigar_code_t *p) {
  if (!p) {
    return 0;
  }

  int len = 0;
  int num_ops = array_list_size(p->ops);

  cigar_op_t *op;

  for (int i = 0; i < num_ops; i++) {
    op = array_list_get(i, p->ops);

    if (op->name == 'M' || op->name == 'I' || op->name == '=') {
      len += op->number;
    }
  }

  return len;
}

//--------------------------------------------------------------------------------------

cigar_code_t *generate_cigar_code(char *query_map, char *ref_map, unsigned int map_len,
				  unsigned int query_start, unsigned int ref_start,
				  unsigned int query_len, unsigned int ref_len,
				  int *distance, int ref_type) {
  cigar_code_t *p = cigar_code_new();

  unsigned int perfect = 0;  
  unsigned int deletions_tot = 0;
  unsigned int insertions_tot = 0;
  unsigned int map_ref_len;
  unsigned int map_seq_len;
  unsigned int last_h, last_h_aux;
  unsigned int number_op = 0;
  short int cigar_soft;
  short int value = 0;
  int dist = 0;
  
  unsigned char status;
  unsigned char transition;
  char operation;


  if (query_start > 0) {
    if (ref_type == FIRST_SW) {
      // Normal Case
      cigar_code_append_op(cigar_op_new(query_start, 'H'), p);
    } else {
      // Middle or last ref
      if (ref_start == 0) {
	      cigar_code_append_op(cigar_op_new(query_start, 'I'), p);
      } else {
        if (ref_start == query_start) {
          cigar_code_append_op(cigar_op_new(query_start, 'M'), p);
        } else {
          if (ref_start > query_start) {
            cigar_code_append_op(cigar_op_new(ref_start - query_start, 'D'), p);
            cigar_code_append_op(cigar_op_new(query_start, 'M'), p);
          } else {
            cigar_code_append_op(cigar_op_new(query_start - ref_start, 'I'), p);
            cigar_code_append_op(cigar_op_new(ref_start, 'M'), p);
          } 
        }
      }
    }
  } else if (ref_start > 0) {
    if (ref_type != FIRST_SW) {
      cigar_code_append_op(cigar_op_new(ref_start, 'D'), p);
    } 
  }
  
  // First Status
  if (query_map[0] != '-' && ref_map[0] != '-') {
    status = CIGAR_MATCH_MISMATCH;

    // Soft clipping
    cigar_soft = 0;

    while ((ref_map[cigar_soft] != '-') && (query_map[cigar_soft] != '-') && 
	         (ref_map[cigar_soft] != query_map[cigar_soft])) {
      cigar_soft++;
      value++;
    }

    if (value > 0) {
      if (ref_start == 0 && query_start == 0) {
	      cigar_code_append_op(cigar_op_new(value, 'S'), p);
      } else {
	      value = 0;
      }
    } 
  } else if (query_map[0] == '-') {
    if (ref_map[0] == '-') {
      status = CIGAR_PADDING;
    } else {
      status = CIGAR_DELETION;
    }
  } else if(ref_map[0] == '-') {
    status = CIGAR_INSERTION;
  }
  
  for (int i = value; i < map_len; i++) {
    // Transition
    if (query_map[i] != '-' && ref_map[i] != '-') {
      transition = CIGAR_MATCH_MISMATCH;

      if (query_map[i] == ref_map[i]) {
        perfect++;
      } else {
         dist++;
      }
    } else if(query_map[i] == '-') {
      if (ref_map[i] == '-') {
        transition = CIGAR_PADDING;
      } else {
        transition = CIGAR_DELETION;
        deletions_tot++;
      }
    } else if (ref_map[i] == '-') {
      transition = CIGAR_INSERTION;
      insertions_tot++;
    }
    
    if (transition != status) {
      // Insert operation in cigar string
      operation = select_op(status);
      cigar_code_append_op(cigar_op_new(number_op, operation), p);
      number_op = 1;
      status = transition;
    } else {
      number_op++;
    }
  }
  
  operation = select_op(status);
  
  // Hard and Soft clipped end
  if (status == CIGAR_MATCH_MISMATCH) {
    cigar_soft = map_len - 1;
    value = 0;

    while (cigar_soft >= 0               && 
	        (ref_map[cigar_soft] != '-')   && 
          (query_map[cigar_soft] != '-') && 
	        (query_map[cigar_soft] != ref_map[cigar_soft])) {
      cigar_soft--;
      value++;
    }
    
    cigar_code_append_op(cigar_op_new(number_op - value, operation), p);
    
    if (value > 0) {
      number_op -= value;
      cigar_code_append_op(cigar_op_new(value, 'S'), p);
    }
  } else {
    cigar_code_append_op(cigar_op_new(number_op, operation), p);
  }

  map_seq_len  = ((map_len - deletions_tot) + query_start);
  map_ref_len  = ((map_len - insertions_tot) + ref_start);

  if (map_seq_len < query_len) {
    last_h = query_len - map_seq_len;
    
    if (ref_type == LAST_SW) {
      // Normal Case
      cigar_code_append_op(cigar_op_new(last_h, 'H'), p);
    } else {
      // Middle or first ref
      if (map_ref_len == ref_len) {
	      cigar_code_append_op(cigar_op_new(last_h, 'I'), p);
      } else {
        last_h_aux = ref_len - map_ref_len;
        
        if (last_h_aux == last_h) {
          cigar_code_append_op(cigar_op_new(last_h, 'M'), p);
        } else {	  
          if (last_h_aux > last_h) {
            cigar_code_append_op(cigar_op_new(last_h_aux - last_h, 'D'), p);
            cigar_code_append_op(cigar_op_new(last_h, 'M'), p);
          } else {
            cigar_code_append_op(cigar_op_new(last_h - last_h_aux, 'I'), p);
            cigar_code_append_op(cigar_op_new(last_h_aux, 'M'), p);
          } 
        }
      }
    }
  } else if (map_ref_len < ref_len) {
    if (ref_type != LAST_SW) {
      cigar_code_append_op(cigar_op_new(ref_len - map_ref_len, 'D'), p);
    }
  }
  
  *distance = dist;
  p->distance = dist;

  return p;
}
