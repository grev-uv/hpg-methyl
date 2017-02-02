#ifndef BREAKPOINT_H
#define BREAKPOINT_H

#include "containers/array_list.h"
#include "containers/linked_list.h"

#include "aligners/bwt/genome.h"
#include "bioformats/bam/alignment.h"
#include "bioformats/fastq/fastq_read.h"
//--------------------------------------------------------------------------------------

#define FIRST_SW 0
#define MIDDLE_SW 1
#define LAST_SW 2

//--------------------------------------------------------------------------------------
//Metaexon constants
#define METAEXON_NORMAL 1
#define METAEXON_LEFT_END 2
#define METAEXON_RIGHT_END 3

//====================================================================================
//  Input structure for CIGAR format
//====================================================================================

typedef struct cigar_op {
  int number;
  char name;
} cigar_op_t;

cigar_op_t *cigar_op_new(int number, char name);
void cigar_op_free(cigar_op_t *cigar_op);

typedef struct cigar_code {
  int distance;
  char *cigar_str;
  array_list_t *ops;
} cigar_code_t;

cigar_code_t *cigar_code_new();
cigar_code_t *cigar_code_new_by_string(char *cigar_str);
void cigar_code_free(cigar_code_t* p);

char *new_cigar_code_string(cigar_code_t* p);
char *cigar_code_get_string(cigar_code_t *p);
int cigar_code_get_num_ops(cigar_code_t *p);
cigar_code_t *cigar_code_merge_sp(cigar_code_t *cc_left,
				  cigar_code_t *cc_middle, 
				  cigar_code_t *cc_right,
				  int l_flank, int r_flank);
  
//-----------------------------------------------------------------------------------

cigar_op_t *cigar_code_get_first_op(cigar_code_t *p);
cigar_op_t *cigar_code_get_op(int index, cigar_code_t *p);
cigar_op_t *cigar_code_get_last_op(cigar_code_t *p);

//-----------------------------------------------------------------------------------

void cigar_code_inc_distance(int distance, cigar_code_t *p);
void cigar_code_add_op(cigar_op_t *op, cigar_code_t *p);
void cigar_code_append_new_op(int value, char name, cigar_code_t *p);
void cigar_code_append_op(cigar_op_t *op, cigar_code_t *p);

int cigar_code_nt_length(cigar_code_t *p);

cigar_code_t *generate_cigar_code(char *query_map, char *ref_map, unsigned int map_len,
				  unsigned int query_start, unsigned int ref_start,
				  unsigned int query_len, unsigned int ref_len,
				  int *distance, int ref_type);
void cigar_code_update(cigar_code_t *p);

//====================================================================================
//  Input structure for breakpoint definition
//====================================================================================

#define LEFT_SIDE  0
#define RIGHT_SIDE 1

typedef struct breakpoint_info {
  int side;
  int num_M;
  int index_M;
  cigar_code_t *cigar;
} breakpoint_info_t;

breakpoint_info_t *breakpoint_info_new(int side, int num_M, int index_M,
				       cigar_code_t *cigar);
void breakpoint_info_free(breakpoint_info_t *p);

//--------------------------------------------------------------------------------------

typedef struct breakpoint {
  int chr_index;
  int position;
  int coverage;
  array_list_t *info_list;
} breakpoint_t;

breakpoint_t *breakpoint_new(int chr_index, int position);
void breakpoint_free(breakpoint_t *p);

breakpoint_t *compute_breakpoint(int chr_index, int alig_position, 
				 cigar_code_t *cigar, array_list_t *list);

void display_breakpoints(array_list_t *list);

//====================================================================================
//  Main structure for metaexons
//====================================================================================

typedef struct metaexon {
  size_t start;
  size_t end;
  unsigned char left_closed; //Is the left extrem metaexon closed? True or False
  unsigned char right_closed; //Is the left extrem metaexon closed? True or False
  array_list_t *left_breaks;
  array_list_t *right_breaks;
} metaexon_t;

//--------------------------------------------------------------------------------------

typedef struct metaexon_pair {
  linked_list_item_t *first;
  linked_list_item_t *last;  
} metaexon_pair_t;

typedef struct metaexons {
  unsigned int       num_chromosomes;
  size_t             *num_chunks;
  pthread_mutex_t    *mutex;
  linked_list_t      **metaexons_list;  
  metaexon_pair_t    **bypass_pointer;
} metaexons_t;

#endif  // SW_SERVER_H
