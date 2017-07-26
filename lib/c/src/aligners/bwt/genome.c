#include "genome.h"

#define NUCLEOTIDES_NUM 5

//---------------------------------------------------------------------------------
//                     GLOBAL VARIABLES

const char NUCLEOTIDES[] = {'A', 'C', 'G', 'N', 'T'};
const unsigned char TOTAL_CODES = NUCLEOTIDES_NUM*NUCLEOTIDES_NUM*NUCLEOTIDES_NUM + NUCLEOTIDES_NUM*NUCLEOTIDES_NUM + NUCLEOTIDES_NUM; 

//------------------------------------------------------------------------------------

genome_t* genome_new(char* sequence_filename, char* directory) {
  const int MAXLINE = 4096;
  genome_t* genome_p = (genome_t*) calloc(1, sizeof(genome_t));
  
  // genome file
  //
  size_t dna_size;
  size_t j, i;
  char path[strlen(directory) + 512];
  
  // read compressed genome file
  sprintf(path, "%s/%s", directory, sequence_filename);




  LOG_DEBUG("Loading Binary DNA File");
  genome_p->X = load_binary_dna(path, &dna_size);
  LOG_DEBUG("Load DNA File Done!");

  genome_p->code_table = load_array_codes();

  // read index file
  sprintf(path, "%s/index", directory);
  //printf("directorio Index: %s\n", path);
  unsigned int chromosome, num_chromosomes = 0;
  unsigned int offset = 0;
  char* p;
  char line[MAXLINE];
  char value[1024];

  FILE *fd = fopen(path, "r");
  if(fd == NULL) { 
    LOG_FATAL_F("FILE: '%s' not found", path);
  }

  int ch, number_of_lines = 0;
  do {
    ch = fgetc(fd);
    if (ch == '\n')
      num_chromosomes++;
  } while (ch != EOF);

  fseek(fd, 0, SEEK_SET);

  genome_p->num_chromosomes = num_chromosomes;
  genome_p->chr_name = (char **) calloc(num_chromosomes, sizeof(char *));
  genome_p->chr_name_length = (size_t *) calloc(num_chromosomes, sizeof(size_t));
  genome_p->chr_size = (size_t *) calloc(num_chromosomes, sizeof(size_t));
  genome_p->chr_offset = (size_t *) calloc(num_chromosomes, sizeof(size_t));

  chromosome = 0;
  while (fgets(line, MAXLINE, fd)) {
    i = 0; j= 1;
    genome_p->chr_name[chromosome] = (char *) calloc(strlen(line) + 10, sizeof(char));
    while(line[j] != ' ' ){ genome_p->chr_name[chromosome][i++] = line[j++]; }
    genome_p->chr_name[chromosome][i] = '\0';
    genome_p->chr_name_length[chromosome] = strlen(genome_p->chr_name[chromosome]);
    
    j++;
    while(line[j] != ' '){j++;}
    
    i=0;j++;

    while(line[j] != '\n'){value[i++] = line[j++];}
    value[i] = '\0';
    
    sscanf(value, "%lu", &genome_p->chr_size[chromosome]);
    genome_p->chr_offset[chromosome] = offset;
    offset += (genome_p->chr_size[chromosome] + 1);
    chromosome++;
  }

  fclose(fd);
  
  return genome_p;
}

//-----------------------------------------------------

void genome_free(genome_t* p) {
  if (p == NULL) {
    return;
  }

  if (p->code_table != NULL){
    for (unsigned int i = 0; i < TOTAL_CODES; i++) {
      free(p->code_table[i]);
    }
    free(p->code_table);
  }

  if (p->chr_name != NULL) {
    for (unsigned int i = 0; i < p->num_chromosomes; i++) {
      free(p->chr_name[i]);
    }
    free(p->chr_name);
  }

  if (p->chr_name_length) free(p->chr_name_length);
  if (p->chr_size) free(p->chr_size);
  if (p->chr_offset) free(p->chr_offset);
  if (p->X) free(p->X);

  free(p);
}

//------------------------------------------------------------------------------------

void genome_read_sequence(char* sequence, unsigned int strand, char* chromosome, 
			  unsigned long int* start_p, unsigned long int* end_p, 
			  genome_t* genome_p) {
  unsigned int i, chr;

  for(i=0 ; i< genome_p->num_chromosomes ; i++) {
    if (strcmp(chromosome, (char *) genome_p->chr_name[i]) == 0) {
      chr = i;
      break;
    }
  }

  genome_read_sequence_by_chr_index(sequence, strand, chr, start_p, end_p, genome_p);
}

//------------------------------------------------------------------------------------

cp_hashtable *load_hasthable_codes() {
  cp_hashtable *t = cp_hashtable_create(400,
					cp_hash_istring,
					(cp_compare_fn)strcmp); 
  
  /*const unsigned char COMBINATORIAL = (NUCLEOTIDES_NUM * NUCLEOTIDES_NUM * NUCLEOTIDES_NUM) +  (NUCLEOTIDES_NUM * NUCLEOTIDES_NUM) +  NUCLEOTIDES_NUM;
  
  unsigned char *id_array = (unsigned char *)malloc(sizeof(unsigned char)*COMBINATORIAL); 
  */
  size_t id = 0;
  char combination[4];
  
  combination[3] = '\0';
  
  for(unsigned int nt_1 = 0; nt_1 < NUCLEOTIDES_NUM; nt_1++){
    combination[0] = NUCLEOTIDES[nt_1];
    for(unsigned int nt_2 = 0; nt_2 < NUCLEOTIDES_NUM; nt_2++){
      combination[1] = NUCLEOTIDES[nt_2];
      for(unsigned int nt_3 = 0; nt_3 < NUCLEOTIDES_NUM; nt_3++){
	combination[2] = NUCLEOTIDES[nt_3];
	cp_hashtable_put(t, strdup(combination), (void *)id);
	id++;
      }
    }
  }
  
  // printf("Table size %d\n", t->table_size);
  
  combination[2] = '\0';  
  for(unsigned int nt_1 = 0; nt_1 < NUCLEOTIDES_NUM; nt_1++){
    combination[0] = NUCLEOTIDES[nt_1];
    for(unsigned int nt_2 = 0; nt_2 < NUCLEOTIDES_NUM; nt_2++){
      combination[1] = NUCLEOTIDES[nt_2];
      cp_hashtable_put(t, strdup(combination), (void *)id);
      id++;	
    }
  }
  
  //printf("Table size %d\n", t->table_size);
  
  combination[1] = '\0';
  for(unsigned int nt = 0; nt < NUCLEOTIDES_NUM; nt++){
    combination[0] = NUCLEOTIDES[nt];
    cp_hashtable_put(t, strdup(combination), (void *)id);
    id++;	
  }

  return t;

}


char** load_array_codes() {
  char **id_array = (char **)malloc(sizeof(char *)*TOTAL_CODES);
  
  unsigned char id = 0;
  char combination[4];
  
  combination[3] = '\0';
  
  for (unsigned int nt_1 = 0; nt_1 < NUCLEOTIDES_NUM; nt_1++) {
    combination[0] = NUCLEOTIDES[nt_1];

    for (unsigned int nt_2 = 0; nt_2 < NUCLEOTIDES_NUM; nt_2++) {
      combination[1] = NUCLEOTIDES[nt_2];

      for (unsigned int nt_3 = 0; nt_3 < NUCLEOTIDES_NUM; nt_3++) {
	      combination[2] = NUCLEOTIDES[nt_3];
	      id_array[id] = strdup(combination);
	      id++;
      }
    }
  }
  
  combination[2] = '\0';
  
  for (unsigned int nt_1 = 0; nt_1 < NUCLEOTIDES_NUM; nt_1++) {
    combination[0] = NUCLEOTIDES[nt_1];

    for (unsigned int nt_2 = 0; nt_2 < NUCLEOTIDES_NUM; nt_2++) {
      combination[1] = NUCLEOTIDES[nt_2];
      id_array[id] = strdup(combination);
      id++;	
    }
  }
   
  combination[1] = '\0';

  for (unsigned int nt = 0; nt < NUCLEOTIDES_NUM; nt++) {
    combination[0] = NUCLEOTIDES[nt];
    id_array[id] = strdup(combination);
    id++;	
  }

  return id_array;
}


void code_binary_file_generator(size_t chunk_size, char *dna_filename, char *dna_binary_filename,  cp_hashtable *t) {
  if (chunk_size <= 0) { 
    chunk_size = 100000000; //100MB 
  }

  FILE *binary_fd, *fd;
  fd = fopen(dna_filename, "r");
  if (fd == NULL) {  LOG_FATAL_F("Error opening file %s", dna_filename); }

  binary_fd = fopen (dna_binary_filename, "wb");
  if (binary_fd == NULL) { LOG_FATAL_F("Error opening file %s", dna_binary_filename); }

  char *dna_chunk = (char *)malloc(sizeof(char)*chunk_size);
  
  size_t codes_allocate = chunk_size;
  unsigned char *code_values = (unsigned char *)malloc(sizeof(unsigned char)*codes_allocate);
  size_t code_pos = 0;

  size_t dna_len;
  char key[4];
  unsigned char max_chunk = 3;
  unsigned char actual_nt = 0;
  
  unsigned char value;
  unsigned char *value_ptr;
  size_t nt = 0;
  unsigned char key_chunk = 3;
  LOG_DEBUG("Process DNA File\n");

  while (!feof(fd)) {
    fgets(dna_chunk, chunk_size, fd);
    if (dna_chunk[0] != '>') {
      dna_len = strlen(dna_chunk);

      for (unsigned int c = 0; c < dna_len; c++) {
	    if (dna_chunk[c] != '\n') {
	      if (dna_chunk[c] == 'a' ||
	          dna_chunk[c] == 'c' ||
	          dna_chunk[c] == 'g' ||
	          dna_chunk[c] == 't' ||
	          dna_chunk[c] == 'n') {
	        dna_chunk[c] = dna_chunk[c] - 32;
	      }

	      key[actual_nt++] = dna_chunk[c];

	      if (actual_nt ==  max_chunk) {
	        key[actual_nt] = '\0';
	        value = (unsigned char)cp_hashtable_get(t, key);
	        code_values[code_pos++] = value;

	        if (code_pos >= codes_allocate) {
	          fwrite(code_values, sizeof(unsigned char), code_pos, binary_fd);
	          code_pos = 0;
	        }

	        actual_nt = 0;
	      }
	    }
      } //End for
    } else {
      LOG_DEBUG_F("Process: %s", &dna_chunk[1]);      
    }//End if strcmp
  }
    
  if(actual_nt > 0) {
    key[actual_nt] = '\0';
    value = (unsigned char)cp_hashtable_get(t, key);
    code_values[code_pos++] = value;	
  }

  if (code_pos >= 0) {
    fwrite(code_values, sizeof(unsigned char), code_pos, binary_fd);
    code_pos = 0;
  }
      
  fclose(fd);
  fclose(binary_fd);
}

unsigned char *load_binary_dna(char *dna_binary_filename, size_t *size){
  FILE *binary_fd = fopen (dna_binary_filename, "rb");
  if (!binary_fd) {
    LOG_FATAL_F("Error to opening '%s' file\n", dna_binary_filename);
  }

  struct stat st;                                                                                                                                                 
  stat(dna_binary_filename, &st);     
  *size = st.st_size;                                                                                                               
  LOG_DEBUG_F("Size File %lu\n", *size);

  unsigned char *dna_encoding = (unsigned char *)malloc(sizeof(unsigned char)*(*size));

  fread(dna_encoding, sizeof(unsigned char), *size, binary_fd);
  
  fclose(binary_fd);

  return dna_encoding;
}

void get_genome_sequence(char *sequence, unsigned int chromosome, unsigned int strand, size_t start, size_t end, char **array_codes, char *dna_encoding){
  size_t group_start = start/3;
  size_t group_end   = end/3;

  size_t nucleotide_start = start%3;
  size_t nucleotide_end   = end%3;

  //size_t total_nt = end - start;
  //size_t total_groups = total_nt/3;
  //size_t rest_groups = total_nt%3;

  unsigned int seq_pos = 0;
  
  char nt_code[4];
  unsigned char id;
  
  if (start >= end)
    return;

  for(unsigned int i = nucleotide_start; i < 3; i++){
    sequence[seq_pos++] = array_codes[dna_encoding[group_start]][i];
  }

  group_start++;
  
  while (group_start < group_end) {
    sequence[seq_pos++] = array_codes[dna_encoding[group_start]][0];
    sequence[seq_pos++] = array_codes[dna_encoding[group_start]][1];
    sequence[seq_pos++] = array_codes[dna_encoding[group_start]][2];
    group_start++;
  }

  if (group_start <= group_end) {
    for(unsigned int i = 0; i <= nucleotide_end; i++){
      sequence[seq_pos++] = array_codes[dna_encoding[group_start]][i];
    }
  }

  sequence[seq_pos] = '\0';
}


void genome_read_sequence_by_chr_index(char* sequence, unsigned int strand, 
				       unsigned int chr, size_t *start_p,
				       size_t *end_p, genome_t* genome_p) {
  size_t s, e;
  
  if (*start_p < 1) (*start_p) = 1;
  if (*start_p > genome_p->chr_size[chr]) (*start_p) = genome_p->chr_size[chr];
  if (*end_p < 1) (*end_p) = 1;
  if (*end_p > genome_p->chr_size[chr]) (*end_p) = genome_p->chr_size[chr];

  s = (*start_p) + genome_p->chr_offset[chr] - 1;
  e = (*end_p) + genome_p->chr_offset[chr] - 1;
  
  size_t group_start = s/3;
  size_t group_end   = e/3;

  size_t nucleotide_start = s%3;
  size_t nucleotide_end   = e%3;

  unsigned int seq_pos = 0;
  
  char nt_code[4];
  unsigned char id;
  
  for(unsigned int i = nucleotide_start; i < 3; i++){
    sequence[seq_pos++] = genome_p->code_table[genome_p->X[group_start]][i];
  }

  group_start++;
  
  while (group_start < group_end) {
    sequence[seq_pos++] = genome_p->code_table[genome_p->X[group_start]][0];
    sequence[seq_pos++] = genome_p->code_table[genome_p->X[group_start]][1];
    sequence[seq_pos++] = genome_p->code_table[genome_p->X[group_start]][2];
    group_start++;
  }
  
  if (group_start <= group_end) {
    for(unsigned int i = 0; i <= nucleotide_end; i++){
      sequence[seq_pos++] = genome_p->code_table[genome_p->X[group_start]][i];
    }
  }

  sequence[seq_pos] = '\0';
}


//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------                                                                                  
void generate_codes(char *dna_binary_filename, char *dna_filename){
  LOG_DEBUG("Loading hashtable Codes ...\n");
  cp_hashtable *t = load_hasthable_codes();
  LOG_DEBUG("Loading done!\n");

  LOG_DEBUG("Genrate Binary Genome File...\n");
  code_binary_file_generator(100000000, dna_filename,
                             dna_binary_filename, t);
  LOG_DEBUG("Generate done! Happy usage!\n");
}

//------------------------------------------------------------------------------------
