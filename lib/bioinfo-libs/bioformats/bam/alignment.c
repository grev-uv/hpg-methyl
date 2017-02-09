#include "alignment.h"


//------------------------------------------------------------------------------------

alignment_t* alignment_new() {
    alignment_t* alignment_p = (alignment_t*) calloc(1, sizeof(alignment_t));
    return alignment_p;
}

//------------------------------------------------------------------------------------

void alignment_init_single_end(char* query_name, char* sequence, char* quality, short int strand, 
                    short int chromosome, int position, char* cigar, short int num_cigar_operations, 
                    int map_quality, short int is_seq_mapped, short int secondary_alignment, 
                    int optional_fields_length, char *optional_fields, alignment_t* alignment_p) {
    alignment_p->query_name = query_name;
    alignment_p->sequence = sequence;
    alignment_p->quality = quality;
    alignment_p->cigar = cigar;

    alignment_p->chromosome = chromosome;
    alignment_p->position = position;
    alignment_p->mate_position = -1; //0;
    alignment_p->mate_chromosome = -1; //;
    alignment_p->template_length = 0; //single end, no template
    alignment_p->map_quality = map_quality;
    alignment_p->num_cigar_operations = num_cigar_operations;

    alignment_p->is_paired_end = 0;
    alignment_p->is_paired_end_mapped = 0;
    alignment_p->is_seq_mapped = is_seq_mapped;
    alignment_p->is_mate_mapped = 0;
    alignment_p->seq_strand = strand;
    alignment_p->mate_strand = 0; //unused value in single end
    alignment_p->pair_num = 0;    //unused value in single end
    alignment_p->secondary_alignment = secondary_alignment;
    alignment_p->fails_quality_check = 0;
    alignment_p->pc_optical_duplicate = 0;
    
    alignment_p->optional_fields = (uint8_t *)optional_fields;
    alignment_p->optional_fields_length = optional_fields_length;
    alignment_p->optional_tags = array_list_new(30, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
}

alignment_t* alignment_new_by_bam(bam1_t* bam_p, int base_quality) {
    //memory allocation for the structure
    alignment_t* alignment_p = (alignment_t*) calloc(1, sizeof(alignment_t));

    //numeric data
    alignment_p->num_cigar_operations = (int) bam_p->core.n_cigar;
    alignment_p->chromosome = bam_p->core.tid;
    alignment_p->position = bam_p->core.pos;
    alignment_p->mate_chromosome = bam_p->core.mtid;
    alignment_p->mate_position = bam_p->core.mpos;
    alignment_p->map_quality = bam_p->core.qual;
    alignment_p->template_length = bam_p->core.isize;

    //memory allocation for inner fields according to indicated sizes
    alignment_p->query_name = (char*) calloc(bam_p->core.l_qname, sizeof(char));
    alignment_p->sequence = (char*) calloc(bam_p->core.l_qseq + 1, sizeof(char));
    alignment_p->quality = (char*) calloc(bam_p->core.l_qseq + 1, sizeof(char));   //same length as sequence
    
    alignment_p->optional_fields = (uint8_t*) calloc(bam_p->l_aux, sizeof(uint8_t));
    alignment_p->optional_fields_length = bam_p->l_aux;

    // Parse the optional tags
    if (alignment_p->optional_tags) {
        array_list_free(alignment_p->optional_tags, NULL);
    }

    alignment_p->optional_tags = bam_string_to_tag(alignment_p->optional_fields, alignment_p->optional_fields_length);

    //copy the data between structures
    strcpy(alignment_p->query_name, bam1_qname(bam_p));
    strcpy(alignment_p->sequence, convert_to_sequence_string(bam1_seq(bam_p), bam_p->core.l_qseq));

    convert_to_quality_string_length(alignment_p->quality, bam1_qual(bam_p), bam_p->core.l_qseq, base_quality);
    alignment_p->cigar = convert_to_cigar_string(bam1_cigar(bam_p), alignment_p->num_cigar_operations);

    memcpy(alignment_p->optional_fields, bam1_aux(bam_p), bam_p->l_aux);

    //flags
    uint32_t flag = (uint32_t) bam_p->core.flag;
    alignment_p->is_paired_end = (flag & BAM_FPAIRED) ? 1 : 0;
    alignment_p->is_paired_end_mapped = (flag & BAM_FPROPER_PAIR) ? 1 : 0;
    alignment_p->is_seq_mapped = (flag & BAM_FUNMAP) ? 0 : 1; //in bam structure is negative flag!!!
    alignment_p->is_mate_mapped = (flag & BAM_FMUNMAP) ? 0 : 1; //in bam structure is negative flag!!!
    alignment_p->seq_strand = (flag & BAM_FREVERSE) ? 1 : 0;
    alignment_p->mate_strand = (flag & BAM_FMREVERSE) ? 1 : 0;

    if (flag & BAM_FREAD1) {
        alignment_p->pair_num = 1;
    } else if (flag & BAM_FREAD2) {
        alignment_p->pair_num = 2;
    } else {
        alignment_p->pair_num = 0;
    }

    alignment_p->secondary_alignment = (flag & BAM_FSECONDARY) ? 1 : 0;
    alignment_p->fails_quality_check = (flag & BAM_FQCFAIL) ? 1 : 0;
    alignment_p->pc_optical_duplicate = (flag & BAM_FDUP) ? 1 : 0;

    return alignment_p;
}

//------------------------------------------------------------------------------------

void alignment_free_tags_cb(void* data) {
    bam_tag_t* tg = (bam_tag_t*)data;
    bam_tag_free(tg);
}

//------------------------------------------------------------------------------------

void alignment_free(alignment_t* p) {
  if (p == NULL) return;

  if (p->query_name != NULL) free(p->query_name);
  if (p->sequence != NULL) free(p->sequence);
  if (p->quality != NULL) free(p->quality);
  if (p->cigar != NULL) free(p->cigar);
  if (p->optional_fields != NULL) free(p->optional_fields);
  if (p->optional_tags != NULL) array_list_free(p->optional_tags, alignment_free_tags_cb);

  free(p);
}

//------------------------------------------------------------------------------------

bam1_t* convert_to_bam(alignment_t* alignment_p, int base_quality) {
    bam1_t* bam_p = bam_init1();
    int data_length = 0, sequence_length = 0, copy_length = 0, index_to_data = 0;
    size_t op_length = 0;
    uint8_t* data = NULL;
    char* op_string = NULL;

    sequence_length = strlen(alignment_p->sequence);

    // Serialize the optional fields
    op_string = bam_tag_list_to_string(alignment_p->optional_tags, &op_length);
    alignment_p->optional_fields_length = op_length;

    //data length is the sum of lengths if five codified fields (cigar, query name, sequence, quality and optional info)
    data_length = (4 * alignment_p->num_cigar_operations) + strlen(alignment_p->query_name) + 1 + ((sequence_length + 1) / 2) + sequence_length + alignment_p->optional_fields_length;

    //memory allocation for data vector from data length
    data = (uint8_t*) calloc(data_length, sizeof(uint8_t));

    //copy query name
    copy_length = strlen(alignment_p->query_name) + 1;
    strcat(alignment_p->query_name, "\0");
    memcpy(&data[index_to_data], alignment_p->query_name, copy_length);
    index_to_data += copy_length;

    //convert cigar to uint32_t format
    convert_to_cigar_uint32_t(&data[index_to_data], alignment_p->cigar, alignment_p->num_cigar_operations);

    copy_length = (4 * alignment_p->num_cigar_operations);
    index_to_data += copy_length;

    //convert sequence to uint8_t format
    convert_to_sequence_uint8_t(&data[index_to_data], alignment_p->sequence, sequence_length);

    copy_length = ((sequence_length + 1) / 2);
    index_to_data += copy_length;

    //convert quality to uint8_t format
    convert_to_quality_uint8_t(&data[index_to_data], alignment_p->quality, sequence_length, base_quality);

    copy_length = sequence_length;
    index_to_data += copy_length;

    //copy optional fields
    memcpy(&data[index_to_data], op_string, op_length);
    free(op_string);

    //finally data is assigned
    bam_p->data = data;

    //filling bam1_t (not core data)
    bam_p->l_aux = alignment_p->optional_fields_length;
    bam_p->data_len = data_length;
    bam_p->m_data = data_length;

    //filling bam1_core_t structure
    bam_p->core.tid = (int32_t) alignment_p->chromosome;
    bam_p->core.pos = (int32_t) alignment_p->position;
    bam_p->core.mtid = (int32_t) alignment_p->mate_chromosome;
    bam_p->core.mpos = (int32_t) alignment_p->mate_position;
    bam_p->core.qual = (uint32_t) alignment_p->map_quality;
    bam_p->core.isize = (int32_t) alignment_p->template_length;
    bam_p->core.l_qname = strlen(alignment_p->query_name) + 1;
    bam_p->core.n_cigar = (uint32_t) alignment_p->num_cigar_operations;
    bam_p->core.l_qseq = (int32_t)(int32_t)bam_cigar2qlen(&bam_p->core, bam1_cigar(bam_p)); //lenght from CIGAR

    //setting flags
    if (alignment_p->is_paired_end)   bam_p->core.flag += BAM_FPAIRED;
    if (alignment_p->is_paired_end_mapped) bam_p->core.flag += BAM_FPROPER_PAIR;
    if (!alignment_p->is_seq_mapped)   bam_p->core.flag += BAM_FUNMAP;   //in bam structure is negative flag!!!
    if ((!alignment_p->is_mate_mapped) && (alignment_p->is_paired_end))   bam_p->core.flag += BAM_FMUNMAP; //in bam structure is negative flag!!!
    if (alignment_p->seq_strand)    bam_p->core.flag += BAM_FREVERSE;
    if (alignment_p->mate_strand)   bam_p->core.flag += BAM_FMREVERSE;

    if (alignment_p->pair_num == 1) {
        bam_p->core.flag += BAM_FREAD1;
    } else if (alignment_p->pair_num == 2) {
        bam_p->core.flag += BAM_FREAD2;
    }

    if (alignment_p->secondary_alignment)    bam_p->core.flag += BAM_FSECONDARY;
    if (alignment_p->fails_quality_check)  bam_p->core.flag += BAM_FQCFAIL;
    if (alignment_p->pc_optical_duplicate) bam_p->core.flag += BAM_FDUP;

    //bin field requieres core
    bam_p->core.bin = bam_reg2bin(alignment_p->position, bam_calend(&bam_p->core, bam1_cigar(bam_p)));

    return bam_p;
}

//------------------------------------------------------------------------------------

bam_header_t* bam_header_new(int specie, int assembly, char* file_path) {
    bamFile bam_header_file;
    bam_header_t* bam_header_p;

    if ((specie == HUMAN) && (assembly == NCBI37)) {
        bam_header_file = bam_open(file_path, "r");
        bam_header_p = bam_header_read(bam_header_file);
    bam_close(bam_header_file);
    }

    return bam_header_p;
}

//------------------------------------------------------------------------------------

void bam_header_free(bam_header_t *header) {
    bam_header_destroy(header);
}

//------------------------------------------------------------------------------------

short int is_secondary_alignment(alignment_t *alignment) {
  return alignment->secondary_alignment;
}

//------------------------------------------------------------------------------------

void set_secondary_alignment(short int set, alignment_t *alignment) {
  alignment->secondary_alignment = set;
} 

/* **********************************************************************
 *      	Functions to manage bam1_t coded fields    		*
 * *********************************************************************/

char* convert_to_sequence_string(uint8_t* sequence_p, int sequence_length) {
  char* sequence_string = sequence_p;

  for (int i = 0; i < sequence_length; i++) {
    switch (bam1_seqi(sequence_p, i)) {
      case 1:
        sequence_string[i] = 'A';
        break;
      case 2:
        sequence_string[i] = 'C';
        break;
      case 4:
        sequence_string[i] = 'G';
        break;
      case 8:
        sequence_string[i] = 'T';
        break;
      case 15:
        sequence_string[i] = 'N';
        break;
    }
  }
  
  sequence_string[sequence_length] = '\0';
  return sequence_string;
}

//------------------------------------------------------------------------------------

char* convert_to_quality_string_length(char* quality_dest, uint8_t* quality_src, int quality_length, int base_quality) {
  for (int i = 0; i < quality_length; i++) {
    quality_dest[i] = quality_src[i] + base_quality;
  }

  quality_dest[quality_length] = '\0';
  return quality_dest;
}

//------------------------------------------------------------------------------------

char* convert_to_cigar_string(uint32_t* cigar_p, int num_cigar_operations) {
    //asumming not more than 3 digits per operation
    char* cigar_string = (char*) calloc(max(MIN_ALLOCATED_SIZE_FOR_CIGAR_STRING, 10 * num_cigar_operations), sizeof(char));
    uint32_t cigar_int;

    for (int i = 0; i < num_cigar_operations; i++) {
        cigar_int = cigar_p[i];

        switch (cigar_int&BAM_CIGAR_MASK) {
            case BAM_CMATCH:  //M: match or mismatch
                sprintf(cigar_string, "%s%u%c", cigar_string, cigar_int >> BAM_CIGAR_SHIFT, 'M');
                break;
            case BAM_CINS:  //I: insertion to the reference
                sprintf(cigar_string, "%s%u%c", cigar_string, cigar_int >> BAM_CIGAR_SHIFT, 'I');
                break;
            case BAM_CDEL:  //D: deletion from the reference
                sprintf(cigar_string, "%s%u%c", cigar_string, cigar_int >> BAM_CIGAR_SHIFT, 'D');
                break;
            case BAM_CREF_SKIP: //N: skip on the reference
                sprintf(cigar_string, "%s%u%c", cigar_string, cigar_int >> BAM_CIGAR_SHIFT, 'N');
                break;
            case BAM_CSOFT_CLIP: //S: clip on the read with clipped sequence
                sprintf(cigar_string, "%s%u%c", cigar_string, cigar_int >> BAM_CIGAR_SHIFT, 'S');
                break;
            case BAM_CHARD_CLIP: //H: clip on the read with clipped sequence trimmed off
                sprintf(cigar_string, "%s%u%c", cigar_string, cigar_int >> BAM_CIGAR_SHIFT, 'H');
                break;
            case BAM_CPAD:  //P: padding
                sprintf(cigar_string, "%s%u%c", cigar_string, cigar_int >> BAM_CIGAR_SHIFT, 'P');
                break;
            case BAM_CEQUAL:  //=: match
                sprintf(cigar_string, "%s%u%c", cigar_string, cigar_int >> BAM_CIGAR_SHIFT, '=');
                break;
            case BAM_CDIFF:  //X: mismatch
                sprintf(cigar_string, "%s%u%c", cigar_string, cigar_int >> BAM_CIGAR_SHIFT, 'X');
                break;
        }
    }

    sprintf(cigar_string, "%s%s", cigar_string, "\0");
    return cigar_string;
}

//------------------------------------------------------------------------------------

void convert_to_cigar_uint32_t(uint8_t* data, char* cigar, int num_cigar_operations) {
    int cigar_string_length = strlen(cigar);
    uint32_t cigar_uint32_position;
    int cigar_position, cigar_operation, cigar_acc_num_operations = 0;

    for (int i = 0; i < cigar_string_length; i++) {
        cigar_position = (int) cigar[i];

        if (cigar_position < 58) {  //numeric
            if (cigar_acc_num_operations == 0) cigar_acc_num_operations += atoi(&cigar[i]);
        } else {    //character: cigar operation
            switch (cigar_position) {
                case 77:  //M: match or mismatch
                    cigar_operation = BAM_CMATCH;
                    break;
                case 73:  //I: insertion to the reference
                    cigar_operation = BAM_CINS;
                    break;
                case 68:  //D: deletion from the reference
                    cigar_operation = BAM_CDEL;
                    break;
                case 78:  //N: skip on the reference
                    cigar_operation = BAM_CREF_SKIP;
                    break;
                case 83:  //S: clip on the read with clipped sequence
                    cigar_operation = BAM_CSOFT_CLIP;
                    break;
                case 72:  //H: clip on the read with clipped sequence trimmed off
                    cigar_operation = BAM_CHARD_CLIP;
                    break;
                case 80:  //P: padding
                    cigar_operation = BAM_CPAD;
                    break;
                case 61:  //=: match
                    cigar_operation = BAM_CEQUAL;
                    break;
                case 88:  //X: mismatch
                    cigar_operation = BAM_CDIFF;
                    break;
            }

            cigar_uint32_position = (uint32_t)((cigar_acc_num_operations << 4) + cigar_operation);

            memcpy(data, &cigar_uint32_position, 4);
            data += 4;

            cigar_acc_num_operations = 0;
        }
    }
}

//------------------------------------------------------------------------------------

void convert_to_sequence_uint8_t(uint8_t* data, char* sequence_p, int sequence_length) {
  uint8_t nts_uint8 = 0;
  
  int nt_int;
  for (int i = 0; i < sequence_length; i++) {
    nt_int = (int) sequence_p[i];
      
    switch (nt_int) {
      case 97:   // 'a'
      case 65:   // 'A'
    nts_uint8 = nts_uint8 + 1;
    break;
      case 116:  // 't'
      case 84:   // 'T'
        nts_uint8 = nts_uint8 + 8;
        break;
      case 99:    // 'c'
      case 67:   // 'C'
        nts_uint8 = nts_uint8 + 2;
        break;
      case 103:   // 'g'
      case 71:   // 'G'
    nts_uint8 = nts_uint8 + 4;
    break;
      case 110:   //'n'
      case 78:   // 'N'
    nts_uint8 = nts_uint8 + 15;
    break;
    }

    if (i & 1) { //even number
      *data++ = nts_uint8;
      nts_uint8 = 0;
    } else {
      nts_uint8 = nts_uint8 << 4;
    }
  }
  
  //if the last position is odd
  if (sequence_length & 1) {
    *data++ = nts_uint8;
  }
}

//------------------------------------------------------------------------------------

void convert_to_quality_uint8_t(uint8_t* data, char* quality_p, int quality_length, int base_quality) {
    for (int i = 0; i < quality_length; i++) {
        *data++ = (uint8_t)(quality_p[i] - base_quality);
    }
}

//------------------------------------------------------------------------------------

char select_op(unsigned char status){
  char operation; 
  switch(status){        
  case CIGAR_MATCH_MISMATCH: 
    operation = 'M';   
    break;             
  case CIGAR_DELETION:       
    operation = 'D';   
    break;             
  case CIGAR_INSERTION:      
    operation = 'I';   
    break;         
  case CIGAR_PERFECT_MATCH:
    operation = '=';
    break;
  case CIGAR_PADDING:
    operation = 'P';
    break;
  }                      
  return operation;       
}                      
