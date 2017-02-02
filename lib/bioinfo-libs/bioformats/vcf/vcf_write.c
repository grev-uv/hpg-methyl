#include "vcf_write.h"

int write_vcf_file(vcf_file_t *file, FILE *fd) {
    assert(file);
    assert(fd);
    
    // Write header: file format, header entries and delimiter
    if (write_vcf_header(file, fd) > 0) {
        return 1;
    }
    
    // Write records: grouped in batches
    vcf_batch_t *batch;
    while ((batch = fetch_vcf_batch(file)) != NULL) {
        if (write_vcf_batch(batch, fd) > 0) {
            return 1;
        }
    }
    
    return 0;
}

int write_vcf_header(vcf_file_t* file, FILE* fd) {
    assert(file);
    assert(fd);
    
    // Write header entries (including fileformat)
    for (int i = 0; i < file->header_entries->size; i++) {
        if (write_vcf_header_entry((vcf_header_entry_t*) array_list_get(i, file->header_entries), fd) > 0) {
            return 1;
        }
    }
    
    // Write delimiter
    if (write_vcf_delimiter(file, fd) > 0) {
        return 1;
    }
    
    return 0;
}


int write_vcf_fileformat(vcf_file_t *file, FILE *fd) {
    assert(file);
    assert(fd);
    
    if (fprintf(fd, "##fileformat=%s\n", file->format) < 0) {
        return 1;
    }
    
    return 0;
}

int write_vcf_header_entry(vcf_header_entry_t *entry, FILE *fd) {
    assert(entry);
    assert(fd);
    
    if (fprintf(fd, "##") < 0) {
        return 1;
    }

    
    // Entries with the form ##value
    if (entry->name == NULL) {
        char *value = array_list_get(0, entry->values); // Just first (and only) value
        if (fprintf(fd, "%s\n", value) < 0) {
            return 1;
        }
        return 0;
    }
    
    // Name for entries with the form ##name=...
    if (fprintf(fd, "%s", entry->name) < 0) {
        return 1;
    }
    
    // Entries with the form ##name=value
    if (entry->values->size == 1) {
        char *value = array_list_get(0, entry->values);
        if (fprintf(fd, "=%s\n", value) < 0) {
            return 1;
        }
    }
    
    // Entries with the form ##name=<field_id=value,field_id=value,...>
    else if (entry->values->size > 1) {
        if (fprintf(fd, "=<") < 0) {
            return 1;
        }
        
        for (int i = 0; i < entry->values->size; i++) {
            char *value = array_list_get(i, entry->values);
            if (fprintf(fd, "%s", value) < 0) {
                return 1;
            }
            if (i < entry->values->size - 1) {
                if (fprintf(fd, ",") < 0) {
                    return 1;
                }
            }
        }
        
        if (fprintf(fd, ">\n") < 0) {
            return 1;
        }
    }
    
    return 0;
}

int write_vcf_delimiter(vcf_file_t *file, FILE *fd) {
    assert(file);
    assert(fd);
    
    if (write_vcf_delimiter_from_samples(file->samples_names->items, file->samples_names->size, fd) > 0) {
        return 1;
    }
//     if (fprintf(fd, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT") < 0) {
//         return 1;
//     }
//     
//     for (int i = 0; i < file->samples_names->size; i++) {
//         if (fprintf(fd, "\t%s", (char*) array_list_get(i, file->samples_names)) < 0) {
//             return 1;
//         }
//     }
//        
//     if (fprintf(fd, "\n") < 0) {
//         return 1;
//     }
    
    return 0;
}

int write_vcf_delimiter_from_samples(char **sample_names, int num_samples, FILE *fd) {
    assert(sample_names);
    assert(fd);
    
    if (fprintf(fd, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT") < 0) {
        return 1;
    }
    
    for (int i = 0; i < num_samples; i++) {
        if (fprintf(fd, "\t%s", sample_names[i]) < 0) {
            return 1;
        }
    }
       
    if (fprintf(fd, "\n") < 0) {
        return 1;
    }
    
    return 0;
}


int write_vcf_batch(vcf_batch_t *batch, FILE *fd) {
    assert(batch);
    assert(fd);
    
    vcf_record_t *record;
    for (int i = 0; i < batch->records->size; i++) {
        record = array_list_get(i, batch->records);
        if (write_vcf_record(record, fd) > 0) {
            return 1;
        }
    }
    
    return 0;
}

int write_vcf_record(vcf_record_t* record, FILE *fd) {
    assert(record);
    assert(fd);
    
    if (fprintf(fd, "%.*s\t%ld\t%.*s\t%.*s\t%.*s\t", record->chromosome_len, record->chromosome, record->position, 
                record->id_len, record->id, record->reference_len, record->reference, record->alternate_len, record->alternate) < 0) {
        return 1;
    }
    
    if (record->quality < 0) {
        if (fprintf(fd, ".\t") < 0) {
            return 1;
        }
    } else {
        if (fprintf(fd, "%.2f\t", record->quality) < 0) {
            return 1;
        }
    }
    if (fprintf(fd, "%.*s\t%.*s\t%.*s", record->filter_len, record->filter, record->info_len, record->info, record->format_len, record->format) < 0) {
        return 1;
    }

    for (int i = 0; i < record->samples->size; i++) {
        if (fprintf(fd, "\t%s", (char*) array_list_get(i, record->samples)) < 0) {
            return 1;
        }
    }

    if (fprintf(fd, "\n") < 0) {
        return 1;
    }
    
    return 0;
}
