#include "bed_read.h"

//====================================================================================
// bed reading functions
//====================================================================================


/* ************ Header management functions **********************/

bed_header_entry_t* bed_header_entry_new() {
    bed_header_entry_t *entry = (bed_header_entry_t*) malloc (sizeof(bed_header_entry_t));
    entry->text = NULL;
    entry->text_len = 0;
    return entry;
}

void set_bed_header_entry_text(char *text, size_t length, bed_header_entry_t *entry) {
    assert(text);
    assert(entry);
    entry->text = strndup(text, length);
    entry->text_len = length;
    LOG_DEBUG_F("set text: %.*s\n", entry->text_len, entry->text);
}


/* ************ Record management functions **********************/

bed_record_t* bed_record_new() {
    return calloc (1, sizeof(bed_record_t));
}

void set_bed_record_chromosome(char* chromosome, size_t length, bed_record_t* bed_record) {
    assert(chromosome);
    assert(bed_record);
    
    if (starts_with_n(chromosome, "chrom", 5)) {
        bed_record->sequence = chromosome + 5;
        bed_record->sequence_len = length - 5;
    } else if (starts_with_n(chromosome, "chr", 3)) {
        bed_record->sequence = chromosome + 3;
        bed_record->sequence_len = length - 3;
    } else {
        bed_record->sequence = chromosome;
        bed_record->sequence_len = length;
    }
    LOG_DEBUG_F("set sequence: %.*s\n", bed_record->sequence_len, bed_record->sequence);
}

void set_bed_record_start(unsigned long start, bed_record_t* bed_record) {
    bed_record->start = start;
    LOG_DEBUG_F("set start: %ld\n", bed_record->start);
}

void set_bed_record_end(unsigned long end, bed_record_t* bed_record) {
    bed_record->end = end;
    LOG_DEBUG_F("set end: %ld\n", bed_record->end);
}

void set_bed_record_name(char* name, size_t length, bed_record_t* bed_record) {
    assert(name);
    assert(bed_record);
    bed_record->name = name;
    bed_record->name_len = length;
    LOG_DEBUG_F("set name: %.*s\n", bed_record->name_len, bed_record->name);
}

void set_bed_record_score(float score, bed_record_t* bed_record) {
    bed_record->score = score;
    LOG_DEBUG_F("set score: %d\n", bed_record->score);
}

void set_bed_record_strand(char strand, bed_record_t* bed_record) {
    bed_record->strand = strand;
    LOG_DEBUG_F("set strand: %c\n", bed_record->strand);
}

void set_bed_record_thickstart(unsigned long thickstart, bed_record_t* bed_record) {
    bed_record->thickstart = thickstart;
    LOG_DEBUG_F("set thick start: %ld\n", bed_record->thickstart);
}

void set_bed_record_thickend(unsigned long thickend, bed_record_t* bed_record) {
    bed_record->thickend = thickend;
    LOG_DEBUG_F("set thick end: %ld\n", bed_record->thickend);
}

void set_bed_record_itemrgb(char* item_rgb, size_t length, bed_record_t* bed_record) {
    assert(item_rgb);
    assert(bed_record);
    bed_record->item_rgb = item_rgb;
    bed_record->item_rgb_len = length;
    LOG_DEBUG_F("set item_rgb: %.*s\n", bed_record->item_rgb_len, bed_record->item_rgb);
}

void set_bed_record_blockcount(int block_count, bed_record_t* bed_record) {
    bed_record->block_count = block_count;
    LOG_DEBUG_F("set start: %ld\n", bed_record->block_count);
}

void set_bed_record_blocksizes(char* block_sizes, size_t length, bed_record_t* bed_record) {
    assert(block_sizes);
    assert(bed_record);
    bed_record->block_sizes = block_sizes;
    bed_record->block_sizes_len = length;
    LOG_DEBUG_F("set block_sizes: %.*s\n", bed_record->block_sizes_len, bed_record->block_sizes);
}

void set_bed_record_blockstarts(char* block_starts, size_t length, bed_record_t* bed_record) {
    assert(block_starts);
    assert(bed_record);
    bed_record->block_starts = block_starts;
    bed_record->block_starts_len = length;
    LOG_DEBUG_F("set block_starts: %.*s\n", bed_record->block_starts_len, bed_record->block_starts);
}
