#include "gff_read.h"

//====================================================================================
// gff reading functions
//====================================================================================


/* ************ Header management functions **********************/

gff_header_entry_t* gff_header_entry_new() {
    gff_header_entry_t *entry = (gff_header_entry_t*) malloc (sizeof(gff_header_entry_t));
    entry->text = NULL;
    entry->text_len = 0;
    return entry;
}

void set_gff_header_entry_text(char *text, size_t length, gff_header_entry_t *entry) {
    assert(text);
    assert(entry);
    entry->text = strndup(text, length);
    entry->text_len = length;
    LOG_DEBUG_F("set text: %.*s\n", entry->text_len, entry->text);
}


/* ************ Record management functions **********************/

gff_record_t* gff_record_new() {
    return calloc (1, sizeof(gff_record_t));
}

void set_gff_record_sequence(char* sequence, size_t length, gff_record_t* gff_record) {
    assert(sequence);
    assert(gff_record);
    
    if (starts_with_n(sequence, "chrom", 5)) {
        gff_record->sequence = sequence + 5;
        gff_record->sequence_len = length - 5;
    } else if (starts_with_n(sequence, "chr", 3)) {
        gff_record->sequence = sequence + 3;
        gff_record->sequence_len = length - 3;
    } else {
        gff_record->sequence = sequence;
        gff_record->sequence_len = length;
    }
    LOG_DEBUG_F("set sequence: %.*s\n", gff_record->sequence_len, gff_record->sequence);
}

void set_gff_record_source(char* source, size_t length, gff_record_t* gff_record) {
    assert(source);
    assert(gff_record);
    gff_record->source = source;
    gff_record->source_len = length;
    LOG_DEBUG_F("set source: %.*s\n", gff_record->source_len, gff_record->source);
}

void set_gff_record_feature(char* feature, size_t length, gff_record_t* gff_record) {
    assert(feature);
    assert(gff_record);
    gff_record->feature = feature;
    gff_record->feature_len = length;
    LOG_DEBUG_F("set feature: %.*s\n", gff_record->feature_len, gff_record->feature);
}

void set_gff_record_start(unsigned long start, gff_record_t* gff_record) {
    gff_record->start = start;
    LOG_DEBUG_F("set start: %ld\n", gff_record->start);
}

void set_gff_record_end(unsigned long end, gff_record_t* gff_record) {
    gff_record->end = end;
    LOG_DEBUG_F("set end: %ld\n", gff_record->end);
}

void set_gff_record_score(float score, gff_record_t* gff_record) {
    gff_record->score = score;
    LOG_DEBUG_F("set score: %d\n", gff_record->score);
}

void set_gff_record_strand(char strand, gff_record_t* gff_record) {
    gff_record->strand = strand;
    LOG_DEBUG_F("set strand: %c\n", gff_record->strand);
}

void set_gff_record_frame(char frame, gff_record_t* gff_record) {
    gff_record->frame = frame;
    LOG_DEBUG_F("set frame: %c\n", gff_record->frame);
}

void set_gff_record_attribute(char* attribute, size_t length, gff_record_t* gff_record) {
    assert(attribute);
    assert(gff_record);
    gff_record->attribute = attribute;
    gff_record->attribute_len = length;
    LOG_DEBUG_F("set attribute: %.*s\n", gff_record->attribute_len, gff_record->attribute);
}
