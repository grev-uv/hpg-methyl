
#line 1 "bed.ragel"
#include "bed_reader.h"

static size_t lines = 1;
static size_t num_records = 0;
static size_t num_batches = 0;

static bed_record_t *current_record;
static bed_header_entry_t *current_header_entry;
static bed_batch_t *current_batch;


#line 15 "bed_reader.c"
static const int bed_start = 35;
static const int bed_first_final = 35;
static const int bed_error = 0;

static const int bed_en_main = 35;


#line 273 "bed.ragel"



int bed_ragel_read(list_t *batches_list, size_t batch_size, bed_file_t *file) {
    int cs;
    char *p = file->data;
    char *pe = p + file->data_len;
    char *eof = pe;
    char *ts, *te;
    int stack[4];
    int top, act;

    current_header_entry = bed_header_entry_new();
    current_batch = bed_batch_new(batch_size);

    
#line 40 "bed_reader.c"
	{
	cs = bed_start;
	}

#line 45 "bed_reader.c"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
case 35:
	switch( (*p) ) {
		case 10: goto tr52;
		case 99: goto tr53;
		case 115: goto tr54;
	}
	if ( 0 <= (*p) )
		goto tr51;
	goto tr3;
tr0:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
	goto st0;
tr3:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 71 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'chrom' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr8:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 85 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'start' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr13:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 99 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'end' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr25:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 141 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'strand' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr27:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr30:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr33:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 195 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockCount' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr35:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 181 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'itemRgb' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr38:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 181 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'itemRgb' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 195 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockCount' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr40:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 169 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickEnd' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 181 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'itemRgb' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 195 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockCount' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr42:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 155 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickStart' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 169 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickEnd' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 181 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'itemRgb' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 195 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockCount' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr44:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 129 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'score' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 155 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickStart' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 169 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickEnd' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 181 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'itemRgb' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 195 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockCount' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr46:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 111 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'name' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr56:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 111 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'name' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 129 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'score' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 141 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'strand' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 155 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickStart' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 169 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickEnd' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 181 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'itemRgb' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 195 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockCount' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr98:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 129 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'score' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 141 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'strand' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 155 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickStart' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 169 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickEnd' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 181 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'itemRgb' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 195 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockCount' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr102:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 141 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'strand' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 155 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickStart' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 169 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickEnd' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 181 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'itemRgb' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 195 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockCount' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr108:
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 111 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'name' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 129 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'score' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 155 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickStart' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 169 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickEnd' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 181 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'itemRgb' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 195 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockCount' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
#line 440 "bed_reader.c"
st0:
cs = 0;
	goto _out;
tr51:
#line 24 "bed.ragel"
	{
        current_header_entry = bed_header_entry_new();
        ts = p;
    }
	goto st1;
st1:
	if ( ++p == pe )
		goto _test_eof1;
case 1:
#line 455 "bed_reader.c"
	if ( (*p) == 10 )
		goto tr2;
	if ( 0 <= (*p) )
		goto st1;
	goto tr0;
tr2:
#line 29 "bed.ragel"
	{
        set_bed_header_entry_text(ts, p-ts, current_header_entry);
        add_bed_header_entry(current_header_entry, file);
    }
#line 19 "bed.ragel"
	{
        lines++;
//        printf("lines read = %d\n", lines);
    }
	goto st36;
tr55:
#line 29 "bed.ragel"
	{
        set_bed_header_entry_text(ts, p-ts, current_header_entry);
        add_bed_header_entry(current_header_entry, file);
    }
#line 19 "bed.ragel"
	{
        lines++;
//        printf("lines read = %d\n", lines);
    }
#line 24 "bed.ragel"
	{
        current_header_entry = bed_header_entry_new();
        ts = p;
    }
	goto st36;
st36:
	if ( ++p == pe )
		goto _test_eof36;
case 36:
#line 494 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr55;
		case 99: goto tr53;
		case 115: goto tr54;
	}
	if ( 0 <= (*p) )
		goto tr51;
	goto tr3;
tr53:
#line 24 "bed.ragel"
	{
        current_header_entry = bed_header_entry_new();
        ts = p;
    }
#line 38 "bed.ragel"
	{
        current_record = bed_record_new();
    }
#line 63 "bed.ragel"
	{
        ts = p;
    }
	goto st2;
tr111:
#line 38 "bed.ragel"
	{
        current_record = bed_record_new();
    }
#line 63 "bed.ragel"
	{
        ts = p;
    }
	goto st2;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
#line 532 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr2;
		case 104: goto st3;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
	switch( (*p) ) {
		case 10: goto tr2;
		case 114: goto st4;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
	switch( (*p) ) {
		case 10: goto tr2;
		case 95: goto st5;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 0 <= (*p) && (*p) <= 47 )
				goto st1;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto st1;
		} else
			goto st5;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto st1;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) )
				goto st1;
		} else
			goto st5;
	} else
		goto st5;
	goto tr3;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	switch( (*p) ) {
		case 9: goto tr7;
		case 10: goto tr2;
		case 95: goto st5;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 0 <= (*p) && (*p) <= 47 )
				goto st1;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto st1;
		} else
			goto st5;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto st1;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) )
				goto st1;
		} else
			goto st5;
	} else
		goto st5;
	goto tr3;
tr7:
#line 67 "bed.ragel"
	{
        set_bed_record_chromosome(ts, p-ts, current_record);
    }
	goto st6;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
#line 620 "bed_reader.c"
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto tr9;
	goto tr8;
tr9:
#line 75 "bed.ragel"
	{
        ts = p;
    }
	goto st7;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
#line 642 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr10;
		case 10: goto tr2;
		case 46: goto st33;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st7;
	goto tr8;
tr10:
#line 79 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_start(atol(field), current_record);
        free(field);
    }
	goto st8;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
#line 669 "bed_reader.c"
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto tr14;
	goto tr13;
tr14:
#line 89 "bed.ragel"
	{
        ts = p;
    }
	goto st9;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
#line 691 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr15;
		case 10: goto tr2;
		case 46: goto st31;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st9;
	goto tr13;
tr15:
#line 93 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_end(atol(field), current_record);
        free(field);
    }
	goto st37;
st37:
	if ( ++p == pe )
		goto _test_eof37;
case 37:
#line 718 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr57;
		case 43: goto tr58;
		case 45: goto tr58;
		case 46: goto tr59;
		case 95: goto tr61;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 0 <= (*p) && (*p) <= 47 )
				goto st1;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto st1;
		} else
			goto tr60;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto st1;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) )
				goto st1;
		} else
			goto tr61;
	} else
		goto tr61;
	goto tr56;
tr63:
#line 29 "bed.ragel"
	{
        set_bed_header_entry_text(ts, p-ts, current_header_entry);
        add_bed_header_entry(current_header_entry, file);
    }
#line 19 "bed.ragel"
	{
        lines++;
//        printf("lines read = %d\n", lines);
    }
	goto st38;
tr57:
#line 29 "bed.ragel"
	{
        set_bed_header_entry_text(ts, p-ts, current_header_entry);
        add_bed_header_entry(current_header_entry, file);
    }
#line 19 "bed.ragel"
	{
        lines++;
//        printf("lines read = %d\n", lines);
    }
#line 42 "bed.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (batch_size > 0 && current_batch->records->size == batch_size) {
            list_item_t *item = list_item_new(num_records, 1, current_batch); 
            list_insert_item(item, batches_list);
//             printf("Batch %d added - %zu records\t", batches, current_batch->records->size);
            current_batch = bed_batch_new(batch_size);
            
            if (p+1) {
                current_batch->text = p+1;
//                 printf("batch text = '%.*s'\n", 50, current_batch->text);
            }
            num_batches++;
            num_records = 0;
        }

        // If not a blank line, add current record to current batch
        add_record_to_bed_batch(current_record, current_batch);
        num_records++;
    }
	goto st38;
tr67:
#line 29 "bed.ragel"
	{
        set_bed_header_entry_text(ts, p-ts, current_header_entry);
        add_bed_header_entry(current_header_entry, file);
    }
#line 19 "bed.ragel"
	{
        lines++;
//        printf("lines read = %d\n", lines);
    }
#line 215 "bed.ragel"
	{
        set_bed_record_blockstarts(ts, p-ts, current_record);
    }
#line 42 "bed.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (batch_size > 0 && current_batch->records->size == batch_size) {
            list_item_t *item = list_item_new(num_records, 1, current_batch); 
            list_insert_item(item, batches_list);
//             printf("Batch %d added - %zu records\t", batches, current_batch->records->size);
            current_batch = bed_batch_new(batch_size);
            
            if (p+1) {
                current_batch->text = p+1;
//                 printf("batch text = '%.*s'\n", 50, current_batch->text);
            }
            num_batches++;
            num_records = 0;
        }

        // If not a blank line, add current record to current batch
        add_record_to_bed_batch(current_record, current_batch);
        num_records++;
    }
	goto st38;
st38:
	if ( ++p == pe )
		goto _test_eof38;
case 38:
#line 833 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr55;
		case 35: goto tr62;
		case 99: goto tr53;
		case 115: goto tr54;
	}
	if ( 0 <= (*p) )
		goto tr51;
	goto tr3;
tr62:
#line 24 "bed.ragel"
	{
        current_header_entry = bed_header_entry_new();
        ts = p;
    }
	goto st10;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
#line 854 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr2;
		case 127: goto st1;
	}
	if ( (*p) > 31 ) {
		if ( 32 <= (*p) && (*p) <= 126 )
			goto st39;
	} else if ( (*p) >= 0 )
		goto st1;
	goto tr0;
st39:
	if ( ++p == pe )
		goto _test_eof39;
case 39:
	switch( (*p) ) {
		case 10: goto tr63;
		case 127: goto st1;
	}
	if ( (*p) > 31 ) {
		if ( 32 <= (*p) && (*p) <= 126 )
			goto st39;
	} else if ( (*p) >= 0 )
		goto st1;
	goto tr0;
tr54:
#line 24 "bed.ragel"
	{
        current_header_entry = bed_header_entry_new();
        ts = p;
    }
#line 38 "bed.ragel"
	{
        current_record = bed_record_new();
    }
#line 63 "bed.ragel"
	{
        ts = p;
    }
	goto st11;
tr112:
#line 38 "bed.ragel"
	{
        current_record = bed_record_new();
    }
#line 63 "bed.ragel"
	{
        ts = p;
    }
	goto st11;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
#line 908 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr2;
		case 99: goto st12;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
	switch( (*p) ) {
		case 10: goto tr2;
		case 97: goto st13;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
	switch( (*p) ) {
		case 10: goto tr2;
		case 102: goto st14;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
	switch( (*p) ) {
		case 10: goto tr2;
		case 102: goto st15;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
	switch( (*p) ) {
		case 10: goto tr2;
		case 111: goto st16;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
	switch( (*p) ) {
		case 10: goto tr2;
		case 108: goto st17;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
	switch( (*p) ) {
		case 10: goto tr2;
		case 100: goto st4;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
tr58:
#line 133 "bed.ragel"
	{
        ts = p;
    }
	goto st18;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
#line 992 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr26;
		case 10: goto tr2;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr25;
tr26:
#line 137 "bed.ragel"
	{
        set_bed_record_strand(*ts, current_record);
    }
	goto st40;
tr104:
#line 137 "bed.ragel"
	{
        set_bed_record_strand(*ts, current_record);
    }
#line 149 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_thickstart(atol(field), current_record);
        free(field);
    }
#line 163 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_thickend(atol(field), current_record);
        free(field);
    }
#line 177 "bed.ragel"
	{
        set_bed_record_itemrgb(ts, p-ts, current_record);
    }
#line 189 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_blockcount(atoi(field), current_record);
        free(field);
    }
#line 203 "bed.ragel"
	{
        set_bed_record_blocksizes(ts, p-ts, current_record);
    }
	goto st40;
st40:
	if ( ++p == pe )
		goto _test_eof40;
case 40:
#line 1042 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr57;
		case 46: goto tr64;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto tr65;
	goto tr42;
tr64:
#line 145 "bed.ragel"
	{
        ts = p;
    }
#line 159 "bed.ragel"
	{
        ts = p;
    }
#line 173 "bed.ragel"
	{
        ts = p;
    }
#line 185 "bed.ragel"
	{
        ts = p;
    }
#line 199 "bed.ragel"
	{
        ts = p;
    }
#line 211 "bed.ragel"
	{
        ts = p;
    }
	goto st41;
st41:
	if ( ++p == pe )
		goto _test_eof41;
case 41:
#line 1086 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr66;
		case 10: goto tr67;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr42;
tr66:
#line 149 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_thickstart(atol(field), current_record);
        free(field);
    }
#line 163 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_thickend(atol(field), current_record);
        free(field);
    }
#line 177 "bed.ragel"
	{
        set_bed_record_itemrgb(ts, p-ts, current_record);
    }
#line 189 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_blockcount(atoi(field), current_record);
        free(field);
    }
#line 203 "bed.ragel"
	{
        set_bed_record_blocksizes(ts, p-ts, current_record);
    }
	goto st42;
st42:
	if ( ++p == pe )
		goto _test_eof42;
case 42:
#line 1126 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr57;
		case 46: goto tr68;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto tr69;
	goto tr40;
tr68:
#line 159 "bed.ragel"
	{
        ts = p;
    }
#line 173 "bed.ragel"
	{
        ts = p;
    }
#line 185 "bed.ragel"
	{
        ts = p;
    }
#line 199 "bed.ragel"
	{
        ts = p;
    }
#line 211 "bed.ragel"
	{
        ts = p;
    }
	goto st43;
st43:
	if ( ++p == pe )
		goto _test_eof43;
case 43:
#line 1166 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr70;
		case 10: goto tr67;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr40;
tr70:
#line 163 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_thickend(atol(field), current_record);
        free(field);
    }
#line 177 "bed.ragel"
	{
        set_bed_record_itemrgb(ts, p-ts, current_record);
    }
#line 189 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_blockcount(atoi(field), current_record);
        free(field);
    }
#line 203 "bed.ragel"
	{
        set_bed_record_blocksizes(ts, p-ts, current_record);
    }
	goto st44;
st44:
	if ( ++p == pe )
		goto _test_eof44;
case 44:
#line 1200 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr57;
		case 46: goto tr71;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto tr72;
	goto tr38;
tr71:
#line 173 "bed.ragel"
	{
        ts = p;
    }
#line 185 "bed.ragel"
	{
        ts = p;
    }
#line 199 "bed.ragel"
	{
        ts = p;
    }
#line 211 "bed.ragel"
	{
        ts = p;
    }
	goto st45;
st45:
	if ( ++p == pe )
		goto _test_eof45;
case 45:
#line 1236 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr73;
		case 10: goto tr67;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr38;
tr73:
#line 177 "bed.ragel"
	{
        set_bed_record_itemrgb(ts, p-ts, current_record);
    }
#line 189 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_blockcount(atoi(field), current_record);
        free(field);
    }
#line 203 "bed.ragel"
	{
        set_bed_record_blocksizes(ts, p-ts, current_record);
    }
	goto st46;
tr91:
#line 177 "bed.ragel"
	{
        set_bed_record_itemrgb(ts, p-ts, current_record);
    }
#line 203 "bed.ragel"
	{
        set_bed_record_blocksizes(ts, p-ts, current_record);
    }
	goto st46;
st46:
	if ( ++p == pe )
		goto _test_eof46;
case 46:
#line 1274 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr57;
		case 46: goto tr74;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto tr75;
	goto tr33;
tr74:
#line 185 "bed.ragel"
	{
        ts = p;
    }
#line 199 "bed.ragel"
	{
        ts = p;
    }
#line 211 "bed.ragel"
	{
        ts = p;
    }
	goto st47;
st47:
	if ( ++p == pe )
		goto _test_eof47;
case 47:
#line 1306 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr76;
		case 10: goto tr67;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr33;
tr76:
#line 189 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_blockcount(atoi(field), current_record);
        free(field);
    }
#line 203 "bed.ragel"
	{
        set_bed_record_blocksizes(ts, p-ts, current_record);
    }
	goto st48;
st48:
	if ( ++p == pe )
		goto _test_eof48;
case 48:
#line 1330 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr57;
		case 46: goto tr77;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto tr78;
	goto tr30;
tr77:
#line 199 "bed.ragel"
	{
        ts = p;
    }
#line 211 "bed.ragel"
	{
        ts = p;
    }
	goto st49;
st49:
	if ( ++p == pe )
		goto _test_eof49;
case 49:
#line 1358 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr79;
		case 10: goto tr67;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr30;
tr79:
#line 203 "bed.ragel"
	{
        set_bed_record_blocksizes(ts, p-ts, current_record);
    }
	goto st50;
st50:
	if ( ++p == pe )
		goto _test_eof50;
case 50:
#line 1376 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr57;
		case 46: goto tr80;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto tr81;
	goto tr27;
tr80:
#line 211 "bed.ragel"
	{
        ts = p;
    }
	goto st51;
st51:
	if ( ++p == pe )
		goto _test_eof51;
case 51:
#line 1400 "bed_reader.c"
	if ( (*p) == 10 )
		goto tr67;
	if ( 0 <= (*p) )
		goto st1;
	goto tr27;
tr81:
#line 211 "bed.ragel"
	{
        ts = p;
    }
	goto st52;
st52:
	if ( ++p == pe )
		goto _test_eof52;
case 52:
#line 1416 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr67;
		case 44: goto st19;
		case 46: goto st20;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st52;
	goto tr27;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st52;
	goto tr27;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st53;
	goto tr27;
st53:
	if ( ++p == pe )
		goto _test_eof53;
case 53:
	switch( (*p) ) {
		case 10: goto tr67;
		case 44: goto st19;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st53;
	goto tr27;
tr78:
#line 199 "bed.ragel"
	{
        ts = p;
    }
#line 211 "bed.ragel"
	{
        ts = p;
    }
	goto st54;
st54:
	if ( ++p == pe )
		goto _test_eof54;
case 54:
#line 1492 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr79;
		case 10: goto tr67;
		case 44: goto st21;
		case 46: goto st22;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st54;
	goto tr30;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st54;
	goto tr30;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st55;
	goto tr30;
st55:
	if ( ++p == pe )
		goto _test_eof55;
case 55:
	switch( (*p) ) {
		case 9: goto tr79;
		case 10: goto tr67;
		case 44: goto st21;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st55;
	goto tr30;
tr75:
#line 185 "bed.ragel"
	{
        ts = p;
    }
#line 199 "bed.ragel"
	{
        ts = p;
    }
#line 211 "bed.ragel"
	{
        ts = p;
    }
	goto st56;
st56:
	if ( ++p == pe )
		goto _test_eof56;
case 56:
#line 1574 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr76;
		case 10: goto tr67;
		case 44: goto st21;
		case 46: goto st23;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st56;
	goto tr33;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st57;
	goto tr33;
st57:
	if ( ++p == pe )
		goto _test_eof57;
case 57:
	switch( (*p) ) {
		case 9: goto tr76;
		case 10: goto tr67;
		case 44: goto st21;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st57;
	goto tr33;
tr72:
#line 173 "bed.ragel"
	{
        ts = p;
    }
#line 185 "bed.ragel"
	{
        ts = p;
    }
#line 199 "bed.ragel"
	{
        ts = p;
    }
#line 211 "bed.ragel"
	{
        ts = p;
    }
	goto st58;
st58:
	if ( ++p == pe )
		goto _test_eof58;
case 58:
#line 1645 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr73;
		case 10: goto tr67;
		case 44: goto st24;
		case 46: goto st26;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st58;
	goto tr38;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st59;
	goto tr35;
st59:
	if ( ++p == pe )
		goto _test_eof59;
case 59:
	switch( (*p) ) {
		case 9: goto tr91;
		case 10: goto tr67;
		case 44: goto st24;
		case 46: goto st25;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st59;
	goto tr35;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st60;
	goto tr35;
st60:
	if ( ++p == pe )
		goto _test_eof60;
case 60:
	switch( (*p) ) {
		case 9: goto tr91;
		case 10: goto tr67;
		case 44: goto st24;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st60;
	goto tr35;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st61;
	goto tr38;
st61:
	if ( ++p == pe )
		goto _test_eof61;
case 61:
	switch( (*p) ) {
		case 9: goto tr73;
		case 10: goto tr67;
		case 44: goto st24;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st61;
	goto tr38;
tr69:
#line 159 "bed.ragel"
	{
        ts = p;
    }
#line 173 "bed.ragel"
	{
        ts = p;
    }
#line 185 "bed.ragel"
	{
        ts = p;
    }
#line 199 "bed.ragel"
	{
        ts = p;
    }
#line 211 "bed.ragel"
	{
        ts = p;
    }
	goto st62;
st62:
	if ( ++p == pe )
		goto _test_eof62;
case 62:
#line 1787 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr70;
		case 10: goto tr67;
		case 44: goto st24;
		case 46: goto st27;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st62;
	goto tr40;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st63;
	goto tr40;
st63:
	if ( ++p == pe )
		goto _test_eof63;
case 63:
	switch( (*p) ) {
		case 9: goto tr70;
		case 10: goto tr67;
		case 44: goto st24;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st63;
	goto tr40;
tr65:
#line 145 "bed.ragel"
	{
        ts = p;
    }
#line 159 "bed.ragel"
	{
        ts = p;
    }
#line 173 "bed.ragel"
	{
        ts = p;
    }
#line 185 "bed.ragel"
	{
        ts = p;
    }
#line 199 "bed.ragel"
	{
        ts = p;
    }
#line 211 "bed.ragel"
	{
        ts = p;
    }
	goto st64;
st64:
	if ( ++p == pe )
		goto _test_eof64;
case 64:
#line 1866 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr66;
		case 10: goto tr67;
		case 44: goto st24;
		case 46: goto st28;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st64;
	goto tr42;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st65;
	goto tr42;
st65:
	if ( ++p == pe )
		goto _test_eof65;
case 65:
	switch( (*p) ) {
		case 9: goto tr66;
		case 10: goto tr67;
		case 44: goto st24;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st65;
	goto tr42;
tr59:
#line 103 "bed.ragel"
	{
        ts = p;
    }
#line 115 "bed.ragel"
	{
        ts = p;
    }
#line 133 "bed.ragel"
	{
        ts = p;
    }
#line 145 "bed.ragel"
	{
        ts = p;
    }
#line 159 "bed.ragel"
	{
        ts = p;
    }
#line 173 "bed.ragel"
	{
        ts = p;
    }
#line 185 "bed.ragel"
	{
        ts = p;
    }
#line 199 "bed.ragel"
	{
        ts = p;
    }
#line 211 "bed.ragel"
	{
        ts = p;
    }
	goto st66;
st66:
	if ( ++p == pe )
		goto _test_eof66;
case 66:
#line 1957 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr97;
		case 10: goto tr67;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr56;
tr47:
#line 107 "bed.ragel"
	{
        set_bed_record_name(ts, p-ts, current_record);
    }
	goto st67;
tr97:
#line 107 "bed.ragel"
	{
        set_bed_record_name(ts, p-ts, current_record);
    }
#line 119 "bed.ragel"
	{
        float score = -1.0f;
        if (strncmp(".", ts, 1) != 0) {
            char *field = strndup(ts, p-ts);
            score = atof(field);
            free(field);
        }
        set_bed_record_score(score, current_record);
    }
#line 137 "bed.ragel"
	{
        set_bed_record_strand(*ts, current_record);
    }
#line 149 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_thickstart(atol(field), current_record);
        free(field);
    }
#line 163 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_thickend(atol(field), current_record);
        free(field);
    }
#line 177 "bed.ragel"
	{
        set_bed_record_itemrgb(ts, p-ts, current_record);
    }
#line 189 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_blockcount(atoi(field), current_record);
        free(field);
    }
#line 203 "bed.ragel"
	{
        set_bed_record_blocksizes(ts, p-ts, current_record);
    }
	goto st67;
tr109:
#line 107 "bed.ragel"
	{
        set_bed_record_name(ts, p-ts, current_record);
    }
#line 119 "bed.ragel"
	{
        float score = -1.0f;
        if (strncmp(".", ts, 1) != 0) {
            char *field = strndup(ts, p-ts);
            score = atof(field);
            free(field);
        }
        set_bed_record_score(score, current_record);
    }
#line 149 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_thickstart(atol(field), current_record);
        free(field);
    }
#line 163 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_thickend(atol(field), current_record);
        free(field);
    }
#line 177 "bed.ragel"
	{
        set_bed_record_itemrgb(ts, p-ts, current_record);
    }
#line 189 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_blockcount(atoi(field), current_record);
        free(field);
    }
#line 203 "bed.ragel"
	{
        set_bed_record_blocksizes(ts, p-ts, current_record);
    }
	goto st67;
st67:
	if ( ++p == pe )
		goto _test_eof67;
case 67:
#line 2063 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr57;
		case 43: goto tr58;
		case 45: goto tr58;
		case 46: goto tr99;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto tr100;
	goto tr98;
tr99:
#line 115 "bed.ragel"
	{
        ts = p;
    }
#line 133 "bed.ragel"
	{
        ts = p;
    }
#line 145 "bed.ragel"
	{
        ts = p;
    }
#line 159 "bed.ragel"
	{
        ts = p;
    }
#line 173 "bed.ragel"
	{
        ts = p;
    }
#line 185 "bed.ragel"
	{
        ts = p;
    }
#line 199 "bed.ragel"
	{
        ts = p;
    }
#line 211 "bed.ragel"
	{
        ts = p;
    }
	goto st68;
st68:
	if ( ++p == pe )
		goto _test_eof68;
case 68:
#line 2117 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr101;
		case 10: goto tr67;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr98;
tr101:
#line 119 "bed.ragel"
	{
        float score = -1.0f;
        if (strncmp(".", ts, 1) != 0) {
            char *field = strndup(ts, p-ts);
            score = atof(field);
            free(field);
        }
        set_bed_record_score(score, current_record);
    }
#line 137 "bed.ragel"
	{
        set_bed_record_strand(*ts, current_record);
    }
#line 149 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_thickstart(atol(field), current_record);
        free(field);
    }
#line 163 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_thickend(atol(field), current_record);
        free(field);
    }
#line 177 "bed.ragel"
	{
        set_bed_record_itemrgb(ts, p-ts, current_record);
    }
#line 189 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_blockcount(atoi(field), current_record);
        free(field);
    }
#line 203 "bed.ragel"
	{
        set_bed_record_blocksizes(ts, p-ts, current_record);
    }
	goto st69;
tr105:
#line 119 "bed.ragel"
	{
        float score = -1.0f;
        if (strncmp(".", ts, 1) != 0) {
            char *field = strndup(ts, p-ts);
            score = atof(field);
            free(field);
        }
        set_bed_record_score(score, current_record);
    }
#line 149 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_thickstart(atol(field), current_record);
        free(field);
    }
#line 163 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_thickend(atol(field), current_record);
        free(field);
    }
#line 177 "bed.ragel"
	{
        set_bed_record_itemrgb(ts, p-ts, current_record);
    }
#line 189 "bed.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_bed_record_blockcount(atoi(field), current_record);
        free(field);
    }
#line 203 "bed.ragel"
	{
        set_bed_record_blocksizes(ts, p-ts, current_record);
    }
	goto st69;
st69:
	if ( ++p == pe )
		goto _test_eof69;
case 69:
#line 2209 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr57;
		case 43: goto tr58;
		case 45: goto tr58;
		case 46: goto tr103;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto tr65;
	goto tr102;
tr103:
#line 133 "bed.ragel"
	{
        ts = p;
    }
#line 145 "bed.ragel"
	{
        ts = p;
    }
#line 159 "bed.ragel"
	{
        ts = p;
    }
#line 173 "bed.ragel"
	{
        ts = p;
    }
#line 185 "bed.ragel"
	{
        ts = p;
    }
#line 199 "bed.ragel"
	{
        ts = p;
    }
#line 211 "bed.ragel"
	{
        ts = p;
    }
	goto st70;
st70:
	if ( ++p == pe )
		goto _test_eof70;
case 70:
#line 2259 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr104;
		case 10: goto tr67;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr102;
tr100:
#line 115 "bed.ragel"
	{
        ts = p;
    }
#line 145 "bed.ragel"
	{
        ts = p;
    }
#line 159 "bed.ragel"
	{
        ts = p;
    }
#line 173 "bed.ragel"
	{
        ts = p;
    }
#line 185 "bed.ragel"
	{
        ts = p;
    }
#line 199 "bed.ragel"
	{
        ts = p;
    }
#line 211 "bed.ragel"
	{
        ts = p;
    }
	goto st71;
st71:
	if ( ++p == pe )
		goto _test_eof71;
case 71:
#line 2301 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr105;
		case 10: goto tr67;
		case 44: goto st24;
		case 46: goto st29;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st71;
	goto tr44;
st29:
	if ( ++p == pe )
		goto _test_eof29;
case 29:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st72;
	goto tr44;
st72:
	if ( ++p == pe )
		goto _test_eof72;
case 72:
	switch( (*p) ) {
		case 9: goto tr105;
		case 10: goto tr67;
		case 44: goto st24;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st72;
	goto tr44;
tr60:
#line 103 "bed.ragel"
	{
        ts = p;
    }
#line 115 "bed.ragel"
	{
        ts = p;
    }
#line 145 "bed.ragel"
	{
        ts = p;
    }
#line 159 "bed.ragel"
	{
        ts = p;
    }
#line 173 "bed.ragel"
	{
        ts = p;
    }
#line 185 "bed.ragel"
	{
        ts = p;
    }
#line 199 "bed.ragel"
	{
        ts = p;
    }
#line 211 "bed.ragel"
	{
        ts = p;
    }
	goto st73;
st73:
	if ( ++p == pe )
		goto _test_eof73;
case 73:
#line 2388 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr109;
		case 10: goto tr67;
		case 44: goto st24;
		case 46: goto st29;
		case 95: goto st30;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 0 <= (*p) && (*p) <= 47 )
				goto st1;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto st1;
		} else
			goto st73;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto st1;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) )
				goto st1;
		} else
			goto st30;
	} else
		goto st30;
	goto tr108;
tr61:
#line 103 "bed.ragel"
	{
        ts = p;
    }
	goto st30;
st30:
	if ( ++p == pe )
		goto _test_eof30;
case 30:
#line 2427 "bed_reader.c"
	switch( (*p) ) {
		case 9: goto tr47;
		case 10: goto tr2;
		case 95: goto st30;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 0 <= (*p) && (*p) <= 47 )
				goto st1;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto st1;
		} else
			goto st30;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto st1;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) )
				goto st1;
		} else
			goto st30;
	} else
		goto st30;
	goto tr46;
st31:
	if ( ++p == pe )
		goto _test_eof31;
case 31:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st32;
	goto tr13;
st32:
	if ( ++p == pe )
		goto _test_eof32;
case 32:
	switch( (*p) ) {
		case 9: goto tr15;
		case 10: goto tr2;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st32;
	goto tr13;
st33:
	if ( ++p == pe )
		goto _test_eof33;
case 33:
	if ( (*p) == 10 )
		goto tr2;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st34;
	goto tr8;
st34:
	if ( ++p == pe )
		goto _test_eof34;
case 34:
	switch( (*p) ) {
		case 9: goto tr10;
		case 10: goto tr2;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st1;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st1;
	} else
		goto st34;
	goto tr8;
tr52:
#line 24 "bed.ragel"
	{
        current_header_entry = bed_header_entry_new();
        ts = p;
    }
	goto st74;
st74:
	if ( ++p == pe )
		goto _test_eof74;
case 74:
#line 2529 "bed_reader.c"
	switch( (*p) ) {
		case 10: goto tr2;
		case 99: goto tr111;
		case 115: goto tr112;
	}
	if ( 0 <= (*p) )
		goto st1;
	goto tr3;
	}
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof36: cs = 36; goto _test_eof; 
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof4: cs = 4; goto _test_eof; 
	_test_eof5: cs = 5; goto _test_eof; 
	_test_eof6: cs = 6; goto _test_eof; 
	_test_eof7: cs = 7; goto _test_eof; 
	_test_eof8: cs = 8; goto _test_eof; 
	_test_eof9: cs = 9; goto _test_eof; 
	_test_eof37: cs = 37; goto _test_eof; 
	_test_eof38: cs = 38; goto _test_eof; 
	_test_eof10: cs = 10; goto _test_eof; 
	_test_eof39: cs = 39; goto _test_eof; 
	_test_eof11: cs = 11; goto _test_eof; 
	_test_eof12: cs = 12; goto _test_eof; 
	_test_eof13: cs = 13; goto _test_eof; 
	_test_eof14: cs = 14; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 
	_test_eof16: cs = 16; goto _test_eof; 
	_test_eof17: cs = 17; goto _test_eof; 
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof40: cs = 40; goto _test_eof; 
	_test_eof41: cs = 41; goto _test_eof; 
	_test_eof42: cs = 42; goto _test_eof; 
	_test_eof43: cs = 43; goto _test_eof; 
	_test_eof44: cs = 44; goto _test_eof; 
	_test_eof45: cs = 45; goto _test_eof; 
	_test_eof46: cs = 46; goto _test_eof; 
	_test_eof47: cs = 47; goto _test_eof; 
	_test_eof48: cs = 48; goto _test_eof; 
	_test_eof49: cs = 49; goto _test_eof; 
	_test_eof50: cs = 50; goto _test_eof; 
	_test_eof51: cs = 51; goto _test_eof; 
	_test_eof52: cs = 52; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 
	_test_eof53: cs = 53; goto _test_eof; 
	_test_eof54: cs = 54; goto _test_eof; 
	_test_eof21: cs = 21; goto _test_eof; 
	_test_eof22: cs = 22; goto _test_eof; 
	_test_eof55: cs = 55; goto _test_eof; 
	_test_eof56: cs = 56; goto _test_eof; 
	_test_eof23: cs = 23; goto _test_eof; 
	_test_eof57: cs = 57; goto _test_eof; 
	_test_eof58: cs = 58; goto _test_eof; 
	_test_eof24: cs = 24; goto _test_eof; 
	_test_eof59: cs = 59; goto _test_eof; 
	_test_eof25: cs = 25; goto _test_eof; 
	_test_eof60: cs = 60; goto _test_eof; 
	_test_eof26: cs = 26; goto _test_eof; 
	_test_eof61: cs = 61; goto _test_eof; 
	_test_eof62: cs = 62; goto _test_eof; 
	_test_eof27: cs = 27; goto _test_eof; 
	_test_eof63: cs = 63; goto _test_eof; 
	_test_eof64: cs = 64; goto _test_eof; 
	_test_eof28: cs = 28; goto _test_eof; 
	_test_eof65: cs = 65; goto _test_eof; 
	_test_eof66: cs = 66; goto _test_eof; 
	_test_eof67: cs = 67; goto _test_eof; 
	_test_eof68: cs = 68; goto _test_eof; 
	_test_eof69: cs = 69; goto _test_eof; 
	_test_eof70: cs = 70; goto _test_eof; 
	_test_eof71: cs = 71; goto _test_eof; 
	_test_eof29: cs = 29; goto _test_eof; 
	_test_eof72: cs = 72; goto _test_eof; 
	_test_eof73: cs = 73; goto _test_eof; 
	_test_eof30: cs = 30; goto _test_eof; 
	_test_eof31: cs = 31; goto _test_eof; 
	_test_eof32: cs = 32; goto _test_eof; 
	_test_eof33: cs = 33; goto _test_eof; 
	_test_eof34: cs = 34; goto _test_eof; 
	_test_eof74: cs = 74; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 1: 
	case 10: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
	break;
	case 37: 
	case 40: 
	case 42: 
	case 44: 
	case 46: 
	case 48: 
	case 50: 
	case 67: 
	case 69: 
#line 42 "bed.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (batch_size > 0 && current_batch->records->size == batch_size) {
            list_item_t *item = list_item_new(num_records, 1, current_batch); 
            list_insert_item(item, batches_list);
//             printf("Batch %d added - %zu records\t", batches, current_batch->records->size);
            current_batch = bed_batch_new(batch_size);
            
            if (p+1) {
                current_batch->text = p+1;
//                 printf("batch text = '%.*s'\n", 50, current_batch->text);
            }
            num_batches++;
            num_records = 0;
        }

        // If not a blank line, add current record to current batch
        add_record_to_bed_batch(current_record, current_batch);
        num_records++;
    }
	break;
	case 2: 
	case 3: 
	case 4: 
	case 5: 
	case 11: 
	case 12: 
	case 13: 
	case 14: 
	case 15: 
	case 16: 
	case 17: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 71 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'chrom' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 6: 
	case 7: 
	case 33: 
	case 34: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 85 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'start' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 8: 
	case 9: 
	case 31: 
	case 32: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 99 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'end' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 30: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 111 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'name' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 18: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 141 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'strand' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 19: 
	case 20: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 41: 
	case 43: 
	case 45: 
	case 47: 
	case 49: 
	case 51: 
	case 52: 
	case 53: 
	case 54: 
	case 55: 
	case 56: 
	case 57: 
	case 58: 
	case 59: 
	case 60: 
	case 61: 
	case 62: 
	case 63: 
	case 64: 
	case 65: 
	case 66: 
	case 68: 
	case 70: 
	case 71: 
	case 72: 
	case 73: 
#line 215 "bed.ragel"
	{
        set_bed_record_blockstarts(ts, p-ts, current_record);
    }
#line 42 "bed.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (batch_size > 0 && current_batch->records->size == batch_size) {
            list_item_t *item = list_item_new(num_records, 1, current_batch); 
            list_insert_item(item, batches_list);
//             printf("Batch %d added - %zu records\t", batches, current_batch->records->size);
            current_batch = bed_batch_new(batch_size);
            
            if (p+1) {
                current_batch->text = p+1;
//                 printf("batch text = '%.*s'\n", 50, current_batch->text);
            }
            num_batches++;
            num_records = 0;
        }

        // If not a blank line, add current record to current batch
        add_record_to_bed_batch(current_record, current_batch);
        num_records++;
    }
	break;
	case 21: 
	case 22: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 24: 
	case 25: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 181 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'itemRgb' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 23: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 195 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockCount' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 26: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 181 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'itemRgb' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 195 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockCount' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 27: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 169 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickEnd' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 181 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'itemRgb' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 195 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockCount' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 28: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 155 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickStart' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 169 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickEnd' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 181 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'itemRgb' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 195 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockCount' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 29: 
#line 34 "bed.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 129 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'score' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 155 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickStart' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 169 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'thickEnd' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 181 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'itemRgb' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 195 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockCount' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 207 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockSizes' field\n", num_batches * batch_size + num_records, file->filename);
    }
#line 219 "bed.ragel"
	{
        printf("Line %zu (%s): Error in 'blockStarts' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
#line 2948 "bed_reader.c"
	}
	}

	_out: {}
	}

#line 291 "bed.ragel"
 

    // Insert the last batch
    if (!bed_batch_is_empty(current_batch)) {
        list_item_t *item = list_item_new(num_records, 1, current_batch); 
        list_insert_item(item, batches_list);
        LOG_DEBUG_F("Batch added - %zu records (last)\n", current_batch->records->size);
    }

    if ( cs < 
#line 2966 "bed_reader.c"
35
#line 300 "bed.ragel"
 ) {
        LOG_INFO_F("Last state is %d, but %d was expected\n", 
                cs, 
#line 2972 "bed_reader.c"
35
#line 302 "bed.ragel"
);
    } 

    LOG_INFO_F("BED records read = %zu\n", num_batches * batch_size + num_records);

    // Free current_xxx pointers if not needed in another module
    //bed_header_entry_free(current_header_entry);

    return cs < 
#line 2984 "bed_reader.c"
35
#line 310 "bed.ragel"
;
}
