
#line 1 "gff.ragel"
#include "gff_reader.h"

static size_t lines = 1;
static size_t num_records = 0;
static size_t num_batches = 0;

static gff_record_t *current_record;
static gff_header_entry_t *current_header_entry;
static gff_batch_t *current_batch;


#line 15 "gff_reader.c"
static const int gff_start = 62;
static const int gff_first_final = 62;
static const int gff_error = 0;

static const int gff_en_main = 62;


#line 225 "gff.ragel"



int gff_ragel_read(list_t *batches_list, size_t batch_size, gff_file_t *file) {
    int cs;
    char *p = file->data;
    char *pe = p + file->data_len;
    char *eof = pe;
    char *ts, *te;
    int stack[4];
    int top, act;

    current_header_entry = gff_header_entry_new();
    current_batch = gff_batch_new(batch_size);

    
#line 40 "gff_reader.c"
	{
	cs = gff_start;
	}

#line 45 "gff_reader.c"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
case 62:
	switch( (*p) ) {
		case 10: goto st63;
		case 35: goto st29;
		case 95: goto tr97;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr97;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr97;
	} else
		goto tr97;
	goto tr0;
tr0:
#line 71 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'seqname' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr3:
#line 83 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'source' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr7:
#line 95 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'feature' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr11:
#line 109 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'start' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr16:
#line 123 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'end' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr21:
#line 141 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'score' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr25:
#line 153 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'strand' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr28:
#line 165 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'frame' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr45:
#line 34 "gff.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
	goto st0;
tr51:
#line 34 "gff.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 71 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'seqname' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr54:
#line 34 "gff.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 83 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'source' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr58:
#line 34 "gff.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 95 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'feature' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr62:
#line 34 "gff.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 109 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'start' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr67:
#line 34 "gff.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 123 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'end' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr72:
#line 34 "gff.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 141 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'score' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr76:
#line 34 "gff.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 153 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'strand' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr79:
#line 34 "gff.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 165 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'frame' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr98:
#line 177 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'attribute' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
tr108:
#line 34 "gff.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 177 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'attribute' field\n", num_batches * batch_size + num_records, file->filename);
    }
	goto st0;
#line 216 "gff_reader.c"
st0:
cs = 0;
	goto _out;
st63:
	if ( ++p == pe )
		goto _test_eof63;
case 63:
	switch( (*p) ) {
		case 10: goto st63;
		case 95: goto tr97;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr97;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr97;
	} else
		goto tr97;
	goto tr0;
tr97:
#line 38 "gff.ragel"
	{
        current_record = gff_record_new();
    }
#line 63 "gff.ragel"
	{
        ts = p;
    }
	goto st1;
st1:
	if ( ++p == pe )
		goto _test_eof1;
case 1:
#line 251 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto tr1;
		case 95: goto st1;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st1;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st1;
	} else
		goto st1;
	goto tr0;
tr1:
#line 67 "gff.ragel"
	{
        set_gff_record_sequence(ts, p-ts, current_record);
    }
	goto st2;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
#line 275 "gff_reader.c"
	switch( (*p) ) {
		case 46: goto tr4;
		case 95: goto tr5;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr5;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr5;
	} else
		goto tr5;
	goto tr3;
tr4:
#line 75 "gff.ragel"
	{
        ts = p;
    }
	goto st3;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
#line 299 "gff_reader.c"
	if ( (*p) == 9 )
		goto tr6;
	goto tr3;
tr6:
#line 79 "gff.ragel"
	{
        set_gff_record_source(ts, p-ts, current_record);
    }
	goto st4;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
#line 313 "gff_reader.c"
	switch( (*p) ) {
		case 46: goto tr8;
		case 95: goto tr9;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr9;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr9;
	} else
		goto tr9;
	goto tr7;
tr8:
#line 87 "gff.ragel"
	{
        ts = p;
    }
	goto st5;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
#line 337 "gff_reader.c"
	if ( (*p) == 9 )
		goto tr10;
	goto tr7;
tr10:
#line 91 "gff.ragel"
	{
        set_gff_record_feature(ts, p-ts, current_record);
    }
	goto st6;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
#line 351 "gff_reader.c"
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr12;
	goto tr11;
tr12:
#line 99 "gff.ragel"
	{
        ts = p;
    }
	goto st7;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
#line 365 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto tr13;
		case 46: goto st25;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st7;
	goto tr11;
tr13:
#line 103 "gff.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_gff_record_start(atol(field), current_record);
        free(field);
    }
	goto st8;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
#line 385 "gff_reader.c"
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr17;
	goto tr16;
tr17:
#line 113 "gff.ragel"
	{
        ts = p;
    }
	goto st9;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
#line 399 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto tr18;
		case 46: goto st23;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st9;
	goto tr16;
tr18:
#line 117 "gff.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_gff_record_end(atol(field), current_record);
        free(field);
    }
	goto st10;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
#line 419 "gff_reader.c"
	if ( (*p) == 46 )
		goto tr22;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr23;
	goto tr21;
tr22:
#line 127 "gff.ragel"
	{
        ts = p;
    }
	goto st11;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
#line 435 "gff_reader.c"
	if ( (*p) == 9 )
		goto tr24;
	goto tr21;
tr24:
#line 131 "gff.ragel"
	{
        float score = -1.0f;
        if (strncmp(".", ts, 1) != 0) {
            char *field = strndup(ts, p-ts);
            score = atof(field);
            free(field);
        }
        set_gff_record_score(score, current_record);
    }
	goto st12;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
#line 455 "gff_reader.c"
	if ( (*p) == 43 )
		goto tr26;
	if ( 45 <= (*p) && (*p) <= 46 )
		goto tr26;
	goto tr25;
tr26:
#line 145 "gff.ragel"
	{
        ts = p;
    }
	goto st13;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
#line 471 "gff_reader.c"
	if ( (*p) == 9 )
		goto tr27;
	goto tr25;
tr27:
#line 149 "gff.ragel"
	{
        set_gff_record_strand(*ts, current_record);
    }
	goto st14;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
#line 485 "gff_reader.c"
	if ( (*p) == 46 )
		goto tr29;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr30;
	goto tr28;
tr29:
#line 157 "gff.ragel"
	{
        ts = p;
    }
	goto st15;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
#line 501 "gff_reader.c"
	if ( (*p) == 9 )
		goto tr31;
	goto tr28;
tr31:
#line 161 "gff.ragel"
	{
        set_gff_record_frame(*ts, current_record);
    }
	goto st64;
st64:
	if ( ++p == pe )
		goto _test_eof64;
case 64:
#line 515 "gff_reader.c"
	if ( (*p) == 10 )
		goto tr99;
	if ( 32 <= (*p) && (*p) <= 126 )
		goto tr100;
	goto tr98;
tr99:
#line 169 "gff.ragel"
	{
        ts = p;
    }
#line 173 "gff.ragel"
	{
        set_gff_record_attribute(ts, p-ts, current_record);
    }
#line 42 "gff.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (batch_size > 0 && current_batch->records->size == batch_size) {
            list_item_t *item = list_item_new(num_records, 1, current_batch); 
            list_insert_item(item, batches_list);
//             printf("Batch %d added - %zu records\t", batches, current_batch->records->size);
            current_batch = gff_batch_new(batch_size);
            
            if (p+1) {
                current_batch->text = p+1;
//                 printf("batch text = '%.*s'\n", 50, current_batch->text);
            }
            num_batches++;
            num_records = 0;
        }

        // If not a blank line, add current record to current batch
        add_record_to_gff_batch(current_record, current_batch);
        num_records++;
    }
#line 19 "gff.ragel"
	{
        lines++;
//        printf("lines read = %d\n", lines);
    }
	goto st65;
tr103:
#line 19 "gff.ragel"
	{
        lines++;
//        printf("lines read = %d\n", lines);
    }
	goto st65;
tr104:
#line 173 "gff.ragel"
	{
        set_gff_record_attribute(ts, p-ts, current_record);
    }
#line 42 "gff.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (batch_size > 0 && current_batch->records->size == batch_size) {
            list_item_t *item = list_item_new(num_records, 1, current_batch); 
            list_insert_item(item, batches_list);
//             printf("Batch %d added - %zu records\t", batches, current_batch->records->size);
            current_batch = gff_batch_new(batch_size);
            
            if (p+1) {
                current_batch->text = p+1;
//                 printf("batch text = '%.*s'\n", 50, current_batch->text);
            }
            num_batches++;
            num_records = 0;
        }

        // If not a blank line, add current record to current batch
        add_record_to_gff_batch(current_record, current_batch);
        num_records++;
    }
#line 19 "gff.ragel"
	{
        lines++;
//        printf("lines read = %d\n", lines);
    }
	goto st65;
st65:
	if ( ++p == pe )
		goto _test_eof65;
case 65:
#line 600 "gff_reader.c"
	switch( (*p) ) {
		case 10: goto st66;
		case 35: goto st16;
		case 95: goto tr97;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr97;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr97;
	} else
		goto tr97;
	goto tr0;
st66:
	if ( ++p == pe )
		goto _test_eof66;
case 66:
	if ( (*p) == 10 )
		goto st66;
	goto st0;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st67;
	goto st0;
st67:
	if ( ++p == pe )
		goto _test_eof67;
case 67:
	if ( (*p) == 10 )
		goto tr103;
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st67;
	goto st0;
tr100:
#line 169 "gff.ragel"
	{
        ts = p;
    }
	goto st68;
st68:
	if ( ++p == pe )
		goto _test_eof68;
case 68:
#line 648 "gff_reader.c"
	if ( (*p) == 10 )
		goto tr104;
	if ( 32 <= (*p) && (*p) <= 126 )
		goto st68;
	goto tr98;
tr30:
#line 157 "gff.ragel"
	{
        ts = p;
    }
	goto st17;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
#line 664 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto tr31;
		case 46: goto st18;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st17;
	goto tr28;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st19;
	goto tr28;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
	if ( (*p) == 9 )
		goto tr31;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st19;
	goto tr28;
tr23:
#line 127 "gff.ragel"
	{
        ts = p;
    }
	goto st20;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
#line 698 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto tr24;
		case 46: goto st21;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st20;
	goto tr21;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st22;
	goto tr21;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
	if ( (*p) == 9 )
		goto tr24;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st22;
	goto tr21;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st24;
	goto tr16;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
	if ( (*p) == 9 )
		goto tr18;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st24;
	goto tr16;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st26;
	goto tr11;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
	if ( (*p) == 9 )
		goto tr13;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st26;
	goto tr11;
tr9:
#line 87 "gff.ragel"
	{
        ts = p;
    }
	goto st27;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
#line 764 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto tr10;
		case 95: goto st27;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st27;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st27;
	} else
		goto st27;
	goto tr7;
tr5:
#line 75 "gff.ragel"
	{
        ts = p;
    }
	goto st28;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
#line 788 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto tr6;
		case 95: goto st28;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st28;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st28;
	} else
		goto st28;
	goto tr3;
st29:
	if ( ++p == pe )
		goto _test_eof29;
case 29:
	if ( (*p) == 35 )
		goto st30;
	goto st0;
st30:
	if ( ++p == pe )
		goto _test_eof30;
case 30:
	if ( 0 <= (*p) )
		goto tr46;
	goto tr45;
tr46:
#line 24 "gff.ragel"
	{
        current_header_entry = gff_header_entry_new();
        ts = p;
    }
	goto st31;
st31:
	if ( ++p == pe )
		goto _test_eof31;
case 31:
#line 827 "gff_reader.c"
	if ( (*p) == 10 )
		goto tr48;
	if ( 0 <= (*p) )
		goto st31;
	goto tr45;
tr48:
#line 29 "gff.ragel"
	{
        set_gff_header_entry_text(ts, p-ts, current_header_entry);
        add_gff_header_entry(current_header_entry, file);
    }
#line 19 "gff.ragel"
	{
        lines++;
//        printf("lines read = %d\n", lines);
    }
	goto st69;
tr50:
#line 24 "gff.ragel"
	{
        current_header_entry = gff_header_entry_new();
        ts = p;
    }
#line 29 "gff.ragel"
	{
        set_gff_header_entry_text(ts, p-ts, current_header_entry);
        add_gff_header_entry(current_header_entry, file);
    }
#line 19 "gff.ragel"
	{
        lines++;
//        printf("lines read = %d\n", lines);
    }
	goto st69;
st69:
	if ( ++p == pe )
		goto _test_eof69;
case 69:
#line 866 "gff_reader.c"
	switch( (*p) ) {
		case 10: goto tr48;
		case 35: goto st32;
		case 95: goto tr107;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 0 <= (*p) && (*p) <= 47 )
				goto st31;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto st31;
		} else
			goto tr107;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto st31;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) )
				goto st31;
		} else
			goto tr107;
	} else
		goto tr107;
	goto tr51;
st32:
	if ( ++p == pe )
		goto _test_eof32;
case 32:
	switch( (*p) ) {
		case 10: goto tr48;
		case 35: goto st33;
	}
	if ( 0 <= (*p) )
		goto st31;
	goto tr45;
st33:
	if ( ++p == pe )
		goto _test_eof33;
case 33:
	if ( (*p) == 10 )
		goto tr50;
	if ( 0 <= (*p) )
		goto tr46;
	goto tr45;
tr107:
#line 38 "gff.ragel"
	{
        current_record = gff_record_new();
    }
#line 63 "gff.ragel"
	{
        ts = p;
    }
	goto st34;
st34:
	if ( ++p == pe )
		goto _test_eof34;
case 34:
#line 927 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto tr52;
		case 10: goto tr48;
		case 95: goto st34;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 0 <= (*p) && (*p) <= 47 )
				goto st31;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto st31;
		} else
			goto st34;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto st31;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) )
				goto st31;
		} else
			goto st34;
	} else
		goto st34;
	goto tr51;
tr52:
#line 67 "gff.ragel"
	{
        set_gff_record_sequence(ts, p-ts, current_record);
    }
	goto st35;
st35:
	if ( ++p == pe )
		goto _test_eof35;
case 35:
#line 964 "gff_reader.c"
	switch( (*p) ) {
		case 10: goto tr48;
		case 46: goto tr55;
		case 95: goto tr56;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 0 <= (*p) && (*p) <= 47 )
				goto st31;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto st31;
		} else
			goto tr56;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto st31;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) )
				goto st31;
		} else
			goto tr56;
	} else
		goto tr56;
	goto tr54;
tr55:
#line 75 "gff.ragel"
	{
        ts = p;
    }
	goto st36;
st36:
	if ( ++p == pe )
		goto _test_eof36;
case 36:
#line 1001 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto tr57;
		case 10: goto tr48;
	}
	if ( 0 <= (*p) )
		goto st31;
	goto tr54;
tr57:
#line 79 "gff.ragel"
	{
        set_gff_record_source(ts, p-ts, current_record);
    }
	goto st37;
st37:
	if ( ++p == pe )
		goto _test_eof37;
case 37:
#line 1019 "gff_reader.c"
	switch( (*p) ) {
		case 10: goto tr48;
		case 46: goto tr59;
		case 95: goto tr60;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 0 <= (*p) && (*p) <= 47 )
				goto st31;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto st31;
		} else
			goto tr60;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto st31;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) )
				goto st31;
		} else
			goto tr60;
	} else
		goto tr60;
	goto tr58;
tr59:
#line 87 "gff.ragel"
	{
        ts = p;
    }
	goto st38;
st38:
	if ( ++p == pe )
		goto _test_eof38;
case 38:
#line 1056 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto tr61;
		case 10: goto tr48;
	}
	if ( 0 <= (*p) )
		goto st31;
	goto tr58;
tr61:
#line 91 "gff.ragel"
	{
        set_gff_record_feature(ts, p-ts, current_record);
    }
	goto st39;
st39:
	if ( ++p == pe )
		goto _test_eof39;
case 39:
#line 1074 "gff_reader.c"
	if ( (*p) == 10 )
		goto tr48;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st31;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st31;
	} else
		goto tr63;
	goto tr62;
tr63:
#line 99 "gff.ragel"
	{
        ts = p;
    }
	goto st40;
st40:
	if ( ++p == pe )
		goto _test_eof40;
case 40:
#line 1096 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto tr64;
		case 10: goto tr48;
		case 46: goto st58;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st31;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st31;
	} else
		goto st40;
	goto tr62;
tr64:
#line 103 "gff.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_gff_record_start(atol(field), current_record);
        free(field);
    }
	goto st41;
st41:
	if ( ++p == pe )
		goto _test_eof41;
case 41:
#line 1123 "gff_reader.c"
	if ( (*p) == 10 )
		goto tr48;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st31;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st31;
	} else
		goto tr68;
	goto tr67;
tr68:
#line 113 "gff.ragel"
	{
        ts = p;
    }
	goto st42;
st42:
	if ( ++p == pe )
		goto _test_eof42;
case 42:
#line 1145 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto tr69;
		case 10: goto tr48;
		case 46: goto st56;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st31;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st31;
	} else
		goto st42;
	goto tr67;
tr69:
#line 117 "gff.ragel"
	{
        char *field = strndup(ts, p-ts);
        set_gff_record_end(atol(field), current_record);
        free(field);
    }
	goto st43;
st43:
	if ( ++p == pe )
		goto _test_eof43;
case 43:
#line 1172 "gff_reader.c"
	switch( (*p) ) {
		case 10: goto tr48;
		case 46: goto tr73;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st31;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st31;
	} else
		goto tr74;
	goto tr72;
tr73:
#line 127 "gff.ragel"
	{
        ts = p;
    }
	goto st44;
st44:
	if ( ++p == pe )
		goto _test_eof44;
case 44:
#line 1196 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto tr75;
		case 10: goto tr48;
	}
	if ( 0 <= (*p) )
		goto st31;
	goto tr72;
tr75:
#line 131 "gff.ragel"
	{
        float score = -1.0f;
        if (strncmp(".", ts, 1) != 0) {
            char *field = strndup(ts, p-ts);
            score = atof(field);
            free(field);
        }
        set_gff_record_score(score, current_record);
    }
	goto st45;
st45:
	if ( ++p == pe )
		goto _test_eof45;
case 45:
#line 1220 "gff_reader.c"
	switch( (*p) ) {
		case 10: goto tr48;
		case 43: goto tr77;
	}
	if ( (*p) < 45 ) {
		if ( 0 <= (*p) && (*p) <= 44 )
			goto st31;
	} else if ( (*p) > 46 ) {
		if ( 47 <= (*p) )
			goto st31;
	} else
		goto tr77;
	goto tr76;
tr77:
#line 145 "gff.ragel"
	{
        ts = p;
    }
	goto st46;
st46:
	if ( ++p == pe )
		goto _test_eof46;
case 46:
#line 1244 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto tr78;
		case 10: goto tr48;
	}
	if ( 0 <= (*p) )
		goto st31;
	goto tr76;
tr78:
#line 149 "gff.ragel"
	{
        set_gff_record_strand(*ts, current_record);
    }
	goto st47;
st47:
	if ( ++p == pe )
		goto _test_eof47;
case 47:
#line 1262 "gff_reader.c"
	switch( (*p) ) {
		case 10: goto tr48;
		case 46: goto tr80;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st31;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st31;
	} else
		goto tr81;
	goto tr79;
tr80:
#line 157 "gff.ragel"
	{
        ts = p;
    }
	goto st48;
st48:
	if ( ++p == pe )
		goto _test_eof48;
case 48:
#line 1286 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto tr82;
		case 10: goto tr48;
	}
	if ( 0 <= (*p) )
		goto st31;
	goto tr79;
tr82:
#line 161 "gff.ragel"
	{
        set_gff_record_frame(*ts, current_record);
    }
	goto st70;
st70:
	if ( ++p == pe )
		goto _test_eof70;
case 70:
#line 1304 "gff_reader.c"
	switch( (*p) ) {
		case 10: goto tr109;
		case 127: goto st31;
	}
	if ( (*p) > 31 ) {
		if ( 32 <= (*p) && (*p) <= 126 )
			goto tr110;
	} else if ( (*p) >= 0 )
		goto st31;
	goto tr108;
tr112:
#line 29 "gff.ragel"
	{
        set_gff_header_entry_text(ts, p-ts, current_header_entry);
        add_gff_header_entry(current_header_entry, file);
    }
#line 19 "gff.ragel"
	{
        lines++;
//        printf("lines read = %d\n", lines);
    }
	goto st71;
tr113:
#line 24 "gff.ragel"
	{
        current_header_entry = gff_header_entry_new();
        ts = p;
    }
#line 29 "gff.ragel"
	{
        set_gff_header_entry_text(ts, p-ts, current_header_entry);
        add_gff_header_entry(current_header_entry, file);
    }
#line 19 "gff.ragel"
	{
        lines++;
//        printf("lines read = %d\n", lines);
    }
	goto st71;
tr109:
#line 29 "gff.ragel"
	{
        set_gff_header_entry_text(ts, p-ts, current_header_entry);
        add_gff_header_entry(current_header_entry, file);
    }
#line 19 "gff.ragel"
	{
        lines++;
//        printf("lines read = %d\n", lines);
    }
#line 169 "gff.ragel"
	{
        ts = p;
    }
#line 173 "gff.ragel"
	{
        set_gff_record_attribute(ts, p-ts, current_record);
    }
#line 42 "gff.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (batch_size > 0 && current_batch->records->size == batch_size) {
            list_item_t *item = list_item_new(num_records, 1, current_batch); 
            list_insert_item(item, batches_list);
//             printf("Batch %d added - %zu records\t", batches, current_batch->records->size);
            current_batch = gff_batch_new(batch_size);
            
            if (p+1) {
                current_batch->text = p+1;
//                 printf("batch text = '%.*s'\n", 50, current_batch->text);
            }
            num_batches++;
            num_records = 0;
        }

        // If not a blank line, add current record to current batch
        add_record_to_gff_batch(current_record, current_batch);
        num_records++;
    }
	goto st71;
tr115:
#line 29 "gff.ragel"
	{
        set_gff_header_entry_text(ts, p-ts, current_header_entry);
        add_gff_header_entry(current_header_entry, file);
    }
#line 19 "gff.ragel"
	{
        lines++;
//        printf("lines read = %d\n", lines);
    }
#line 173 "gff.ragel"
	{
        set_gff_record_attribute(ts, p-ts, current_record);
    }
#line 42 "gff.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (batch_size > 0 && current_batch->records->size == batch_size) {
            list_item_t *item = list_item_new(num_records, 1, current_batch); 
            list_insert_item(item, batches_list);
//             printf("Batch %d added - %zu records\t", batches, current_batch->records->size);
            current_batch = gff_batch_new(batch_size);
            
            if (p+1) {
                current_batch->text = p+1;
//                 printf("batch text = '%.*s'\n", 50, current_batch->text);
            }
            num_batches++;
            num_records = 0;
        }

        // If not a blank line, add current record to current batch
        add_record_to_gff_batch(current_record, current_batch);
        num_records++;
    }
	goto st71;
st71:
	if ( ++p == pe )
		goto _test_eof71;
case 71:
#line 1426 "gff_reader.c"
	switch( (*p) ) {
		case 10: goto tr48;
		case 35: goto st49;
		case 95: goto tr107;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 0 <= (*p) && (*p) <= 47 )
				goto st31;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto st31;
		} else
			goto tr107;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto st31;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) )
				goto st31;
		} else
			goto tr107;
	} else
		goto tr107;
	goto tr51;
st49:
	if ( ++p == pe )
		goto _test_eof49;
case 49:
	switch( (*p) ) {
		case 10: goto tr48;
		case 35: goto st73;
		case 127: goto st31;
	}
	if ( (*p) > 31 ) {
		if ( 32 <= (*p) && (*p) <= 126 )
			goto st72;
	} else if ( (*p) >= 0 )
		goto st31;
	goto tr45;
tr114:
#line 24 "gff.ragel"
	{
        current_header_entry = gff_header_entry_new();
        ts = p;
    }
	goto st72;
st72:
	if ( ++p == pe )
		goto _test_eof72;
case 72:
#line 1479 "gff_reader.c"
	switch( (*p) ) {
		case 10: goto tr112;
		case 127: goto st31;
	}
	if ( (*p) > 31 ) {
		if ( 32 <= (*p) && (*p) <= 126 )
			goto st72;
	} else if ( (*p) >= 0 )
		goto st31;
	goto tr45;
st73:
	if ( ++p == pe )
		goto _test_eof73;
case 73:
	switch( (*p) ) {
		case 10: goto tr113;
		case 127: goto tr46;
	}
	if ( (*p) > 31 ) {
		if ( 32 <= (*p) && (*p) <= 126 )
			goto tr114;
	} else if ( (*p) >= 0 )
		goto tr46;
	goto tr45;
tr110:
#line 169 "gff.ragel"
	{
        ts = p;
    }
	goto st74;
st74:
	if ( ++p == pe )
		goto _test_eof74;
case 74:
#line 1514 "gff_reader.c"
	switch( (*p) ) {
		case 10: goto tr115;
		case 127: goto st31;
	}
	if ( (*p) > 31 ) {
		if ( 32 <= (*p) && (*p) <= 126 )
			goto st74;
	} else if ( (*p) >= 0 )
		goto st31;
	goto tr108;
tr81:
#line 157 "gff.ragel"
	{
        ts = p;
    }
	goto st50;
st50:
	if ( ++p == pe )
		goto _test_eof50;
case 50:
#line 1535 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto tr82;
		case 10: goto tr48;
		case 46: goto st51;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st31;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st31;
	} else
		goto st50;
	goto tr79;
st51:
	if ( ++p == pe )
		goto _test_eof51;
case 51:
	if ( (*p) == 10 )
		goto tr48;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st31;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st31;
	} else
		goto st52;
	goto tr79;
st52:
	if ( ++p == pe )
		goto _test_eof52;
case 52:
	switch( (*p) ) {
		case 9: goto tr82;
		case 10: goto tr48;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st31;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st31;
	} else
		goto st52;
	goto tr79;
tr74:
#line 127 "gff.ragel"
	{
        ts = p;
    }
	goto st53;
st53:
	if ( ++p == pe )
		goto _test_eof53;
case 53:
#line 1592 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto tr75;
		case 10: goto tr48;
		case 46: goto st54;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st31;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st31;
	} else
		goto st53;
	goto tr72;
st54:
	if ( ++p == pe )
		goto _test_eof54;
case 54:
	if ( (*p) == 10 )
		goto tr48;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st31;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st31;
	} else
		goto st55;
	goto tr72;
st55:
	if ( ++p == pe )
		goto _test_eof55;
case 55:
	switch( (*p) ) {
		case 9: goto tr75;
		case 10: goto tr48;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st31;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st31;
	} else
		goto st55;
	goto tr72;
st56:
	if ( ++p == pe )
		goto _test_eof56;
case 56:
	if ( (*p) == 10 )
		goto tr48;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st31;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st31;
	} else
		goto st57;
	goto tr67;
st57:
	if ( ++p == pe )
		goto _test_eof57;
case 57:
	switch( (*p) ) {
		case 9: goto tr69;
		case 10: goto tr48;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st31;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st31;
	} else
		goto st57;
	goto tr67;
st58:
	if ( ++p == pe )
		goto _test_eof58;
case 58:
	if ( (*p) == 10 )
		goto tr48;
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st31;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st31;
	} else
		goto st59;
	goto tr62;
st59:
	if ( ++p == pe )
		goto _test_eof59;
case 59:
	switch( (*p) ) {
		case 9: goto tr64;
		case 10: goto tr48;
	}
	if ( (*p) < 48 ) {
		if ( 0 <= (*p) && (*p) <= 47 )
			goto st31;
	} else if ( (*p) > 57 ) {
		if ( 58 <= (*p) )
			goto st31;
	} else
		goto st59;
	goto tr62;
tr60:
#line 87 "gff.ragel"
	{
        ts = p;
    }
	goto st60;
st60:
	if ( ++p == pe )
		goto _test_eof60;
case 60:
#line 1713 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto tr61;
		case 10: goto tr48;
		case 95: goto st60;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 0 <= (*p) && (*p) <= 47 )
				goto st31;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto st31;
		} else
			goto st60;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto st31;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) )
				goto st31;
		} else
			goto st60;
	} else
		goto st60;
	goto tr58;
tr56:
#line 75 "gff.ragel"
	{
        ts = p;
    }
	goto st61;
st61:
	if ( ++p == pe )
		goto _test_eof61;
case 61:
#line 1750 "gff_reader.c"
	switch( (*p) ) {
		case 9: goto tr57;
		case 10: goto tr48;
		case 95: goto st61;
	}
	if ( (*p) < 65 ) {
		if ( (*p) < 48 ) {
			if ( 0 <= (*p) && (*p) <= 47 )
				goto st31;
		} else if ( (*p) > 57 ) {
			if ( 58 <= (*p) && (*p) <= 64 )
				goto st31;
		} else
			goto st61;
	} else if ( (*p) > 90 ) {
		if ( (*p) < 97 ) {
			if ( 91 <= (*p) && (*p) <= 96 )
				goto st31;
		} else if ( (*p) > 122 ) {
			if ( 123 <= (*p) )
				goto st31;
		} else
			goto st61;
	} else
		goto st61;
	goto tr54;
	}
	_test_eof63: cs = 63; goto _test_eof; 
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof4: cs = 4; goto _test_eof; 
	_test_eof5: cs = 5; goto _test_eof; 
	_test_eof6: cs = 6; goto _test_eof; 
	_test_eof7: cs = 7; goto _test_eof; 
	_test_eof8: cs = 8; goto _test_eof; 
	_test_eof9: cs = 9; goto _test_eof; 
	_test_eof10: cs = 10; goto _test_eof; 
	_test_eof11: cs = 11; goto _test_eof; 
	_test_eof12: cs = 12; goto _test_eof; 
	_test_eof13: cs = 13; goto _test_eof; 
	_test_eof14: cs = 14; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 
	_test_eof64: cs = 64; goto _test_eof; 
	_test_eof65: cs = 65; goto _test_eof; 
	_test_eof66: cs = 66; goto _test_eof; 
	_test_eof16: cs = 16; goto _test_eof; 
	_test_eof67: cs = 67; goto _test_eof; 
	_test_eof68: cs = 68; goto _test_eof; 
	_test_eof17: cs = 17; goto _test_eof; 
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 
	_test_eof21: cs = 21; goto _test_eof; 
	_test_eof22: cs = 22; goto _test_eof; 
	_test_eof23: cs = 23; goto _test_eof; 
	_test_eof24: cs = 24; goto _test_eof; 
	_test_eof25: cs = 25; goto _test_eof; 
	_test_eof26: cs = 26; goto _test_eof; 
	_test_eof27: cs = 27; goto _test_eof; 
	_test_eof28: cs = 28; goto _test_eof; 
	_test_eof29: cs = 29; goto _test_eof; 
	_test_eof30: cs = 30; goto _test_eof; 
	_test_eof31: cs = 31; goto _test_eof; 
	_test_eof69: cs = 69; goto _test_eof; 
	_test_eof32: cs = 32; goto _test_eof; 
	_test_eof33: cs = 33; goto _test_eof; 
	_test_eof34: cs = 34; goto _test_eof; 
	_test_eof35: cs = 35; goto _test_eof; 
	_test_eof36: cs = 36; goto _test_eof; 
	_test_eof37: cs = 37; goto _test_eof; 
	_test_eof38: cs = 38; goto _test_eof; 
	_test_eof39: cs = 39; goto _test_eof; 
	_test_eof40: cs = 40; goto _test_eof; 
	_test_eof41: cs = 41; goto _test_eof; 
	_test_eof42: cs = 42; goto _test_eof; 
	_test_eof43: cs = 43; goto _test_eof; 
	_test_eof44: cs = 44; goto _test_eof; 
	_test_eof45: cs = 45; goto _test_eof; 
	_test_eof46: cs = 46; goto _test_eof; 
	_test_eof47: cs = 47; goto _test_eof; 
	_test_eof48: cs = 48; goto _test_eof; 
	_test_eof70: cs = 70; goto _test_eof; 
	_test_eof71: cs = 71; goto _test_eof; 
	_test_eof49: cs = 49; goto _test_eof; 
	_test_eof72: cs = 72; goto _test_eof; 
	_test_eof73: cs = 73; goto _test_eof; 
	_test_eof74: cs = 74; goto _test_eof; 
	_test_eof50: cs = 50; goto _test_eof; 
	_test_eof51: cs = 51; goto _test_eof; 
	_test_eof52: cs = 52; goto _test_eof; 
	_test_eof53: cs = 53; goto _test_eof; 
	_test_eof54: cs = 54; goto _test_eof; 
	_test_eof55: cs = 55; goto _test_eof; 
	_test_eof56: cs = 56; goto _test_eof; 
	_test_eof57: cs = 57; goto _test_eof; 
	_test_eof58: cs = 58; goto _test_eof; 
	_test_eof59: cs = 59; goto _test_eof; 
	_test_eof60: cs = 60; goto _test_eof; 
	_test_eof61: cs = 61; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 30: 
	case 31: 
	case 32: 
	case 33: 
	case 49: 
#line 34 "gff.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
	break;
	case 1: 
#line 71 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'seqname' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 2: 
	case 3: 
	case 28: 
#line 83 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'source' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 4: 
	case 5: 
	case 27: 
#line 95 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'feature' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 6: 
	case 7: 
	case 25: 
	case 26: 
#line 109 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'start' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 8: 
	case 9: 
	case 23: 
	case 24: 
#line 123 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'end' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 10: 
	case 11: 
	case 20: 
	case 21: 
	case 22: 
#line 141 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'score' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 12: 
	case 13: 
#line 153 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'strand' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 14: 
	case 15: 
	case 17: 
	case 18: 
	case 19: 
#line 165 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'frame' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 34: 
#line 34 "gff.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 71 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'seqname' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 35: 
	case 36: 
	case 61: 
#line 34 "gff.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 83 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'source' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 37: 
	case 38: 
	case 60: 
#line 34 "gff.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 95 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'feature' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 39: 
	case 40: 
	case 58: 
	case 59: 
#line 34 "gff.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 109 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'start' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 41: 
	case 42: 
	case 56: 
	case 57: 
#line 34 "gff.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 123 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'end' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 43: 
	case 44: 
	case 53: 
	case 54: 
	case 55: 
#line 34 "gff.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 141 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'score' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 45: 
	case 46: 
#line 34 "gff.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 153 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'strand' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 47: 
	case 48: 
	case 50: 
	case 51: 
	case 52: 
#line 34 "gff.ragel"
	{
        printf("Line %zu (%s): Error in header\n", lines, file->filename);
    }
#line 165 "gff.ragel"
	{
        printf("Line %zu (%s): Error in 'frame' field\n", num_batches * batch_size + num_records, file->filename);
    }
	break;
	case 68: 
	case 74: 
#line 173 "gff.ragel"
	{
        set_gff_record_attribute(ts, p-ts, current_record);
    }
#line 42 "gff.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (batch_size > 0 && current_batch->records->size == batch_size) {
            list_item_t *item = list_item_new(num_records, 1, current_batch); 
            list_insert_item(item, batches_list);
//             printf("Batch %d added - %zu records\t", batches, current_batch->records->size);
            current_batch = gff_batch_new(batch_size);
            
            if (p+1) {
                current_batch->text = p+1;
//                 printf("batch text = '%.*s'\n", 50, current_batch->text);
            }
            num_batches++;
            num_records = 0;
        }

        // If not a blank line, add current record to current batch
        add_record_to_gff_batch(current_record, current_batch);
        num_records++;
    }
	break;
	case 64: 
	case 70: 
#line 169 "gff.ragel"
	{
        ts = p;
    }
#line 173 "gff.ragel"
	{
        set_gff_record_attribute(ts, p-ts, current_record);
    }
#line 42 "gff.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (batch_size > 0 && current_batch->records->size == batch_size) {
            list_item_t *item = list_item_new(num_records, 1, current_batch); 
            list_insert_item(item, batches_list);
//             printf("Batch %d added - %zu records\t", batches, current_batch->records->size);
            current_batch = gff_batch_new(batch_size);
            
            if (p+1) {
                current_batch->text = p+1;
//                 printf("batch text = '%.*s'\n", 50, current_batch->text);
            }
            num_batches++;
            num_records = 0;
        }

        // If not a blank line, add current record to current batch
        add_record_to_gff_batch(current_record, current_batch);
        num_records++;
    }
	break;
#line 2092 "gff_reader.c"
	}
	}

	_out: {}
	}

#line 243 "gff.ragel"
 

    // Insert the last batch
    if (!gff_batch_is_empty(current_batch)) {
        list_item_t *item = list_item_new(num_records, 1, current_batch); 
        list_insert_item(item, batches_list);
        LOG_DEBUG_F("Batch added - %zu records (last)\n", current_batch->records->size);
    }

    if ( cs < 
#line 2110 "gff_reader.c"
62
#line 252 "gff.ragel"
 ) {
        LOG_INFO_F("Last state is %d, but %d was expected\n", 
                cs, 
#line 2116 "gff_reader.c"
62
#line 254 "gff.ragel"
);
    } 

    LOG_INFO_F("GFF records read = %zu\n", num_batches * batch_size + num_records);

    // Free current_xxx pointers if not needed in another module
    //gff_header_entry_free(current_header_entry);

    return cs < 
#line 2128 "gff_reader.c"
62
#line 262 "gff.ragel"
;
}
