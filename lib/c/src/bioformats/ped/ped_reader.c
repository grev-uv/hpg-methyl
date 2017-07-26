
#line 1 "ped.ragel"
#include "ped_reader.h"

static size_t lines = 0;
static size_t num_records = 0;
static size_t genotype = 0;

static ped_record_t *current_record;
static ped_batch_t *current_batch;


#line 14 "ped_reader.c"
static const int ped_start = 21;
static const int ped_first_final = 21;
static const int ped_error = 0;

static const int ped_en_main = 21;


#line 201 "ped.ragel"




int ped_ragel_read(list_t *batches_list, size_t batch_size, ped_file_t *file)
{
    int cs;
    char *p = file->data;
    char *pe = p + file->data_len;
    char *eof = pe;
    char *ts;
    int custom_field_count = 0;

    current_batch = ped_batch_new(batch_size);

    
#line 41 "ped_reader.c"
	{
	cs = ped_start;
	}

#line 46 "ped_reader.c"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
case 21:
	switch( (*p) ) {
		case 10: goto st22;
		case 35: goto st16;
	}
	if ( 33 <= (*p) && (*p) <= 126 )
		goto tr36;
	goto tr0;
tr0:
#line 53 "ped.ragel"
	{
        LOG_ERROR_F("Line %zu (%s): Error in 'family' field\n", lines + 1, file->filename);
    }
	goto st0;
tr3:
#line 65 "ped.ragel"
	{
        LOG_ERROR_F("Line %zu (%s): Error in 'individual' field\n", lines + 1, file->filename);
    }
	goto st0;
tr7:
#line 77 "ped.ragel"
	{
        LOG_ERROR_F("Line %zu (%s): Error in 'father' field\n", lines + 1, file->filename);
    }
	goto st0;
tr11:
#line 89 "ped.ragel"
	{
        LOG_ERROR_F("Line %zu (%s): Error in 'mother' field\n", lines + 1, file->filename);
    }
	goto st0;
tr15:
#line 109 "ped.ragel"
	{
        LOG_ERROR_F("Line %zu (%s): Error in 'sex' field\n", lines + 1, file->filename);
    }
	goto st0;
tr19:
#line 124 "ped.ragel"
	{
        LOG_ERROR_F("Line %zu (%s): Error in 'phenotype' field\n", lines + 1, file->filename);
    }
	goto st0;
tr26:
#line 141 "ped.ragel"
	{
        LOG_ERROR_F("Line %zu (%s): Error in 'header' field\n", lines + 1, file->filename);
    }
	goto st0;
tr44:
#line 161 "ped.ragel"
	{
        LOG_ERROR_F("Line %zu (%s): Error in '%s' field\n", lines + 1, file->filename, current_record->custom_field);
    }
	goto st0;
#line 108 "ped_reader.c"
st0:
cs = 0;
	goto _out;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
	if ( (*p) == 10 )
		goto st22;
	goto st0;
tr36:
#line 22 "ped.ragel"
	{
        current_record = create_ped_record();
        genotype = 0;
    }
#line 45 "ped.ragel"
	{
        ts = p;
    }
	goto st1;
st1:
	if ( ++p == pe )
		goto _test_eof1;
case 1:
#line 134 "ped_reader.c"
	if ( (*p) == 9 )
		goto tr1;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto st1;
	goto tr0;
tr1:
#line 49 "ped.ragel"
	{
        set_ped_record_family_id(strndup(ts, p-ts), current_record);
    }
	goto st2;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
#line 150 "ped_reader.c"
	if ( (*p) == 95 )
		goto tr4;
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr4;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr4;
	} else
		goto tr4;
	goto tr3;
tr4:
#line 57 "ped.ragel"
	{
        ts = p;
    }
	goto st3;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
#line 172 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr5;
		case 95: goto st3;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st3;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st3;
	} else
		goto st3;
	goto tr3;
tr5:
#line 61 "ped.ragel"
	{
        set_ped_record_individual_id(strndup(ts, p-ts), current_record);
    }
	goto st4;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
#line 196 "ped_reader.c"
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
#line 69 "ped.ragel"
	{
        ts = p;
    }
	goto st5;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
#line 220 "ped_reader.c"
	if ( (*p) == 9 )
		goto tr10;
	goto tr7;
tr10:
#line 73 "ped.ragel"
	{
        set_ped_record_father_id(strndup(ts, p-ts), current_record);
    }
	goto st6;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
#line 234 "ped_reader.c"
	switch( (*p) ) {
		case 46: goto tr12;
		case 95: goto tr13;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr13;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr13;
	} else
		goto tr13;
	goto tr11;
tr12:
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st7;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
#line 258 "ped_reader.c"
	if ( (*p) == 9 )
		goto tr14;
	goto tr11;
tr14:
#line 85 "ped.ragel"
	{
        set_ped_record_mother_id(strndup(ts, p-ts), current_record);
    }
	goto st8;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
#line 272 "ped_reader.c"
	if ( (*p) == 46 )
		goto tr16;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr17;
	goto tr15;
tr16:
#line 93 "ped.ragel"
	{
        ts = p;
    }
	goto st9;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
#line 288 "ped_reader.c"
	if ( (*p) == 9 )
		goto tr18;
	goto tr15;
tr18:
#line 97 "ped.ragel"
	{
        char *field = strndup(ts, p-ts);
        enum Sex sex = UNKNOWN_SEX;
        if (atoi(field) == 1) {
            sex = MALE;
        } else if (atoi(field) == 2) {
            sex = FEMALE;
        }
        set_ped_record_sex(sex, current_record);
        free(field);    // Not set as ped_record_t variable -> not freed later
    }
	goto st10;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
#line 310 "ped_reader.c"
	switch( (*p) ) {
		case 32: goto tr20;
		case 95: goto tr20;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr20;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr20;
	} else
		goto tr20;
	goto tr19;
tr20:
#line 113 "ped.ragel"
	{
        ts = p;
    }
	goto st23;
tr42:
#line 117 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            set_ped_record_phenotype(field, current_record, file);
        }
    }
#line 145 "ped.ragel"
	{
        custom_field_count = 6;
    }
	goto st23;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
#line 347 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr39;
		case 10: goto tr40;
		case 32: goto tr42;
		case 95: goto st23;
	}
	if ( (*p) < 48 ) {
		if ( 11 <= (*p) && (*p) <= 13 )
			goto tr41;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 90 ) {
			if ( 97 <= (*p) && (*p) <= 122 )
				goto st23;
		} else if ( (*p) >= 65 )
			goto st23;
	} else
		goto st23;
	goto tr19;
tr39:
#line 117 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            set_ped_record_phenotype(field, current_record, file);
        }
    }
#line 145 "ped.ragel"
	{
        custom_field_count = 6;
    }
	goto st24;
tr49:
#line 153 "ped.ragel"
	{
        char* field_name = strndup(ts, p-ts);
        custom_field_count++;
        if (custom_field_count == file->num_field) {
            set_ped_record_custom_field(field_name, current_record, file);
        }
    }
	goto st24;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
#line 393 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto st24;
		case 10: goto tr46;
	}
	if ( (*p) > 13 ) {
		if ( 32 <= (*p) && (*p) <= 126 )
			goto tr48;
	} else if ( (*p) >= 11 )
		goto st26;
	goto tr44;
tr40:
#line 117 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            set_ped_record_phenotype(field, current_record, file);
        }
    }
#line 145 "ped.ragel"
	{
        custom_field_count = 6;
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        if (current_record) {
            add_record_to_ped_batch(current_record, current_batch);
            num_records++;
        }
        current_record = NULL;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st25;
tr46:
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        if (current_record) {
            add_record_to_ped_batch(current_record, current_batch);
            num_records++;
        }
        current_record = NULL;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st25;
tr50:
#line 153 "ped.ragel"
	{
        char* field_name = strndup(ts, p-ts);
        custom_field_count++;
        if (custom_field_count == file->num_field) {
            set_ped_record_custom_field(field_name, current_record, file);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        if (current_record) {
            add_record_to_ped_batch(current_record, current_batch);
            num_records++;
        }
        current_record = NULL;
    }
#line 18 "ped.ragel"
	{
        lines++;
    }
	goto st25;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
#line 499 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr46;
		case 32: goto st26;
	}
	if ( (*p) < 33 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto st26;
	} else if ( (*p) > 34 ) {
		if ( 36 <= (*p) && (*p) <= 126 )
			goto tr36;
	} else
		goto tr36;
	goto tr0;
tr41:
#line 117 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            set_ped_record_phenotype(field, current_record, file);
        }
    }
#line 145 "ped.ragel"
	{
        custom_field_count = 6;
    }
	goto st26;
tr51:
#line 153 "ped.ragel"
	{
        char* field_name = strndup(ts, p-ts);
        custom_field_count++;
        if (custom_field_count == file->num_field) {
            set_ped_record_custom_field(field_name, current_record, file);
        }
    }
	goto st26;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
#line 540 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr46;
		case 32: goto st26;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st26;
	goto st0;
tr48:
#line 149 "ped.ragel"
	{
        ts = p;
    }
	goto st27;
tr52:
#line 153 "ped.ragel"
	{
        char* field_name = strndup(ts, p-ts);
        custom_field_count++;
        if (custom_field_count == file->num_field) {
            set_ped_record_custom_field(field_name, current_record, file);
        }
    }
	goto st27;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
#line 568 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr49;
		case 10: goto tr50;
		case 32: goto tr52;
	}
	if ( (*p) > 13 ) {
		if ( 33 <= (*p) && (*p) <= 126 )
			goto st27;
	} else if ( (*p) >= 11 )
		goto tr51;
	goto tr44;
tr17:
#line 93 "ped.ragel"
	{
        ts = p;
    }
	goto st11;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
#line 590 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr18;
		case 46: goto st12;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st11;
	goto tr15;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st13;
	goto tr15;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
	if ( (*p) == 9 )
		goto tr18;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st13;
	goto tr15;
tr13:
#line 81 "ped.ragel"
	{
        ts = p;
    }
	goto st14;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
#line 624 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr14;
		case 95: goto st14;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st14;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st14;
	} else
		goto st14;
	goto tr11;
tr9:
#line 69 "ped.ragel"
	{
        ts = p;
    }
	goto st15;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
#line 648 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr10;
		case 95: goto st15;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st15;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st15;
	} else
		goto st15;
	goto tr7;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
	switch( (*p) ) {
		case 9: goto st17;
		case 32: goto tr28;
		case 95: goto tr29;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr29;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr29;
	} else
		goto tr29;
	goto tr26;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
	switch( (*p) ) {
		case 32: goto tr29;
		case 95: goto tr29;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr29;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr29;
	} else
		goto tr29;
	goto tr26;
tr29:
#line 128 "ped.ragel"
	{
        ts = p;
    }
	goto st18;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
#line 707 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr30;
		case 10: goto tr31;
		case 32: goto st18;
		case 95: goto st18;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto st18;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st18;
	} else
		goto st18;
	goto tr26;
tr30:
#line 132 "ped.ragel"
	{
        char* field_name = strndup(ts, p-ts);
        custom_field_count++;
        if (file->variable_field && !strcmp(field_name, file->variable_field)) {
            file->num_field = custom_field_count;
        }
        free(field_name);
    }
	goto st19;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
#line 738 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto st19;
		case 10: goto st28;
		case 32: goto tr29;
		case 95: goto tr29;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr29;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr29;
	} else
		goto tr29;
	goto tr26;
tr31:
#line 132 "ped.ragel"
	{
        char* field_name = strndup(ts, p-ts);
        custom_field_count++;
        if (file->variable_field && !strcmp(field_name, file->variable_field)) {
            file->num_field = custom_field_count;
        }
        free(field_name);
    }
	goto st28;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
#line 769 "ped_reader.c"
	if ( (*p) == 10 )
		goto st22;
	if ( (*p) > 34 ) {
		if ( 36 <= (*p) && (*p) <= 126 )
			goto tr36;
	} else if ( (*p) >= 33 )
		goto tr36;
	goto tr0;
tr28:
#line 128 "ped.ragel"
	{
        ts = p;
    }
	goto st20;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
#line 788 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr30;
		case 10: goto tr31;
		case 32: goto tr29;
		case 95: goto tr29;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr29;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr29;
	} else
		goto tr29;
	goto tr26;
	}
	_test_eof22: cs = 22; goto _test_eof; 
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
	_test_eof23: cs = 23; goto _test_eof; 
	_test_eof24: cs = 24; goto _test_eof; 
	_test_eof25: cs = 25; goto _test_eof; 
	_test_eof26: cs = 26; goto _test_eof; 
	_test_eof27: cs = 27; goto _test_eof; 
	_test_eof11: cs = 11; goto _test_eof; 
	_test_eof12: cs = 12; goto _test_eof; 
	_test_eof13: cs = 13; goto _test_eof; 
	_test_eof14: cs = 14; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 
	_test_eof16: cs = 16; goto _test_eof; 
	_test_eof17: cs = 17; goto _test_eof; 
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof28: cs = 28; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 24: 
	case 25: 
	case 26: 
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        if (current_record) {
            add_record_to_ped_batch(current_record, current_batch);
            num_records++;
        }
        current_record = NULL;
    }
	break;
	case 1: 
#line 53 "ped.ragel"
	{
        LOG_ERROR_F("Line %zu (%s): Error in 'family' field\n", lines + 1, file->filename);
    }
	break;
	case 2: 
	case 3: 
#line 65 "ped.ragel"
	{
        LOG_ERROR_F("Line %zu (%s): Error in 'individual' field\n", lines + 1, file->filename);
    }
	break;
	case 4: 
	case 5: 
	case 15: 
#line 77 "ped.ragel"
	{
        LOG_ERROR_F("Line %zu (%s): Error in 'father' field\n", lines + 1, file->filename);
    }
	break;
	case 6: 
	case 7: 
	case 14: 
#line 89 "ped.ragel"
	{
        LOG_ERROR_F("Line %zu (%s): Error in 'mother' field\n", lines + 1, file->filename);
    }
	break;
	case 8: 
	case 9: 
	case 11: 
	case 12: 
	case 13: 
#line 109 "ped.ragel"
	{
        LOG_ERROR_F("Line %zu (%s): Error in 'sex' field\n", lines + 1, file->filename);
    }
	break;
	case 10: 
#line 124 "ped.ragel"
	{
        LOG_ERROR_F("Line %zu (%s): Error in 'phenotype' field\n", lines + 1, file->filename);
    }
	break;
	case 16: 
	case 17: 
	case 18: 
	case 19: 
	case 20: 
#line 141 "ped.ragel"
	{
        LOG_ERROR_F("Line %zu (%s): Error in 'header' field\n", lines + 1, file->filename);
    }
	break;
	case 27: 
#line 153 "ped.ragel"
	{
        char* field_name = strndup(ts, p-ts);
        custom_field_count++;
        if (custom_field_count == file->num_field) {
            set_ped_record_custom_field(field_name, current_record, file);
        }
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        if (current_record) {
            add_record_to_ped_batch(current_record, current_batch);
            num_records++;
        }
        current_record = NULL;
    }
	break;
	case 23: 
#line 117 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            set_ped_record_phenotype(field, current_record, file);
        }
    }
#line 145 "ped.ragel"
	{
        custom_field_count = 6;
    }
#line 27 "ped.ragel"
	{
        // If batch is full, add to the list of batches and create a new, empty one
        if (ped_batch_is_full(current_batch))
        {
            list_item_t *item = list_item_new(num_records, 1, current_batch);
            list_insert_item(item, batches_list);
            LOG_DEBUG_F("Batch added - %zu records\n", current_batch->length);
            current_batch = ped_batch_new(batch_size);
        }

        // Add current record to current batch
        if (current_record) {
            add_record_to_ped_batch(current_record, current_batch);
            num_records++;
        }
        current_record = NULL;
    }
	break;
#line 973 "ped_reader.c"
	}
	}

	_out: {}
	}

#line 221 "ped.ragel"
 

    // Insert the last batch
    if (!ped_batch_is_empty(current_batch))
    {
        list_item_t *item = list_item_new(num_records, 1, current_batch); 
        list_insert_item(item, batches_list);
        LOG_DEBUG_F("Batch added - %zu records (last)\n", current_batch->length);
    }

    if ( cs < 
#line 992 "ped_reader.c"
21
#line 231 "ped.ragel"
 ) 
    {
        LOG_ERROR("The file was not successfully read\n");
        LOG_INFO_F("Last state is %d, but %d was expected\n", 
                cs, 
#line 1000 "ped_reader.c"
21
#line 235 "ped.ragel"
);
    } 

    LOG_INFO_F("PED records read = %zu\n", num_records);

    return cs < 
#line 1009 "ped_reader.c"
21
#line 240 "ped.ragel"
;
}
