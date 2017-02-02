
#line 1 "ped.ragel"
#include "ped_reader.h"

static size_t lines = 0;
static size_t num_records = 0;
static size_t genotype = 0;

static ped_record_t *current_record;
static ped_batch_t *current_batch;


#line 14 "ped_reader.c"
static const int ped_start = 16;
static const int ped_first_final = 16;
static const int ped_error = 0;

static const int ped_en_main = 16;


#line 165 "ped.ragel"




int ped_ragel_read(list_t *batches_list, size_t batch_size, ped_file_t *file)
{
    int cs;
    char *p = file->data;
    char *pe = p + file->data_len;
    char *eof = pe;
    char *ts, *te;
    int stack[4];
    int top, act;

    current_batch = ped_batch_new(batch_size);

    
#line 40 "ped_reader.c"
	{
	cs = ped_start;
	}

#line 45 "ped_reader.c"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
case 16:
	if ( (*p) == 10 )
		goto st17;
	if ( 33 <= (*p) && (*p) <= 126 )
		goto tr27;
	goto tr0;
tr0:
#line 53 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
	goto st0;
tr3:
#line 65 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
	goto st0;
tr7:
#line 79 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
	goto st0;
tr11:
#line 93 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	goto st0;
tr15:
#line 113 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	goto st0;
tr19:
#line 130 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	goto st0;
#line 93 "ped_reader.c"
st0:
cs = 0;
	goto _out;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
	if ( (*p) == 10 )
		goto st17;
	goto st0;
tr27:
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
#line 119 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr1;
		case 32: goto tr1;
	}
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
#line 137 "ped_reader.c"
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
#line 159 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr5;
		case 32: goto tr5;
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
#line 184 "ped_reader.c"
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
#line 208 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr10;
		case 32: goto tr10;
	}
	goto tr7;
tr10:
#line 73 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_father_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st6;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
#line 226 "ped_reader.c"
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
#line 83 "ped.ragel"
	{
        ts = p;
    }
	goto st7;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
#line 250 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr14;
		case 32: goto tr14;
	}
	goto tr11;
tr14:
#line 87 "ped.ragel"
	{
        if (strncmp("0", ts, 1)) {
            set_ped_record_mother_id(strndup(ts, p-ts), current_record);
        }
    }
	goto st8;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
#line 268 "ped_reader.c"
	if ( (*p) == 46 )
		goto tr16;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr17;
	goto tr15;
tr16:
#line 97 "ped.ragel"
	{
        ts = p;
    }
	goto st9;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
#line 284 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr18;
		case 32: goto tr18;
	}
	goto tr15;
tr18:
#line 101 "ped.ragel"
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
#line 308 "ped_reader.c"
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr20;
	goto tr19;
tr20:
#line 117 "ped.ragel"
	{
        ts = p;
    }
	goto st18;
tr32:
#line 121 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st18;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
#line 333 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr30;
		case 46: goto tr31;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr32;
	goto tr29;
tr29:
#line 121 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st19;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
#line 356 "ped_reader.c"
	if ( (*p) == 10 )
		goto tr34;
	goto st19;
tr30:
#line 121 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
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
	goto st20;
tr34:
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
	goto st20;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
#line 421 "ped_reader.c"
	switch( (*p) ) {
		case 10: goto tr34;
		case 32: goto st21;
	}
	if ( (*p) > 13 ) {
		if ( 33 <= (*p) && (*p) <= 126 )
			goto tr27;
	} else if ( (*p) >= 9 )
		goto st21;
	goto tr0;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
	switch( (*p) ) {
		case 10: goto tr34;
		case 32: goto st21;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st21;
	goto st0;
tr31:
#line 121 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st22;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
#line 458 "ped_reader.c"
	if ( (*p) == 10 )
		goto tr34;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st23;
	goto st19;
tr37:
#line 121 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
        }
    }
	goto st23;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
#line 479 "ped_reader.c"
	if ( (*p) == 10 )
		goto tr30;
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr37;
	goto tr29;
tr17:
#line 97 "ped.ragel"
	{
        ts = p;
    }
	goto st11;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
#line 495 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr18;
		case 32: goto tr18;
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
	switch( (*p) ) {
		case 9: goto tr18;
		case 32: goto tr18;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st13;
	goto tr15;
tr13:
#line 83 "ped.ragel"
	{
        ts = p;
    }
	goto st14;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
#line 532 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr14;
		case 32: goto tr14;
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
#line 557 "ped_reader.c"
	switch( (*p) ) {
		case 9: goto tr10;
		case 32: goto tr10;
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
	}
	_test_eof17: cs = 17; goto _test_eof; 
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
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 
	_test_eof21: cs = 21; goto _test_eof; 
	_test_eof22: cs = 22; goto _test_eof; 
	_test_eof23: cs = 23; goto _test_eof; 
	_test_eof11: cs = 11; goto _test_eof; 
	_test_eof12: cs = 12; goto _test_eof; 
	_test_eof13: cs = 13; goto _test_eof; 
	_test_eof14: cs = 14; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 19: 
	case 20: 
	case 21: 
	case 22: 
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
        printf("Line %zu (%s): Error in 'family' field\n", lines, file->filename);
    }
	break;
	case 2: 
	case 3: 
#line 65 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'individual' field\n", lines, file->filename);
    }
	break;
	case 4: 
	case 5: 
	case 15: 
#line 79 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'father' field\n", lines, file->filename);
    }
	break;
	case 6: 
	case 7: 
	case 14: 
#line 93 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'mother' field\n", lines, file->filename);
    }
	break;
	case 8: 
	case 9: 
	case 11: 
	case 12: 
	case 13: 
#line 113 "ped.ragel"
	{
        printf("Line %zu (%s): Error in 'sex' field\n", lines, file->filename);
    }
	break;
	case 10: 
#line 130 "ped.ragel"
	{
        printf("Line %zu: Error in 'phenotype' field\n", lines);
    }
	break;
	case 18: 
	case 23: 
#line 121 "ped.ragel"
	{
        if (strncmp(".", ts, 1)) {
            char *field = strndup(ts, p-ts);
            float phenotype = atof(field);
            set_ped_record_phenotype(phenotype, current_record);
            free(field);
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
#line 698 "ped_reader.c"
	}
	}

	_out: {}
	}

#line 184 "ped.ragel"
 

    // Insert the last batch
    if (!ped_batch_is_empty(current_batch))
    {
        list_item_t *item = list_item_new(num_records, 1, current_batch); 
        list_insert_item(item, batches_list);
        LOG_DEBUG_F("Batch added - %zu records (last)\n", current_batch->length);
    }

    if ( cs < 
#line 717 "ped_reader.c"
16
#line 194 "ped.ragel"
 ) 
    {
        LOG_ERROR("The file was not successfully read\n");
        LOG_INFO_F("Last state is %d, but %d was expected\n", 
                cs, 
#line 725 "ped_reader.c"
16
#line 198 "ped.ragel"
);
    } 

    LOG_INFO_F("PED records read = %zu\n", num_records);

    return cs < 
#line 734 "ped_reader.c"
16
#line 203 "ped.ragel"
;
}
