#ifndef RESULT_H
#define RESULT_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <libxml/encoding.h>
#include <libxml/xmlwriter.h>
#include <libxml/parser.h>
#include <libxml/tree.h>


// int META_ITEMS_SIZE	= 10;
// int INPUT_ITEMS_SIZE	= 10;
// int OUTPUT_ITEMS_SIZE	= 20;

enum ITEM_TYPE {META_TYPE, INPUT_TYPE, OUTPUT_TYPE};

#define RESULT_FILE_DATA_TYPE			"FILE"
#define RESULT_IMAGE_DATA_TYPE			"IMAGE"
#define RESULT_LINK_DATA_TYPE			"LINK"
#define RESULT_TEXT_DATA_TYPE			"TEXT"
#define RESULT_MESSAGE_DATA_TYPE		"MESSAGE"
#define RESULT_HTML_DATA_TYPE			"HTML"
#define RESULT_DATA_DATA_TYPE			"DATA"

// enum DATA_TYPE {FILE_T, IMAGE_T, LINK_T, TEXT_T, MESSAGE_T, HTML_T, DATA_T};
// enum DATA_TYPE data_type;

typedef struct result_item {
	char *name;
	char *value;
	char *title;
	char *data_type;
	char *tags;
// 	char *style;
	char *group;
	char *context;
	
	enum ITEM_TYPE item_type;
	struct result_item *next_p;
} result_item_t;

// typedef struct result_item_list {
// 	result_item_t *item;
// 	
// 	struct result_item_list *next;
// } result_item_list_t;

typedef struct result_file {
	char *version;
	char *filename;
	
// 	FILE *result_file;
// 	uint16_t num_meta_items;
// 	uint16_t num_input_items;
// 	uint16_t num_output_items;
	
	result_item_t *result_item_first_p;
	result_item_t *result_item_last_p;
	
// 	result_item_t *input_items;
// 	result_item_t *output_items;
} result_file_t;



result_item_t *result_item_new(char *name, char *value, char *title, char *data_type, char *tags, char *group, char *context);

void result_item_free(result_item_t * result_item_p);

void result_item_print(result_item_t * result_item_p);



result_file_t *result_file_new(char *version, char *filename);

void result_file_free(result_file_t *result_file_p);

int result_file_read(result_file_t *result_file);

int result_file_write(char *filename, result_file_t *result_file);

int result_json_file_write(char *filename, result_file_t *result_file);

int result_add_meta_item(result_item_t* item_p, result_file_t *result_file);

int result_add_input_item(result_item_t* item_p, result_file_t *result_file);

int result_add_output_item(result_item_t* item_p, result_file_t *result_file);

void result_file_print(result_file_t *result_file);

#endif
