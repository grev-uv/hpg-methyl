#include "result.h"

void result_json_write_items(result_file_t *result_file, enum ITEM_TYPE item_type, FILE *file);

result_file_t *result_file_open(char *filename);

void result_file_close(result_file_t *result_file);

int result_add_item(result_item_t *item, result_file_t *result_file);

void result_write_items(result_file_t *result_file, enum ITEM_TYPE item_type, xmlNodePtr node);


/**
 * RESULT ITEMS FUNCTIONS
 */
result_item_t *result_item_new(char *name, char *value, char *title, char *data_type, char *tags, char *group, char *context) {
	result_item_t *result_item_p = (result_item_t *)malloc(sizeof(result_item_t));

	result_item_p->name = name;
	result_item_p->value = value;
	result_item_p->title = title;
	result_item_p->data_type = data_type;
	result_item_p->tags = tags;
	result_item_p->group = group;
	result_item_p->context = context;

	//   result_item_p->item_type = item_type;
	result_item_p->next_p = NULL;

	return result_item_p;
}

void result_item_free(result_item_t *result_item_p) {
	free(result_item_p);
}

void result_item_print(result_item_t *result_item_p) {
	printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%i\n", result_item_p->name, result_item_p->value,
			result_item_p->title, result_item_p->data_type, result_item_p->tags, result_item_p->group, result_item_p->context, result_item_p->item_type);
}



/**
 * RESULT FILE FUNCTIONS
 */

result_file_t *result_file_new(char *version, char *filename) {
	result_file_t *result_file_p = (result_file_t *) malloc(sizeof(result_file_t));

	result_file_p->version = version;
	result_file_p->filename = filename;

	result_file_p->result_item_first_p = NULL;
	result_file_p->result_item_last_p = NULL;

	return result_file_p;
}

void result_file_free(result_file_t *result_file_p) {
	result_item_t *result_item_aux_p = result_file_p->result_item_first_p;
	while(result_item_aux_p != NULL) {
		result_file_p->result_item_first_p = result_item_aux_p->next_p;
		result_item_free(result_item_aux_p);
		result_item_aux_p = result_file_p->result_item_first_p;
	}
	free(result_file_p);
}


int result_file_read(result_file_t *result_file) {
	return 1;
}

int result_file_write(char *filename, result_file_t *result_file) {
	xmlDocPtr doc = NULL;       /* document pointer */
	xmlNodePtr root_node = NULL;

	/*
	 * Creates a new document, a node and set it as a root node
	 */
	doc = xmlNewDoc(BAD_CAST "1.0");
	root_node = xmlNewNode(NULL, BAD_CAST "result");
	xmlDocSetRootElement(doc, root_node);

	xmlNodePtr metadata_node = xmlNewChild(root_node, NULL, BAD_CAST "metadata", BAD_CAST NULL);
	result_write_items(result_file, META_TYPE, metadata_node);

	xmlNodePtr input_node = xmlNewChild(root_node, NULL, BAD_CAST "input", BAD_CAST NULL);
	result_write_items(result_file, INPUT_TYPE, input_node);

	xmlNodePtr output_node = xmlNewChild(root_node, NULL, BAD_CAST "output", BAD_CAST NULL);
	result_write_items(result_file, OUTPUT_TYPE, output_node);

	/*
	 * Dumping document to stdio or file
	 */
	xmlSaveFormatFileEnc(filename, doc, "UTF-8", 1);

	/*free the document */
	xmlFreeDoc(doc);

	/*
	 *Free the global variables that may
	 *have been allocated by the parser.
	 */
	xmlCleanupParser();

	/*
	 * this is to debug memory for regression tests
	 */
	xmlMemoryDump();
	return(0);
}

int result_json_file_write(char *filename, result_file_t *result_file) {
	FILE *file = fopen(filename, "w");
	if (file == NULL) {
		printf("I couldn't open %s for writing.\n", filename);
		exit(0);
	}
	fprintf(file, "{\n");
	fprintf(file, "\t\"version\": \"%s\",\n", result_file->version);
	result_json_write_items(result_file, META_TYPE, file);
	fprintf(file, ",\n");
	result_json_write_items(result_file, INPUT_TYPE, file);
	fprintf(file, ",\n");
	result_json_write_items(result_file, OUTPUT_TYPE, file);
	fprintf(file, "\n}\n");
	fclose(file);
	return(0);
}

int result_add_meta_item(result_item_t *meta_item, result_file_t *result_file) {
	meta_item->item_type = META_TYPE;
	return result_add_item(meta_item, result_file);
}

int result_add_input_item(result_item_t* input_item, result_file_t *result_file) {
	input_item->item_type = INPUT_TYPE;
	return result_add_item(input_item, result_file);
}

int result_add_output_item(result_item_t* output_item, result_file_t *result_file) {
	output_item->item_type = OUTPUT_TYPE;
	return result_add_item(output_item, result_file);
}

void result_file_print(result_file_t *result_file) {
	result_item_t *result_item_aux_p = result_file->result_item_first_p;
	while(result_item_aux_p != NULL) {
		result_item_print(result_item_aux_p);
		result_item_aux_p = result_item_aux_p->next_p;
	}
}



/**
 * PRIVATE FUNCTIONS
 */
int result_add_item(result_item_t *item, result_file_t *result_file) {
	item->next_p = NULL;
	if(result_file->result_item_first_p == NULL) {
		result_file->result_item_first_p = item;
	}else {
		result_file->result_item_last_p->next_p = item;
	}
	result_file->result_item_last_p = item;
	return 1;
}


void result_write_items(result_file_t *result_file, enum ITEM_TYPE item_type, xmlNodePtr node) {
	result_item_t *result_item_aux_p = result_file->result_item_first_p;
	while(result_item_aux_p != NULL) {
		if(result_item_aux_p->item_type == item_type) {
			xmlNodePtr node1 = xmlNewChild(node, NULL, BAD_CAST "item", BAD_CAST result_item_aux_p->value);
			xmlNewProp(node1, BAD_CAST "name", BAD_CAST result_item_aux_p->name);
			xmlNewProp(node1, BAD_CAST "title", BAD_CAST result_item_aux_p->title);
			xmlNewProp(node1, BAD_CAST "type", BAD_CAST result_item_aux_p->data_type);
			xmlNewProp(node1, BAD_CAST "tags", BAD_CAST result_item_aux_p->tags);
			xmlNewProp(node1, BAD_CAST "group", BAD_CAST result_item_aux_p->group);
			xmlNewProp(node1, BAD_CAST "context", BAD_CAST result_item_aux_p->context);
		}
		result_item_aux_p = result_item_aux_p->next_p;
	}
}

void result_json_write_items(result_file_t *result_file, enum ITEM_TYPE item_type, FILE *file) {
	if(item_type == 0) {
		fprintf(file, "\t\"%s\": [\n", "metadata");
	}else {
		if(item_type == 1) {
			fprintf(file, "\t\"%s\": [\n", "input");
		}else {
			fprintf(file, "\t\"%s\": [\n", "output");
		}
	}
	result_item_t *result_item_aux_p = result_file->result_item_first_p;
	int first = 0;
	while(result_item_aux_p != NULL) {
		if(result_item_aux_p->item_type == item_type) {
			if(first) {
				fprintf(file, ",\n");
			}
			fprintf(file, "\t\t{\n");
			fprintf(file, "\t\t\"name\": \"%s\",\n", result_item_aux_p->name);
			fprintf(file, "\t\t\"title\": \"%s\",\n", result_item_aux_p->title);
			fprintf(file, "\t\t\"type\": \"%s\",\n", result_item_aux_p->data_type);
			fprintf(file, "\t\t\"tags\": \"%s\",\n", result_item_aux_p->tags);
			fprintf(file, "\t\t\"group\": \"%s\",\n", result_item_aux_p->group);
			fprintf(file, "\t\t\"context\": \"%s\"", result_item_aux_p->context);
			fprintf(file, "\n\t\t}");
			first = 1;
		}
		result_item_aux_p = result_item_aux_p->next_p;
	}
	fprintf(file, "\n\t]");
}
