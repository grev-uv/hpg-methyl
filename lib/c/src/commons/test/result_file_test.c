#include <stdio.h>
#include <stdlib.h>
#include <check.h>

#include "../result.h" 

Suite *create_test_suite();
// list_t *read_test_datasuite(vcf_file_t *file);
// void free_test_datasuite(list_t *datasuite, vcf_file_t *file);

// result_item_t *result_item_p;
// result_file_t *result_file_p;


/* ******************************
 *      Unchecked fixtures      *
 * ******************************/
void setup_result_item(void) {
	printf("Begin Result Item tests\n");
}

void teardown_result_item(void) {
	printf("Finished Result Item tests\n\n");
}

void setup_result_file(void) {
	printf("Begin Result File tests\n");
}

void teardown_result_file(void) {
	printf("Finished Result File tests\n\n");
}
/* ******************************
 *          Unit tests         *
 * ******************************/

START_TEST (result_item_test) {
	// char *tags = (char*) malloc(5*sizeof(char));
	// strcpy(tags, "a,b");
	result_item_t *result_item_p = result_item_new("name", "value", "title", RESULT_MESSAGE_DATA_TYPE, "a,b", "group", "context");
	result_item_print(result_item_p);
	fail_unless(strcmp(result_item_p->name, "name") == 0, "Error in 'name' test");
	result_item_free(result_item_p);
}
END_TEST


START_TEST (result_file_test) {
	result_file_t *result_file_p = result_file_new("v0.7", "/tmp/result.xml");

	result_item_t *result_meta_item_p = result_item_new("tool", "variant", "Variant summary", RESULT_MESSAGE_DATA_TYPE, "", "Summary", "context");
	result_item_t *result_input_item_p = result_item_new("file", "test1.vcf", "VCF file", RESULT_FILE_DATA_TYPE, "", "Input", "context");
	result_item_t *result_output_item_p = result_item_new("outdir", "/tmp/variant.txt", "Result file", RESULT_FILE_DATA_TYPE, "TABLE", "Results", "context");
	result_item_t *result_output_item2_p = result_item_new("outdir", "/tmp/log.txt", "Result file", RESULT_FILE_DATA_TYPE, "TABLE", "Results", "context");

	result_add_meta_item(result_meta_item_p, result_file_p);
	result_add_input_item(result_input_item_p, result_file_p);
	result_add_output_item(result_output_item_p, result_file_p);
	result_add_output_item(result_output_item2_p, result_file_p);

	result_file_print(result_file_p);
	//   fail_unless(result_add_meta_item(result_item_p, result_file_p) == 1, "no va1");
	//   fail_unless(result_file_p->num_meta_items == 1, "no va2");

	//   result_item_free(result_input_item_p);
	//   result_item_free(result_output_item_p);

	result_file_write("/tmp/result.xml", result_file_p);

	result_json_file_write("/tmp/result.json", result_file_p);

	result_file_free(result_file_p);
}
END_TEST


Suite *create_test_suite() {
	// SNP filter (include and exclude)
	TCase *tc_result_item = tcase_create("Result Item test");
	tcase_add_unchecked_fixture(tc_result_item, setup_result_item, teardown_result_item);
	tcase_add_test(tc_result_item, result_item_test);

	TCase *tc_result_file = tcase_create("Result File test");
	tcase_add_unchecked_fixture(tc_result_file, setup_result_file, teardown_result_file);
	tcase_add_test(tc_result_file, result_file_test);

	// Add test cases to a test suite
	Suite *fs = suite_create("Job Result invocation");
	suite_add_tcase(fs, tc_result_item);
	suite_add_tcase(fs, tc_result_file);

	return fs;
}


/* ******************************
 *      Main entry point        *
 * ******************************/

int main (int argc, char **argv) {
	//     vcf_file_t *file = vcf_open("CEU.exon.2010_03.genotypes__head400.vcf");
	//     list_t *batches = read_test_datasuite(file);
	//     datasuite = batches->first_p->data_p;

	Suite *fs = create_test_suite();
	SRunner *fs_runner = srunner_create(fs);
	srunner_run_all(fs_runner, CK_NORMAL);
	int number_failed = srunner_ntests_failed (fs_runner);
	srunner_free (fs_runner);

	//     free_test_datasuite(datasuite, file);   // TODO exceeds check timeout
	//     vcf_close(file);

	return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
