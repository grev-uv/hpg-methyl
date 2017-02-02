#include <stdlib.h>
#include <check.h>

#include "region_table.h"


Suite *create_test_suite();


static region_table_t *table;

static region_t *reg_2, *reg_3, *reg_4_1, *reg_4_2, *reg_X_1, *reg_X_2, *reg_Y;


/* ******************************
 * 	Unchecked fixtures	*
 * ******************************/

void setup_region_table(void)
{
	table = create_table(NULL);
}

void teardown_region_table(void)
{
	cp_hashtable_destroy(table->storage);
}

void setup_regions(void)
{
	reg_2 = (region_t*) malloc (sizeof(region_t));
	reg_2->chromosome = "2";
	reg_2->start_position = 10000;
	reg_2->end_position = 20000;
	
	reg_3 = (region_t*) malloc (sizeof(region_t));
	reg_3->chromosome = "3";
	reg_3->start_position = 30000;
	reg_3->end_position = 40000;
	
	reg_4_1 = (region_t*) malloc (sizeof(region_t));
	reg_4_1->chromosome = "4";
	reg_4_1->start_position = 60000;
	reg_4_1->end_position = 80000;
	
	reg_4_2 = (region_t*) malloc (sizeof(region_t));
	reg_4_2->chromosome = "4";
	reg_4_2->start_position = 200000;
	reg_4_2->end_position = 250000;
	
	reg_X_1 = (region_t*) malloc (sizeof(region_t));
	reg_X_1->chromosome = "X";
	reg_X_1->start_position = 1000000;
	reg_X_1->end_position = 1500000;
	
	reg_X_2 = (region_t*) malloc (sizeof(region_t));
	reg_X_2->chromosome = "X";
	reg_X_2->start_position = 1000000; // Same start_position
	reg_X_2->end_position = 1800000;   // Different end_position
	
	reg_Y = (region_t*) malloc (sizeof(region_t));
	reg_Y->chromosome = "Y";
	reg_Y->start_position = 2000000;
	reg_Y->end_position = 3000000;
}

void teardown_regions(void) { }


/* ******************************
 * 	     Unit tests	        *
 * ******************************/

/* Data structure initialization */

START_TEST(table_structure)
{
	fail_unless(table->max_chromosomes == 25, "The default number of chromosomes for HSA is 25");
	fail_if(cp_hashtable_count(table->storage) > 0, "There must be no elements in the table after its creation");
}
END_TEST


/* Data structure manipulation */

START_TEST(insert_region_and_chromosome)
{
	// Insert a chromosome not via a region
	fail_if(insert_chromosome("1", table) == NULL, "Insertion of chromosome 1 must be successfully performed");
	fail_unless(cp_hashtable_count(table->storage) == 1, "There must be one element in the chromosome table");
	
	// Insert a region in chromosome 2
	fail_unless(insert_region(reg_2, table) == 0, "Insertion of region 2:10000-20000 must be successfully performed");
	fail_unless(contains_region(reg_2, table) == 1, "Region 2:10000-20000 must be found after its insertion");
	fail_unless(cp_hashtable_count(table->storage) == 2, "There must be 2 elements in the chromosome table");
	
	// Insert the same chromosome the region defined
	fail_if(insert_chromosome("2", table) == NULL, "Insertion of chromosome 2 must be successfully performed");
	fail_unless(cp_hashtable_count(table->storage) == 2, "There must be 2 elements in the chromosome table");
	fail_unless(contains_region(reg_2, table) == 1, "Region 2:10000-20000 must not be deleted after trying to re-insert chr2");
}
END_TEST

START_TEST(insert_several_regions)
{
	// Insert regions in different chromosomes, in the same chromosome and in the same start position
	fail_unless(insert_region(reg_2, table) == 0, "Insertion of region in chr2 must be successfully performed");
	fail_unless(insert_region(reg_3, table) == 0, "Insertion of region in chr3 must be successfully performed");
	fail_unless(insert_region(reg_4_1, table) == 0, "Insertion of region in chr4, position 60K must be successfully performed");
	fail_unless(insert_region(reg_4_2, table) == 0, "Insertion of region in chr4, position 200-250K must be successfully performed");
	fail_unless(insert_region(reg_X_1, table) == 0, "Insertion of region in chrX, position 1M-1.5M must be successfully performed");
	fail_unless(insert_region(reg_X_2, table) == 0, "Insertion of region in chrX, position 1M-1.8M must be successfully performed");
	fail_unless(insert_region(reg_Y, table) == 0, "Insertion of region in chrY must be successfully performed");
	
// 	void **keys = cp_hashtable_get_keys(table->storage);
// 	for (int i = 0; i < cp_hashtable_count(table->storage); i++)
// 	{
// 		printf("%d = %s\n", i, (char*) keys[i]);
// 	}
	// Check number of elements inserted in the table and each tree
	fail_unless(cp_hashtable_count(table->storage) == 5, "There must be 5 elements in the chromosome table");
	
	fail_unless(cp_avltree_count(get_chromosome("2", table)) == 1, "There must be 1 element(s) in chr2");
	fail_unless(cp_avltree_count(get_chromosome("3", table)) == 1, "There must be 1 element(s) in chr3");
	fail_unless(cp_avltree_count(get_chromosome("4", table)) == 2, "There must be 2 element(s) in chr4");
	fail_unless(cp_avltree_count(get_chromosome("X", table)) == 1, "There must be 1 element(s) in chrX");
	fail_unless(cp_avltree_count(get_chromosome("Y", table)) == 1, "There must be 1 element(s) in chrY");
	
	fail_unless(cp_vector_size(cp_avltree_get(get_chromosome("X", table), &(reg_X_2->start_position))) == 2, "There must be 2 element(s) in X:1000000");
}
END_TEST

START_TEST(remove_regions_and_chromosomes)
{
	// Insert regions
	fail_unless(insert_region(reg_2, table) == 0, "Insertion of region in chr2 must be successfully performed");
	fail_unless(insert_region(reg_3, table) == 0, "Insertion of region in chr3 must be successfully performed");
	fail_unless(insert_region(reg_4_1, table) == 0, "Insertion of region in chr4, position 60K must be successfully performed");
	fail_unless(insert_region(reg_4_2, table) == 0, "Insertion of region in chr4, position 200-250K must be successfully performed");
	fail_unless(insert_region(reg_X_1, table) == 0, "Insertion of region in chrX, position 1M-1.5M must be successfully performed");
	fail_unless(insert_region(reg_X_2, table) == 0, "Insertion of region in chrX, position 1M-1.8M must be successfully performed");
	fail_unless(insert_region(reg_Y, table) == 0, "Insertion of region in chrY must be successfully performed");
	
	fail_unless(cp_hashtable_count(table->storage) == 5, "There must be 5 elements in the chromosome table");
	
	// Delete regions
	fail_if(remove_region(reg_3, table) == NULL, "Deletion of region in chr3 must be successfully performed");
	fail_unless(cp_avltree_count(get_chromosome("3", table)) == 0, "There must be 0 element(s) in chr3");
	
	fail_if(remove_region(reg_4_2, table) == NULL, "Deletion of region in chr4, position 200-250K must be successfully performed");
	fail_unless(cp_avltree_count(get_chromosome("4", table)) == 1, "There must be 1 element(s) in chr4");
	
	fail_if(remove_region(reg_X_1, table) == NULL, "Deletion of region in chrX, position 1M-1.5M must be successfully performed");
	fail_unless(cp_avltree_count(get_chromosome("X", table)) == 1, "There must be 1 element(s) in chrX");
	
	fail_if(remove_region(reg_X_2, table) == NULL, "Deletion of region in chrX, position 1M-1.8M must be successfully performed");
	fail_unless(cp_avltree_count(get_chromosome("X", table)) == 0, "There must be 0 element(s) in chrX");
	
	fail_unless(cp_hashtable_count(table->storage) == 5, "There must be 5 elements in the chromosome table");
	
	// Delete chromosomes
	fail_if(remove_chromosome("3", table) == NULL, "Deletion of chr3 must be successfully performed");
	fail_unless(cp_hashtable_count(table->storage) == 4, "chr3: There must be 4 elements in the chromosome table");
	
	fail_if(remove_chromosome("4", table) == NULL, "Deletion of chr4 must be successfully performed");
	fail_unless(cp_hashtable_count(table->storage) == 3, "chr4: There must be 3 elements in the chromosome table");
}
END_TEST


/* Data structure search */

START_TEST(search_region)
{
	// Insert regions
	fail_unless(insert_region(reg_2, table) == 0, "Insertion of region in chr2 must be successfully performed");
	fail_unless(insert_region(reg_3, table) == 0, "Insertion of region in chr3 must be successfully performed");
	fail_unless(insert_region(reg_4_1, table) == 0, "Insertion of region in chr4, position 60-80K must be successfully performed");
	fail_unless(insert_region(reg_4_2, table) == 0, "Insertion of region in chr4, position 200-250K must be successfully performed");
	fail_unless(insert_region(reg_X_1, table) == 0, "Insertion of region in chrX, position 1M-1.5M must be successfully performed");
	fail_unless(insert_region(reg_X_2, table) == 0, "Insertion of region in chrX, position 1M-1.8M must be successfully performed");
	fail_unless(insert_region(reg_Y, table) == 0, "Insertion of region in chrY must be successfully performed");
	
	fail_unless(cp_hashtable_count(table->storage) == 5, "There must be 5 elements in the chromosome table");
	
	// TODO Find regions
	
	// Region in the limits of chr3 tree
	region_t *f_reg_3 = (region_t*) malloc (sizeof(region_t));
	f_reg_3->chromosome = "3";
	f_reg_3->start_position = 30000;
	f_reg_3->end_position = 40000;
	fail_if(find_region(f_reg_3, table) == 0, "Region 3:30000-40000 must be found");
	
	// Region contained in chrY tree
	region_t *f_reg_Y = (region_t*) malloc (sizeof(region_t));
	f_reg_Y->chromosome = "Y";
	f_reg_Y->start_position = 2100000;
	f_reg_Y->end_position = 2200000;
	fail_if(find_region(f_reg_Y, table) == 0, "Region Y:2100000-2200000 must be found");
	
	// Region before, in between and after the ones in chr4
	region_t *f_reg_4 = (region_t*) malloc (sizeof(region_t));
	f_reg_4->chromosome = "4";
	f_reg_4->start_position = 30000;
	f_reg_4->end_position = 70000;
	fail_unless(find_region(f_reg_4, table) == 0, "Region 4:30000-70000 must not be found");
	
	f_reg_4->start_position = 30000;
	f_reg_4->end_position = 90000;
	fail_unless(find_region(f_reg_4, table) == 0, "Region 4:30000-90000 must not be found");
	
	f_reg_4->start_position = 100000;
	f_reg_4->end_position = 190000;
	fail_unless(find_region(f_reg_4, table) == 0, "Region 4:100000-190000 must not be found");
	
	f_reg_4->start_position = 160000;
	f_reg_4->end_position = 210000;
	fail_unless(find_region(f_reg_4, table) == 0, "Region 4:100000-210000 must not be found");
	
	f_reg_4->start_position = 240000;
	f_reg_4->end_position = 270000;
	fail_unless(find_region(f_reg_4, table) == 0, "Region 4:240000-270000 must not be found");
	
	f_reg_4->start_position = 280000;
	f_reg_4->end_position = 300000;
	fail_unless(find_region(f_reg_4, table) == 0, "Region 4:280000-300000 must not be found");
	
	
	free(f_reg_3);
	free(f_reg_4);
	free(f_reg_Y);
}
END_TEST


/* ******************************
 * 	Main entry point	*
 * ******************************/


int main(int argc, char *argv[])
{
	Suite *fs = create_test_suite();
	SRunner *fs_runner = srunner_create(fs);
	srunner_run_all(fs_runner, CK_NORMAL);
	int number_failed = srunner_ntests_failed(fs_runner);
	srunner_free (fs_runner);
	
	return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Suite *create_test_suite()
{
	// Creation of the table for storing regions
	TCase *tc_table = tcase_create("Table creation");
	tcase_add_unchecked_fixture(tc_table, setup_region_table, teardown_region_table);
	tcase_add_test(tc_table, table_structure);
	
	// Manipulation of the data structure
	TCase *tc_manipulation = tcase_create("Data structure manipulation");
	tcase_add_unchecked_fixture(tc_manipulation, setup_region_table, teardown_region_table);
	tcase_add_checked_fixture(tc_manipulation, setup_regions, teardown_regions);
	tcase_add_test(tc_manipulation, insert_region_and_chromosome);
	tcase_add_test(tc_manipulation, insert_several_regions);
	tcase_add_test(tc_manipulation, remove_regions_and_chromosomes);
	
	// Region searching
	TCase *tc_searching = tcase_create("Searching in region data structure");
	tcase_add_unchecked_fixture(tc_searching, setup_region_table, teardown_region_table);
	tcase_add_checked_fixture(tc_searching, setup_regions, teardown_regions);
	tcase_add_test(tc_searching, search_region);
	
	// Add test cases to a test suite
	Suite *fs = suite_create("Region searching table");
	suite_add_tcase(fs, tc_table);
	suite_add_tcase(fs, tc_manipulation);
	suite_add_tcase(fs, tc_searching);
	
	return fs;
}
