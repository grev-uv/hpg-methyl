#include <stdio.h>
#include <stdlib.h>
#include <check.h>

#include <bioformats/features/region/region.h>
#include <bioformats/features/region/region_table.h>
#include <bioformats/features/region/region_table_utils.h>


Suite *create_test_suite();

region_table_t *region_table;

char *chromosomes[6] = { "1", "2", "20", "22", "3", "X" };

/* ******************************
 *      Unchecked fixtures      *
 * ******************************/

void setup_gff(void) {
    region_table = parse_regions_from_gff_file(strdup("test_parse_regions_01.gff"), "http://ws.bioinfo.cipf.es", "hsa", "latest");
}

void teardown_gff(void) {
    free_table(region_table);
}


/* ******************************
 *          Unit tests         *
 * ******************************/

START_TEST (gff_file_chromosomes)
{
    // Check number of chromosomes
    fail_unless(cp_hashtable_count(region_table->storage) == 6, "There must be 6 chromosomes");
    
    // Check identity and # of regions of each chromosome
    for (int i = 0; i < 6; i++) {
        fail_if(!get_chromosome(chromosomes[i], region_table), "A chromosome has not been inserted");
    }
    
    // Check number of regions per chromosome
    fail_unless(count_regions_in_chromosome("1", region_table) == 2, "Chromosome 1 must have 2 regions");
    fail_unless(count_regions_in_chromosome("2", region_table) == 3, "Chromosome 2 must have 3 regions");
    fail_unless(count_regions_in_chromosome("3", region_table) == 2, "Chromosome 3 must have 2 regions");
    fail_unless(count_regions_in_chromosome("20", region_table) == 1, "Chromosome 20 must have 1 region");
    fail_unless(count_regions_in_chromosome("22", region_table) == 2, "Chromosome 11 must have 2 regions");
    fail_unless(count_regions_in_chromosome("X", region_table) == 1, "Chromosome X must have 1 regions");
}
END_TEST

START_TEST (gff_file_regions)
{
    region_t region;
    
    // Check regions in chr1
    region.chromosome = "1";
    region.start_position = 10000000;
    region.end_position = 10001000;
    fail_if(!find_region(&region, region_table), "Region 1:10000000-10001000 must have been inserted");
    
    region.start_position = 10010000;
    region.end_position = 10010100;
    fail_if(!find_region(&region, region_table), "Region 1:10010000-10010100 must have been inserted");
    
    // Check regions in chr2
    region.chromosome = "2";
    region.start_position = 20020000;
    region.end_position = 20025000;
    fail_if(!find_region(&region, region_table), "Region 2:20020000-20025000 must have been inserted");
    
    region.start_position = 20000000;
    region.end_position = 20001000;
    fail_if(!find_region(&region, region_table), "Region 2:20000000-20001000 must have been inserted");
    
    region.start_position = 20010000;
    region.end_position = 20010100;
    fail_if(!find_region(&region, region_table), "Region 2:20010000-20010100 must have been inserted");
    
    // Check regions in chr3
    region.chromosome = "3";
    region.start_position = 30020000;
    region.end_position = 30025000;
    fail_if(!find_region(&region, region_table), "Region 3:30020000-30025000 must have been inserted");
    
    region.start_position = 30000000;
    region.end_position = 30001000;
    fail_if(!find_region(&region, region_table), "Region 3:30000000-30001000 must have been inserted");
    
    // Check regions in chr20
    region.chromosome = "20";
    region.start_position = 20010000;
    region.end_position = 20010100;
    fail_if(!find_region(&region, region_table), "Region 20:20010000-20010100 must have been inserted");
    
    // Check regions in chr22
    region.chromosome = "22";
    region.start_position = 22000000;
    region.end_position = 22001000;
    fail_if(!find_region(&region, region_table), "Region 22:22000000-22001000 must have been inserted");
    
    region.start_position = 22010000;
    region.end_position = 22010100;
    fail_if(!find_region(&region, region_table), "Region 22:22010000-22010100 must have been inserted");
    
    // Check regions in chrX
    region.chromosome = "X";
    region.start_position = 95020000;
    region.end_position = 95025000;
    fail_if(!find_region(&region, region_table), "Region X:95020000-95025000 must have been inserted");
    
}
END_TEST



/* ******************************
 *        Main entry point      *
 * ******************************/

int main (int argc, char *argv) {    
    Suite *fs = create_test_suite();
    SRunner *fs_runner = srunner_create(fs);
    srunner_run_all(fs_runner, CK_NORMAL);
    int number_failed = srunner_ntests_failed (fs_runner);
    srunner_free (fs_runner);
    
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}


Suite *create_test_suite() {
    TCase *tc_gff_file = tcase_create("Regions from GFF file");
    tcase_add_unchecked_fixture(tc_gff_file, setup_gff, teardown_gff);
    tcase_add_test(tc_gff_file, gff_file_chromosomes);
    tcase_add_test(tc_gff_file, gff_file_regions);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("Regions parsing");
    suite_add_tcase(fs, tc_gff_file);
    
    return fs;
}

