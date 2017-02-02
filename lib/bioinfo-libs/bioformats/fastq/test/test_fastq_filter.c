/*
 * test_fastq_filter.c
 *
 *  Created on: Sep 29, 2012
 *      Author: imedina
 */

#include <stdio.h>
#include <stdlib.h>
#include <check.h>

#include "containers/array_list.h"

//#include "vcf_batch.h"
#include "fastq_file.h"
#include "fastq_filter.h"


/*
 * How to get the constants defined below
 * cat CEU.exon.2010_03.genotypes__head400.vcf | grep "#" | wc -l
 * Result: 11 lines, so there are 389 records
 * tail -n 389 CEU.exon.2010_03.genotypes__head400.vcf | cut -f 3 | grep "rs" | wc -l
 * tail -n 389 CEU.exon.2010_03.genotypes__head400.vcf | cut -f 3 | grep "\." | wc -l
 */
#define MAX_RECORDS	    389
#define SNPS_IN_FILE	266


Suite *create_test_suite();


array_list_t *passed, *failed;


/* ******************************
 *       Unchecked fixtures     *
 * ******************************/

void setup_snp(void)
{
	printf("Begin SNP filter testing\n");
}

void teardown_snp(void)
{
	printf("Finished SNP filter testing\n");
}


/* ******************************
 *       Checked fixtures       *
 * ******************************/

void create_passed_failed(void)
{

}

void free_passed_failed(void)
{
	//     list_free_deep(passed, vcf_record_free);
	//     list_free_deep(failed, vcf_record_free);
	array_list_free(passed, fastq_read_free);
	array_list_free(failed, fastq_read_free);
	free(passed);
	free(failed);
}


/* ******************************
 *           Unit tests         *
 * ******************************/

START_TEST (fastq_filter_test)
{

}
END_TEST




/* ******************************
 * 	Main entry point	*
 * ******************************/

int main (int argc, char *argv)
{
	fastq_file_t *file = fastq_fopen("2M_reads_se.fastq");

	Suite *fs = create_test_suite();
	SRunner *fs_runner = srunner_create(fs);
	srunner_run_all(fs_runner, CK_NORMAL);
	int number_failed = srunner_ntests_failed (fs_runner);
	srunner_free (fs_runner);

	fastq_fclose(file);

	return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}


Suite *create_test_suite()
{
	// SNP filter (include and exclude)
	TCase *tc_snp = tcase_create("SNP filters");
	tcase_add_unchecked_fixture(tc_snp, setup_snp, teardown_snp);
	tcase_add_checked_fixture(tc_snp, create_passed_failed, free_passed_failed);
	tcase_add_test(tc_snp, fastq_filter_test);


	// Add test cases to a test suite
	Suite *fs = suite_create("VCF filters");
	suite_add_tcase(fs, tc_snp);

	return fs;
}


