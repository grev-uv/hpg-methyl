#include <stdio.h>
#include <stdlib.h>
#include <check.h>

#include "containers/array_list.h" 
#include "containers/list.h" 

#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_file_structure.h>
#include <bioformats/vcf/vcf_filters.h>
#include <bioformats/vcf/vcf_util.h>


/*
 * How to get the constants defined below
 * cat CEU.exon.2010_03.genotypes__head400.vcf | grep "#" | wc -l
 * Result: 11 lines, so there are 389 records
 * tail -n 389 CEU.exon.2010_03.genotypes__head400.vcf | cut -f 3 | grep "rs" | wc -l
 * tail -n 389 CEU.exon.2010_03.genotypes__head400.vcf | cut -f 3 | grep "\." | wc -l
 */
#define MAX_RECORDS	    389
#define SNPS_IN_FILE	266

static char *url = "http://ws.bioinfo.cipf.es";
static char *species = "hsa";
static char *version = "latest";


Suite *create_test_suite();
void read_test_datasuite(vcf_file_t *file);
//variant_stats_t ** <NOMBRE_FUNCION_TO_MOLONA>(array_list_t* datasuite);

array_list_t *datasuite, *missing_values_datasuite, *num_alleles_datasuite, *quality_datasuite;
array_list_t *passed, *failed;
variant_stats_t **input_stats;

filter_t *coverage_f, *quality_f, *missing_values_f, *num_alleles_f, *region_f, *snp_f;
filter_chain *chain;



/* ******************************
 *       Unchecked fixtures     *
 * ******************************/

void setup_snp(void)
{
	printf("Begin SNP filter testing\n");
	snp_f = snp_filter_new(1);
}

void teardown_snp(void)
{
	printf("Finished SNP filter testing\n");
	snp_f->free_func(snp_f);
}

void setup_region(void)
{
	printf("Begin region filter testing\n");
}

void teardown_region(void)
{
	printf("Finished region filter testing\n");
}

void setup_quality(void)
{
    printf("Begin quality filter testing\n");
    quality_f = quality_filter_new(30);
}

void teardown_quality(void)
{
    printf("Finished quality filter testing\n");
    quality_f->free_func(quality_f);
}

void setup_coverage(void)
{
    printf("Begin coverage filter testing\n");
    coverage_f = coverage_filter_new(50);
}

void teardown_coverage(void)
{
    printf("Finished coverage filter testing\n");
    coverage_f->free_func(coverage_f);
}

void setup_missing_values(void)
{
    printf("Begin missing values filter testing\n");
    missing_values_f = missing_values_filter_new(0.2);
}

void teardown_missing_values(void)
{
    printf("Finished missing values filter testing\n");
    missing_values_f->free_func(missing_values_f);
}

void setup_num_alleles(void)
{
    printf("Begin allele count filter testing\n");
    num_alleles_f = num_alleles_filter_new(2);
}

void teardown_num_alleles(void)
{
    printf("Finished allele count filter testing\n");
    num_alleles_f->free_func(num_alleles_f);
}

void setup_snp_region(void) {
	printf("Begin SNP + region filter chain testing\n");
	
	snp_f = snp_filter_new(1);
	char input[] = "1:1000000-5000000";
	region_f = region_filter_new(input, 0, NULL, url, species, version);
	
                            
	chain = add_to_filter_chain(region_f, chain);
	chain = add_to_filter_chain(snp_f, chain);
}

void teardown_snp_region(void) {
	printf("Finished SNP + region filter chain testing\n");
	
	snp_f->free_func(snp_f);
	region_f->free_func(region_f);
//     free_filter_chain(chain);
}

/* ******************************
 *       Checked fixtures       *
 * ******************************/

void create_passed_failed(void) {
    failed = array_list_new(MAX_RECORDS, 1, COLLECTION_MODE_ASYNCHRONIZED);
}

void free_passed_failed(void) {
    array_list_free(passed, NULL);
    array_list_free(failed, NULL);
}


/* ******************************
 *           Unit tests         *
 * ******************************/

START_TEST (coverage_basic) {
    ((coverage_filter_args*) coverage_f->args)->min_coverage = 3000;
    passed = coverage_f->filter_func(datasuite, failed, NULL, coverage_f->name, coverage_f->args);
    
    fail_unless(passed->size == 364, "C3000: The number of coverages found is not correct");
    // size(accepted + rejected) = size(whole input)
    fail_unless(passed->size + failed->size == datasuite->size,
            "C3000: The sum of the number of accepted and rejected records must be the same as the input records");
    
    // no accepted ID < min_qual
    for (int i = 0; i < passed->size; i++) {
        vcf_record_t *record = passed->items[i];
        char *aux_info = strndup(record->info, record->info_len);
        char *dp = get_field_value_in_info("DP", aux_info);
        fail_if(dp != NULL && atoi(dp) < 3000, "C3000: An accepted record can't have less coverage than specified");
    }
    
    // no rejected ID >= min_cov
    for (int i = 0; i < failed->size; i++) {
        vcf_record_t *record = failed->items[i];
        char *aux_info = strndup(record->info, record->info_len);
        char *dp = get_field_value_in_info("DP", aux_info);
        fail_unless(dp != NULL && atoi(dp) < 3000, "C3000: A rejected record can't have greater or equal coverage than specified");
    }
}
END_TEST

START_TEST (coverage_all_excluded) {
    ((coverage_filter_args*) coverage_f->args)->min_coverage = 20000;
    passed = coverage_f->filter_func(datasuite, failed,NULL, coverage_f->name, coverage_f->args);
    
    fail_unless(passed->size == 0, "C20000: The number of coverages found is not correct");
    // size(accepted + rejected) = size(whole input)
    fail_unless(passed->size + failed->size == datasuite->size,
            "C20000: The sum of the number of accepted and rejected records must be the same as the input records");
    
    // no rejected ID >= min_cov
    for (int i = 0; i < failed->size; i++) {
        vcf_record_t *record = failed->items[i];
        char *aux_info = strndup(record->info, record->info_len);
        char *dp = get_field_value_in_info("DP", aux_info);
        fail_unless(dp != NULL && atoi(dp) < 20000, "C20000: A rejected record can't have greater or equal coverage than specified");
    }
}
END_TEST

START_TEST (coverage_all_included) {
    ((coverage_filter_args*) coverage_f->args)->min_coverage = 60;
    passed = coverage_f->filter_func(datasuite, failed,NULL, coverage_f->name, coverage_f->args);
    
    fail_unless(passed->size == 389, "C60: The number of coverages found is not correct");
    // size(accepted + rejected) = size(whole input)
    fail_unless(passed->size + failed->size == datasuite->size,
            "C60: The sum of the number of accepted and rejected records must be the same as the input records");
    
    // no accepted ID < min_cov
    for (int i = 0; i < passed->size; i++) {
        vcf_record_t *record = passed->items[i];
        char *aux_info = strndup(record->info, record->info_len);
        char *dp = get_field_value_in_info("DP", aux_info);
        fail_if(dp != NULL && atoi(dp) < 60, "C60: An accepted record can't have less coverage than specified");
    }
}
END_TEST

START_TEST (missing_values) {
    ((missing_values_filter_args*) missing_values_f->args)->max_missing = 0.15;
	
    list_t *output_list = (list_t*)malloc(sizeof(list_t));
    list_init("list",1,10000,output_list);
    file_stats_t * file_s = file_stats_new();
    get_variants_stats(((vcf_record_t**)missing_values_datasuite->items), missing_values_datasuite->size, NULL, NULL,0, output_list, file_s);
    fail_if(output_list->length == 0, "There must be one element processed");
    variant_stats_t **input_stats_array = (variant_stats_t**) list_to_array(output_list);
    
    
    passed = missing_values_f->filter_func(missing_values_datasuite, failed,input_stats_array, missing_values_f->name, missing_values_f->args);

    fail_unless(passed->size == 16, "Missing0.2: The number of missing values found is not correct");
    // size(accepted + rejected) = size(whole input)
    fail_unless(passed->size + failed->size == missing_values_datasuite->size,
            "Missing0.2: The sum of the number of accepted and rejected records must be the same as the input records");
    
    // no rejected ID > max_miss
    for (int i = 0; i < failed->size; i++) {
        vcf_record_t *record = failed->items[i];
        fail_if(record->position != 1105366 && record->position != 1105411 && record->position != 3537996,
                "Failed records must be in chromosome 1, positions (1105366, 1105411, 11633148)");
    }
    list_free_deep(output_list, ((void*(*)(void*)) variant_stats_free));
    file_stats_free(file_s);
    free(input_stats_array);
}
END_TEST

START_TEST (num_alelles_basic) {
    // TODO check for different number of alleles in a file which includes a variety of them
    ((num_alleles_filter_args*) num_alleles_f->args)->num_alleles = 2;
    
    list_t *output_list = (list_t*)malloc(sizeof(list_t));
    list_init("list",1,10000,output_list);
    file_stats_t * file_s = file_stats_new();
    get_variants_stats(((vcf_record_t**)num_alleles_datasuite->items), num_alleles_datasuite->size, NULL, NULL,0, output_list, file_s);
    fail_if(output_list->length == 0, "There must be one element processed");
    variant_stats_t **input_stats_array = (variant_stats_t**) list_to_array(output_list);
    
    passed = num_alleles_f->filter_func(num_alleles_datasuite, failed,input_stats_array, num_alleles_f->name, num_alleles_f->args);
    
    fail_unless(passed->size == 37, "#alleles basic: The number of occurrences found is not correct");
    // size(accepted + rejected) = size(whole input)
    fail_unless(passed->size + failed->size == num_alleles_datasuite->size,
            "#alleles basic: The sum of the number of accepted and rejected records must be the same as the input records");
    
    int num_alternates;
    // no accepted ID != num alleles
    for (int i = 0; i < passed->size; i++) {
        vcf_record_t *record = passed->items[i];
        char **alternates = split(strndup(record->alternate, record->alternate_len), ",", &num_alternates);
        fail_if(num_alternates != 1, "#alleles basic: An accepted record must have the specified number of alleles");
        free(alternates);
    }
    
    int num_multiallelic = 0;
    char **alternates = NULL;
    // no rejected ID == num_alleles
    for (int i = 0; i < failed->size; i++) {
        vcf_record_t *record = failed->items[i];
        fail_if(!strncmp(".", record->alternate, record->alternate_len), "No-allele counts as an alternate");
        
        alternates = split(strndup(record->alternate, record->alternate_len), ",", &num_alternates);
        if (num_alternates > 1) { 
            num_multiallelic++; 
        }
        free(alternates);
    }
    
    fail_if(num_multiallelic != 3, "#alleles basic: The number of multiallelic records found is not correct");
    list_free_deep(output_list, ((void*(*)(void*)) variant_stats_free));
    file_stats_free(file_s);
    free(input_stats_array);
}
END_TEST

START_TEST (num_alelles_all_included) {
    // TODO check for biallelic variants in a file which contains biallelics only
    ((num_alleles_filter_args*) num_alleles_f->args)->num_alleles = 2;
        
    list_t *output_list = (list_t*)malloc(sizeof(list_t));
    list_init("list",1,10000,output_list);
    file_stats_t * file_s = file_stats_new();
    get_variants_stats(((vcf_record_t**)datasuite->items), datasuite->size, NULL, NULL,0, output_list, file_s);
    fail_if(output_list->length == 0, "There must be one element processed");
    variant_stats_t **input_stats_array = (variant_stats_t**) list_to_array(output_list);
    
    passed = num_alleles_f->filter_func(datasuite, failed, input_stats_array, num_alleles_f->name, num_alleles_f->args);
    
    fail_unless(passed->size == 389, "All biallelic: The number of occurrences found is not correct");
    // size(accepted + rejected) = size(whole input)
    fail_unless(passed->size + failed->size == datasuite->size,
            "All biallelic: The sum of the number of accepted and rejected records must be the same as the input records");
    
    int num_alternates;
    // no accepted ID != num_alleles
    for (int i = 0; i < passed->size; i++) {
        vcf_record_t *record = passed->items[i];
        char **alternates = split(record->alternate, ",", &num_alternates);
        fail_if(num_alternates != 1, "All biallelic: An accepted record can't have a distinct number of alleles than specified");
        free(alternates);
    }
    list_free_deep(output_list, ((void*(*)(void*)) variant_stats_free));
    file_stats_free(file_s);
    free(input_stats_array);
}
END_TEST

START_TEST (num_alelles_all_excluded) {
    // TODO check for multiallelic variants in a file which contains biallelics only
    ((num_alleles_filter_args*) num_alleles_f->args)->num_alleles = 3;
        
    list_t *output_list = (list_t*)malloc(sizeof(list_t));
    list_init("list",1,10000,output_list);
    file_stats_t * file_s = file_stats_new();
    get_variants_stats(((vcf_record_t**)datasuite->items), datasuite->size, NULL, NULL,0, output_list, file_s);
    fail_if(output_list->length == 0, "There must be one element processed");
    variant_stats_t **input_stats_array = (variant_stats_t**) list_to_array(output_list);
    
    passed = num_alleles_f->filter_func(datasuite, failed, input_stats_array, num_alleles_f->name, num_alleles_f->args);
    
    fail_unless(passed->size == 0, "None multiallelic: The number of occurrences found is not correct");
    // size(accepted + rejected) = size(whole input)
    fail_unless(passed->size + failed->size == datasuite->size,
            "None multiallelic: The sum of the number of accepted and rejected records must be the same as the input records");
    
    int num_alternates;
    // no accepted ID < min_qual
    for (int i = 0; i < passed->size; i++) {
        vcf_record_t *record = passed->items[i];
        char **alternates = split(record->alternate, ",", &num_alternates);
        fail_if(num_alternates == 3, "None multiallelic: An accepted record can't have the same number of alleles as specified");
        free(alternates);
    }
    list_free_deep(output_list, ((void*(*)(void*)) variant_stats_free));
    file_stats_free(file_s);
    free(input_stats_array);
}
END_TEST


/*
START_TEST (snp_include) {
	((snp_filter_args*) snp_f->args)->include_snps = 1;
        
    list_t *output_list = (list_t*)malloc(sizeof(list_t));
    list_init("list",1,10000,output_list);
    file_stats_t * file_s = file_stats_new();
    get_variants_stats(((vcf_record_t**)datasuite->items), datasuite->size, NULL, NULL,0, output_list, file_s);
    fail_if(output_list->length == 0, "There must be one element processed");
    variant_stats_t **input_stats_array = (variant_stats_t**) list_to_array(output_list);
    
	passed = snp_f->filter_func(datasuite, failed, input_stats_array, snp_f->name, snp_f->args);
	
	// size(accepted) = SNPS_IN_FILE
	fail_unless(passed->size == SNPS_IN_FILE, "The number of SNP recognized is not correct");
	// size(accepted + rejected) = size(whole input)
	fail_unless(passed->size + failed->size == datasuite->size,
		    "The sum of the number of accepted and rejected records must be the same as the input records");
	
	// no accepted ID = '.'
    for (int i = 0; i < passed->size; i++) {
        vcf_record_t *record = passed->items[i];
		fail_if(!strncmp(".", record->id, record->id_len), "A known SNP must have an ID");
	}
	
	// no rejected ID != '.'
    for (int i = 0; i < failed->size; i++) {
        vcf_record_t *record = failed->items[i];
		fail_if(strncmp(".", record->id, record->id_len), "An unknown SNP can't have an ID defined");
	}
    list_free_deep(output_list, ((void*(*)(void*)) variant_stats_free));
    file_stats_free(file_s);
    free(input_stats_array);
}
END_TEST
*/
/*
START_TEST (snp_exclude) {
	((snp_filter_args*) snp_f->args)->include_snps = 0;
        
    list_t *output_list = (list_t*)malloc(sizeof(list_t));
    list_init("list",1,10000,output_list);
    file_stats_t * file_s = file_stats_new();
    get_variants_stats(((vcf_record_t**)datasuite->items), datasuite->size, NULL, NULL,0, output_list, file_s);
    fail_if(output_list->length == 0, "There must be one element processed");
    variant_stats_t **input_stats_array = (variant_stats_t**) list_to_array(output_list);
    
	passed = snp_f->filter_func(datasuite, failed,input_stats_array, snp_f->name, snp_f->args);
	
	// size(failed) = SNPS_IN_FILE
	fail_unless(failed->size == SNPS_IN_FILE, "The number of SNP recognized is not correct");
	// size(accepted + rejected) = size(whole input)
	fail_unless(passed->size + failed->size == datasuite->size,
		    "The sum of the number of accepted and rejected records must be the same as the input records");
	
	// no accepted ID != '.'
    for (int i = 0; i < passed->size; i++) {
        vcf_record_t *record = passed->items[i];
		fail_if(strncmp(".", record->id, record->id_len), "An unknown SNP can't have an ID defined");
	}
	
	// no rejected ID = '.'
    for (int i = 0; i < failed->size; i++) {
        vcf_record_t *record = failed->items[i];
		fail_if(!strncmp(".", record->id, record->id_len), "A known SNP must have an ID");
	}
    list_free_deep(output_list, ((void*(*)(void*)) variant_stats_free));
    file_stats_free(file_s);
    free(input_stats_array);
}
END_TEST
*/
/*
START_TEST (region_chrom_1) {
	// create filter for just one chromosome
	char input[] = "1";
	region_f = region_filter_new(input, 0, "1", url, species, version);
        
    list_t *output_list = (list_t*)malloc(sizeof(list_t));
    list_init("list",1,10000,output_list);
    file_stats_t * file_s = file_stats_new();
    get_variants_stats(((vcf_record_t**)datasuite->items), datasuite->size, NULL, NULL,0, output_list, file_s);
    fail_if(output_list->length == 0, "There must be one element processed");
    variant_stats_t **input_stats_array = (variant_stats_t**) list_to_array(output_list);
    
	passed = region_f->filter_func(datasuite, failed,input_stats_array, region_f->name, region_f->args);
	
	// size(accepted + rejected) = size(whole input)
	fail_unless(passed->size + failed->size == datasuite->size,
		    "The sum of the number of accepted and rejected records must be the same as the input records");
	
	// no accepted chromosome != 1
    for (int i = 0; i < passed->size; i++) {
        vcf_record_t *record = passed->items[i];
		fail_if(strncmp("1", record->chromosome, record->chromosome_len), 
			"The record must be in chromosome 1");
	}
	
	// no rejected chromosome == 1
    for (int i = 0; i < failed->size; i++) {
        vcf_record_t *record = failed->items[i];
		fail_if(!strncmp("1", record->chromosome, record->chromosome_len), 
			"The record must not be in chromosome 1");
	}
	
	region_f->free_func(region_f);
}
END_TEST


START_TEST (region_chrom_1_2) {	
	// create filter for both chromosomes in test file
	char input[] = "1,2";
	region_f = region_filter_new(input, 0, "1", url, species, version);

    list_t *output_list = (list_t*)malloc(sizeof(list_t));
    list_init("list",1,10000,output_list);
    file_stats_t * file_s = file_stats_new();
    get_variants_stats(((vcf_record_t**)datasuite->items), datasuite->size, NULL, NULL,0, output_list, file_s);
    fail_if(output_list->length == 0, "There must be one element processed");
    variant_stats_t **input_stats_array = (variant_stats_t**) list_to_array(output_list);
    printf("Hoo\n");
	passed = region_f->filter_func(datasuite, failed, input_stats_array, region_f->name, region_f->args);
	printf("Hoo\n");
	// all records pass
	fail_if(passed->size < datasuite->size, "All records must pass the test");
	fail_if(failed->size > 0, "There must not be rejected records");
	
	// size(accepted + rejected) = size(whole input)
	fail_unless(passed->size + failed->size == datasuite->size,
		    "The sum of the number of accepted and rejected records must be the same as the input records");
	
	// all records must be in chrom 1/2
    for (int i = 0; i < passed->size; i++) {
        vcf_record_t *record = passed->items[i];
		fail_unless(!strncmp("1", record->chromosome, record->chromosome_len) || !strncmp("2", record->chromosome, record->chromosome_len), 
			"The record must be in chromosome 1 or 2");
	}
	
	region_f->free_func(region_f);
}
END_TEST
*/


/* ******************************
 * 	Main entry point	*
 * ******************************/

int main (int argc, char *argv) {
    vcf_file_t *file = vcf_open("CEU.exon.2010_03.genotypes__head400.vcf", 10, VCF_FILE_VCF);
    vcf_file_t *quality_file = vcf_open("qualities.vcf", 10, VCF_FILE_VCF);
    vcf_file_t *num_alleles_file = vcf_open("num_alleles_test.vcf", 10, VCF_FILE_VCF);
    vcf_file_t *missing_values_file = vcf_open("missing_values.vcf", 10, VCF_FILE_VCF);
    
    read_test_datasuite(file);
    datasuite = ((vcf_batch_t*) file->record_batches->first_p->data_p)->records;
    
    read_test_datasuite(quality_file);
    quality_datasuite = ((vcf_batch_t*) quality_file->record_batches->first_p->data_p)->records;
    
    read_test_datasuite(num_alleles_file);
    num_alleles_datasuite = ((vcf_batch_t*) num_alleles_file->record_batches->first_p->data_p)->records;
	
    read_test_datasuite(missing_values_file);
    missing_values_datasuite = ((vcf_batch_t*) missing_values_file->record_batches->first_p->data_p)->records;
    
    Suite *fs = create_test_suite();
    SRunner *fs_runner = srunner_create(fs);
    srunner_run_all(fs_runner, CK_NORMAL);
    int number_failed = srunner_ntests_failed (fs_runner);
    srunner_free (fs_runner);

    vcf_close(file);
    vcf_close(quality_file);
    vcf_close(num_alleles_file);
	
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}


Suite *create_test_suite() {
    // Coverage filter
    TCase *tc_coverage = tcase_create("Coverage filters");
    tcase_add_unchecked_fixture(tc_coverage, setup_coverage, teardown_coverage);
    tcase_add_checked_fixture(tc_coverage, create_passed_failed, free_passed_failed);
    tcase_add_test(tc_coverage, coverage_basic);
    tcase_add_test(tc_coverage, coverage_all_included);
    tcase_add_test(tc_coverage, coverage_all_excluded);
    
    // Missing values filter
    TCase *tc_missing_values = tcase_create("Missing values filters");
    tcase_add_unchecked_fixture(tc_missing_values, setup_missing_values, teardown_missing_values);
    tcase_add_checked_fixture(tc_missing_values, create_passed_failed, free_passed_failed);
    tcase_add_test(tc_missing_values, missing_values);
    
    // Allele number filter
    TCase *tc_num_alleles = tcase_create("Allele count filters");
    tcase_add_unchecked_fixture(tc_num_alleles, setup_num_alleles, teardown_num_alleles);
    tcase_add_checked_fixture(tc_num_alleles, create_passed_failed, free_passed_failed);
    tcase_add_test(tc_num_alleles, num_alelles_basic);
    tcase_add_test(tc_num_alleles, num_alelles_all_included);
    tcase_add_test(tc_num_alleles, num_alelles_all_excluded);
    tcase_set_timeout(tc_num_alleles, 0);
    
    // Quality filter
    TCase *tc_quality = tcase_create("Quality filters");
    tcase_add_unchecked_fixture(tc_quality, setup_quality, teardown_quality);
    tcase_add_checked_fixture(tc_quality, create_passed_failed, free_passed_failed);
    //tcase_add_test(tc_quality, quality_limit_bound);
    //tcase_add_test(tc_quality, quality_limit_over_bound);
    //tcase_add_test(tc_quality, quality_all_included);
    //tcase_add_test(tc_quality, quality_all_excluded);
    
    // Region filter (chromosome, chrom+start position, chrom+start+end position)
    TCase *tc_region = tcase_create("Region filters");
    tcase_add_unchecked_fixture(tc_region, setup_region, teardown_region);
    tcase_add_checked_fixture(tc_region, create_passed_failed, free_passed_failed);
    //tcase_add_test(tc_region, region_chrom_1);
    //tcase_add_test(tc_region, region_chrom_1_2);
    //tcase_add_test(tc_region, region_chrom_start);
    //tcase_add_test(tc_region, region_chrom_start_end);
	
    // SNP filter (include and exclude)
    TCase *tc_snp = tcase_create("SNP filters");
    tcase_add_unchecked_fixture(tc_snp, setup_snp, teardown_snp);
    tcase_add_checked_fixture(tc_snp, create_passed_failed, free_passed_failed);
    //tcase_add_test(tc_snp, snp_include);
    //tcase_add_test(tc_snp, snp_exclude);
    
    // Chains of filter (SNP+region...)
    TCase *tc_filterchain = tcase_create("Filter chains");
    //tcase_add_unchecked_fixture(tc_filterchain, setup_snp_region, teardown_snp_region);
    //tcase_add_checked_fixture(tc_filterchain, create_passed_failed, free_passed_failed);
    //tcase_add_test(tc_filterchain, snpinclude_regionchromstartend_chain);

    // Add test cases to a test suite
    Suite *fs = suite_create("VCF filters");
    suite_add_tcase(fs, tc_coverage);
    suite_add_tcase(fs, tc_missing_values);
    suite_add_tcase(fs, tc_num_alleles);
    suite_add_tcase(fs, tc_quality);
    suite_add_tcase(fs, tc_region);
    suite_add_tcase(fs, tc_snp);
    suite_add_tcase(fs, tc_filterchain);

    return fs;
}

void read_test_datasuite(vcf_file_t *file) {
    if (vcf_parse_batches(400, file)) {
            fprintf(stderr, "Error reading file\n");
    }

    printf("Read %zu/%zu batches\n", file->record_batches->length, file->record_batches->max_length);
    printf("Batch contains %zu records\n", ((vcf_batch_t*) file->record_batches->first_p->data_p)->records->size);
    
}
