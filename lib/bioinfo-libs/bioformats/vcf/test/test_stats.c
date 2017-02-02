#include <stdio.h>
#include <stdlib.h>

#include <check.h>

#include <bioformats/vcf/vcf_file_structure.h>
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_stats.h>


static vcf_record_t *record;
static list_t *output_list;
static list_item_t *record_item;

static file_stats_t *file_stats;


Suite *create_test_suite();
void read_test_datasuite(vcf_file_t *file);


/* ******************************
 *      Unchecked fixtures      *
 * ******************************/

void setup_stats(void) {
    record = vcf_record_new();
    set_vcf_record_chromosome("1", 1, record);
    set_vcf_record_position(1234567, record);
    set_vcf_record_id("rs12", 4, record);
    set_vcf_record_reference("G", 1, record);
    set_vcf_record_filter("PASS", 4, record);
    set_vcf_record_info("NS=3", 4, record);
    
    output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", 1, 8, output_list);
    
    file_stats = file_stats_new();
}

void teardown_stats(void) {
    free(output_list);
    file_stats_free(file_stats);
}

/* ******************************
 *          Unit tests         *
 * ******************************/

START_TEST (biallelic) {
    set_vcf_record_alternate("T", 1, record);
    set_vcf_record_format("GC:GT", 5, record);
    
    size_t sample_idx = 0;
    char *sample = strdup("1:0/0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("2:1/0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("1:0/1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("3:0/0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("1:1/1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("1:./1"); add_vcf_record_sample(sample, strlen(sample), record);
    
    get_variants_stats(&record, 1, output_list, file_stats);
    fail_if(output_list->length == 0, "There must be one element processed");
    
    variant_stats_t *result = (variant_stats_t*) output_list->first_p->data_p;
    
    fail_unless(strcmp(result->ref_allele, "G") == 0, "The reference allele should be G");
    fail_unless(result->num_alleles == 2, "There should be 2 alleles");
    fail_unless(strcmp(result->alternates[0], "T") == 0, "The alternate allele should be T");
    fail_unless(result->alleles_count[0] == 6, "There should be 6 reference alleles read");
    fail_unless(result->alleles_count[1] == 5, "There should be 5 alternate alleles read");
    
    fail_unless(result->genotypes_count[0] == 2, "There should be 2 0/0 alleles read");
    fail_unless(result->genotypes_count[1] == 1, "There should be 1 0/1 alleles read");
    fail_unless(result->genotypes_count[2] == 1, "There should be 1 1/0 alleles read");
    fail_unless(result->genotypes_count[3] == 1, "There should be 1 1/1 alleles read");
    
    fail_unless(result->missing_alleles == 1, "There should be 1 missing allele");
    fail_unless(result->missing_genotypes == 1, "There should be 1 missing genotype");
}
END_TEST


START_TEST (multiallelic) {
    set_vcf_record_alternate("T,GT", 4, record);
    set_vcf_record_format("GT:GC", 5, record);
    
    size_t sample_idx = 0;
    char *sample = strdup("0/0:1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("1/0:2"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/1:1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("1/2:3"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("1/1:1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("1/.:1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/2:1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("2/1:1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/1:1"); add_vcf_record_sample(sample, strlen(sample), record);
    
    get_variants_stats(&record, 1, output_list, file_stats);
    fail_if(output_list->length == 0, "There must be one element processed");
    
    variant_stats_t *result = (variant_stats_t*) output_list->first_p->data_p;
    
    fail_unless(strcmp(result->ref_allele, "G") == 0, "The reference allele should be G");
    fail_unless(result->num_alleles == 3, "There should be 3 alleles");
    fail_unless(strcmp(result->alternates[0], "T") == 0, "The #1 alternate allele should be T");
    fail_unless(strcmp(result->alternates[1], "GT") == 0, "The #2 alternate allele should be GT");
    
    fail_unless(result->alleles_count[0] == 6, "There should be 6 reference alleles read");
    fail_unless(result->alleles_count[1] == 8, "There should be 8 alternate (#1) alleles read");
    fail_unless(result->alleles_count[2] == 3, "There should be 3 alternate (#2) alleles read");
    
    fail_unless(result->genotypes_count[0] == 1, "There should be 1 0/0 alleles read");
    fail_unless(result->genotypes_count[1] == 2, "There should be 1 0/1 alleles read");
    fail_unless(result->genotypes_count[2] == 1, "There should be 1 0/2 alleles read");
    fail_unless(result->genotypes_count[3] == 1, "There should be 1 1/0 alleles read");
    fail_unless(result->genotypes_count[4] == 1, "There should be 1 1/1 alleles read");
    fail_unless(result->genotypes_count[5] == 1, "There should be 1 1/2 alleles read");
    fail_unless(result->genotypes_count[6] == 0, "There should be 1 2/0 alleles read");
    fail_unless(result->genotypes_count[7] == 1, "There should be 1 2/1 alleles read");
    fail_unless(result->genotypes_count[8] == 0, "There should be 0 2/2 alleles read");
    
    fail_unless(result->missing_alleles == 3, "There should be 1 missing allele");
    fail_unless(result->missing_genotypes == 2, "There should be 1 missing genotype");
}
END_TEST


START_TEST (homozygous) {
    set_vcf_record_alternate(".", 1, record);
    set_vcf_record_format("GT:GQ:DP:HQ", 11, record);
    
    size_t sample_idx = 0;
    char *sample = strdup("0/0:1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/0:2"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/.:1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./0:3"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/0:1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:1"); add_vcf_record_sample(sample, strlen(sample), record);
    
    get_variants_stats(&record, 1, output_list, file_stats);
    fail_if(output_list->length == 0, "There must be one element processed");
    
    variant_stats_t *result = (variant_stats_t*) output_list->first_p->data_p;
    
    fail_unless(strcmp(result->ref_allele, "G") == 0, "The reference allele should be G");
    fail_unless(result->num_alleles == 2, "There should be 2 alleles");
    fail_unless(strcmp(result->alternates[0], ".") == 0, "The alternate allele should be .");
    fail_unless(result->alleles_count[0] == 8, "There should be 8 reference alleles read");
    
    fail_unless(result->genotypes_count[0] == 3, "There should be 3 0/0 alleles read");
    
    fail_unless(result->missing_alleles == 4, "There should be 4 missing allele");
    fail_unless(result->missing_genotypes == 3, "There should be 3 missing genotype");
}
END_TEST


START_TEST (from_CEU_exon) {
    set_vcf_record_reference("A", 1, record);
    set_vcf_record_alternate("G", 1, record);
    set_vcf_record_format("GT:DP", 5, record);
    
    size_t sample_idx = 0;
    char *sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("1/1:31"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/1:2"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("1/1:12"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/1:20"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/1:13"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:3"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/1:11"); add_vcf_record_sample(sample, strlen(sample), record);
    
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/1:39"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/1:26"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("1/1:29"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:3"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:1"); add_vcf_record_sample(sample, strlen(sample), record);
    
    sample = strdup("./.:2"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("1/1:36"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/1:19"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/1:27"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:2"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/1:5"); add_vcf_record_sample(sample, strlen(sample), record);
    
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("1/1:29"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:3"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("1/1:19"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("1/1:5"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/1:3"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:3"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/1:4"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:6"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/1:4"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/1:3"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    
    sample = strdup("./.:1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/1:3"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:0"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:1"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/0:44"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("./.:3"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("1/1:11"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/1:21"); add_vcf_record_sample(sample, strlen(sample), record);
    sample = strdup("0/0:53"); add_vcf_record_sample(sample, strlen(sample), record);
    
    get_variants_stats(&record, 1, output_list, file_stats);
    fail_if(output_list->length == 0, "There must be one element processed");
    
    variant_stats_t *result = (variant_stats_t*) output_list->first_p->data_p;
    
    fail_unless(strcmp(result->ref_allele, "A") == 0, "The reference allele should be A");
    fail_unless(result->num_alleles == 2, "There should be 2 alleles");
    fail_unless(strcmp(result->alternates[0], "G") == 0, "The alternate allele should be G");
    fail_unless(result->alleles_count[0] == 19, "There should be 19 reference alleles read");
    fail_unless(result->alleles_count[1] == 31, "There should be 31 alternate alleles read");
    
    fail_unless(result->genotypes_count[0] == 2, "There should be 2 0/0 alleles read");
    fail_unless(result->genotypes_count[1] == 15, "There should be 15 0/1 alleles read");
    fail_unless(result->genotypes_count[2] == 0, "There should be 0 1/0 alleles read");
    fail_unless(result->genotypes_count[3] == 8, "There should be 8 1/1 alleles read");
    
    fail_unless(result->missing_alleles == 130, "There should be 130 missing allele");
    fail_unless(result->missing_genotypes == 65, "There should be 65 missing genotype");
}
END_TEST


START_TEST(mendelian_errors) {
    vcf_file_t *vcf_file = vcf_open("sample_stats.vcf", 10);
    read_test_datasuite(vcf_file);
    array_list_t *datasuite = ((vcf_batch_t*) vcf_file->record_batches->first_p->data_p)->records;
    
    ped_file_t *ped_file = ped_open("sample_stats.ped");
    ped_read(ped_file);
    
    khash_t(ids) *sample_ids = associate_samples_and_positions(vcf_file);
    individual_t **individuals = sort_individuals(vcf_file, ped_file);
    
    for (int i = 0; i < get_num_vcf_samples(vcf_file); i++) {
        fail_unless(individuals[i], "All individuals must exist");
        if (i % 3 == 0) {
            fail_unless(individuals[i]->father, "Children must have a father");
            fail_unless(individuals[i]->mother, "Children must have a mother");
        }
    }
    
    sample_stats_t **sample_stats = malloc (get_num_vcf_samples(vcf_file) * sizeof(sample_stats_t*));
    for (int j = 0; j < get_num_vcf_samples(vcf_file); j++) {
        sample_stats[j] = sample_stats_new(array_list_get(j, vcf_file->samples_names));
    }
    
    get_sample_stats((vcf_record_t**) datasuite->items, datasuite->size, individuals, sample_ids, sample_stats, file_stats);
    
    sample_stats_t *cur_stats = sample_stats[0];  // 3376
    fail_if(cur_stats->missing_genotypes != 2, "3376: 2 missing GTs expected");
    fail_if(cur_stats->mendelian_errors != 4, "3376: 4 mendelian errors expected");
    
    cur_stats = sample_stats[1];  // 3378
    fail_if(cur_stats->missing_genotypes != 1, "3378: 1 missing GTs expected");
    fail_if(cur_stats->mendelian_errors != 0, "3378: 0 mendelian errors expected");
    
    cur_stats = sample_stats[2];  // 3377
    fail_if(cur_stats->missing_genotypes != 1, "3377: 1 missing GTs expected");
    fail_if(cur_stats->mendelian_errors != 0, "3377: 0 mendelian errors expected");
    
    cur_stats = sample_stats[3];  // 3385
    fail_if(cur_stats->missing_genotypes != 1, "3385: 1 missing GTs expected");
    fail_if(cur_stats->mendelian_errors != 3, "3385: 3 mendelian errors expected");
    
    cur_stats = sample_stats[4];  // 3386
    fail_if(cur_stats->missing_genotypes != 2, "3386: 2 missing GTs expected");
    fail_if(cur_stats->mendelian_errors != 0, "3386: 0 mendelian errors expected");
    
    cur_stats = sample_stats[5];  // 3387
    fail_if(cur_stats->missing_genotypes != 1, "3387: 1 missing GTs expected");
    fail_if(cur_stats->mendelian_errors != 0, "3387: 0 mendelian errors expected");
}
END_TEST


/* ******************************
 *      Main entry point        *
 * ******************************/

int main (int argc, char *argv) {
    Suite *fs = create_test_suite();
    SRunner *fs_runner = srunner_create(fs);
    srunner_run_all(fs_runner, CK_NORMAL);
    int number_failed = srunner_ntests_failed (fs_runner);
    srunner_free (fs_runner);
    
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}


Suite *create_test_suite(void) {
    TCase *tc_variant_stats = tcase_create("Variant stats");
    tcase_add_unchecked_fixture(tc_variant_stats, setup_stats, teardown_stats);
    tcase_add_test(tc_variant_stats, biallelic);
    tcase_add_test(tc_variant_stats, multiallelic);
    tcase_add_test(tc_variant_stats, homozygous);
    tcase_add_test(tc_variant_stats, from_CEU_exon);
    
    TCase *tc_sample_stats = tcase_create("Sample stats");
    tcase_add_test(tc_sample_stats, mendelian_errors);
            
    // Add test cases to a test suite
    Suite *fs = suite_create("VCF stats");
    suite_add_tcase(fs, tc_variant_stats);
    suite_add_tcase(fs, tc_sample_stats);
    
    return fs;
}

void read_test_datasuite(vcf_file_t *file) {
    if (vcf_parse_batches(400, file)) {
        fprintf(stderr, "Error reading file\n");
    }

    printf("Read %zu/%zu batches\n", file->record_batches->length, file->record_batches->max_length);
    printf("Batch contains %zu records\n", ((vcf_batch_t*) file->record_batches->first_p->data_p)->records->size);
    
}
