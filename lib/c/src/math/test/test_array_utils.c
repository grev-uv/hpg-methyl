#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <check.h>

#include "array_utils.h"
#include "p_adjust.h"

Suite *create_test_suite();

unsigned int length = 20;
double *data1, *data2;


/* ******************************
 *        Checked fixtures      *
 * ******************************/

void setup_data1(void) {
    data1 = (double*) malloc (length * sizeof(double));
    for (int i = 0; i < length; i++) {
        data1[i] = i*0.001;
    }
}

void setup_data2(void) {
    data2 = (double*) malloc (length * sizeof(double));
    for (int i = 0; i < length; i++) {
        data2[i] = 2 + i;
    }
}

void teardown_data1(void) {
    free(data1);
}

void teardown_data2(void) {
    free(data2);
}


/* ******************************
 *          Unit tests          *
 * ******************************/

/* Array-scalar operations */
START_TEST(test_array_scalar_multiply) {
    double scalar = 2.0;
    fail_unless(array_scalar_multiply(NULL, length, scalar) == -1, "The 1st arg in array_scalar_multiply is not valid");
    fail_unless(array_scalar_multiply(data1, 0, scalar) == -2, "The 2nd arg in array_scalar_multiply is not valid");
    
    int ret = array_scalar_multiply(data1, length, scalar);
    fail_unless(ret == 0, "An error in array_scalar_multiply has occurred");
    for (int i = 0; i < length; i++) {
        fail_if(data1[i] - (i * 2.0) > 0, "The input array elements must have been multiplied by 2.0");
    }
}
END_TEST

START_TEST(test_array_scalar_sum) {
    double scalar = 5.0;
    fail_unless(array_scalar_sum(NULL, length, scalar) == -1, "The 1st arg in array_scalar_sum is not valid");
    fail_unless(array_scalar_sum(data1, 0, scalar) == -2, "The 2nd arg in array_scalar_sum is not valid");
    
    int ret = array_scalar_sum(data1, length, scalar);
    
    fail_unless(ret == 0, "An error in array_scalar_sum has occurred");
    for (int i = 0; i < length; i++) {
        fail_if(data1[i] - (i + 5.0) > 0, "The input array elements must have been added with 5.0");
    }
}
END_TEST


/* Array-array operations */
START_TEST(test_array_sum) {
    fail_unless(array_sum(NULL, data2, length) == -1, "The 1st arg in array_sum is not valid");
    fail_unless(array_sum(data1, NULL, length) == -2, "The 2nd arg in array_sum is not valid");
    fail_unless(array_sum(data1, data2, 0) == -3, "The 3rd arg in array_sum is not valid");
    
    double *aux = (double*) malloc (length * sizeof(double));
    memcpy(aux, data1, length * sizeof(double));
    int ret = array_sum(aux, data2, length);
    
    fail_unless(ret == 0, "An error in array_sum has occurred");
    for (int i = 0; i < length; i++) {
        fail_if(aux[i] - (data1[i] + data2[i]) > 0, "The input array elements must have been added");
    }
    
    free(aux);
}
END_TEST

START_TEST(test_array_substract) {
    fail_unless(array_substract(NULL, data2, length) == -1, "The 1st arg in array_substract is not valid");
    fail_unless(array_substract(data1, NULL, length) == -2, "The 2nd arg in array_substract is not valid");
    fail_unless(array_substract(data1, data2, 0) == -3, "The 3rd arg in array_substract is not valid");
    
    double *aux = (double*) malloc (length * sizeof(double));
    memcpy(aux, data1, length * sizeof(double));
    int ret = array_substract(aux, data2, length);
    
    fail_unless(ret == 0, "An error in array_substract has occurred");
    for (int i = 0; i < length; i++) {
//         printf("aux[%d] = %e, data1[i] + data2[i] = %e\n", aux[i], (data1[i] - data2[i]));
        fail_if(aux[i] - (data1[i] - data2[i]) > 0, "The input array elements must have been substracted");
    }
    
    free(aux);
}
END_TEST

START_TEST(test_array_dotproduct) {
    double res, expect_res = 0.0;
    
    fail_unless(array_dotproduct(NULL, data2, length, &res) == -1, "The 1st arg in array_log is not valid");
    fail_unless(array_dotproduct(data1, NULL, length, &res) == -2, "The 2nd arg in array_log is not valid");
    fail_unless(array_dotproduct(data1, data2, 0, &res) == -3, "The 3rd arg in array_log is not valid");
    fail_unless(array_dotproduct(data1, data2, length, NULL) == -4, "The 4th arg in array_log is not valid");
    
    for (int i = 0; i < length; i++) {
        expect_res += data1[i] * data2[i];
    }
    int ret = array_dotproduct(data1, data2, length, &res);
    fail_unless(ret == 0, "An error in array_sum has occurred");
    fail_unless(res == expect_res, "The dot product has not been correctly calculated");
}
END_TEST


/* Logarithms of array elements */
START_TEST(test_array_log) {
    fail_unless(array_log(NULL, length) == -1, "The 1st arg in array_log is not valid");
    fail_unless(array_log(data1, 0) == -2, "The 2nd arg in array_log is not valid");
    
    double *aux = (double*) malloc (length * sizeof(double));
    memcpy(aux, data1, length * sizeof(double));
    int ret = array_log(data1, length);
    
    fail_unless(ret == 0, "An error in array_log has occurred");
    for (int i = 0; i < length; i++) {
        fail_if(data1[i] - log(aux[i]) > 0, "The input array elements must have been transformed in their ln");
    }
}
END_TEST

START_TEST(test_array_log10) {
    fail_unless(array_log10(NULL, length) == -1, "The 1st arg in array_log10 is not valid");
    fail_unless(array_log10(data1, 0) == -2, "The 2nd arg in array_log10 is not valid");
    
    double *aux = (double*) malloc (length * sizeof(double));
    memcpy(aux, data1, length * sizeof(double));
    int ret = array_log10(data1, length);
    
    fail_unless(ret == 0, "An error in array_log10 has occurred");
    for (int i = 0; i < length; i++) {
        fail_if(data1[i] - log10(aux[i]) > 0, "The input array elements must have been transformed in their ln");
    }
}
END_TEST

START_TEST(test_array_log_base) {
    int base = 2;
    fail_unless(array_log_base(NULL, length, base) == -1, "The 1st arg in array_log_base is not valid");
    fail_unless(array_log_base(data1, 0, base) == -2, "The 2nd arg in array_log_base is not valid");
    
    double *aux = (double*) malloc (length * sizeof(double));
    memcpy(aux, data1, length * sizeof(double));
    int ret = array_log_base(data1, length, base);
    
    fail_unless(ret == 0, "An error in array_log_base has occurred");
    for (int i = 0; i < length; i++) {
        fail_if(data1[i] - (log(aux[i])/log(base)) > 0, "The input array elements must have been transformed in their ln");
    }
}
END_TEST


/* Array elements accumulation */
START_TEST(test_array_accum) {
    double res, expected_res = 0.0;
    
    fail_unless(array_accum(NULL, length, &res) == -1, "The 1st arg in array_accum is not valid");
    fail_unless(array_accum(data1, 0, &res) == -2, "The 2nd arg in array_accum is not valid");
    fail_unless(array_accum(data1, length, NULL) == -3, "The 3rd arg in array_accum is not valid");
    
    for (int i = 0; i < length; i++) {
        expected_res += data1[i];
    }
    int ret = array_accum(data1, length, &res);
    
    fail_unless(ret == 0, "An error in array_accum has occurred");
    for (int i = 0; i < length; i++) {
        fail_unless(res - expected_res == 0, "The accumulated sum has not been correctly calculated");
    }
}
END_TEST

START_TEST(test_array_accum_range) {
    double res, expected_res = 0.0;
    size_t init_idx = 5;
    
    fail_unless(array_accum_range(NULL, init_idx, length, &res) == -1, "The 1st arg in array_accum_range is not valid");
    fail_unless(array_accum_range(data1, length, init_idx, &res) == -2, "The 2nd arg in array_accum_range is not valid (1)");
    fail_unless(array_accum_range(data1, init_idx, init_idx, &res) == -2, "The 2nd arg in array_accum_range is not valid (2)");
    fail_unless(array_accum_range(data1, init_idx, length, NULL) == -4, "The 3rd arg in array_accum_range is not valid");
    
    for (int i = init_idx; i < length; i++) {
        expected_res += data1[i];
    }
    int ret = array_accum_range(data1, init_idx, length, &res);
    
    fail_unless(ret == 0, "An error in array_accum_range has occurred");
    for (int i = 0; i < length; i++) {
        fail_unless(res - expected_res == 0, "The accum_rangeulated sum has not been correctly calculated");
    }
}
END_TEST

START_TEST(test_array_order) {
    printf("Ordering data1 array\n");
    size_t *indices = (size_t*)malloc(length*sizeof(size_t));
    
    printf("Nominal p-values:\n");
    array_printf(data1, length, "%f ");
    printf("\n");
    BH_correction(data1, length);
    printf("BH p-adjust::\n");
    array_printf(data1, length, "%f ");
    printf("\n\n");
    
    data1[8] = INFINITY;
    data1[4] = NAN;
    array_order(data1, length, 0, indices);
    
    double *double_ordered = (double*)malloc(length*sizeof(double));
    array_ordered(data1, length, indices, double_ordered);
    printf("Original with NaN and INFINITY:\n");
    array_printf(data1, length, "%f ");
    printf("\n");

    printf("Ordered desc:\n");
    array_printf(double_ordered, length, "%f ");
    printf("\n");

    printf("Ordered asc:\n");
    array_order(data1, length, 1, indices);
    array_ordered(data1, length, indices, double_ordered);
    array_printf(double_ordered, length, "%f ");
    printf("\n\n");
    
    free(indices);
    free(double_ordered);
   
}
END_TEST

START_TEST(test_array_io) {
    printf("Printing data1 array\n");
    size_t chars_printed = array_printf(data1, length, "%f ");
    fail_unless(chars_printed > 0, "An error printing array has occurred");
   
    printf("Saving data1 array to '/tmp/array_fprintf.txt'\n");
    FILE *file = fopen("/tmp/array_fprintf.txt", "w");
    chars_printed = array_fprintf(data1, length, "%f\n", file);
    fclose(file);
    fail_unless(chars_printed > 0, "An error saving array has occurred");
    
    printf("Reading data1 array from '/tmp/array_fprintf.txt'\n");
    file = fopen("/tmp/array_fprintf.txt", "r");
    double *data_file = (double*)malloc(length*sizeof(double));
    array_fread(file, data_file, length);
    fclose(file);
    array_printf(data_file, length, "%f ");
    printf("\n");
}
END_TEST
/* ******************************
 *        Main entry point      *
 * ******************************/


int main(int argc, char *argv[]) {
    Suite *fs = create_test_suite();
    SRunner *fs_runner = srunner_create(fs);
    srunner_run_all(fs_runner, CK_NORMAL);
    int number_failed = srunner_ntests_failed(fs_runner);
    srunner_free (fs_runner);
    
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Suite *create_test_suite() {
    TCase *tc_array_scalar = tcase_create("Array-scalar operations");
    tcase_add_checked_fixture(tc_array_scalar, setup_data1, teardown_data1);
    tcase_add_test(tc_array_scalar, test_array_scalar_multiply);
    tcase_add_test(tc_array_scalar, test_array_scalar_sum);
    
    TCase *tc_array_array = tcase_create("Array-array operations");
    tcase_add_checked_fixture(tc_array_array, setup_data1, teardown_data1);
    tcase_add_checked_fixture(tc_array_array, setup_data2, teardown_data2);
    tcase_add_test(tc_array_array, test_array_sum);
    tcase_add_test(tc_array_array, test_array_substract);
    tcase_add_test(tc_array_array, test_array_dotproduct);
    
    TCase *tc_log = tcase_create("Logarithms of array elements");
    tcase_add_checked_fixture(tc_log, setup_data1, teardown_data1);
    tcase_add_test(tc_log, test_array_log);
    tcase_add_test(tc_log, test_array_log10);
    tcase_add_test(tc_log, test_array_log_base);
    
    TCase *tc_accum = tcase_create("Array elements accumulation");
    tcase_add_checked_fixture(tc_accum, setup_data1, teardown_data1);
    tcase_add_test(tc_accum, test_array_accum);
    tcase_add_test(tc_accum, test_array_accum_range);
    
    TCase *tc_order = tcase_create("Array elements order");
    tcase_add_checked_fixture(tc_order, setup_data1, teardown_data1);
    tcase_add_test(tc_order, test_array_order);
    
    TCase *tc_io = tcase_create("Array elements IO");
    tcase_add_checked_fixture(tc_io, setup_data1, teardown_data1);
    tcase_add_test(tc_io, test_array_io);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("Array utilities");
    suite_add_tcase(fs, tc_array_scalar);
    suite_add_tcase(fs, tc_array_array);
    suite_add_tcase(fs, tc_log);
    suite_add_tcase(fs, tc_accum);
    suite_add_tcase(fs, tc_order);
    suite_add_tcase(fs, tc_io);
    
    return fs;
}
