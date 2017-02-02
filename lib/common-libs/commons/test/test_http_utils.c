#include <stdio.h>
#include <stdlib.h>
#include <check.h>

#include "../http_utils.h"

Suite *create_test_suite();
static size_t callback(char *contents, size_t size, size_t nmemb, void *userdata);


/* ******************************
 *  Unchecked fixtures  *
 * ******************************/

void setup_http_get(void)
{
    printf("Begin HTTP GET request\n");
}

void teardown_http_get(void)
{
    printf("Finished HTTP GET request\n");
}

void setup_http_post(void)
{
    printf("Begin HTTP POST request\n");
}

void teardown_http_post(void)
{
    printf("Finished HTTP POST request\n");
}


/* ******************************
 *       Unit tests         *
 * ******************************/

START_TEST (http_get_url_noparams)
{
   int ret_code = http_get("www.google.com", NULL, NULL, 0, callback);
   fail_unless(ret_code == 0, "There must be no error when requesting www.google.com");
   
   ret_code = http_get("www.madeupurl.com", NULL, NULL, 0, callback);
   fail_if(ret_code == 0, "There must be errors when requesting www.madeupurl.com");
}
END_TEST


START_TEST (http_get_url_params)
{
   char *params[2] = { "answer", "id" };
   char *params_values[2] = { "1046380353", "1044780608" };
   int ret_code = http_get("faq.cprogramming.com/cgi-bin/smartfaq.cgi", params, params_values, 0, callback);
   fail_unless(ret_code == 0, "There must be no error when requesting faq.cprogramming.com/cgi-bin/smartfaq.cgi?answer=1046380353&id=1044780608");
}
END_TEST


START_TEST (http_post_url_noparams)
{   
   int ret_code = http_post("www.google.com", NULL, NULL, 0, callback);
   fail_unless(ret_code == 0, "There must be no error when requesting www.google.com");
   
   ret_code = http_post("www.madeupurl.com", NULL, NULL, 0, callback);
   fail_if(ret_code == 0, "There must be errors when requesting www.madeupurl.com");
}
END_TEST


START_TEST (http_post_url_params)
{   
   char *params[2] = { "answer", "id" };
   char *params_values[2] = { "1046380353", "1044780608" };
   int ret_code = http_post("faq.cprogramming.com/cgi-bin/smartfaq.cgi", params, params_values, 0, callback);
   fail_unless(ret_code == 0, "There must be no error when requesting faq.cprogramming.com/cgi-bin/smartfaq.cgi?answer=1046380353&id=1044780608");
}
END_TEST


/* ******************************
 *  Main entry point    *
 * ******************************/

int main (int argc, char *argv[])
{    
    Suite *fs = create_test_suite();
    SRunner *fs_runner = srunner_create(fs);
    srunner_run_all(fs_runner, CK_NORMAL);
    int number_failed = srunner_ntests_failed (fs_runner);
    srunner_free (fs_runner);
    
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}


Suite *create_test_suite()
{
    // HTTP GET requests
    TCase *tc_get = tcase_create("HTTP GET method");
    tcase_add_unchecked_fixture(tc_get, setup_http_get, teardown_http_get);
    tcase_add_test(tc_get, http_get_url_noparams);
    tcase_add_test(tc_get, http_get_url_params);
    
    // HTTP POST requests
    TCase *tc_post = tcase_create("HTTP POST method");
    tcase_add_unchecked_fixture(tc_post, setup_http_post, teardown_http_post);
    tcase_add_test(tc_post, http_post_url_noparams);
    tcase_add_test(tc_post, http_post_url_params);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("HTTP utilities");
    suite_add_tcase(fs, tc_get);
    suite_add_tcase(fs, tc_post);
    
    return fs;
}

static size_t callback(char *contents, size_t size, size_t nmemb, void *userdata) {
    size_t realsize = size * nmemb;
    printf("Response of %zu bytes\n", realsize);
    return realsize;
}
