#include <stdio.h>

#include <check.h>

#include "containers/linked_list.h"


linked_list_t* list;

Suite *create_test_suite(void);


/* **************************
 *      Checked fixtures    *
 * **************************/

void create_linked_list() {
    list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    for (int i = 5; i >= 0; i--) {
        linked_list_insert((void *)i, list);
    }
}

void free_linked_list() {
    linked_list_free(list, NULL);
}



START_TEST(test_empty_list) {
    linked_list_t *empty_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    
    fail_if(linked_list_get_first(empty_list), "First: An empty list has no elements");
    fail_if(linked_list_get_last(empty_list), "Last: An empty list has no elements");
    
    fail_if(linked_list_remove_first(empty_list), "Remove first: An empty list has no elements to remove");
    fail_if(linked_list_remove_last(empty_list), "Remove last: An empty list has no elements to remove");
    fail_if(linked_list_remove_at(2, empty_list), "Remove at: An empty list has no elements to remove");
}
END_TEST


START_TEST(test_insert_and_remove) {
    // [0->1->2->3->4->5->]
    fail_if((int) linked_list_get(3, list) != 3, "Get: Position #3 must contain number 3");
    fail_if((int) linked_list_get_first(list) != 0, "Get: Position #0 (first) must contain number 0");
    fail_if((int) linked_list_get_last(list) != 5, "Get: Position #5 (last) must contain number 5");
    
    fail_if((int) linked_list_remove_first(list) != 0, "Remove: Position #0 (first) must contain number 0");
    fail_if((int) linked_list_remove_last(list) != 5, "Remove: Position #5 (last) must contain number 5");
    fail_if((int) linked_list_remove_at(2, list) != 3, "Remove: Position #2 must contain number 3");
    
    // [1->2->4->]
    fail_if(linked_list_size(list) != 3, "The list must contain 3 elements");
    
    linked_list_insert_at(2, (void *)3, list);
    // [1->2->3->4->]
    fail_if(linked_list_size(list) != 4, "The list must contain 4 elements");
    for (int i = 0; i < linked_list_size(list); i++) {
        fail_if((int) linked_list_get(i, list) != i + 1, "Each position must contain an element of the same value + 1");
    }
    
    linked_list_insert_first((void *)0, list);
    linked_list_insert_last((void *)5, list);
    // [0->1->2->3->4->5->]
    fail_if(linked_list_size(list) != 6, "The list must contain 6 elements");
    for (int i = 0; i < linked_list_size(list); i++) {
        fail_if((int) linked_list_get(i, list) != i, "Each position must contain an element of the same value");
    }
}
END_TEST


START_TEST(test_iterators) {
    // [0->1->2->3->4->5->]
    linked_list_iterator_t* iterator = linked_list_iterator_new(list);
    fail_if(linked_list_iterator_curr(iterator) != 0, "A new iterator must be at position #0");
    
    fail_if(linked_list_iterator_next(iterator) != 1, "Iterator must be at position #1");

    fail_if(linked_list_iterator_prev(iterator) != 0, "Iterator must be at position #0");

    fail_if(linked_list_iterator_prev(iterator) != NULL, "Position #0 must have no previous position");

    // Iterator is in NULL position, inserts at the end
    linked_list_iterator_insert((void*) 8, iterator);
    // [0->1->2->3->4->5->8->]
    fail_if(linked_list_size(list) != 7, "The list must contain 7 elements");
    fail_if((int) linked_list_get_last(list) != 8, "List tail must be 8");
    for (int i = 0; i < linked_list_size(list) - 1; i++) {
        fail_if((int) linked_list_get(i, list) != i, "Each position must contain an element of the same value");
    }

    linked_list_iterator_first(iterator);
    linked_list_iterator_next(iterator);
    linked_list_iterator_next(iterator);
    linked_list_iterator_insert((void *)8, iterator);
    // [0->1->8->2->3->4->5->8->]
    fail_if(linked_list_size(list) != 8, "The list must contain 8 elements");
    fail_if((int) linked_list_get(2, list) != 8, "Position #2 must contain value 8");

    linked_list_iterator_last(iterator);
    linked_list_iterator_next(iterator);
    linked_list_iterator_insert((void *)8, iterator);
    // [0->1->8->2->3->4->5->8->8->]
    fail_if(linked_list_size(list) != 9, "The list must contain 9 elements");
    fail_if((int) linked_list_get(7, list) != 8, "Position #7 must contain value 8");
    fail_if((int) linked_list_get(8, list) != 8, "Position #8 must contain value 8");

    void* item = linked_list_iterator_remove(iterator);
    // [0->1->8->2->3->4->5->8->]
    fail_if(linked_list_size(list) != 8, "The list must contain 8 elements");
    fail_if((int) item != 8, "The removed item must be 8");
    fail_if((int) linked_list_get_last(list) != 8, "Position #7 must contain value 8");

    linked_list_iterator_first(iterator);
    item = linked_list_iterator_remove(iterator);
    // [1->8->2->3->4->5->8->]
    fail_if(linked_list_size(list) != 7, "The list must contain 7 elements");
    fail_if((int) item != 0, "The removed item must be 0");
    fail_if((int) linked_list_get_first(list) != 1, "Position #0 must contain value 1");
    
    linked_list_iterator_next(iterator);
    linked_list_iterator_next(iterator);
    item = linked_list_iterator_remove(iterator);
    // [1->8->2->3->4->5->8->]
    fail_if(linked_list_size(list) != 6, "The list must contain 6 elements");
    fail_if((int) item != 2, "The removed item must be 2");
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
    TCase *tc_insert_remove = tcase_create("Insert and remove");
    tcase_add_checked_fixture(tc_insert_remove, create_linked_list, free_linked_list);
    tcase_add_test(tc_insert_remove, test_empty_list);
    tcase_add_test(tc_insert_remove, test_insert_and_remove);
    
    TCase *tc_iterators = tcase_create("Iterators");
    tcase_add_checked_fixture(tc_iterators, create_linked_list, free_linked_list);
    tcase_add_test(tc_iterators, test_iterators);
    
    // Add test cases to a test suite
    Suite *fs = suite_create("Linked list");
    suite_add_tcase(fs, tc_insert_remove);
    suite_add_tcase(fs, tc_iterators);
    
    return fs;
}

