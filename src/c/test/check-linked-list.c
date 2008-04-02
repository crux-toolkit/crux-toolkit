#include "check-parameter.h"
#include "../linked_list.h"
#include "../utils.h"

// declare things to set up
int *data1, *data2, *data3;
LINKED_LIST_T *list1, *list2;

void list_setup(){
  data1 = (int*)malloc(sizeof(int));
  data2 = (int*)malloc(sizeof(int));
  data3 = (int*)malloc(sizeof(int));
  *data1 = 7;
  *data2 = 8;
  *data3 = 78;
  list1 = new_list(data1);
  list2 = new_list(data3);
}

void list_teardown(){
  if( data1 ){ free(data1); }
  if( data2 ){ free(data2); }
  if( data3 ){ free(data3); }
  if( list1 ){ delete_linked_list(list1); }
  if( list2 ){ delete_linked_list(list2); }
}

START_TEST(test_create){
  fail_unless( list1 != NULL, "Failed to create a new list." );
  fail_unless( get_data_linked_list(list1) == data1, 
               "List head does not point to correct data." );
}
END_TEST

START_TEST(test_add){
  // add to list
  LINKED_LIST_T* node_added = add_linked_list( list1, data2 );
  fail_unless( get_data_linked_list(node_added) == data2, 
               "Failed to correctly add data to a new list node");
  fail_unless( get_next_linked_list(list1) == node_added, 
               "Failed to add a node to list1");

  // add again to head
  LINKED_LIST_T* end = node_added;
  node_added = add_linked_list( list1, data1 );
  fail_unless( get_data_linked_list(node_added) == data1, 
               "Failed to correctly add data to a new list node");
  fail_unless( get_next_linked_list(end) == node_added, 
               "Failed to add a node to list1");

  // add to end
  end = node_added;
  node_added = add_linked_list( end, data3);
  fail_unless( get_next_linked_list(end) == node_added, 
               "Failed to add a node to end of list");

  // can you get the new end from the head?
  LINKED_LIST_T* cur_node = list1;
  while( has_next_linked_list(cur_node) == TRUE ) {
    cur_node = get_next_linked_list(cur_node);
  }
  fail_unless( cur_node == node_added, "Did not reach the end from the head");
}
END_TEST

START_TEST(test_combine){
  // lengthen list1
  add_linked_list(list1, data1);
  add_linked_list(list1, data1);
  LINKED_LIST_T* end = add_linked_list(list1, data1);

  // add on list2
  combine_lists(list1, list2);

  // check that the end of list1 points to list2
  fail_unless( get_next_linked_list(end) == list2,
               "Failed to combine lists." );
}
END_TEST

START_TEST(test_copy){
  // copy a list with one element
  LINKED_LIST_T* copy = copy_list(list1);
  fail_unless( get_data_linked_list(copy) == data1,
               "Failed to copy list1");
  delete_linked_list(copy);

  // copy a list with three elements
  add_linked_list(list1, data2);
  add_linked_list(list1, data3);
  copy = copy_list(list1);
  fail_unless( get_data_linked_list(copy) == data1,
               "Failed to copy first element of 3 element list");
  LINKED_LIST_T* next = get_next_linked_list(copy);
  fail_unless( get_data_linked_list(next) == data2,
               "Failed to copy second element of 3 element list");
  next = get_next_linked_list(next);
  fail_unless( get_data_linked_list(next) == data3,
               "Failed to copy third element of 3 element list");
  fail_unless( get_next_linked_list(next) == NULL, 
               "End of copied list does not point to NULL");
}
END_TEST

START_TEST(test_delete){
  // create a list where each element points to same data
  LINKED_LIST_T* head = new_list(data1);
  int i=0;
  for(i=0; i<6; i++){
    add_linked_list(head, data1);
  }
  // delete the list and make sure that data remains unchanged
  delete_linked_list(head);
  fail_unless( *data1 == 7, "Data was corrupted after deleting list.");
}
END_TEST

/* Test Case for boundry conditions */
START_TEST(test_null){
  LINKED_LIST_T* null_list = NULL;
  // get data from empty list, return NULL
  fail_unless( get_data_linked_list( null_list ) == NULL, 
               "Failed to return null from get data on null list");
  // get next from empty list, return NULL
  fail_unless( get_next_linked_list( null_list) == NULL, 
               "Failed to get NULL from next of empty list");
  // has next from empty list, return FALSE
  fail_unless( has_next_linked_list( null_list) == FALSE, 
               "Failed on has next for empty list");
  // add to empty list, check that a new node is returned
  fail_unless( add_linked_list(null_list, NULL) != NULL,
               "Failed to add to an empty list");
  // copy a null list and get NULL
  fail_unless( copy_list(null_list) == NULL,
               "Failed to return null from copying an empty list");
  // combine to empty list, get second
  // combine list with empty, end points to NULL
  // delete null list
}
END_TEST

Suite* list_suite(){
  Suite* s = suite_create("Linked-list");
  // Test basic features
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  tcase_add_test(tc_core, test_add);
  tcase_add_test(tc_core, test_combine);
  tcase_add_test(tc_core, test_copy);
  tcase_add_test(tc_core, test_delete);
  tcase_add_checked_fixture(tc_core, list_setup, list_teardown);

  // Test boundry conditions
  TCase *tc_limits = tcase_create("Limits");
  suite_add_tcase(s, tc_limits);
  tcase_add_test(tc_limits, test_null);
  tcase_add_checked_fixture(tc_limits, list_setup, list_teardown);

  return s;
}













