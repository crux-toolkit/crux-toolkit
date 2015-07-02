#include "check-parameter.h"
#include "linked_list.h"
#include "utils.h"

// declare things to set up
int *data1, *data2, *data3;
LINKED_LIST_T *list1, *list2, *list3;
LIST_POINTER_T *nodea, *nodeb;

void list_setup(){
  data1 = (int*)malloc(sizeof(int));
  data2 = (int*)malloc(sizeof(int));
  data3 = (int*)malloc(sizeof(int));
  *data1 = 7;
  *data2 = 8;
  *data3 = 78;
  list1 = new_list(data1);
  list2 = new_list(data3);
  list3 = new_empty_list();
  nodea = NULL;
  nodeb = NULL;
}

void list_teardown(){
  if( data1 ){ free(data1); }
  if( data2 ){ free(data2); }
  if( data3 ){ free(data3); }
  if( list1 ){ delete_linked_list(list1); }
  if( list2 ){ delete_linked_list(list2); }
}

START_TEST(test_create){
  // test empty list
  fail_unless( list3 != NULL, "Failed to create a new empty list.");
  LIST_POINTER_T* first_node = get_first_linked_list(list3);
  fail_unless( first_node == NULL, "Empty list should not have a first node");
  fail_unless( is_empty_linked_list( list3 ) == true, 
               "List3 should be empty on creation.");

  // test one-item ilst
  fail_unless( list1 != NULL, "Failed to create a new list." );
  fail_unless( is_empty_linked_list(list1) == false,
               "List 1 should not be empty on creation.");
  first_node = get_first_linked_list( list1 );
  fail_unless( first_node != NULL,
               "List1 as created should have first item" );
  fail_unless( get_data_linked_list( first_node ) == data1, 
               "List head does not point to correct data." );
  fail_unless( has_next_linked_list( first_node ) == false,
               "First node of one-item list should not have next.");
  fail_unless( get_next_linked_list( first_node ) == NULL,
               "Next of newly created list is not NULL.");
}
END_TEST

START_TEST(test_push_front){
  // add to empty list
  LIST_POINTER_T* new_node = push_front_linked_list(list3, data1);
  fail_unless( is_empty_linked_list(list3) == false, 
               "Previously empty list should not be after push front.");
  fail_unless( new_node != NULL, 
               "Node added to empty list should not be NULL");
  fail_unless( new_node == get_first_linked_list(list3),
               "Head of empty list not set correctly after push front.");
  fail_unless( data1 == get_data_linked_list(new_node),
                "Data not correctly set to node added to empty list.");

  // add to list with one item
  // first find the existing node
  LIST_POINTER_T* was_first_node = get_first_linked_list( list1 );
  LIST_POINTER_T* end_node = was_first_node; // one item in list 

  LIST_POINTER_T* node_added = push_front_linked_list( list1, data2 );
  fail_unless( node_added != NULL, "Push_front returned a NULL pointer."); 
  // check list front
  fail_unless( get_first_linked_list(list1) == node_added,
               "Head of list does not point to added node after push front");
  // check data
  fail_unless( get_data_linked_list(node_added) == data2, 
               "Failed to correctly add data2 to a new list node");
  // check next
  fail_unless( get_next_linked_list( node_added ) == was_first_node, 
               "Node added to head does not have correct next.");

  // add again to same list
  was_first_node = node_added;
  node_added = push_front_linked_list( list1, data3 );
  fail_unless( node_added != NULL, "Second push front returned NULL pointer");
  fail_unless( get_first_linked_list(list1) == node_added,
               "Second push_front did not set list head to added node.");
  fail_unless( get_data_linked_list(node_added) == data3, 
               "Failed to correctly add data3 to a new list node");
  fail_unless( get_next_linked_list( node_added ) == was_first_node, 
               "Second push front did not set head->next correctly.");

  // can you get the new end from the head?
  LIST_POINTER_T* cur_node = get_first_linked_list(list1);
  while( has_next_linked_list(cur_node) == true ) {
    cur_node = get_next_linked_list(cur_node);
  }
  fail_unless( cur_node == end_node, 
               "After two push_back, did not reach the end from the head");

}
END_TEST

START_TEST(test_push_back){
  // add to empty list
  LIST_POINTER_T* new_node = push_back_linked_list(list3, data1);
  fail_unless( is_empty_linked_list(list3) == false, 
               "Previously empty list should not be after push back.");
  fail_unless( new_node != NULL, 
               "Node added to empty list should not be NULL");
  fail_unless( new_node == get_first_linked_list(list3),
               "Head of empty list not set correctly after push back.");
  fail_unless( data1 == get_data_linked_list(new_node),
                "Data not correctly set for node added to empty list.");

  // add to end
  LIST_POINTER_T* head = get_first_linked_list(list1);
  LIST_POINTER_T* was_end = head; // one item in list 

  LIST_POINTER_T* node_added = push_back_linked_list(list1, data2);
  // check returned pointer
  fail_unless(node_added != NULL, "Push back returned null pointer");
  // chech list head
  fail_unless( head == get_first_linked_list(list1),
               "Push back changed the head of the list");
  // check data
  fail_unless( get_data_linked_list(node_added) == data2, 
               "Data pushed back is not correct.");
  // check next
  fail_unless( get_next_linked_list(was_end) == node_added, 
               "Push back did not add correctly to the list end");
  fail_unless( has_next_linked_list(node_added) == false, 
               "Push back did not set end node next to NULL.");

  // repeat for a second data item
  was_end = node_added; 
  node_added = push_back_linked_list(list1, data3);
  fail_unless(node_added != NULL, "Push back returned null pointer");
  fail_unless( head == get_first_linked_list(list1),
               "Push back changed the head of the list");
  fail_unless( get_data_linked_list(node_added) == data3, 
               "Data pushed back is not correct.");
  fail_unless( get_next_linked_list(was_end) == node_added, 
               "Push back did not add correctly to the list end");
  fail_unless( has_next_linked_list(node_added) == false, 
               "Push back did not set end node next to NULL.");

  // can you get the new end from the head?
  LIST_POINTER_T* cur_node = get_first_linked_list(list1);
  while( has_next_linked_list(cur_node) == true ) {
    cur_node = get_next_linked_list(cur_node);
  }
  fail_unless( cur_node == node_added, 
               "After two push_back, did not reach the end from the head");

}
END_TEST

START_TEST(test_combine){
  // lengthen list1
  push_back_linked_list(list1, data1);
  push_back_linked_list(list1, data1);
  LIST_POINTER_T* end1 = push_back_linked_list(list1, data1);
  LIST_POINTER_T* begin2 = get_first_linked_list(list2);


  // add on list2
  combine_lists(list1, list2);

  // check that the end of list1 points to list2
  fail_unless( get_next_linked_list(end1) == get_first_linked_list(list2),
               "Failed to combine lists." );
  // check that list2 is unchanged
  fail_unless( get_data_linked_list( begin2 ) == data3,
               "List2 was changed after being added to list1");


  // just get rid of the head and the elements remain in list1
  free(list2);
  // set list2 to NULL so we don't try to delete it in clean up
  list2 = NULL;
  fail_unless( get_next_linked_list(end1) == begin2, 
               "Elements of list2 were lost when head was deleted" );

  // add an empty list to list1
  combine_lists(list1, list3);

  // end of list1 should be same
  fail_unless( has_next_linked_list( begin2 ) == false,
               "No new nodes should be on the end of list1");

  // add list1 to an empty list
  fail_unless( list3 != NULL, "List3 should not be null");
  fail_unless( is_empty_linked_list(list3) == true,
               "List3 should still be empty");
  fail_unless( get_first_linked_list(list3) == NULL, 
               "first element of list3 is not null");
  combine_lists(list3, list1);
  fail_unless( is_empty_linked_list(list3) == false,
               "List3 should no longer be empty");
}
END_TEST

START_TEST(test_copy){
  // copy a list with one element
  LINKED_LIST_T* copy = copy_list(list1);
  fail_unless( copy != NULL, "Copy returned a null pointer");
  fail_unless( get_data_linked_list( get_first_linked_list(copy) ) == data1,
               "Failed to copy list1");
  delete_linked_list(copy);

  // copy a list with three elements
  push_back_linked_list(list1, data2);
  push_back_linked_list(list1, data3);
  copy = copy_list(list1);
  fail_unless( copy != NULL, "Copy returned a null pointer");
  fail_unless( copy != list1, 
               "Copy should not point to list1 after creating a copy");
  LIST_POINTER_T* next = get_first_linked_list(copy);
  fail_unless( get_data_linked_list(next) == data1,
               "Failed to copy first element of 3 element list");

  // delete original list and check that copy is still there
  delete_linked_list(list1);
  list1 = NULL;

  fail_unless( get_data_linked_list(next) == data1,
               "Failed to copy first element of 3 element list after delete.");
  next = get_next_linked_list(next);
  fail_unless( get_data_linked_list(next) == data2,
               "Failed to copy second element of 3 element list");
  next = get_next_linked_list(next);
  fail_unless( get_data_linked_list(next) == data3,
               "Failed to copy third element of 3 element list");
  fail_unless( get_next_linked_list(next) == NULL, 
               "End of copied list does not point to NULL");

  // copy an empty list
  delete_linked_list(copy);
  copy = NULL;
  copy = copy_list( list3 );
  fail_unless( copy != NULL, 
               "Copied list should not be null after copying empty list.");
  fail_unless( is_empty_linked_list(copy) == true, 
               "Copied list should still be empty.");
  fail_unless( get_first_linked_list(copy) == NULL,
               "Should be no elements to get from empty list.");
}
END_TEST

START_TEST(test_clear){

  // clear one-element list
  clear_list(list1);
  fail_unless( is_empty_linked_list(list1) == true,
               "List should be empty after being cleared.");
  fail_unless( *data1 == 7, "Data was corrupted after clear");

  // clear empty list
  clear_list(list1);
  fail_unless( is_empty_linked_list(list1) == true,
               "Empty list should be empty after being cleared.");

  // add to cleared list
  int i = 0;
  for(i=0; i< 5; i++){
    push_back_linked_list(list1, data1);
  }
  fail_unless( !is_empty_linked_list(list1) == true,
               "List should not be empty.");

  // clear multi-element list
  clear_list(list1);
  fail_unless( is_empty_linked_list(list1) == true,
               "List should be empty after being cleared.");
  fail_unless( *data1 == 7, "Data was corrupted after multi-element clear");
}
END_TEST

START_TEST(test_delete){
  // create a list where each element points to same data
  LINKED_LIST_T* head = new_list(data1);
  int i=0;
  for(i=0; i<6; i++){
    //add_linked_list(head, data1);
    push_back_linked_list(head, data1);
  }
  // delete the list and make sure that data remains unchanged
  delete_linked_list(head);
  fail_unless( *data1 == 7, "Data was corrupted after deleting list.");
}
END_TEST

START_TEST(test_delete_next){
  // create a list 
  LINKED_LIST_T* head = new_list(data1);// node 1
  push_back_linked_list(head, data1); // node 2
  push_back_linked_list(head, data2); // node 3
  push_back_linked_list(head, data2); // node 4
  push_back_linked_list(head, data3); // node 5

  // to delete node 4, get a pointer to 3
  nodea = get_first_linked_list(head); // 1
  nodea = get_next_linked_list(nodea); // 2
  nodea = get_next_linked_list(nodea); // 3
  // also get node 4
  nodeb = get_next_linked_list(nodea);

  // confirm data
  fail_unless(data2 == (int*)get_data_linked_list(nodea),
              "Node 3 should point to %i but instead points to %i.",
              *data2, *(int*)get_data_linked_list(nodea));
  fail_unless(data2 == (int*)get_data_linked_list(nodeb),
              "Node 4 should point to %i but instead points to %i.",
              *data2, *(int*)get_data_linked_list(nodeb));

  delete_next_list_node(nodea); // now 3 points to 5
  // 3 still has data2
  fail_unless(data2 == (int*)get_data_linked_list(nodea),
              "Node 3 should point to %i but instead points to %i.",
              *data2, *(int*)get_data_linked_list(nodea));

  // 3 now points to node 5 which has data3
  nodeb = get_next_linked_list(nodea);
  fail_unless(data3 == (int*)get_data_linked_list(nodeb),
              "Node after 3 should point to %i but instead points to %i.",
              *data3, *(int*)get_data_linked_list(nodeb));

  // if we delete the end, node 3 should point to null
  delete_next_list_node(nodea); // now 3 points to null
  nodeb = get_next_linked_list(nodea);
  fail_unless(nodeb == NULL, 
              "After deleting two nodes, node 3 should be at the end.");

}
END_TEST

START_TEST(test_pop){
  // pop front/back emptylist
  fail_unless( pop_back_linked_list( list3 ) == NULL,
               "Data from back of empty list should be NULL");
  fail_unless( pop_front_linked_list( list3 ) == NULL,
               "Data from front of empty list should be NULL");

  // pop front/back list with one element
  fail_unless( pop_front_linked_list(list1) == data1,
               "Data1 should have been popped from front of list.");
  fail_unless( is_empty_linked_list(list1) == true,
               "List1 should be empty after only element popped.");
  fail_unless( pop_back_linked_list(list1) == NULL,
               "Now empty list should return NULL from pop_back.");

  fail_unless( pop_back_linked_list(list2) == data3,
               "Data1 should have been popped from front of list.");
  fail_unless( is_empty_linked_list(list2) == true,
               "List2 should be empty after only element popped.");
  fail_unless( pop_front_linked_list(list2) == NULL,
               "Now empty list should return NULL from pop_front.");

  // pop front/back of list with three elements
  push_front_linked_list(list1, data1);
  push_front_linked_list(list1, data2);
  push_front_linked_list(list1, data3);

  fail_unless( pop_front_linked_list(list1) == data3,
               "Data3 should have been popped from front of list.");
  fail_unless( pop_front_linked_list(list1) == data2,
               "Data2 should have been popped from front of list.");
  fail_unless( pop_front_linked_list(list1) == data1,
               "Data1 should have been popped from front of list.");

  push_front_linked_list(list1, data1);
  push_front_linked_list(list1, data2);
  push_front_linked_list(list1, data3);

  fail_unless( pop_back_linked_list(list1) == data1,
               "Data1 should have been popped from front of list.");
  fail_unless( pop_back_linked_list(list1) == data2,
               "Data2 should have been popped from front of list.");
  fail_unless( pop_back_linked_list(list1) == data3,
               "Data3 should have been popped from front of list.");
}
END_TEST

/* Test Case for boundry conditions */
START_TEST(test_check_null){
  LINKED_LIST_T* null_list = NULL;
  LIST_POINTER_T* null_list_p = NULL;

  fail_unless( is_empty_linked_list(null_list) == true, 
               "NULL list should be empty");
  // has next from empty list, return false
  fail_unless( has_next_linked_list( null_list_p) == false, 
               "Failed on has next for empty list");

}
END_TEST

START_TEST(test_get_null){
  LINKED_LIST_T* null_list = NULL;
  LIST_POINTER_T* null_list_p = NULL;
  fail_unless( get_first_linked_list(null_list) == NULL,
               "First of null list should be NULL.");
  fail_unless( get_data_linked_list( null_list_p ) == NULL, 
               "Failed to return null from get data on null list");
  // get next from empty list, return NULL
  fail_unless( get_next_linked_list( null_list_p ) == NULL, 
               "Failed to get NULL from next of empty list");
}
END_TEST

START_TEST(test_push_null){
  LINKED_LIST_T* null_list = NULL;

  // add data to null list
  fail_unless( push_back_linked_list( null_list, data1 ) == NULL,
               "Incorrect return for adding data to null list");
  fail_unless( push_front_linked_list( null_list, data1 ) == NULL,
               "Incorrect return for adding data to null list");
}
END_TEST

START_TEST(test_pop_null){
  LINKED_LIST_T* null_list = NULL;
  fail_unless( pop_front_linked_list(null_list) == NULL,
               "Pop front should return null from null list.");
  fail_unless( pop_back_linked_list(null_list) == NULL,
               "Pop back should return null from null list.");
}
END_TEST

START_TEST(test_null_data){
  // add null data to real list
  LIST_POINTER_T* added_node = push_front_linked_list(list1, NULL);
  // added node not NULL
  fail_unless( added_node != NULL,
               "Push front null data should return non-null node.");
  // get first is correct
  fail_unless( added_node == get_first_linked_list(list1),
               "Adding null data did not set head correctly");
  // data is null
  fail_unless( get_data_linked_list(added_node) == NULL,
               "Data should be NULL");

  added_node = push_back_linked_list(list1, NULL);
  // added node not NULL
  fail_unless( added_node != NULL,
               "Push back null data should return non-null node.");
  // data is null
  fail_unless( get_data_linked_list(added_node) == NULL,
               "Data should be NULL");

  // list1 should now have null data at beginning and end
  fail_unless( pop_front_linked_list(list1) == NULL,
               "Failed to pop NULL data off front of list.");
  fail_unless( pop_back_linked_list(list1) == NULL,
               "Failed to pop NULL data off end of list.");

  // create a list with NULL data
  LINKED_LIST_T* list = new_list(NULL);
  fail_unless( list != NULL,
               "New list with null data should not be null");
  LIST_POINTER_T* list_p = get_first_linked_list(list);
  fail_unless( list_p != NULL, 
               "Node with null data should not itself be NULL");
  fail_unless( get_data_linked_list( list_p ) == NULL,
               "NULL data is not null.");
  delete_linked_list( list );
}
END_TEST

START_TEST(test_do_null){
  LINKED_LIST_T* null_list = NULL;
  LIST_POINTER_T* null_list_p = NULL;
  // copy a null list and get NULL
  fail_unless( copy_list(null_list) == NULL,
               "Failed to return null from copying an empty list");
  // combine to empty list, get second
  fail_unless( combine_lists(null_list, list1) == list1,
               "Combining to a null list should return the second");
  // combine list with empty, end points to NULL
  LINKED_LIST_T* combined = combine_lists(list2, null_list);
  fail_unless( get_first_linked_list(combined) == get_first_linked_list(list2),
               "first and last element of combined and list2 should be same");

  // delete null list
  delete_linked_list(null_list);
  delete_list_node(null_list_p);
}
END_TEST

Suite* list_suite(){
  Suite* s = suite_create("Linked-list");
  // Test basic features
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  tcase_add_test(tc_core, test_push_front);
  tcase_add_test(tc_core, test_push_back);
  tcase_add_test(tc_core, test_combine);
  tcase_add_test(tc_core, test_copy);
  tcase_add_test(tc_core, test_clear);
  tcase_add_test(tc_core, test_delete);
  tcase_add_test(tc_core, test_pop);
  tcase_add_test(tc_core, test_delete_next);
  tcase_add_checked_fixture(tc_core, list_setup, list_teardown);

  // Test boundry conditions
  TCase *tc_limits = tcase_create("Limits");
  suite_add_tcase(s, tc_limits);
  //  tcase_add_test(tc_limits, test_null);
  tcase_add_test(tc_limits, test_check_null);
  tcase_add_test(tc_limits, test_get_null);
  tcase_add_test(tc_limits, test_push_null);
  tcase_add_test(tc_limits, test_pop_null);
  tcase_add_test(tc_limits, test_do_null);
  tcase_add_test(tc_limits, test_null_data);
  tcase_add_checked_fixture(tc_limits, list_setup, list_teardown);

  return s;
}













