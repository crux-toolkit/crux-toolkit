/**
 * \file linked_list.cpp
 * \brief Method definitions for a linked list class.
 */

#include "linked_list.h"

/**
 * \struct _linked_list_node
 * \brief A node of a linked list
 */
struct _linked_list_node{
  void* data;
  struct _linked_list_node *next;
};
//typedef struct _linked_list_node LIST_POINTER_T*;

/**
 * \struct _linked_list_head
 * \brief The head of a linked list
 */
struct _linked_list_head{
  // TODO change next to head or first
  struct _linked_list_node *next;
  struct _linked_list_node *last;
};
//typedef struct _linked_list_head LINKED_LIST_T;

/**
 * \brief Allocate the first element of a new list.
 * \returns A pointer to the new list.
 */
LINKED_LIST_T* new_empty_list(){
  LINKED_LIST_T* head = (LINKED_LIST_T*)mymalloc(sizeof(LINKED_LIST_T));
  head->next = NULL;
  head->last = NULL;
  return head;
}

/**
 * \brief Private function for adding a node to the list.  
 * Allocates the node, sets data field to new_item, sets next to NULL.
 * \returns A pointer to the newly created node.
 */
LIST_POINTER_T* new_node(void* new_item){
  LIST_POINTER_T* the_node = (LIST_POINTER_T*)mymalloc(sizeof(LIST_POINTER_T));
  the_node->data = new_item;
  the_node->next = NULL;
  return the_node;
}

/**
 * \brief Allocate the first element of a new list.
 * \returns A pointer to the new list.
 */
LINKED_LIST_T* new_list(void* first_item){
  LINKED_LIST_T* head = (LINKED_LIST_T*)mymalloc(sizeof(LINKED_LIST_T));
  head->next = new_node(first_item);
  head->last = head->next;
  return head;
}

/**
 * \brief Add an item to the end of this list.  Traverses to the end
 * of the list, creates a new list node, and adds it to the end of the list.
 * \returns A pointer to the last (newly created) element in the list.
 */
LIST_POINTER_T* push_back_linked_list(LINKED_LIST_T* list_head, void* add_me){

  if( list_head == NULL ){
    carp(CARP_ERROR, "Cannot add to a NULL list.");
    return NULL;
  }

  // adding to an empty list
  if( list_head->next == NULL ){
    list_head->next = new_node(add_me);
    list_head->last = list_head->next;
    return list_head->next;
  }

  // find the end of the list
  LIST_POINTER_T* last_node = list_head->last;
  assert( (list_head->last)->next == NULL);

  last_node->next = new_node(add_me);
  list_head->last = last_node->next;
  return last_node->next;
}

/**
 * \brief Add an item to the beginning of this list.  Creates a new
 * list node, and insersts it between the head and the first node.
 * \returns A pointer to the first (newly created) element in the list.
 */
LIST_POINTER_T* push_front_linked_list(LINKED_LIST_T* list_head, void* add_me){

  if( list_head == NULL ){
    carp(CARP_ERROR, "Cannot add to null list\n");
    return NULL;
  }

  LIST_POINTER_T* add_this_node = new_node(add_me);
  add_this_node->next = list_head->next;
  list_head->next = add_this_node;
  if( list_head->last == NULL ){
    list_head->last = add_this_node;
  }
  return add_this_node;
}

/**
 * \brief Remove the element from the end of the list (farthest from
 * the head) and return a pointer to its data. 
 * Deletes the last list node leaving the data and the rest of the
 * list intact. 
 * \returns A pointer to the data in the last item in the list, NULL
 * if list is empty.
 */
void* pop_back_linked_list(LINKED_LIST_T* list){

  if( list == NULL || list->next == NULL ){
    return NULL;
  }

  // find the end of the list
  LIST_POINTER_T* last_node = list->last;
  void* data = last_node->data;

  // special case of one item in list
  if( last_node == list->next ){
    list->last = NULL;
    list->next = NULL;
  }else{
    // find the node pointing to last
    LIST_POINTER_T* second_to_last_node = list->next;

    while( second_to_last_node->next != last_node ){
      second_to_last_node = second_to_last_node->next;
    }

    list->last = second_to_last_node;
    second_to_last_node->next = NULL;
  }
  free(last_node);
  return data;
}

/**
 * \brief Remove the first element from the list, returning a pointer
 * to the data.
 * Deletes the first list node leaving the data and the
 * head of the list intact. Head now points to what was the second
 * element, if any.
 * \returns A pointer to the data in the first item in the list, NULL
 * if list is empty.
 */
void* pop_front_linked_list(LINKED_LIST_T* list){

  if( list == NULL || list->next == NULL ){
    return NULL;
  }

  LIST_POINTER_T* first_node = list->next;
  void* data = first_node->data;
  list->next = first_node->next;
  // if only one node, change last to NULL
  if( list->next == NULL ){
    list->last = NULL;
  }

  free(first_node);
  return data;
}

/**
 * \brief Get a pointer to the first element in the list.  Can be used
 * to start a traverse of the list.
 * \returns A pointer to the first list node.  NULL if empty list.
 */
LIST_POINTER_T* get_first_linked_list(LINKED_LIST_T* head){
  if( head == NULL ){
    return NULL;
  }
  return head->next;
}

/**
 * \brief Get the next element in the list.  Can be used for
 * traversing a list.
 * \returns A pointer to the next list node.  NULL if at end of list.
 */
LIST_POINTER_T* get_next_linked_list(LIST_POINTER_T* node){
  if( node == NULL ){
    return NULL;
  }
  return node->next;
}

/**
 * \brief Get the item pointed to by this node in the list.
 * \returns Returns a void pointer to the item held at this node.
 */
void* get_data_linked_list(LIST_POINTER_T* node){
  if( node == NULL ){
    return NULL;
  }
  return node->data;
}

/**
 * \brief Is this list empty?
 * \returns true if node->next is NULL, else false
 */
bool is_empty_linked_list(LINKED_LIST_T* head){
  if( head == NULL || head->next == NULL ){
    return true;
  }
  return false;
}

/**
 * \brief Is this list element at the end of the list?
 * \returns false if node->next is NULL, else true
 */
bool has_next_linked_list(LIST_POINTER_T* node){
  if( node == NULL || node->next == NULL ){
    return false;
  }
  return true;
}

/**
 * \brief Combines two lists by having the end of the first list point
 * to the first element of the second.  If the first list is null, returns
 * a pointer to the second list. Does not change the second list.
 * \returns Returns a pointer to the beginning of the combined lists.
 */
LINKED_LIST_T* combine_lists(LINKED_LIST_T* first, LINKED_LIST_T* second){
  if( first == NULL ){
    return second;
  }
  if( second == NULL ){
    return first;
  }

  // check for empty list
  if( first->next == NULL ){
    first->next = second->next;
    first->last = second->last;
  }else{
    // find the end of the first list
    LIST_POINTER_T* last_node = first->last;
    
    // set to the beginning of the second list
    last_node->next = second->next;
    
    // set the end of the first to the end of the last  
    first->last = second->last;
  }
  return first;
}

/**
 * \brief Create a copy of each node in the original list, linked in
 * the same order.
 * \return A pointer to the head of the new list.
 */
LINKED_LIST_T* copy_list(LINKED_LIST_T* original){
  if( original == NULL ){
    return NULL;
  }

  LINKED_LIST_T* copy = new_empty_list();
  // check for empty list
  if( is_empty_linked_list( original ) ){
    return copy;
  }

  //initialize first element
  LIST_POINTER_T* cur_original = get_first_linked_list(original);
  copy->next = new_node( get_data_linked_list(cur_original) );

  LIST_POINTER_T* cur_copy = copy->next;
  cur_original = cur_original->next;
  
  while( cur_original != NULL ){
    cur_copy->next = new_node( cur_original->data );
    cur_copy = cur_copy->next;
    cur_original = cur_original->next;
  }
  
  copy->last = cur_copy->next;
  return copy;
}

// TODO? delete_only_this_list_node?
/**
 * \brief Deletes the given list node and all list nodes after this
 * one WITHOUT deleting the data pointed to.
 * \returns void
 */
void delete_list_node(LIST_POINTER_T* list){
  if( list == NULL ){
    return;
  }

  if( list->next != NULL ){
    delete_list_node(list->next);
  }// else next is null and this is the end of the list
  
  free(list);
  return;
}

/**
 * \brief Removes all nodes from a list while leaving the data intact.
 */
void clear_list(LINKED_LIST_T* list){
  if( list == NULL ){
    return;
  }
  LIST_POINTER_T* first = list->next;
  delete_list_node(first);
  list->next = NULL;
}



/**
 * \brief Deletes the given list WITHOUT deleting the data pointed to.
 * \returns void
 */
void delete_linked_list(LINKED_LIST_T* list){
  if( list ){
    if( list->next != NULL ){
      delete_list_node(list->next);
    }// else next is null and this is the end of the list
    
    free(list);
  }
  return;
}

/**
 * \brief Deletes only the node AFTER this one.  The given node then
 * points to what the next node pointed to before. The data remains
 * untouched. 
 * \returns A pointer to the new list element that follows the given one.
 */
LIST_POINTER_T* delete_next_list_node(LIST_POINTER_T* pre_node){

  if( pre_node == NULL ){
    return NULL;
  }
  // before: head ....pre-> delete-> next...NULL
  // after : head ....pre-> next...NULL

  LIST_POINTER_T* delete_node = pre_node->next;
  if( delete_node == NULL ){
    return NULL;
  }
  LIST_POINTER_T* next_node = delete_node->next;

  pre_node->next = next_node;
  free(delete_node);

  return next_node;
}


/*
// add node after this but before what follows
void insert_linked_list(LIST_POINTER_T* node){
}
*/
