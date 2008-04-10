/**
 * \file linked_list.c
 * Method definitions for a linked list class.
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


/**
 * \brief Allocate the first element of a new list.
 * \returns A pointer to the new list.
 */
LINKED_LIST_T* new_list(void* first_item){
  LINKED_LIST_T* head = (LINKED_LIST_T*)mymalloc(sizeof(LINKED_LIST_T));
  head->data = first_item;
  head->next = NULL;
  return head;
}

/**
 * \brief Add an item to the end of this list.  Traverses to the end
 * of the list, creates a new list node, and adds it to the end of the list.
 * \returns A pointer to the last (newly created) element in the list.
 */
LINKED_LIST_T* add_linked_list(LINKED_LIST_T* list, void* add_me){

  if( list == NULL ){
    return new_list(add_me);
  }

  // find the end of the list
  while(list->next != NULL ){
    list = list->next;
  }
  list->next = new_list(add_me);
  return list->next;
}

/**
 * \brief Get the item pointed to by this node in the list.
 * \returns Returns a void pointer to the item held at this node.
 */
void* get_data_linked_list(LINKED_LIST_T* list){
  if( list == NULL ){
    return NULL;
  }
  return list->data;
}

/**
 * \brief Get the next element in the list.
 * \returns A pointer to the next list node.  NULL if at end of list.
 */
LINKED_LIST_T* get_next_linked_list(LINKED_LIST_T* node){
  if( node == NULL ){
    return NULL;
  }
  return node->next;
}

/**
 * \brief Is this list element at the end of the list?
 * \returns FALSE if node->next is NULL, else TRUE
 */
BOOLEAN_T has_next_linked_list(LINKED_LIST_T* node){
  if( node == NULL || node->next == NULL ){
    return FALSE;
  }
  return TRUE;
}

/**
 * \brief Combines two lists by having the end of the first list point
 * to the beginning of the second.  If the first list is null, returns
 * a pointer to the second list.
 * \returns Returns a pointer to the beginning of the combined lists.
 */
LINKED_LIST_T* combine_lists(LINKED_LIST_T* first, LINKED_LIST_T* second){
  if( first == NULL ){
    return second;
  }
  while( first->next != NULL ){
    first = first->next;
  }
  first->next = second;

  return first;
}

/**
 * \brief Create a copy of each node in the original list, linked in
 * the same order.
 * \return A pointer to the head of the new list.
 */
LINKED_LIST_T* copy_list(LINKED_LIST_T* original){
  LINKED_LIST_T* copy = NULL;

  if( original != NULL ){
    copy = new_list( original->data );
    copy->next = copy_list(original->next);
  }

  return copy;
}

/**
 * \brief Deletes the given list node and all list nodes after this
 * one WITHOUT deleting the data pointed to.
 * \returns void
 */
void delete_linked_list(LINKED_LIST_T* list){
  if( list ){
    if( list->next != NULL ){
      delete_linked_list(list->next);
    }// else next is null and this is the end of the list
    
    free(list);
  }
  return;
}

/**
 * \brief Deletes the given list node leaving the data intact and
 * returns a pointer to the next item in the list, NULL if this node
 * is the last.
 */
LINKED_LIST_T* get_next_free_this_linked_list(LINKED_LIST_T* list){

  if( list == NULL ){
    return NULL;
  }

  LINKED_LIST_T* next_node = list->next;

  free(list);
  return next_node;
}

