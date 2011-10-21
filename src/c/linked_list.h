/**
 * \file linked_list.h
 * Header file for a linked list class.
 */
#ifndef LINKED_LIST_H
#define LINKED_LIST_H

#include <assert.h>
#include "utils.h"
#include "objects.h"
#include "carp.h"

/**
 * \brief Allocate the first element of a new list.
 * \returns A pointer to the new list.
 */
LINKED_LIST_T* new_empty_list();

/**
 * \brief Allocate the first element of a new list.
 * \returns A pointer to the new list.
 */
LINKED_LIST_T* new_list(void* first_item);

/**
 * \brief Add an item to the end of this list.  Traverses to the end
 * of the list, creates a new list node, and adds it to the end of the list.
 * \returns A pointer to the last (newly created) element in the list.
 */
LIST_POINTER_T* push_back_linked_list(LINKED_LIST_T* list_head, void* add_me);

/**
 * \brief Add an item to the beginning of this list.  Creates a new
 * list node, and insersts it between the head and the first node.
 * \returns A pointer to the first (newly created) element in the list.
 */
LIST_POINTER_T* push_front_linked_list(LINKED_LIST_T* list_head, void* add_me);

/**
 * \brief Remove the element from the end of the list (farthest from
 * the head) and return a pointer to its data. 
 * Deletes the last list node leaving the data and the rest of the
 * list intact. 
 * \returns A pointer to the data in the last item in the list, NULL
 * if list is empty.
 */
void* pop_back_linked_list(LINKED_LIST_T* list);

/**
 * \brief Remove the first element from the list, returning a pointer
 * to the data.
 * Deletes the first list node leaving the data and the
 * head of the list intact. Head now points to what was the second
 * element, if any.
 * \returns A pointer to the data in the first item in the list, NULL
 * if list is empty.
 */
void* pop_front_linked_list(LINKED_LIST_T* list);

/**
 * \brief Get a pointer to the first element in the list.  Can be used
 * to start a traverse of the list.
 * \returns A pointer to the first list node.  NULL if empty list.
 */
LIST_POINTER_T* get_first_linked_list(LINKED_LIST_T* head);

/**
 * \brief Get the next element in the list.  Can be used for
 * traversing a list.
 * \returns A pointer to the next list node.  NULL if at end of list.
 */
LIST_POINTER_T* get_next_linked_list(LIST_POINTER_T* node);

/**
 * \brief Get the item pointed to by this node in the list.
 * \returns Returns a void pointer to the item held at this node.
 */
void* get_data_linked_list(LIST_POINTER_T* list);

/**
 * \brief Is this list empty?
 * \returns TRUE if node->next is NULL, else FALSE
 */
bool is_empty_linked_list(LINKED_LIST_T* head);

/**
 * \brief Is this list element at the end of the list?
 * \returns FALSE if node->next is NULL, else TRUE
 */
bool has_next_linked_list(LIST_POINTER_T* node);

/**
 * \brief Combines two lists by having the end of the first list point
 * to the beginning of the second.  If the first list is null, returns
 * a pointer to the second list.
 * \returns Returns a pointer to the beginning of the combined lists.
 */
LINKED_LIST_T* combine_lists(LINKED_LIST_T* beginning, LINKED_LIST_T* add_me);

/**
 * \brief Create a copy of each node in the original list, linked in
 * the same order.
 * \return A pointer to the head of the new list.
 */
LINKED_LIST_T* copy_list(LINKED_LIST_T* original);

/**
 * \brief Removes all nodes from a list while leaving the data intact.
 */
void clear_list(LINKED_LIST_T* list);

/**
 * \brief Deletes the given list node and all list nodes after this
 * one WITHOUT deleting the data pointed to.
 * \returns void
 */
void delete_linked_list(LINKED_LIST_T* list);

/**
 * \brief Deletes the given list node and all list nodes after this
 * one WITHOUT deleting the data pointed to.
 * \returns void
 */
void delete_list_node(LIST_POINTER_T* list);

/**
 * \brief Deletes only the node AFTER this one.  The given node then
 * points to what the next node pointed to before. 
 * \returns A pointer to the new list element that follows the given one.
 */
LIST_POINTER_T* delete_next_list_node(LIST_POINTER_T* pre_node);













#endif // LINKED_LIST_H
