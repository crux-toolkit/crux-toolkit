/**
 * \file linked_list.h
 * Header file for a linked list class.
 */

#include "utils.h"
#include "objects.h"

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
LINKED_LIST_T* add_linked_list(LINKED_LIST_T* list, void* add_me);

/**
 * \brief Get the item pointed to by this node in the list.
 * \returns Returns a void pointer to the item held at this node.
 */
void* get_data_linked_list(LINKED_LIST_T* list);

/**
 * \brief Get the next element in the list.
 * \returns A pointer to the next list node.  NULL if at end of list.
 */
LINKED_LIST_T* get_next_linked_list(LINKED_LIST_T* node);

/**
 * \brief Is this list element at the end of the list?
 * \returns FALSE if node->next is NULL, else TRUE
 */
BOOLEAN_T has_next_linked_list(LINKED_LIST_T* node);

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
 * \brief Deletes the given list node and all list nodes after this
 * one WITHOUT deleting the data pointed to.
 * \returns void
 */
void delete_linked_list(LINKED_LIST_T* list);

/**
 * \brief Deletes the given list node leaving the data intact and
 * returns a pointer to the next item in the list, NULL if this node
 * is the last.
 */
LINKED_LIST_T* get_next_free_this_linked_list(LINKED_LIST_T* list);



















