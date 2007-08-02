/**
 * \file hash.h 
 * $Revision: 1.1 $
 * \brief Object for hashing.
 */
#include "crux-utils.h"
#include "objects.h"

#ifndef __HASH__
#define __HASH__

/**
 * Create new hashtable with capacity.
 *\returns the hash table
 */
HASH_T* new_hash(
  unsigned int capacity ///< The capacity of the new hash table
  );

/**
 * free hashtable
 */
void free_hash(
  HASH_T* h ///< Hash object to free -in
  );

/**
 * add key and value to hash table.
 * Must add a heap allocated key, value may be NULL
 *\returns TRUE if successfully adds to new record, else FALSE
 */
BOOLEAN_T add_hash(
  HASH_T* h, ///< Hash object to add -in/out
  char *key, ///< key of the record to add -in
  void *value ///< value to add to be hashed if needed -in
  );

/**
 * Get the value of the record for the hash key
 *\return the value for the hash record, returns NULL if can't find key
 */
void* get_hash_value(
  HASH_T* h, ///< working hash object -in
  char *key  ///< the key of the record to retrieve -in
  );

/**
 * Get the count of the record for the hash key
 *\return the count for the hash record, returns NULL if can't find key
 */
int get_hash_count(
  HASH_T* h, ///< working hash object -in
  char *key  ///< the key of the record to retrieve -in
  );

/**
 * Remove key from table, returning value.
 *\returns the value of removing record, returns NULL if can't find key
 */
void* remove_hash(
  HASH_T* h, ///< working hash object -in
  char *key  ///< the key of the record to remove -in
  );

/**
 *\returns total number of keys in the hashtable.
 */  
unsigned int hash_size(
  HASH_T* h ///< working hash object -in
  );

#endif


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
