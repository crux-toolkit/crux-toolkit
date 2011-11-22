/*************************************************************************//**
 * \file hash.cpp
 * AUTHOR: David Crawshaw, Chris Park
 * \brief Object for hashing.
 ****************************************************************************/
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hash.h"
//#include "carp.h"
#include "parse_arguments.h"
#include "crux-utils.h"
#include "objects.h"
#include "parameter.h"

// TODO why does hash include parameter and not the other way around?

// Table is sized by primes to minimise clustering.
static const unsigned int sizes[] = {
    53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317,
    196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843,
    50331653, 100663319, 201326611, 402653189, 805306457, 1610612741
};

static const unsigned int sizes_count = sizeof(sizes) / sizeof(sizes[0]);
static const FLOAT_T load_factor = 0.65;

/**
 * \struct record
 * \brief record for each value/key pair
 */
struct record {
  unsigned int hash; ///< hash algorithm code
  char* key; ///< the key for the record
  void* value; ///< value the record contains
  int count; ///< count for duplicate adds(if adding two of the same key count 1 is added)
};

/**
 * \struct hash
 * \brief hash table, contains the records
 */
struct our_hash {
  RECORD_T* records; ///< record holds key & values
  unsigned int records_count; ///< number of records
  unsigned int size_index; ///< index into the size array, thus can get the size of the hash table
};

/**
 *\struct hash_iterator
 *\brief An object that iterates over the keys in a hash 
 */
struct hash_iterator {
  HASH_T* hash; ///< the hash to iterate over
  int hash_idx;   ///< current hash key to return
  int hash_total; ///< total number of hash keys
};

// Function definition, description found below 
bool add_hash_when_grow(
  HASH_T* h, ///< Hash object to add -in/out
  const char *key, ///< key of the record to add -in
  const void *value, ///< value to add to be hashed if needed -in
  int count  ///< the count of the record to grow -in
  );

/**
 * Increase the size of hash table
 *\returns the true if successfully inceased hash table size, else false
 */
static int hash_grow(
  HASH_T* h ///< Hash object to grow -out
  )
{
  unsigned int i;
  RECORD_T* old_recs;
  unsigned int old_recs_length;
  
  old_recs_length = sizes[h->size_index];
  old_recs = h->records;
  
  if (h->size_index == sizes_count - 1){
    return false;
  }
  // increase larger hash table
  if ((h->records = (RECORD_T*)mycalloc(sizes[++h->size_index],
                             sizeof(RECORD_T))) == NULL) {
    h->records = old_recs;
    return false;
  }
  
  h->records_count = 0;
  
  // rehash table
  for (i=0; i < old_recs_length; i++){
    if (old_recs[i].hash && old_recs[i].key){
      add_hash_when_grow(h, old_recs[i].key, old_recs[i].value, old_recs[i].count);
    }
  }
  
  free(old_recs);
  
  return true;
}

/**
 * hash algorithm 
 * \returns the array slot
 */
static unsigned int strhash(
  const char *str ///< string into the hash function -in
  )
{
  int c;
  int hash = 5381;
  while ((c = *str++)){
    hash = hash * 33 + c;
  }
  return hash == 0 ? 1 : hash;
}

/**
 * Create new hashtable with capacity.
 *\returns the hash table
 */
HASH_T* new_hash(
  unsigned int capacity ///< The capacity of the new hash table
  ) 
{
  HASH_T* h;
  unsigned int i, sind = 0;
  
  capacity = (unsigned int)((FLOAT_T)capacity / load_factor);
  
  for (i=0; i < sizes_count; i++) 
    if (sizes[i] > capacity) { sind = i; break; }
  
  if ((h = (HASH_T*)malloc(sizeof(HASH_T))) == NULL) return NULL;
  if ((h->records = (RECORD_T*)mycalloc(sizes[sind], sizeof(RECORD_T))) == NULL) {
    free(h);
    return NULL;
  }
  
  h->records_count = 0;
  h->size_index = sind;
  
  return h;
}

/**
 * free hashtable
 */
void free_hash(
  HASH_T* h ///< Hash object to free -in
  )
{
  unsigned int idx = 0;
  unsigned int size = sizes[h->size_index];
  // free up all key & values strings
  for(; idx < size; ++idx){
    if(h->records[idx].key != NULL){
      free(h->records[idx].key);
    }
    if(h->records[idx].value != NULL){
      free(h->records[idx].value);
    }    
  }
  
  free(h->records);
  free(h);
}

/**
 * add key and value to hash table.
 * If key exists, free current value and allocate and set new one
 * If key not found, allocate key, and allocate and set value
 * Does not copy value (for use with void pointers).
 *\returns true if successfully adds to new record, else false
 */
bool add_or_update_hash(
  HASH_T* h, ///< Hash object to add to -in/out
  const char *key, ///< key of the record to add or update -in
  const void *value ///< value to associate with the key -in
  )
{
  RECORD_T* recs;  
  int rc;
  unsigned int off, ind, size, code;

  if (key == NULL || *key == '\0') return false;
  
  code = strhash(key);
  recs = h->records;
  size = sizes[h->size_index];
  
  ind = code % size;
  off = 0;
  
  // probe down until reaching open slot
  // Quadratic probing used
  while (recs[ind].key){     
    // if find duplicate key, thus identical item
    if ((code == recs[ind].hash) && recs[ind].key &&
        strcmp(key, recs[ind].key) == 0){
      // free existing value
      free(recs[ind].value); 
      // set new value
      recs[ind].value = (void*)my_copy_string((char*)value);
      return true;
    }
    else{
      // continue to search
      ind = (code + (int)pow(++off,2.0)) % size;
    }
  }

  // key not found, add it
  // first check size
  if (h->records_count > sizes[h->size_index] * load_factor) {
    rc = hash_grow(h);
    if (rc) return false;
  }


  recs[ind].hash = code;
  recs[ind].key = my_copy_string(key);
  recs[ind].value = (void*)my_copy_string((char*)value);
  recs[ind].count = 1;
  
  h->records_count++;
  
  return true;
}




/**
 * add key and value to hash table.
 * Must add a heap allocated key, value may be NULL
 * If finds duplicate key, just increase count by 1
 *\returns true if successfully adds to new record, else false
 */
bool add_hash(
  HASH_T* h, ///< Hash object to add -in/out
  const char *key, ///< key of the record to add -in
  const void *value ///< value to add to be hashed if needed -in
  )
{
    RECORD_T* recs;
    int rc;
    unsigned int off, ind, size, code;

    if (key == NULL || *key == '\0') return false;
    if (h->records_count > sizes[h->size_index] * load_factor) {
        rc = hash_grow(h);
        if (rc) return false;
    }

    code = strhash(key);
    recs = h->records;
    size = sizes[h->size_index];

    ind = code % size;
    off = 0;

    // probe down until reach open slot
    // Quadratic probing used
    while (recs[ind].key){     
      // if find duplicate key, thus identical item
      if ((code == recs[ind].hash) && recs[ind].key &&
          strcmp(key, recs[ind].key) == 0){
        // increment count
        ++recs[ind].count;        
        return true;
      }
      else{
        // continue to search
        ind = (code + (int)pow(++off,2.0)) % size;
      }
    }
    
    recs[ind].hash = code;
    recs[ind].key = my_copy_string(key);
    recs[ind].value = (void*)my_copy_string((char*)value);
    recs[ind].count = 1;

    h->records_count++;
    
    return true;
}

/**
 * add key and value to hash table.
 *\returns true if successfully adds to new record, else false
 */
bool add_hash_when_grow(
  HASH_T* h, ///< Hash object to add -in/out
  const char *key, ///< key of the record to add -in
  const void *value, ///< value to add to be hashed if needed -in
  int count  ///< the count of the record to grow -in
  )
{
    RECORD_T* recs;
    int rc;
    unsigned int off, ind, size, code;

    if (key == NULL || *key == '\0') return false;
    if (h->records_count > sizes[h->size_index] * load_factor) {
        rc = hash_grow(h);
        if (rc) return false;
    }

    code = strhash(key);
    recs = h->records;
    size = sizes[h->size_index];

    ind = code % size;
    off = 0;

    // probe down until reach open slot
    while (recs[ind].key){
      ind = (code + (int)pow(++off,2.0)) % size;
    }
    
    recs[ind].hash = code;
    recs[ind].key = my_copy_string(key);
    recs[ind].value = (void*)my_copy_string((char*)value);
    recs[ind].count = count;
    
    h->records_count++;
    
    return true;
}

/**
 * Updates the value for the key
 * Must already have a existing value for the key
 * Copies the value, thus no need to pass in a heap allocated value
 *\returns true if successfully updates hash value, else false
 */
bool update_hash_value(
  HASH_T* h, ///< Hash object to add -in/out
  const char *key, ///< key of the record to update -in
  const void *value ///< value to add to be hash -in
  )
{
  RECORD_T* recs;  
  unsigned int off, ind, size, code;
  
  if (key == NULL || *key == '\0') return false;
  
  code = strhash(key);
  recs = h->records;
  size = sizes[h->size_index];
  
  ind = code % size;
  off = 0;
  
  // probe down until reach open slot
  // Quadratic probing used
  while (recs[ind].key){     
    // if find duplicate key, thus identical item
    if ((code == recs[ind].hash) && recs[ind].key &&
        strcmp(key, recs[ind].key) == 0){
      // free existing value
      free(recs[ind].value); 
      // set new value
      recs[ind].value = my_copy_string((char*)value); 
      return true;
    }
    else{
      // continue to search
      ind = (code + (int)pow(++off,2.0)) % size;
    }
  }
  //carp(CARP_ERROR, "Failed to find key %s in hash table", key);
  return false;
}

/**
 * Get the value of the record for the hash key
 *\return the value for the hash record, returns NULL if can't find key
 */
void* get_hash_value(
  HASH_T* h, ///< working hash object -in
  const char *key  ///< the key of the record to retrieve -in
  )
{
  RECORD_T* recs;
  unsigned int off, ind, size;
  unsigned int code = strhash(key);
  
  recs = h->records;
  size = sizes[h->size_index];
  ind = code % size;
  off = 0;
  
  // search on hash which remains even if a record has been removed,
  // so remove_hash() does not need to move any collision records
  while (recs[ind].hash) {
    if ((code == recs[ind].hash) && recs[ind].key &&
        strcmp(key, recs[ind].key) == 0)
      return recs[ind].value;
    ind = (code + (int)pow(++off,2.0)) % size;
  }
  
  return NULL;
}

/**
 * Get a pointer to the variable pointing to the value of the record 
 * for the hash key
 * BAD!! This is not the right thing to do.  It is so parse_arguments can
 * directly insert values into the hash without using the key.  Soon I should
 * change parse_arguments so that it takes a reference to the hash and uses
 * the option name to insert the value into the hash
 *\return a pointer to variable pointing to the value for the hash record, returns NULL if can't find key
 */
void** get_hash_value_ref(
  HASH_T* h, ///< working hash object -in
  char *key  ///< the key of the record to retrieve -in
  )
{
  RECORD_T* recs;
  unsigned int off, ind, size;
  unsigned int code = strhash(key);
  
  recs = h->records;
  size = sizes[h->size_index];
  ind = code % size;
  off = 0;
  
  // search on hash which remains even if a record has been removed,
  // so remove_hash() does not need to move any collision records
  while (recs[ind].hash) {
    if ((code == recs[ind].hash) && recs[ind].key &&
        strcmp(key, recs[ind].key) == 0)
      return &recs[ind].value;
    ind = (code + (int)pow(++off,2.0)) % size;
  }
  
  return NULL;
}

/**
 * Get the count of the record for the hash key
 *\return the count for the hash record, returns -1 if can't find key
 */
int get_hash_count(
  HASH_T* h, ///< working hash object -in
  const char *key  ///< the key of the record to retrieve -in
  )
{
  RECORD_T* recs;
  unsigned int off, ind, size;
  unsigned int code = strhash(key);
  
  recs = h->records;
  size = sizes[h->size_index];
  ind = code % size;
  off = 0;
  
  // search on hash which remains even if a record has been removed,
  // so remove_hash() does not need to move any collision records
  while (recs[ind].hash) {
    if ((code == recs[ind].hash) && recs[ind].key &&
        strcmp(key, recs[ind].key) == 0)
      return recs[ind].count;
    ind = (code + (int)pow(++off,2.0)) % size;
  }
  
  return -1;
}

/**
 * Remove key from table, returning value.
 *\returns the value of removing record, returns NULL if can't find key
 */
void* remove_hash(
  HASH_T* h, ///< working hash object -in
  char *key  ///< the key of the record to remove -in
  )
{
  unsigned int code = strhash(key);
  RECORD_T* recs;
  void * value;
  unsigned int off, ind, size;
  
  recs = h->records;
  size = sizes[h->size_index];
  ind = code % size;
  off = 0;
  
  while (recs[ind].hash) {
    if ((code == recs[ind].hash) && recs[ind].key &&
        strcmp(key, recs[ind].key) == 0) {
      // do not erase hash, so probes for collisions succeed
      value = recs[ind].value;
      // free key
      free(recs[ind].key);
      recs[ind].key = 0;
      recs[ind].value = 0;
      h->records_count--;
      return value;
    }
    ind = (code + (int)pow(++off, 2.0)) % size;
  }
  
  return NULL;
}

/**
 *\returns total number of keys in the hashtable.
 */  
unsigned int hash_size(
  HASH_T* h ///< working hash object -in
  )
{
  return h->records_count;
}


/**
 * hash_iterator routines!
 */

/**
 *\returns a new memory allocated hash iterator
 */
HASH_ITERATOR_T* new_hash_iterator(
  HASH_T* hash ///< the hash collection to iterate -out
  ){
  if (hash == NULL){
    carp(CARP_FATAL, "Null hash collection passed to hash iterator");
  }
  
  // allocate a new hash iterator
  HASH_ITERATOR_T* hash_iterator = 
    (HASH_ITERATOR_T*)mycalloc(1, sizeof(HASH_ITERATOR_T));
  
  // set items
  hash_iterator->hash = hash;
  hash_iterator->hash_idx = 0;
  hash_iterator->hash_total = sizes[hash->size_index];

  return hash_iterator;
}

/**
 * Does the hash_iterator have another hash object to return?
 * \returns true, if hash iterator has a next hash, else false
 */
bool hash_iterator_has_next(
  HASH_ITERATOR_T* hash_iterator ///< the working  hash iterator -in
  )
{
  HASH_T* hash = hash_iterator->hash;
  while (hash_iterator->hash_idx < hash_iterator->hash_total && 
         hash->records[hash_iterator->hash_idx].key == NULL){
    hash_iterator->hash_idx++;
  }
  return (hash_iterator->hash_idx < hash_iterator->hash_total);
}

/**
 * \returns the next the hash struct
 */
char* hash_iterator_next(
  HASH_ITERATOR_T* hash_iterator ///< the working hash iterator -in
  )
{
  HASH_T* hash = hash_iterator->hash;
  return hash->records[hash_iterator->hash_idx++].key;
}

/**
 * free the memory allocated iterator
 */
void free_hash_iterator(
  HASH_ITERATOR_T* hash_iterator ///< the hash iterator to free
  )
{
  if (hash_iterator != NULL){
    free(hash_iterator);
  }
}



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
