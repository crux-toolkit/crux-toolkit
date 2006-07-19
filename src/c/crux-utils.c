#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include "utils.h"

/**
 * returns a heap allocated copy of the src string
 */
char* copy_string(char* src){
  int length = strlen(src) +1; //+\0
  char* copy = 
    (char *)mymalloc(sizeof(char)*length);
  return strncpy(copy, src, length);  
}
