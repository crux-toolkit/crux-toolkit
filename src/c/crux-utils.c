#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "utils.h"

/**
 * returns a heap allocated copy of the src string
 */
char* my_copy_string(char* src){
  int length = strlen(src) +1; //+\0
  char* copy = 
    (char *)mymalloc(sizeof(char)*length);
  return strncpy(copy, src, length);  
}
