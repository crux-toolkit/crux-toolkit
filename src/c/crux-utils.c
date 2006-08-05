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
  int length = strlen(src) +1; // +\0
  char* copy = 
    (char *)mymalloc(sizeof(char)*length);
  return strncpy(copy, src, length);  
}

/**
 * returns copy of the src string upto the specified length
 * includes a null terminating \\0 character
 * the string is heap allocated thus, user must free
 */
char* copy_string_part(char* src, int length){
  char* copy = (char*)mycalloc(length+1, sizeof(char));
  strncpy(copy, src, length);
  copy[length] = '\0';
  return copy;
}

/**
 * \returns the 0 if equal, 1 if float_a is larger, -1 if float_b is larger
 * compare the absolute value of the difference of two numbers with an 
 * appropriate epsilon to get relations.
 * Multiplying the epsilon by sum of the comparands adjusts the comparison 
 * to the range of the numbers, allowing a single epsilon to be used for many, 
 * or perhaps all compares.
 */
inline int compare_float(float float_a, float float_b){
  float EPSILON = 0.0000005;
  float sum = float_a + float_b;
  //a == b
  if( fabsf(float_a - float_b) <= fabsf(sum)* EPSILON ){
    return 0;
  }
  //a > b
  else if((float_a - float_b) > fabsf(sum)* EPSILON){
    return 1;
  }
  // a < b
  else{
    return -1;
  }
}

/**
 * parses the file path  
 * returns NULL if only a filename was passed in
 * ex) ../../file_name => returns ../../
 *     file_name => returns NULL
 *\returns A heap allocated path to the location of the file
 */
char* parse_file_path(char* file){
  int len = strlen(file);
  int end_idx = len;
  int end_path = -1;  //index of where the last "/" is located
  char* path = NULL;

  for(; end_idx > 0; --end_idx){
    if(strcmp(file[end_idx - 1], "/") == 0){
      end_path = end_idx;
      break;
    }
  }
  //there is a "/" in the file
  if(end_path != -1){
    path = copy_string_part(file, end_path);
  }
  return path;
}
