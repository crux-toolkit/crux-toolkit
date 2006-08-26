#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "utils.h"
#include "objects.h"

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
 * parses the filename and path  
 * returns an array A, with A[0] the filename and A[1] the path to the filename
 * returns A[1] NULL if only a filename was passed in
 * ex) ../../file_name => returns filename , ../../
 *     file_name => returns filename, NULL
 *\returns A heap allocated array of both filename and path
 */
char** parse_filename_path(char* file){
  int len = strlen(file);
  int end_idx = len;
  int end_path = -1;  //index of where the last "/" is located
  char* path = NULL;
  char* filename = NULL;
  char** result = (char**)mycalloc(2, sizeof(char*));

  for(; end_idx > 0; --end_idx){
    if(strcmp(&file[end_idx - 1], "/") == 0){
      end_path = end_idx;
      break;
    }
  }
  //copy path, if there is a "/" in the file
  if(end_path != -1){
    path = copy_string_part(file, end_path);
  }
  //copy filename
  filename = copy_string_part(&file[end_idx], len); 
  
  //set result with filename and path
  result[0] = filename;
  result[1] = path;
  
  return result;
}

/**
 * convert the integer into a string
 * \returns a heap allocated string
 */
char* int_to_char(unsigned int i){
  unsigned int digit = i / 10;
  char* int_string = (char*)mycalloc(digit+2, sizeof(char));
  sprintf(int_string, "%d", i);
  return int_string;
}
 

/**
 *prints the peptide type given it's enum value
 */
void print_peptide_type(PEPTIDE_TYPE_T peptide_type, FILE* file){
  if(peptide_type == TRYPTIC){
    fprintf(file, "%s", "TRYPTIC");
  }
  else if(peptide_type == PARTIALLY_TRYPTIC){
    fprintf(file, "%s", "PARTIALLY_TRYPTIC");
  }
  else if(peptide_type == NOT_TRYPTIC){
    fprintf(file, "%s", "NOT_TRYPTIC");
  }
  else if(peptide_type == ANY_TRYPTIC){
    fprintf(file, "%s", "ANY_TRYPTIC");
  }
}

/**
 * given two strings return a concatenated third string
 * \returns a heap allocated string that concatenates the two inputs
 */
char* cat_string(char* string_one, char* string_two){
  int len_one = strlen(string_one);
  int len_two = strlen(string_two);
  
  char* result = (char*)mycalloc(len_one + len_two + 1, sizeof(char));
  strncpy(result, string_one, len_one);
  strncpy(&result[len_one], string_two, len_two);
  return result;
}
