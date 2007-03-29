#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>
#include <time.h>
#include "carp.h"
#include "utils.h"
#include "objects.h"


/**
 * PRECISION, determines the precision of the compare float, users
 * should lower the number if need more precision
 */
#define PRECISION 0.000000005 

/**
 * the maximum error in terms of Units in the Last Place. 
 * This specifies how big an error we are willing to accept in terms of the value of the least significant 
 * digit of the floating point numbers representation. 
 * MAX_ULPS can also be interpreted in terms of how many representable floats 
 * we are willing to accept between A and B. This function will allow MAX_ULPS-1 floats between A and B.
 */
#define MAX_ULPS 2

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
  float EPSILON = PRECISION;
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
 *\returns TRUE if float_a is between the interaval of min and max, else FALSE
 */
inline BOOLEAN_T compare_float_three(float float_a, float min, float max){
  if(compare_float(float_a, min) == -1 ||
     compare_float(float_a, max) ==  1){
    return FALSE;
  }
  return TRUE;
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
    if(strncmp(&file[end_idx - 1], "/", 1) == 0){
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
 * convert the integer into a string
 * \returns a heap allocated string
 */
char* signed_int_to_char(int i){
  int digit = abs(i)/ 10;
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

/**
 * given the path and the filename return a file with path
 * "path/filename"
 * \returns a heap allocated string, "path/filename"
 */
char* get_full_filename(char* path, char* filename){
  char* ready_path = cat_string(path, "/");
  char* result = cat_string(ready_path, filename);
  free(ready_path);
  return result;
}


/**
 * returns the file size of the given filename
 */
long get_filesize(char *FileName){
    struct stat file;
    //return file size
    if(!stat(FileName,&file)){
      return file.st_size;
    }
    return 0;
}

/**
 * deletes a given directory and it's files inside.
 * assumes that there's no sub directories, only files
 * \returns TRUE if successfully deleted directory
 */
BOOLEAN_T delete_dir(char* dir) {
  struct dirent **namelist =NULL;
  int num_file =0;
  int result;

  //does the directory to remove exist?, if so move into it..
  if(chdir(dir) == -1){
    return FALSE;
  }

  //collect all files in dir
  num_file = scandir(".", &namelist, 0, alphasort);

  //delete all files in temp dir
  while(num_file--){
    remove(namelist[num_file]->d_name);
    free(namelist[num_file]);
  }
  free(namelist);

  chdir("..");
  result = rmdir(dir);
  if(result == FALSE){
    return FALSE;
  }
  
  return TRUE;
}

/**
 * given a fasta_file name it returns a name with the name_tag add to the end
 * format: myfasta_nameTag
 * \returns A heap allocated file name of the given fasta file
 */
char* generate_name(
  char* fasta_filename,
  char* name_tag
  )
{
  int len = strlen(fasta_filename);
  int end_idx = len;
  int end_path = len;  //index of where the "." is located in the file
  char* name = NULL;
  
  //cut off the ".fasta" if needed
  for(; end_idx > 0; --end_idx){
    if(strcmp(&fasta_filename[end_idx - 1], ".fasta") == 0){
      end_path = end_idx - 1;
      break;
    }
  }
  
  name = (char*)mycalloc(end_path + strlen(name_tag) + 1, sizeof(char));
  strncpy(name, fasta_filename, end_path);
  strcat(name, name_tag);
  return name;
}

/***
 *
 */



/**
 * checks if each AA is an AA
 *\returns TRUE if sequence is valid else, FALSE
 */
BOOLEAN_T valid_peptide_sequence( char* sequence){
  //iterate over all AA and check if with in range
  while(sequence[0] != '\0'){
    if(sequence[0] < 65 || sequence[0] > 90 ){
      return FALSE;
    }
    ++sequence;
  }
  return TRUE;
}


void swap_quick(
  float* a,
  int idx,
  int jdx
  )
{
  float temp = 0;
  temp = a[idx];
  a[idx] = a[jdx];
  a[jdx] = temp;
}
 
int Random(int i, int j) {
  return i + rand() % (j-i+1);
}

void quick_sort(float a[], int left, int right) {
  int last = left, i;

  if (left >= right) return;
  
  swap_quick(a,left,Random(left,right));
  for (i = left + 1; i <= right; i++)
    if (a[i] > a[left]) ///CHECK THIS!!
      swap_quick(a,++last,i);
  swap_quick(a,left,last);
  quick_sort(a,left,last-1);
  quick_sort(a,last+1,right);
}

void quicksort(float a[], int array_size){
  quick_sort(a, 0, array_size-1);
}

