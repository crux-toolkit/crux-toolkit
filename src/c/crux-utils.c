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
#include "crux-utils.h"
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
 * parses the filename
 * ex) ../../file_name => returns filename
 *\returns A heap allocated array of filename
 */
char* parse_filename(char* file){
  int len = strlen(file);
  int end_idx = len;
  int end_path = -1;  //index of where the last "/" is located
  char* filename = NULL;
  
  for(; end_idx > 0; --end_idx){
    if(strncmp(&file[end_idx - 1], "/", 1) == 0){
      end_path = end_idx;
      break;
    }
  }
  
  //copy filename
  filename = copy_string_part(&file[end_idx], len); 
  
  return filename;
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
  else if(peptide_type == N_TRYPTIC){
    fprintf(file, "%s", "N_TRYPTIC");
  }
  else if(peptide_type == C_TRYPTIC){
    fprintf(file, "%s", "C_TRYPTIC");
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
 * check if the string has the correct suffix
 * \returns TRUE, if the string starts with the suffix, else FALSE
 */
BOOLEAN_T suffix_compare(
  char* string, ///< The string suffix to compare
  char* suffix  ///< The suffix to compare
  )
{
  int len = strlen(string);
  int len_suffix = strlen(suffix);

  if(len_suffix > len){
    return FALSE;
  }
  
  if(strncmp(string, suffix, len_suffix) == 0){
    return TRUE;
  }
  
  return FALSE;
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
 * Suffix may be NULL
 * format: suffix_myfasta_nameTag
 * \returns A heap allocated file name of the given fasta file
 */
char* generate_name(
  char* fasta_filename,
  char* name_tag,
  char* file_extension,
  char* suffix
  )
{
  int len = strlen(fasta_filename);
  int end_idx = len;
  int end_path = len;  //index of where the "." is located in the file
  char* name = NULL;
  char* after_suffix = NULL;
  int suffix_length = 0;
  char** file_n_path = NULL;
  int length = 0;

  //cut off the file extension if needed
  for(; end_idx > 0; --end_idx){
    if(strcmp(&fasta_filename[end_idx - 1], file_extension) == 0){
      end_path = end_idx - 1;
      break;
    }
  }
  
  //check suffix
  if(suffix != NULL){
    suffix_length = strlen(suffix);
    file_n_path = parse_filename_path(fasta_filename);
  }

  name = (char*)mycalloc(suffix_length + end_path + strlen(name_tag) + 1, sizeof(char));
  after_suffix = name;
  
  //if suffix exit add to top
  if(suffix_length != 0){
    length = strlen(file_n_path[1]);
    if(file_n_path[1] != NULL){
      strncpy(name, file_n_path[1], length);
      after_suffix = &name[length];
    }    
    strncpy(after_suffix, suffix, suffix_length);
    after_suffix = &after_suffix[suffix_length];
    
    length = strlen(file_n_path[0]);
    
    strncpy(after_suffix, file_n_path[0], (length - (len-end_path)));
    after_suffix = &after_suffix[(length - (len-end_path))];
    
    free(file_n_path[0]);
    free(file_n_path[1]);
    free(file_n_path);
  }
  else{
    strncpy(after_suffix, fasta_filename, end_path);
  }
  
  strcat(after_suffix, name_tag);
  return name;
}

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

/**
 * Open and create a file handle of a file that is named 
 * and located in user specified location
 * Assumes the directory exists
 *\returns a file handle of a file that is named and located in user specified location
 */
FILE* create_file_in_path(
  char* filename,  ///< the filename to create & open -in
  char* directory  ///< the directory to open the file in -in
  )
{
  char* file_full_path = get_full_filename(directory, filename);
  FILE* file = fopen(file_full_path, "a+");
  
  if(file == NULL){
    carp(CARP_ERROR, "failed to create and open file: %s", file_full_path);
  }
  
  free(file_full_path);
  
  return file;
}

/**
 *\returns a heap allocated feature name array for the algorithm type
 */
char** generate_feature_name_array(
  ALGORITHM_TYPE_T algorithm ///< the algorithm's feature name to produce -in
)
{
  char** name_array = NULL;

  switch(algorithm){
  case PERCOLATOR:
  case CZAR:
  case ALL:
    name_array = (char**)mycalloc(20, sizeof(20));
    name_array[0] =  my_copy_string("XCorr");
    name_array[1] =  my_copy_string("DeltCN");
    name_array[2] =  my_copy_string("DeltLCN");
    name_array[3] =  my_copy_string("Sp");
    name_array[4] =  my_copy_string("lnrSp");
    name_array[5] =  my_copy_string("dM");
    name_array[6] =  my_copy_string("absdM");
    name_array[7] =  my_copy_string("Mass");
    name_array[8] =  my_copy_string("ionFrac");
    name_array[9] =  my_copy_string("lnSM");
    name_array[10] =  my_copy_string("enzN");
    name_array[11] =  my_copy_string("enzC");
    name_array[12] =  my_copy_string("enzInt");
    name_array[13] =  my_copy_string("pepLen");
    name_array[14] =  my_copy_string("charge1");
    name_array[15] =  my_copy_string("charge2");
    name_array[16] =  my_copy_string("charge3");
    name_array[17] =  my_copy_string("numPep");
    name_array[18] =  my_copy_string("numProt");
    name_array[19] =  my_copy_string("pepSite");
  }
  
  return name_array;
}

/**
 * User define our upper and our lower bounds.
 * The random number will always be 
 * between low and high, inclusive.
 * There is no seeding in this function, user must do it for themselves
 *\returns a random number between the interval user provides
 */
int get_random_number_interval(
  int low, ///< the number for lower bound -in
  int high ///< the number for higher bound -in
  )
{  
  return (rand() % (high - low + 1) + low);
}

/**
 *\returns the number of digits in the number
 */
int get_number_digits(
  int number ///< the number to count digits
  )
{
  int idx = 0;
  for(; number >= 10; ++idx){
    number = number/10;    
  }

  return ++idx;
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

