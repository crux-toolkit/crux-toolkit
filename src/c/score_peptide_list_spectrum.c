/*************************************************************************//**
 * \file score_list_spectrum
 * AUTHOR: Chris Park
 * CREATE DATE: 10/13 2006
 * DESCRIPTION: Object for given a list and a spectrum, generate a
 * perliminary scores with the given score function (ex, Sp) 
 *
 * REVISION: 
 ****************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <unistd.h>
#include "carp.h"
#include "parse_arguments.h"
#include "spectrum.h"
#include "spectrum_collection.h"
#include "ion.h"
#include "ion_series.h"
#include "crux-utils.h"
#include "objects.h"
#include "parameter.h"
#include "scorer.h"

/**
 * when wrong command is seen carp, and exit
 */
void wrong_command(char* arg, char* comment){

  if(comment == NULL){
    carp(CARP_FATAL, "incorrect argument: %s", arg);
  }
  else {
    // print comment if given
    carp(CARP_FATAL, "incorrect argument: %\n%s", arg, comment);
  }

}

int main(int argc, char** argv){

  // required variables
  char* ms2_file = NULL;
  int scan_num = 0;
  char* list_file = NULL;
  
  // optional variables
  char* type = "SP";
  char* parameter_file = "crux.params";
  
  // parsing variables
  int result = 0;
  const char * error_message;


 /* Define optional command line arguments */  
  parse_arguments_set_opt(
    "score type", 
    "The type of scoring function to use. sp",
    (void *) &type, 
    STRING_ARG);

  parse_arguments_set_opt(
    "parameter_file",
    "The crux parameter file to parse parameter from.",
    (void *) &parameter_file,
    STRING_ARG);

 /* Define required command line arguments */
  parse_arguments_set_req(
    "list file", 
    "The file that lists the peptide sequence and charge",
    (void *) &list_file, 
    STRING_ARG);

  parse_arguments_set_req(
    "scan number", 
    "The scan number for the MS-MS spectrum to extract from the ms2 file. This is an integer in the range [1, 100000], and uniquely identifies a particular MS-MS spectrum within an .ms2 file.",
    (void *) &scan_num, INT_ARG);

  parse_arguments_set_req(
    "ms2 file name", 
    "A file containing multiple MS-MS spectra in .ms2 format.",
    (void *) &ms2_file,
    STRING_ARG);

 /* Parse the command line */
 if (parse_arguments(argc, argv, 0)) {
   // parsed arguments
   char peptide_sequence[52] = ""; 
   int peptide_charge = 1;
   SCORER_TYPE_T score_type = SP; 

   SPECTRUM_T* spectrum = NULL;
   SPECTRUM_COLLECTION_T * collection = NULL;
   ION_SERIES_T* ion_series = NULL;
   SCORER_T* scorer = NULL;
   FLOAT_T score = 0;
   int  verbosity = CARP_INFO;

   // list file parsing
   FILE* file = NULL;
   char* new_line = NULL;
   int line_length;
   size_t buf_length = 0;
   
   // set verbosity
   set_verbosity_level(verbosity);

   // score type
   if(strcmp(type, "SP")== 0){
     score_type = SP;
   }
   else{
     wrong_command(type, "The type of scoring function to use. sp");
   }
   
   // parse paramter file
   parse_parameter_file(parameter_file);
   
   // set ion constraint to sequest settings
   ION_CONSTRAINT_T* ion_constraint = new_ion_constraint_sequest_sp();  
   
   // read ms2 file
   collection = new_spectrum_collection(ms2_file);
   spectrum = allocate_spectrum();
   
   // search for spectrum with correct scan number
   if(!get_spectrum_collection_spectrum(collection, scan_num, spectrum)){
     carp(CARP_FATAL, "failed to find spectrum with scan_num: %d", scan_num);
   }
   
   // create new scorer
   scorer = new_scorer(SP);  
   
   // open file and 
   file = fopen(list_file, "r");
   
   // check if succesfully opened file
   if(file == NULL){
     carp(CARP_FATAL, "failed to open list file");
   }
   
   // parse each line for peptides
   while((line_length =  getline(&new_line, &buf_length, file)) != -1){
     if(line_length >= 3 && new_line[0] != ' '){ 
       if(sscanf(new_line,"%s %d",// test format:peak line has more than 2 fields
                 peptide_sequence, &peptide_charge) > 0){
         
         // check peptide charge
         if(peptide_charge <1 || peptide_charge >3){
           free(new_line);
           wrong_command("peptide_charge", "The peptide charge. 1|2|3");
         }

         
         // check peptide sequence
         if(!valid_peptide_sequence(peptide_sequence)){           
           wrong_command(peptide_sequence, "not a valid peptide sequence");
         }

         // create new ion series
         ion_series = new_ion_series(peptide_sequence, peptide_charge, ion_constraint);
         
         // now predict ions
         predict_ions(ion_series);

         // calculates the Sp score
         score = score_spectrum_v_ion_series(scorer, spectrum, ion_series);
         
         // print the Sp score
         printf(/*"Sp score is:*/ "%.2f\n", score);
         
         free_ion_series(ion_series);
       }   
     }
   }
   
   fclose(file);
   free(new_line);
   
   // free heap
   free_scorer(scorer);
   free_ion_constraint(ion_constraint);
   free_spectrum_collection(collection);
   free_spectrum(spectrum);
 }
 else{
   char* usage = parse_arguments_get_usage("score_peptide_list_spectrum");
   result = parse_arguments_get_error(&error_message);
   carp(
    CARP_FATAL,
    "Error in command line. Error # %d\n%s\n%s", 
    result,
    error_message,
    usage
   );
 }
 exit(0);
}
