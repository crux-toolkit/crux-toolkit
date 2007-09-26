/*****************************************************************************
 * \file predict_peptide_ion
 * AUTHOR: Chris Park
 * CREATE DATE: 10/05 2006
 * DESCRIPTION: Object for given a peptide and a charge state, predict the ions
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
#include "ion.h"
#include "ion_series.h"
#include "crux-utils.h"
#include "objects.h"

/**
 * when wrong command is seen carp, and exit
 */
void wrong_command(char* arg, char* comment){
  char* usage = parse_arguments_get_usage("predict_peptide_ions");
  carp(CARP_FATAL, "incorrect argument: %s", arg);

  // print comment if given
  if(comment != NULL){
    carp(CARP_FATAL, "%s", comment);
  }

  // FIXME uncomment this print if want to print usage whenever error message is printed
  // fprintf(stderr, "%s", usage);
  free(usage);
  exit(1);
}

int main(int argc, char** argv){

  char* primary_ions = "by";
  char* precursor_ions = "F";
  char* neutral_losses = "all";
  int isotope_count = 0;
  char* flanking = "F";
  char* max_ion_charge = "peptide";
  int nh3_count = 0;
  int h2o_count = 0;
  char* peptide_sequence = "";
  int charge_state = 0;
  
  // parsing variables
  int result = 0;
  char * error_message;

 /* Define optional command line arguments */ 
 parse_arguments_set_opt(
    "primary_ions", 
    "Predict the specified primary ion series. 'b' indicates b-ions only, 'y' indicates y-ions only, 'by' indicates both. b|y|by", 
    (void *) &primary_ions, 
    STRING_ARG);

 parse_arguments_set_opt(
    "precursor_ions", 
    "Predict the precursor ions, and all associated ions (neutral-losses, multiple charge states) consistent with the other specified options. T|F",
    (void *) &precursor_ions, 
    STRING_ARG);

 parse_arguments_set_opt(
    "neutral-losses", 
    "Predict the the specified neutral loss peaks ('none' indicates no neutral losses, 'h2o' indicates water loss only, 'nh3' indicates ammonia losses only, 'all' indicates all kinds of losses).none|h2o|nh3|all",
    (void *) &neutral_losses, 
    STRING_ARG);

 parse_arguments_set_opt(
    "isotope",                          
    "Predict the given number of isotope peaks.0|1|2",
    (void *) &isotope_count,
    INT_ARG);

 parse_arguments_set_opt(
    "flanking",                          
    "Predict flanking peaks for b- and y-ions.T|F",
    (void *) &flanking, 
    STRING_ARG);

 parse_arguments_set_opt(
    "max-ion-charge",                          
    "Predict ions up to this charge-state. The integer options ('1','2' and '3') specify a fixed maximum charge-state. The 'peptide' option indicates that the ions should range up to the maximum charge-state of the peptide itself (thus, a +2 charge state peptide would have ions of +1 and +2).1|2|3|peptide",
    (void *) &max_ion_charge, 
    STRING_ARG);

 parse_arguments_set_opt(
    "nh3",                          
    "Predict peaks with the following max nh3 modification.",
    (void *) &nh3_count,
    INT_ARG);
 parse_arguments_set_opt(
    "h2o",                          
    "Predict peaks with the following max h2o modification.",
    (void *) &h2o_count,
    INT_ARG);


 /* Define required command line arguments */
 parse_arguments_set_req(
    "peptide sequence", 
    "The literal peptide sequence (e.g. EAMAPK) that is used to predict the ions.",
    (void *) &peptide_sequence, STRING_ARG);

 parse_arguments_set_req(
    "charge state", 
    "The charge state of the peptide used to predict the ions.",
    (void *) &charge_state, INT_ARG);
  

 /* Parse the command line */
 if (parse_arguments(argc, argv, 0)) {
   ION_TYPE_T ion_type = BY_ION;
   BOOLEAN_T use_precursor_ions = FALSE;
   int neutral_loss_count[MAX_MODIFICATIONS];
   BOOLEAN_T is_modification = TRUE;
   BOOLEAN_T is_modification_h2o = TRUE;
   BOOLEAN_T is_modification_nh3 = TRUE;
   int max_charge = charge_state;

   // check if charge >= 0
   if(charge_state < 0){
     wrong_command("peptide charge cannot be bellow 0", NULL);
   }

   // check peptide sequence
   if(!valid_peptide_sequence(peptide_sequence)){
     wrong_command(peptide_sequence, "not a valid peptide sequence");
   }

   // primary_ions
   if(strcmp(primary_ions, "b") == 0){
     ion_type = B_ION;
   }
   else if(strcmp(primary_ions, "y")== 0){
     ion_type = Y_ION;
   }
   else if(strcmp(primary_ions, "by")== 0){
     ion_type = BY_ION;
   }
   else{
     wrong_command(primary_ions, "primary_ions are b|y|by");
   }

   // precursor_ions
   if(strcmp(precursor_ions, "F")== 0){
     use_precursor_ions = FALSE;
   }
   else if(strcmp(precursor_ions, "T")== 0){
     use_precursor_ions = TRUE;
   }
   else{
     wrong_command(precursor_ions, "must pick between T|F");
   }

   // neutral_losses
   // initialize
   int modification_idx = 0;
   for(; modification_idx < MAX_MODIFICATIONS; ++modification_idx){
     neutral_loss_count[modification_idx] = 0;
   }

   if(strcmp(neutral_losses, "none")== 0){
     is_modification = FALSE;
     is_modification_h2o = FALSE;
     is_modification_nh3 = FALSE;
   }
   else if(strcmp(neutral_losses, "h2o")== 0){
     is_modification = TRUE;
     is_modification_h2o = TRUE;
     is_modification_nh3 = FALSE;
   }
   else if(strcmp(neutral_losses, "nh3")== 0){
     is_modification = TRUE;
     is_modification_h2o = FALSE;
     is_modification_nh3 = TRUE;
   }
   else if(strcmp(neutral_losses, "all")!= 0){
     wrong_command(neutral_losses, "neutral loss can be none|h2o|nh3|all");
   }
   
   // isotope
   if(isotope_count < 0 || isotope_count > 3){
     // out of bounds...error
     wrong_command("isotope not an integer between 0 and 3", NULL);
   } 
   else{
     neutral_loss_count[ISOTOPE] = isotope_count;
   }

   // flanking
   if(strcmp(flanking, "T")== 0){
     is_modification = TRUE;
     neutral_loss_count[FLANK] = 1;
   }
   else if(strcmp(flanking, "F")== 0){
     neutral_loss_count[FLANK] = 0;
   }
   else{
     wrong_command(flanking, "must pick between T|F");
   }

   // max_ion_charge
   if(strcmp(max_ion_charge, "1")== 0){
     max_charge = 1;
   }
   else if(strcmp(max_ion_charge, "2")== 0){
     max_charge = 2;
   }
   else if(strcmp(max_ion_charge, "3")== 0){
     max_charge = 3;
   }
   else if(strcmp(max_ion_charge, "peptide")== 0){
     max_charge = charge_state;
   }
   else{
     wrong_command(max_ion_charge, "max_ion_charge must be 1|2|3|peptide");
   }
   
   // set nh3, h2o
   if(is_modification_nh3){
     neutral_loss_count[NH3] = nh3_count;
   }
   if(is_modification_h2o){
     neutral_loss_count[H2O] = h2o_count;
   }

   // creat ion_constraint
   ION_CONSTRAINT_T* ion_constraint = 
     new_ion_constraint(MONO, max_charge, ion_type, use_precursor_ions);

   
   // set ion_constraint3 modification counts, if modifications should occur
   if(is_modification){
     set_ion_constraint_modification( ion_constraint, NH3, neutral_loss_count[NH3]);
     set_ion_constraint_modification( ion_constraint, H2O, neutral_loss_count[H2O]);
     set_ion_constraint_modification( ion_constraint, ISOTOPE, neutral_loss_count[ISOTOPE]);
     set_ion_constraint_modification( ion_constraint, FLANK, neutral_loss_count[FLANK]);
   }

   // create ion_series
   ION_SERIES_T* ion_series = new_ion_series(peptide_sequence, charge_state, ion_constraint);
   
   // now predict ions
   predict_ions(ion_series);
   
   // print settings
   printf("# PEPTIDE: %s\n",peptide_sequence);
   printf("# CHARGE: %d\n", charge_state);
   printf("# MAX-ION-CHRAGE: %s\n", max_ion_charge);
   printf("# NH3 modification: %d\n", neutral_loss_count[NH3]);
   printf("# H2O modification: %d\n", neutral_loss_count[H2O] );
   printf("# ISOTOPE modification: %d\n", neutral_loss_count[ISOTOPE] );
   printf("# FLANK modification: %d\n", neutral_loss_count[FLANK]);

   // print ions
   print_ion_series(ion_series, stdout);


   // free
   free_ion_constraint(ion_constraint);
   free_ion_series(ion_series);
 }
 else {
   char* usage = parse_arguments_get_usage("predict_peptide_ions");
   result = parse_arguments_get_error(&error_message);
   fprintf(stderr, "Error in command line. Error # %d\n", result);
   fprintf(stderr, "%s\n", error_message);
   fprintf(stderr, "%s", usage);
   free(usage);
 }
 exit(0);
}
