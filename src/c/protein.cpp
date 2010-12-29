/*************************************************************************//**
 * \file protein.cpp
 * \brief Object for representing a single protein.
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <vector>
#include "utils.h"
#include "parameter.h"
#include "objects.h"
#include "peptide.h"
#include "protein.h"
#include "peptide_src.h"
#include "database.h"
#include "carp.h"
#include "peptide_constraint.h"

using namespace std;

/**
 * Constants
 */
static const int PROTEIN_ID_LENGTH = 100;
static const int PROTEIN_SEQUENCE_LENGTH = 40000;
static const int PROTEIN_ANNOTATION_LENGTH = 100;
static const int LONGEST_LINE  = PROTEIN_ID_LENGTH + PROTEIN_ID_LENGTH;
static const int FASTA_LINE = 50;
static const int SMALLEST_MASS = 57;
static const int LARGEST_MASS = 190;

/**
 * \struct protein 
 * \brief A protein sequence.
 */
struct protein{
  DATABASE_T*  database; ///< Which database is this protein part of
  unsigned long int offset; ///< The file location in the database source file
  unsigned int protein_idx; ///< The index of the protein in it's database.
  BOOLEAN_T    is_light; ///< is the protein a light protein?
  BOOLEAN_T    is_memmap; ///< is the protein produced from memory mapped file
  char*              id; ///< The protein sequence id.
  char*        sequence; ///< The protein sequence.
  unsigned int   length; ///< The length of the protein sequence.
  char*      annotation; ///< Optional protein annotation. 
};    

/**
 * \struct protein_peptide_iterator
 * \brief Object to iterate over the peptides within a protein in an
 * unspecified order. The peptides should satisfy the constraints specified
 * in the peptide_constraint object.
 */
struct protein_peptide_iterator {
  PROTEIN_T* protein; ///< The protein whose peptides to iterate over. 
  unsigned short int cur_start; ///< Start in protein of the current peptide.
  unsigned short int cur_length; ///< The length of the current peptide.
  unsigned int peptide_idx; ///< The index of the current peptide.
  PEPTIDE_CONSTRAINT_T* peptide_constraint; ///< peptide type to iterate over.
  double* mass_array; ///< stores all the peptides' masses
  vector<int>* nterm_cleavage_positions; ///< nterm cleavages that satisfy 
                                        ///< constraint. 1st aa is 1.
  vector<int>* peptide_lengths; ///< all the lengths of valid peptides
  vector<FLOAT_T>* peptide_masses; ///< all the masses of valid peptides
  vector<int>* cumulative_cleavages; ///< cumulative number of cleavages so far
  int current_cleavage_idx; /// where are we in the cleavage positions?
  int num_cleavages; /// how many cleavage positions?

  BOOLEAN_T has_next; ///< is there a next? 
  int num_mis_cleavage; ///< The maximum mis cleavage of the peptide
};

// def below
static BOOLEAN_T read_title_line
  (FILE* fasta_file,
   char* name,
   char* description,
   PROTEIN_T* protein);

// def below
static BOOLEAN_T read_raw_sequence
  (FILE* fasta_file,   // Input Fasta file.
   char* name,         // Sequence ID (used in error messages).
   unsigned int   max_chars,    // Maximum number of allowed characters.
   char* raw_sequence, // Pre-allocated sequence.
   unsigned int* sequence_length // the sequence length -chris added
   );

/**
 * \returns An (empty) protein object.
 */
PROTEIN_T* allocate_protein(void){
  PROTEIN_T* protein = (PROTEIN_T*)mycalloc(1, sizeof(PROTEIN_T));
  return protein;
}

/**
 * \returns A new protein object(heavy).
 * The protein is does not constain a database, users must provide one.
 */
PROTEIN_T* new_protein(
  const char*         id, ///< The protein sequence id. -in
  const char*   sequence, ///< The protein sequence. -in
  unsigned int length, ///< The length of the protein sequence. -in
  const char* annotation,  ///< Optional protein annotation.  -in
  unsigned long int offset, ///< The file location in the source file in the database -in
  unsigned int protein_idx, ///< The index of the protein in it's database.-in  
  DATABASE_T* database ///< the database of its origin
  )
{
  PROTEIN_T* protein = allocate_protein();
  set_protein_id(protein, id);
  set_protein_sequence(protein, sequence);
  set_protein_length(protein, length);
  set_protein_annotation(protein, annotation);
  set_protein_offset(protein, offset);
  set_protein_protein_idx(protein, protein_idx);
  set_protein_is_light(protein, FALSE);
  protein->database = copy_database_ptr(database); 
  protein->is_memmap = FALSE;
  return protein;
}         

/**
 * \returns A new light protein object.
 */
PROTEIN_T* new_light_protein(
  unsigned long int offset, ///< The file location in the source file in the database -in
  unsigned int protein_idx ///< The index of the protein in it's database. -in
  )
{
  PROTEIN_T* protein = allocate_protein();
  set_protein_is_light(protein, TRUE);
  set_protein_offset(protein, offset);
  set_protein_protein_idx(protein, protein_idx);
  return protein;
}


/**
 * convert light protein to heavy, by parsing all the sequence from fasta file
 * \returns TRUE if successfully converts the protein to heavy 
 */
BOOLEAN_T protein_to_heavy(
  PROTEIN_T* protein ///< protein to convert to heavy -in 
  )
{
  // protein already heavy
  if(!protein->is_light){
    return TRUE;
  }
  
  FILE* file = get_database_file(protein->database);
  
  // rewind to the begining of the protein to include ">" line
  fseek(file, protein->offset, SEEK_SET);

  // failed to parse the protein from fasta file
  // protein offset is set in the parse_protein_fasta_file method
  if(!parse_protein_fasta_file(protein ,file)){
    carp(CARP_ERROR, 
         "failed convert protein to heavy, cannot parse fasta file");
    return FALSE;
  }
      
  protein->is_light = FALSE;
  
  return TRUE;
}                            

/**
 * covert heavy protein back to light
 * \returns TRUE if successfully converts the protein to light
 */
BOOLEAN_T protein_to_light(
  PROTEIN_T* protein ///< protein to convert back to light -in 
  )
{
  // protein already light
  if(protein->is_light){
    return TRUE;
  }
  // free all char* in protein object
  free(protein->sequence);
  protein->sequence = NULL;
  free(protein->annotation);
  protein->annotation = NULL;
  free(protein->id);
  protein->id = NULL;
  
  return (protein->is_light = TRUE);
}                            

/**
 * Frees an allocated protein object.
 */
void free_protein(
  PROTEIN_T* protein ///< object to free -in
  )
{
  // FIXME what is the point of checking this?
  if(!protein->is_memmap && !protein->is_light){ 
    if (protein->id != NULL){
      free(protein->id);
    }
    if (protein->sequence != NULL){
      free(protein->sequence);
    }
    if (protein->annotation != NULL){
      free(protein->annotation);
    }
  }
  free(protein);
}

/**
 * Prints a protein object to file.
 * if light protein coverts it to heavy protein
 */
void print_protein(
  PROTEIN_T* protein, ///< protein to print -in
  FILE* file ///< output stream -out
  )
{
  // covnert to heavy protein
  if(protein->is_light){
    protein_to_heavy(protein);
  }
  int   sequence_index;
  int   sequence_length = get_protein_length(protein);
  char* sequence = get_protein_sequence(protein);
  char* id = get_protein_id(protein);
  char* annotation = get_protein_annotation(protein);
  
  fprintf(file, ">%s %s\n", id, annotation);

  sequence_index = 0;
  while (sequence_length - sequence_index > FASTA_LINE) {
    fprintf(file, "%.*s\n", FASTA_LINE, &(sequence[sequence_index]));
    sequence_index += FASTA_LINE;
  }
  fprintf(file, "%s\n\n", &(sequence[sequence_index]));

  free(sequence);
  free(id);
  free(annotation);
}


/**
 * prints a binary representation of the protein
 * 
 * FORMAT
 * <int: id length><char: id><int: annotation length><char: annotation><int: sequence length><char: sequence>
 *
 * make sure when rading the binary data, add one to the length so that it will read in the terminating char as well
 */
void serialize_protein(
  PROTEIN_T* protein, ///< protein to print as binary -in
  FILE* file ///< output stream -out
  )
{
  // covnert to heavy protein
  if(protein->is_light){
    protein_to_heavy(protein);
  }
  
  int id_length = strlen(protein->id);
  int annotation_length = strlen(protein->annotation);

  // write the protein id length
  fwrite(&id_length, sizeof(int), 1, file);
  
  // write the protein id 
 // include "/0"
  fwrite(protein->id, sizeof(char), id_length+1, file);

  // write the protein annotation length
  fwrite(&annotation_length, sizeof(int), 1, file);

  // write the protein annotation
  // include "/0"
  fwrite(protein->annotation, sizeof(char), annotation_length+1, file);
  
  // write the protein sequence length
  fwrite(&protein->length, sizeof(unsigned int), 1, file);
  
  // write the protein sequence
  // include "/0"
  fwrite(protein->sequence, sizeof(char), protein->length+1, file);
}


/**
 * Copies protein object src to dest.
 * assumes that the protein is heavy
 * dest must be a heap allocated object 
 */
void copy_protein(
  PROTEIN_T* src,///< protein to copy -in
  PROTEIN_T* dest ///< protein to copy to -out
  )
{
  char* id = get_protein_id(src);
  char* sequence = get_protein_sequence(src);
  char* annotation = get_protein_annotation(src);
  
  set_protein_id(dest, id);
  set_protein_sequence(dest, sequence);
  set_protein_length(dest, get_protein_length(src));
  set_protein_annotation(dest, annotation);
  set_protein_offset(dest, src->offset);
  set_protein_protein_idx(dest, src->protein_idx);
  set_protein_is_light(dest, src->is_light);
  dest->database = src->database;
  
  free(id);
  free(sequence);
  free(annotation);
}


/**
 * Parses a protein from an memory mapped binary fasta file
 * the protein_idx field of the protein must be added before or after
 * you parse the protein 
 * protein must be a heap allocated
 * 
 * Assume memmap pointer is set at beginning of protein
 * Assume protein binary format
 * <int: id length><char: id><int: annotation length><char: annotation><int: sequence length><char: sequence>
 *
 * modifies the *memmap pointer!
 * \returns TRUE if success. FALSE is failure.
 */
BOOLEAN_T parse_protein_binary_memmap(
  PROTEIN_T* protein, ///< protein object to fill in -out
  char** memmap ///< a pointer to a pointer to the memory mapped binary fasta file -in
  )
{
  int id_length = 0;
  int annotation_length = 0;
  int sequence_length = 0;

  /* FIXME, maybe use this to check if still within file
  if(*memmap_as_char[0] == EOF){
    carp(CARP_ERROR, "end of file");
  }
  */

  /***set protein ID***/

  // read id length
  id_length = *((int *) *memmap);

  // reset pointer to start of id
  *memmap += sizeof(int);

  // set protein id to mem mapped id
  protein->id = *memmap;

  // reset pointer to move to annotation_length
  *memmap += (id_length + 1);


  /***set protein annotation***/

  // read annotation length
  annotation_length = *((int *) *memmap);

  // reset pointer to start of annotation
  *memmap += sizeof(int);

  // set protein annotation to mem mapped annotation
  protein->annotation = *memmap;

  // reset pointer to move to sequence_length
  *memmap += (annotation_length + 1);


  /***set protein sequence***/
  
  // read sequence length
  sequence_length = *((int *) *memmap);
  protein->length = sequence_length;

  // reset pointer to start of sequence
  *memmap += sizeof(int);

  // set protein annotation to mem mapped sequence
  protein->sequence = *memmap;

  // reset pointer to move to start of next protein
  *memmap += (sequence_length + 1);
  
  // now this protein has been created from memory mapped!
  protein->is_memmap = TRUE;

  return TRUE;
}

// FIXME ID line and annotation might need to be fixed
VERBOSE_T verbosity = NORMAL_VERBOSE;
/**
 * Parses a protein from an open (FASTA) file.
 * the protein_idx field of the protein must be added before or after
 * you parse the protein  
 * \returns TRUE if success. FALSE is failure.
 * protein must be a heap allocated
 */
BOOLEAN_T parse_protein_fasta_file(
  PROTEIN_T* protein, ///< protein object to fill in -out
  FILE* file ///< fasta file -in
  )
{
  static char name[LONGEST_LINE];    ///< Just the sequence ID.
  static char desc[LONGEST_LINE];    ///< Just the comment field.
  static char buffer[PROTEIN_SEQUENCE_LENGTH];///> The sequence to read in.
  static unsigned int sequence_length; // the sequence length

  // Read the title line.
  if (!read_title_line(file, name, desc, protein)) {
    return(FALSE);
  }
  
  buffer[0] = '\0';

  // Read the sequence.
  if (!read_raw_sequence(file, name, PROTEIN_SEQUENCE_LENGTH, buffer, &sequence_length)) {
    carp(CARP_FATAL, "Sequence %s is too long.\n", name);
    exit(1);
  }
    
  // update the protein object.
  set_protein_length(protein, sequence_length);
  set_protein_id(protein, name);
  set_protein_sequence(protein, buffer);
  set_protein_annotation(protein, desc);

  return(TRUE);

}




/**************************************************/

/**
 * FASTA file parsing code
 * AUTHOR: William Stafford Noble
 * modified by Chris Park
 */

/**
 * Find the beginning of the next sequence, and read the sequence ID
 * and the comment.
 */
static BOOLEAN_T read_title_line
  (FILE* fasta_file,
   char* name,
   char* description,
   PROTEIN_T* protein)
{
  static char id_line[LONGEST_LINE];  // Line containing the ID and comment.
  int a_char;                         // The most recently read character.

  // Read until the first occurrence of ">".
  while ((a_char = getc(fasta_file)) != '>') {
    // If we hit the end of the file, return FALSE.
    if (a_char == EOF) {
      return(FALSE);
    }  
  }
  // set protein offset                   FIXME: might not need to "-1" -CHRIS
  protein->offset = ftell(fasta_file) - 1;

  /**
   * chris edited, added this block to make sure all of comment line
   * is read although might not be stored, to ensure the file* is at
   * start of the sequence
   */

  {
    char* new_line = NULL;
    int line_length;
    size_t buf_length = 0;

    if((line_length =  getline(&new_line, &buf_length, fasta_file)) == -1){
      carp(CARP_FATAL, "Error reading Fasta file.\n");
    }
    strncpy(id_line, new_line, LONGEST_LINE-1);
    free(new_line);
  }

  // Remove EOL.
  id_line[strlen(id_line) - 1] = '\0';

  // Extract the ID from the beginning of the line.
  if (sscanf(id_line, "%s", name) != 1) {
    carp(CARP_FATAL, "Error reading sequence ID.\n%s\n", id_line);
  }

  // Store the rest of the line as the comment.
  strcpy(description, &(id_line[strlen(name)+1]));

  return(TRUE);
}


/****************************************************************************
 * Read raw sequence until a '>' is encountered or too many letters
 * are read.  The new sequence is appended to the end of the given
 * sequence.
 *
 * Return: Was the sequence read completely?
 ****************************************************************************/
static BOOLEAN_T read_raw_sequence
  (FILE* fasta_file,   // Input Fasta file.
   char* name,         // Sequence ID (used in error messages).
   unsigned int   max_chars,    // Maximum number of allowed characters.
   char* raw_sequence, // Pre-allocated sequence.
   unsigned int* sequence_length // the sequence length -chris added
   )
{
  int a_char;
  unsigned int i_seq;
  BOOLEAN_T return_value = TRUE;

  // Start at the end of the given sequence.
  i_seq = strlen(raw_sequence);
  assert((unsigned int)strlen(raw_sequence) < max_chars);

  // Read character by character.
  while ((a_char = getc(fasta_file)) != EOF) {

    // Check for the beginning of the next sequence.
    if (a_char == '>') {
      // Put the ">" back onto the stream for the next call to find.
      ungetc(a_char, fasta_file);
      break;
    }

    // Skip non-alphabetic characters.
    if (!isalpha((int)a_char)) {
      if ((a_char != ' ') && (a_char != '\t') && (a_char != '\n') && (a_char != '\r')) {
        carp(CARP_WARNING,"Skipping character %c in sequence %s.",
             a_char, name);
      }

    } else {

      // Convert invalid characters to X.
      a_char = toupper((int)a_char);

      /**
       * Check the ASCII code.  If the char is above or below the
       * A(65)~Z(90)range, convert the character to an 'X'.
       */
      if ( (int)a_char < 65 || (int)a_char  > 90 ) {
        carp(CARP_WARNING, "Converting illegal character %c to X ",
             a_char);
        carp(CARP_WARNING, "in sequence %s.", name);
        a_char = 'X';
      }
      
      raw_sequence[i_seq] = a_char;
      i_seq++;
    }
    if (i_seq >= max_chars) {
      return_value = FALSE;
      break;
    }
  }
  raw_sequence[i_seq] = '\0';
  *sequence_length = i_seq; // chris added

  return(return_value);
}


/**
 * end of FASTA parsing
 * Thanks Bill!
 */


/** 
 * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
 */

/**
 * Additional get and set methods
 */

/**
 *\returns the id of the protein
 * returns a heap allocated new copy of the id
 * user must free the return id
 * assumes that the protein is heavy
 */
char* get_protein_id(
  PROTEIN_T* protein ///< the query protein -in 
  )
{
  if( protein == NULL ){
    carp(CARP_ERROR, "Cannot get protein id from NULL protein.\n");
    return NULL;
  } 

  if(protein->is_light){
    carp(CARP_FATAL, "Cannot get ID from light protein.");
  }
  
  int id_length = strlen(protein->id) +1; // +\0
  char* copy_id = 
    (char *)mymalloc(sizeof(char)*id_length);
  
  strncpy(copy_id, protein->id, id_length); 

  return copy_id;
}

/**
 *\returns a pointer to the id of the protein
 * assumes that the protein is heavy
 */
char* get_protein_id_pointer(
  PROTEIN_T* protein ///< the query protein -in 
  )
{
  if(protein->is_light){
    carp(CARP_FATAL, "Cannot get ID pointer from light protein.");
  }
  return protein->id; 
}

/**
 * sets the id of the protein
 */
void set_protein_id(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  const char* id ///< the sequence to add -in
  )
{
  free(protein->id);
  int id_length = strlen(id) +1; // +\0
  char* copy_id = 
    (char *)mymalloc(sizeof(char)*id_length);
  protein->id =
    strncpy(copy_id, id, id_length);  
}

/**
 *\returns the sequence of the protein
 * returns a heap allocated new copy of the sequence
 * user must free the return sequence 
 * assumes that the protein is heavy
 */
char* get_protein_sequence(
  PROTEIN_T* protein ///< the query protein -in 
  )
{
  if(protein->is_light){
    carp(CARP_FATAL, "Cannot get sequence from light protein.");
  }
  unsigned int sequence_length = strlen(protein->sequence) +1; // +\0
  char * copy_sequence = 
    (char *)mymalloc(sizeof(char)*sequence_length);
  return strncpy(copy_sequence, protein->sequence, sequence_length);  
}

/**
 *\returns a pointer to the sequence of the protein
 * assumes that the protein is heavy
 */
char* get_protein_sequence_pointer(
  PROTEIN_T* protein ///< the query protein -in 
  )
{
  if(protein->is_light){
    carp(CARP_FATAL, "Cannot get sequence pointer from light protein.");
  }
  return protein->sequence;
}

/**
 * sets the sequence of the protein
 */
void set_protein_sequence(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  const char* sequence ///< the sequence to add -in
  )
{
  free(protein->sequence);
  unsigned int sequence_length = strlen(sequence) +1; // +\0
  char * copy_sequence = 
    (char *)mymalloc(sizeof(char)*sequence_length);

  protein->sequence =
    strncpy(copy_sequence, sequence, sequence_length);  
}

/**
 *\returns the length of the protein
 * assumes that the protein is heavy
 */
unsigned int get_protein_length(
  PROTEIN_T* protein ///< the query protein -in 
  )
{
  return protein->length;
}

/**
 * sets the id of the protein
 */
void set_protein_length(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  unsigned int length ///< the length to add -in
  )
{
  protein->length = length;
}

/**
 *\returns the annotation of the protein
 * returns a heap allocated new copy of the annotation
 * user must free the return annotation
 * assumes that the protein is heavy
 */
char* get_protein_annotation(
  PROTEIN_T* protein ///< the query protein -in 
  )
{
  if(protein->is_light){
    carp(CARP_FATAL, "Cannot get annotation from light protein.");
  }
  int annotation_length = strlen(protein->annotation) +1; // +\0
  char * copy_annotation = 
    (char *)mymalloc(sizeof(char)*annotation_length);
  return strncpy(copy_annotation, protein->annotation, annotation_length);  
}

/**
 * sets the annotation of the protein
 */
void set_protein_annotation(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  const char* annotation ///< the sequence to add -in
  )
{
  if( annotation == NULL ){
    return;
  }

  if(!protein->is_light){
    free(protein->annotation);
  }
  int annotation_length = strlen(annotation) +1; // +\0
  char * copy_annotation = 
    (char *)mymalloc(sizeof(char)*annotation_length);
  protein->annotation =
    strncpy(copy_annotation, annotation, annotation_length);  
}

/**
 * sets the offset of the protein in the fasta file
 */
void set_protein_offset(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  unsigned long int offset ///< The file location in the source file in the database -in
  )
{
  protein->offset = offset;
}

/**
 *\returns the offset the protein
 */
unsigned long int get_protein_offset(
  PROTEIN_T* protein ///< the query protein -in 
  )
{
  return protein->offset;
}

/**
 * sets the protein_idx (if, idx=n, nth protein in the fasta file)
 */
void set_protein_protein_idx(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  unsigned int protein_idx ///< The index of the protein in it's database. -in
  )
{
  // carp(CARP_DETAILED_DEBUG, "set protein idx = %i", protein_idx);
  protein->protein_idx = protein_idx;
}

/**
 *\returns the protein_idx field
 */
unsigned int get_protein_protein_idx(
  PROTEIN_T* protein ///< the query protein -in 
  )
{
  return protein->protein_idx;
}

/**
 * sets the is_light field (is the protein a light protein?)
 */
void set_protein_is_light(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  BOOLEAN_T is_light ///< is the protein a light protein? -in
  )
{
  protein->is_light = is_light;
}

/**
 *\returns TRUE if the protein is light protein
 */
BOOLEAN_T get_protein_is_light(
  PROTEIN_T* protein ///< the query protein -in 
  )
{
  return protein->is_light;
}

/**
 * sets the database for protein
 */
void set_protein_database(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  DATABASE_T*  database ///< Which database is this protein part of -in
  )
{
  protein->database = copy_database_ptr(database);
}

/**
 *\returns Which database is this protein part of
 */
DATABASE_T* get_protein_database(
  PROTEIN_T* protein ///< the query protein -in 
  )
{
  return protein->database;
}

/**
 * Iterator
 * iterates over the peptides given a partent protein and constraints
 */

/*
 * Takes a cumulative distribution of peptide masses (the mass_array) and
 * the start index and end index and returns a peptide mass
 */
FLOAT_T calculate_subsequence_mass (
    double* mass_array,
    int start_idx,
    int cur_length
  ){

  FLOAT_T mass_h2o = MASS_H2O_AVERAGE;
  if(get_mass_type_parameter("isotopic-mass") == MONO){
    mass_h2o = MASS_H2O_MONO;
  }

  // carp(CARP_DETAILED_DEBUG, "mass start = %i", start_idx);
  int end_idx = start_idx + cur_length;
  // carp(CARP_DETAILED_DEBUG, "mass end = %i", end_idx);
  FLOAT_T peptide_mass = mass_array[end_idx] - mass_array[start_idx] + mass_h2o;

  return peptide_mass;
}

/**
 * \brief Decide if a residue is in an inclusion list or is not in an
 * exclusion list. 
 *
 * For use with the user-specified enzyme digestion.  Takes an amino
 * acid, a list of amino acids, and a flag for if it is an inclusion
 * list or an exclusion list.  A cleavage can happen before/after the
 * given residue if it is either in the inclusion list or is not in
 * the exculsion list.
 * \returns TRUE if the residue is in the inclusion list or not in the
 * exclusion list.
 */
BOOLEAN_T is_residue_legal(char aa, 
                           char* aa_list, 
                           int list_size, 
                           BOOLEAN_T for_inclusion){

  // The logic for returning for_inclusion:
  // For an inclusion list (TRUE), once we find the aa it passes (TRUE)
  // For an exclusion list (FALSE), once we find the aa, it fails (FALSE)
  int idx=0;
  for(idx=0; idx < list_size; idx++){
    if( aa == aa_list[idx] ){ return for_inclusion; }
  }
  // or if we got to the end of the list and didn't find a match
  // for inclusion, it fails (!TRUE)
  // for exclusion, it passes (!FALSE)
  return ! for_inclusion;
}

/**
 * Compares the first and second amino acids in the given sequence to
 * see if they conform to the cleavage rules of the given enzyme.  For
 * NO_ENZYME, always returns TRUE.
 *
 * \returns TRUE if this is a valid cleavage position for the given enzyme.
 */
BOOLEAN_T valid_cleavage_position(
   char* sequence,
   //   PEPTIDE_TYPE_T cleavage
   ENZYME_T enzyme
){

  switch(enzyme){

  case TRYPSIN:
    if ((sequence[0] == 'K' || sequence[0] == 'R') && (sequence[1] != 'P')){
      return TRUE;
    } else {
      return FALSE;
    }
    break;
    
  case CHYMOTRYPSIN:
    if ((sequence[0] == 'F' || sequence[0] == 'W' || sequence[0] == 'Y') 
        && (sequence[1] != 'P')){
      return TRUE;
    } else {
      return FALSE;
    }
    break;
    break;

  case ELASTASE:
    if ((sequence[0] == 'A' || sequence[0] == 'L' ||
         sequence[0] == 'I' || sequence[0] == 'V') 
        && (sequence[1] != 'P')){
      return TRUE;
    } else {
      return FALSE;
    }
    break;

  case CLOSTRIPAIN:
    if (sequence[0] == 'R'){
      return TRUE;
    } else {
      return FALSE;
    }
    break;

  case CYANOGEN_BROMIDE:
    if (sequence[0] == 'M'){
      return TRUE;
    } else {
      return FALSE;
    }
    break;

  case IODOSOBENZOATE:
    if (sequence[0] == 'W'){
      return TRUE;
    } else {
      return FALSE;
    }
    break;

  case PROLINE_ENDOPEPTIDASE:
    if (sequence[0] == 'P'){
      return TRUE;
    } else {
      return FALSE;
    }
    break;

  case STAPH_PROTEASE:
    if (sequence[0] == 'E'){
      return TRUE;
    } else {
      return FALSE;
    }
    break;

  case ASPN:
    if (sequence[1] == 'D'){
      return TRUE;
    } else {
      return FALSE;
    }
    break;

  case MODIFIED_CHYMOTRYPSIN:
    if ((sequence[0] == 'F' || sequence[0] == 'L' ||
         sequence[0] == 'W' || sequence[0] == 'Y') 
        && (sequence[1] != 'P')){
      return TRUE;
    } else {
      return FALSE;
    }
    break;

  case ELASTASE_TRYPSIN_CHYMOTRYPSIN:
    if ((sequence[0] == 'A' || sequence[0] == 'L' ||
         sequence[0] == 'I' || sequence[0] == 'V' ||
         sequence[0] == 'K' || sequence[0] == 'R' ||
         sequence[0] == 'W' || sequence[0] == 'F' ||
         sequence[0] == 'Y' ) 
        && (sequence[1] != 'P')){
      return TRUE;
    } else {
      return FALSE;
    }
    break;

  case CUSTOM_ENZYME:
    //carp(CARP_FATAL, "The custom enzyme is not yet implmented.");

    return ( is_residue_legal(sequence[0], 
                              pre_cleavage_list,
                              pre_list_size, 
                              pre_for_inclusion)
             && 
             is_residue_legal(sequence[1], 
                              post_cleavage_list,
                              post_list_size, 
                              post_for_inclusion) );
    break;

  case NO_ENZYME:
    return TRUE;
    break;

  case INVALID_ENZYME:
  case NUMBER_ENZYME_TYPES:
    carp(CARP_FATAL, "Cannot generate peptides with invalid enzyme.");
    break;

  }// end switch

  return FALSE;
}

/**
 * \brief Adds cleavages to the protein peptide iterator that obey iterator
 * constraint.
 *
 * Uses the allowed cleavages arrays, and whether skipped cleavages
 * are allowed. 
 * A small inconsistency: 
 *  Allowed cleavages start at 0, while the output cleavages start at 1.
 */
void iterator_add_cleavages(
    PROTEIN_PEPTIDE_ITERATOR_T* iterator, 
    int* nterm_allowed_cleavages, 
    int  nterm_num_cleavages, 
    int* cterm_allowed_cleavages, 
    int  cterm_num_cleavages, 
    BOOLEAN_T skip_cleavage_locations){

  // to avoid checking a lot of C-term before our current N-term cleavage
  int previous_cterm_cleavage_start= 0; 

  PEPTIDE_CONSTRAINT_T* constraint = iterator->peptide_constraint;
  int nterm_idx, cterm_idx;

  // iterate over possible nterm and cterm cleavage locations
  for (nterm_idx=0; nterm_idx < nterm_num_cleavages; nterm_idx++){
    
    int next_cterm_cleavage_start = previous_cterm_cleavage_start;
    BOOLEAN_T no_new_cterm_cleavage_start = TRUE;
    for (cterm_idx = previous_cterm_cleavage_start; 
         cterm_idx < cterm_num_cleavages; cterm_idx++){

      // if we have skipped a cleavage location, break to next nterm
      if(
         (skip_cleavage_locations == FALSE)
         &&
         ((*iterator->cumulative_cleavages)[nterm_allowed_cleavages[nterm_idx]] 
          < 
          (*iterator->cumulative_cleavages)[cterm_allowed_cleavages[cterm_idx]-1])
         ){
        break;
      }
      
      if (cterm_allowed_cleavages[cterm_idx] 
            <= nterm_allowed_cleavages[nterm_idx]){
        continue;
      }

      // check our length constraint
      int length = 
       cterm_allowed_cleavages[cterm_idx] - nterm_allowed_cleavages[nterm_idx];

      if (length < get_peptide_constraint_min_length(constraint)){
        continue;
      } else if (length > get_peptide_constraint_max_length(constraint)){
        break;
      } else if (no_new_cterm_cleavage_start){
        next_cterm_cleavage_start = cterm_idx;
        no_new_cterm_cleavage_start = FALSE;
      }
     
      // check our mass constraint
      FLOAT_T peptide_mass = calculate_subsequence_mass(iterator->mass_array, 
          nterm_allowed_cleavages[nterm_idx], length);

      if ((get_peptide_constraint_min_mass(constraint) <= peptide_mass) && 
          (peptide_mass <= get_peptide_constraint_max_mass(constraint))){ 

        // we have found a peptide
        iterator->nterm_cleavage_positions->push_back(nterm_allowed_cleavages[nterm_idx] + 1);

        iterator->peptide_lengths->push_back(length);
        iterator->peptide_masses->push_back(peptide_mass);
        carp(CARP_DETAILED_DEBUG, 
            "New pep: %i (%i)", nterm_allowed_cleavages[nterm_idx], length);

        iterator->num_cleavages++;
      }
    }
    previous_cterm_cleavage_start = next_cterm_cleavage_start;

  }
}

/**
 * Creates the data structures in the protein_peptide_iterator object necessary
 * for creating peptide objects.
 * - mass_array - cumulative distribution of masses. used to determine 
 *     the mass of any peptide subsequence.
 * - nterm_cleavage_positions - the nterm cleavage positions of the 
 *     peptides that satisfy the protein_peptide_iterator contraints
 * - peptide_lengths - the lengths of the peptides that satisfy the constraints
 * - cumulative_cleavages - cumulative distribution of cleavage positions
 *    used to determine if a cleavage location has been skipped
 */
void prepare_protein_peptide_iterator(
    PROTEIN_PEPTIDE_ITERATOR_T* iterator
  )
{
  prepare_protein_peptide_iterator_mc(iterator,  
    get_boolean_parameter("missed-cleavages"));
}

void prepare_protein_peptide_iterator_mc(
    PROTEIN_PEPTIDE_ITERATOR_T* iterator,
    BOOLEAN_T missed_cleavages)
{
  PROTEIN_T* protein = iterator->protein;
  MASS_TYPE_T mass_type = get_peptide_constraint_mass_type(
      iterator->peptide_constraint);
  double* mass_array = (double*)mycalloc(protein->length+1, sizeof(double));

  //  PEPTIDE_TYPE_T pep_type = get_peptide_type_parameter("cleavages");
  ENZYME_T enzyme = get_peptide_constraint_enzyme(iterator->peptide_constraint);
  FLOAT_T mass_h2o = MASS_H2O_AVERAGE;

  // set correct H2O mass
  if(mass_type == MONO){
    mass_h2o = MASS_H2O_MONO;
  }
  
  // initialize mass matrix and enzyme cleavage positions
  int* cleavage_positions = (int*) mycalloc(protein->length+1, sizeof(int));
  int* non_cleavage_positions = (int*)mycalloc(protein->length+1, sizeof(int));
  int* all_positions = (int*) mycalloc(protein->length+1, sizeof(int));

  // initialize first value in all array except non_cleavage_positions
  unsigned int start_idx = 0;
  mass_array[start_idx] = 0.0;
  int cleavage_position_idx = 0;
  int non_cleavage_position_idx = 0;
  cleavage_positions[cleavage_position_idx++] = 0;

  // calculate our cleavage positions and masses
  for(start_idx = 1; start_idx < protein->length+1; start_idx++){
    int sequence_idx = start_idx - 1; 
    mass_array[start_idx] = mass_array[start_idx-1] + 
      get_mass_amino_acid(protein->sequence[sequence_idx], mass_type);

    // increment cumulative cleavages before we check if current position
    // is a cleavage site because cleavages come *after* the current amino acid
    iterator->cumulative_cleavages->push_back(cleavage_position_idx);

    //if (valid_cleavage_position(protein->sequence + sequence_idx)){ 
    if (valid_cleavage_position(protein->sequence + sequence_idx, enzyme)){ 
      cleavage_positions[cleavage_position_idx++] = sequence_idx + 1;
    } else {
      non_cleavage_positions[non_cleavage_position_idx++] = sequence_idx + 1;
    }

    all_positions[sequence_idx] = sequence_idx;
  }

  // put in the implicit cleavage at end of protein
  if (cleavage_positions[cleavage_position_idx-1] != (int)protein->length){
    cleavage_positions[cleavage_position_idx++] = protein->length; 
  }

  all_positions[protein->length] = (int)protein->length;

  int num_cleavage_positions = cleavage_position_idx;
  int num_non_cleavage_positions = non_cleavage_position_idx;
  iterator->mass_array = mass_array;

  carp(CARP_DETAILED_DEBUG, "num_cleavage_positions = %i", num_cleavage_positions);

  // now determine the cleavage positions that actually match our constraints

  DIGEST_T digestion = 
    get_peptide_constraint_digest(iterator->peptide_constraint);

  switch (digestion){

  case FULL_DIGEST:
      iterator_add_cleavages(iterator,
        cleavage_positions, num_cleavage_positions-1,
        cleavage_positions+1, num_cleavage_positions-1, 
        missed_cleavages);

      break;

  case PARTIAL_DIGEST:
      // add the C-term tryptic cleavage positions.
      iterator_add_cleavages(iterator,
        all_positions, protein->length,
        cleavage_positions+1, num_cleavage_positions-1, 
        missed_cleavages);

      // add the N-term tryptic cleavage positions.
      // no +1 below for non_cleavage_positions below 
      // because it does not include sequence beginning. it is *special*
      iterator_add_cleavages(iterator,
        cleavage_positions, num_cleavage_positions-1,
        non_cleavage_positions, num_non_cleavage_positions-1,
        missed_cleavages);

      break;

  case NON_SPECIFIC_DIGEST:
      iterator_add_cleavages(iterator,
        all_positions, protein->length,
        all_positions+1, protein->length, // len-1?
        TRUE); // for unspecific ends, allow internal cleavage sites
      break;

  case INVALID_DIGEST:
  case NUMBER_DIGEST_TYPES:
    carp(CARP_FATAL, "Invalid digestion type in protein peptide iterator.");
  }

/*
  int idx;
  for (idx=0; idx < iterator->num_cleavages; idx++){
    carp(CARP_DETAILED_DEBUG, "%i->%i", 
         iterator->nterm_cleavage_positions[idx], 
         //iterator->peptide_lengths[idx], 
         //iterator->peptide_lengths[idx], 
       iterator->protein->sequence[iterator->nterm_cleavage_positions[idx]-1]);
  }
*/
  if (iterator->num_cleavages > 0){
    iterator->has_next = TRUE;
  } else { 
    iterator->has_next = FALSE;
  }

  free(cleavage_positions);
  free(non_cleavage_positions);
  free(all_positions);

}

/**
 * \brief Estimate the maximum number of peptides a protein can
 * produce.  Counts the number of subsequences of length
 * min_seq_length, min_seq_length + 1, ..., max_seq_length that can be
 * formed from a protein of the given length.  No enzyme specificity
 * assumed.  
 */
unsigned int count_max_peptides(
 unsigned int protein_length,   ///< length of protein
 unsigned int min_seq_length,   ///< min peptide length
 unsigned int max_seq_length)  ///< max peptide length
{
  if( max_seq_length > protein_length ){
    max_seq_length = protein_length;
  }

  unsigned int total_peptides = 0;
  for(unsigned int len = min_seq_length; len <= max_seq_length; len++){
    total_peptides += protein_length + 1 - len;
  }
  return total_peptides;
}

/**
 * Instantiates a new peptide_iterator from a protein.
 * \returns a PROTEIN_PEPTIDE_ITERATOR_T object.
 * assumes that the protein is heavy
 */
PROTEIN_PEPTIDE_ITERATOR_T* new_protein_peptide_iterator(
  PROTEIN_T* protein, ///< the protein's peptide to iterate -in
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraints -in
  )
{
  PROTEIN_PEPTIDE_ITERATOR_T* iterator = (PROTEIN_PEPTIDE_ITERATOR_T*)
    mycalloc(1, sizeof(PROTEIN_PEPTIDE_ITERATOR_T));

  // initialize iterator
  iterator->peptide_idx = 0;
  iterator->peptide_constraint 
    = copy_peptide_constraint_ptr(peptide_constraint);
  iterator->cur_start = 0; 
  iterator->cur_length = 1;  
  iterator->num_mis_cleavage 
    = get_peptide_constraint_num_mis_cleavage(iterator->peptide_constraint);
  iterator->protein = protein;

  iterator->nterm_cleavage_positions = new vector<int>();
  iterator->peptide_lengths = new vector<int>();
  iterator->peptide_masses = new vector<FLOAT_T>();
  iterator->cumulative_cleavages = new vector<int>();

  // estimate array size and reserve space to avoid resizing vector
  int max_peptides = count_max_peptides(protein->length, 
                                        get_int_parameter("min-length"),
                                        get_int_parameter("max-length"));
  iterator->nterm_cleavage_positions->reserve(max_peptides); 
  iterator->peptide_lengths->reserve(max_peptides);
  iterator->peptide_masses->reserve(max_peptides);
  iterator->cumulative_cleavages->reserve(max_peptides);
  iterator->num_cleavages = 0;

  // prepare the iterator data structures
  prepare_protein_peptide_iterator(iterator);
  return iterator;
}


/**
 * Frees an allocated peptide_iterator object.
 */
void free_protein_peptide_iterator(
  PROTEIN_PEPTIDE_ITERATOR_T* protein_peptide_iterator 
    ///< the iterator to free -in
  )
{
  free_peptide_constraint(protein_peptide_iterator->peptide_constraint);
  free(protein_peptide_iterator->mass_array); 
  delete protein_peptide_iterator->nterm_cleavage_positions; 
  delete protein_peptide_iterator->peptide_lengths; 
  delete protein_peptide_iterator->peptide_masses; 
  delete protein_peptide_iterator->cumulative_cleavages; 
  free(protein_peptide_iterator);
}

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides, FALSE if not.
 */
BOOLEAN_T protein_peptide_iterator_has_next(
  PROTEIN_PEPTIDE_ITERATOR_T* protein_peptide_iterator 
    ///< the iterator of interest -in
  )
{
  return protein_peptide_iterator->has_next;
}

/**
 * \returns The next peptide in the protein, in an unspecified order
 * the Peptide is new heap allocated object, user must free it
 */
PEPTIDE_T* protein_peptide_iterator_next(
  PROTEIN_PEPTIDE_ITERATOR_T* iterator
  )
{
  if( !iterator->has_next){
    return NULL;
  }

  int cleavage_idx = iterator->current_cleavage_idx;
  int current_start = (*iterator->nterm_cleavage_positions)[cleavage_idx];
  int current_length = (*iterator->peptide_lengths)[cleavage_idx];
  FLOAT_T peptide_mass = (*iterator->peptide_masses)[cleavage_idx];

  // create new peptide
  PEPTIDE_T* peptide = new_peptide(current_length, peptide_mass, 
                                   iterator->protein, current_start);//, peptide_type);
  
  // update position of iterator
  ++iterator->current_cleavage_idx;

  // update has_next field
  if (iterator->current_cleavage_idx == iterator->num_cleavages){
    iterator->has_next = FALSE;
  } else {
    iterator->has_next = TRUE;
  }

  return peptide;
}

/**
 *\returns the protein that the iterator was created on
 */
PROTEIN_T* get_protein_peptide_iterator_portein(
  PROTEIN_PEPTIDE_ITERATOR_T* protein_peptide_iterator ///< working protein_peptide_iterator -in
  )
{
  return protein_peptide_iterator->protein;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

