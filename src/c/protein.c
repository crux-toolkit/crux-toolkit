/*****************************************************************************
 * \file protein.c
 * $Revision: 1.39 $
 * \brief: Object for representing a single protein.
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "utils.h"
#include "alphabet.h"
#include "objects.h"
#include "peptide.h"
#include "protein.h"
#include "peptide_src.h"
#include "database.h"
#include "carp.h"
#include "peptide_constraint.h"

/**
 * Constants
 */
#define PROTEIN_ID_LENGTH 100
#define PROTEIN_SEQUENCE_LENGTH 40000
#define PROTEIN_ANNOTATION_LENGTH 100
#define LONGEST_LINE PROTEIN_ID_LENGTH + PROTEIN_ID_LENGTH
#define FASTA_LINE 50
#define SMALLEST_MASS 57
#define LARGEST_MASS 190

/**
 * \struct protein 
 * \brief A protein sequence.
 */
struct protein{
  DATABASE_T*  database; ///< Which database is this protein part of
  unsigned long int offset; ///< The file location in the source file in the database
  unsigned int protein_idx; ///< The index of the protein in it's database.
  BOOLEAN_T    is_light; ///< is the protein a light protein?
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
 * 
 */
struct protein_peptide_iterator {
  PROTEIN_T* protein; ///< The protein whose peptides to iterate over. 
  unsigned short int cur_start; ///< Start in protein of the current peptide.
  unsigned short int cur_length; ///< The length of the current peptide.
  unsigned int peptide_idx; ///< The index of the current peptide.
  PEPTIDE_CONSTRAINT_T* peptide_constraint; ///< The type of peptide to iterate over.
  float** mass_matrix; ///< stores all the peptide's mass
  BOOLEAN_T has_next; ///< is there a next? 
  int num_mis_cleavage; ///< The maximum mis cleavage of the peptide
  unsigned int* seq_marker; ///< The array that marks all the 'K | R | P'
  unsigned int kr_idx; //idx for the closest to cur_start K | R is located
  unsigned int first_kr_idx; //idx for the first K | R is located
  BOOLEAN_T is_kr; //is there a K|R found in this sequence
};

//def bellow
static BOOLEAN_T read_title_line
  (FILE* fasta_file,
   char* name,
   char* description,
   PROTEIN_T* protein);

//def bellow
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
 * \returns A new protein object.
 * The protein is does not constain a database, users must provide one.
 */
PROTEIN_T* new_protein(
  char*         id, ///< The protein sequence id. -in
  char*   sequence, ///< The protein sequence. -in
  unsigned int length, ///< The length of the protein sequence. -in
  char* annotation,  ///< Optional protein annotation.  -in
  unsigned long int offset, ///< The file location in the source file in the database -in
  unsigned int protein_idx ///< The index of the protein in it's database.-in
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
  protein->is_light = FALSE;
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
  //protein already heavy
  if(!protein->is_light){
    return TRUE;
  }
  
  FILE* file = get_database_file(protein->database);
  
  //rewind to the begining of the protein to include ">" line
  fseek(file, protein->offset, SEEK_SET);

  //failed to parse the protein from fasta file
  //protein offset is set in the parse_protein_fasta_file method
  if(!parse_protein_fasta_file(protein ,file)){
    carp(CARP_ERROR, "failed convert protein to heavy, cannot parse fasta file");
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
  //protein already light
  if(protein->is_light){
    return TRUE;
  }
  //free all char* in protein object
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
  if(!protein->is_light){
    free(protein->id);
    free(protein->sequence);
    free(protein->annotation);
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
  //covnert to heavy protein
  if(protein->is_light){
    protein_to_heavy(protein);
  }
  int   sequence_index;
  unsigned int   sequence_length = get_protein_length(protein);
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


//FIXME ID line and annotation might need to be fixed
VERBOSE_T verbosity = NORMAL_VERBOSE;
/**
 * Parses a protein from an open (FASTA) file.
 * the protein_idx field of the protein must be added before or after you parse the protein
 * \returns TRUE if success. FALSE is failure.
 * protein must be a heap allocated
 */
BOOLEAN_T parse_protein_fasta_file(
  PROTEIN_T* protein, ///< protein object to fill in -out
  FILE* file ///< fasta file -in
  )
{
  static char name[LONGEST_LINE];     // Just the sequence ID.
  static char desc[LONGEST_LINE];     // Just the comment field.
  static char buffer[PROTEIN_SEQUENCE_LENGTH];        // The sequence, as it's read in.
  static unsigned int sequence_length; //the sequence length

  // Read the title line.
  if (!read_title_line(file, name, desc, protein)) {
    return(FALSE);
  }
  
  //need this line to initialize alphabet to set for protein instead of DNA
  set_alphabet(verbosity, "ACDEFGHIKLMNPQRSTVWY"); 
  buffer[0] = '\0';

  // Read the sequence.
  if (!read_raw_sequence(file, name, PROTEIN_SEQUENCE_LENGTH, buffer, &sequence_length)) {
    die("Sequence %s is too long.\n", name);
  }
    
  //update the protein object.
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
  //set protein offset                   FIXME: might not need to "-1" -CHRIS
  protein->offset = ftell(fasta_file) - 1;

  //chris edited, added this block to make sure all of comment line is read 
  //although might not be stored, to ensure the file* is at start of the sequence
  {
    char* new_line = NULL;
    int line_length;
    size_t buf_length = 0;

    if((line_length =  getline(&new_line, &buf_length, fasta_file)) == -1){
      die("Error reading Fasta file.\n");
    }
    strncpy(id_line, new_line, LONGEST_LINE-1);
    free(new_line);
  }

  //this is Bill's old code
  /*
  // Read the ID and comment line.
  if (fgets(id_line, LONGEST_LINE-1, fasta_file) == NULL) {
    die("Error reading Fasta file.\n");
  }
  */

  // Remove EOL.
  id_line[strlen(id_line) - 1] = '\0';

  // Extract the ID from the beginning of the line.
  if (sscanf(id_line, "%s", name) != 1) {
    die("Error reading sequence ID.\n%s\n", id_line);
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
  // char a_char;
  // tlb; change a_char to integer so it will compile on SGI
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
	carp(CARP_WARNING,"Skipping character %c in sequence %s.\n",
		a_char, name);
      }

    } else {

      // Convert invalid characters to X.
      a_char = toupper((int)a_char);

      /**
       * this code check the character against to what verbosity you set
       * bill's code, char_in_string can be found in utils.c
       * very slow!!
      if (!char_in_string(get_alphabet(TRUE), a_char)) {
	carp(CARP_WARNING, "Converting illegal character %c to X ",
		a_char);
	carp(CARP_WARNING, "in sequence %s.", name);
	a_char = 'X';
      }
      */
      
      /**
       * To speed up the process, checks the ASCII code, 
       * if the char is above or bellow the A(65)~Z(90)range,
       * converts the character to a 'X'
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
  
  if(protein->is_light){
    die("protein is light, must be heavy");
  }
  
  int id_length = strlen(protein->id) +1; //+\0
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
    die("protein is light, must be heavy");
  }
  return protein->id; 
}

/**
 * sets the id of the protein
 */
void set_protein_id(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  char* id ///< the sequence to add -in
  )
{
  free(protein->id);
  int id_length = strlen(id) +1; //+\0
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
    die("protein is light, must be heavy");
  }
  unsigned int sequence_length = strlen(protein->sequence) +1; //+\0
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
    die("protein is light, must be heavy");
  }
  return protein->sequence;
}

/**
 * sets the sequence of the protein
 */
void set_protein_sequence(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  char* sequence ///< the sequence to add -in
  )
{
  free(protein->sequence);
  unsigned int sequence_length = strlen(sequence) +1; //+\0
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
    die("protein is light, must be heavy");
  }
  int annotation_length = strlen(protein->annotation) +1; //+\0
  char * copy_annotation = 
    (char *)mymalloc(sizeof(char)*annotation_length);
  return strncpy(copy_annotation, protein->annotation, annotation_length);  
}

/**
 * sets the annotation of the protein
 */
void set_protein_annotation(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  char* annotation ///< the sequence to add -in
  )
{
  if(!protein->is_light){
    free(protein->annotation);
  }
  int annotation_length = strlen(annotation) +1; //+\0
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
  protein->database = database;
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

/**
 * find if there's a K or R in the middle of peptide, skips K|R with P after them
 * returns TRUE if no missed cleavage sites discovered
 */
BOOLEAN_T find_krp(
  unsigned int* seq_marker,
  unsigned int kr_idx,
  unsigned int start_idx, ///< the start index of peptide, 1 is the first residue -in 
  unsigned int end_idx
  ){
  if(kr_idx < end_idx && 
     kr_idx >= start_idx){
    //is there a P after the K|R 
    if(seq_marker[kr_idx] == 1){ //check boundary
      return find_krp(seq_marker, seq_marker[kr_idx-1], start_idx, end_idx);
    }
    return FALSE;
  }
  
  else if(kr_idx >= end_idx || seq_marker[kr_idx-1] >= end_idx){
    return TRUE;
  }
  else if(seq_marker[kr_idx-1] < end_idx && 
          seq_marker[kr_idx-1] >= start_idx){
    //is there a P after the K|R 
    if(seq_marker[seq_marker[kr_idx-1]] == 1){ //check boundary
      return find_krp(seq_marker, seq_marker[kr_idx-1], start_idx, end_idx);
    }
    return FALSE;
  }
  return TRUE;
}

//FIXME only examines if there is a mis-cleavage or not
// eventually would like to implement so that it will return the total number of mis-cleavage
/**
 * examines the peptide if it contains miscleavage sites within it's sequence
 * \returns 0 if no miscleavage sites, 1 if there exist at least 1 mis cleavage sites
 */
int examine_peptide_cleavage(
  PROTEIN_PEPTIDE_ITERATOR_T* iterator, ///< working iterator -in
  unsigned int start_idx, ///< the start index of peptide, 1 is the first residue -in 
  unsigned int end_idx ///< the end index of peptide -in
  )
{
  //the K|R|P index array
  unsigned int* seq_marker = iterator->seq_marker; 
  unsigned int kr_idx = iterator->kr_idx;
  
  if(kr_idx == iterator->protein->length + 1 || 
     find_krp(seq_marker, kr_idx, start_idx, end_idx)){
    return 0;
  }
  return 1;
}

/**
 * examines the peptide with context of it's parent protein to determine it's type
 * \returns the peptide type
 */
PEPTIDE_TYPE_T examine_peptide_type(
  PROTEIN_PEPTIDE_ITERATOR_T* iterator, ///< working iterator -in
  unsigned int start_idx, ///< the start index of peptide, 1 is the first residue -in 
  unsigned int end_idx ///< the end index of peptide -in
  )
{
  //the K|R|P index array
  unsigned int* seq_marker = iterator->seq_marker; 
  BOOLEAN_T front = FALSE;
  BOOLEAN_T back = FALSE;
  
  //examine front cleavage site
  if(start_idx == 1 || ((seq_marker[start_idx-2] > 1) && seq_marker[start_idx-1] != 1)){
    front = TRUE;
  }
  
  //exaimine end cleavage site
  if(end_idx == iterator->protein->length || 
     ((seq_marker[end_idx-1] > 1) && seq_marker[end_idx] != 1)){
    back = TRUE;
  }

  if(front && back){
    return TRYPTIC;
  }
  else if(front || back){
    return PARTIALLY_TRYPTIC;
  }
  else{
    return NOT_TRYPTIC;
  }
}

/**
 * Finds the next peptide that fits the constraints
 * \returns TRUE if there is a next peptide. FALSE if not.
 */
BOOLEAN_T iterator_state_help(
  PROTEIN_PEPTIDE_ITERATOR_T* iterator, 
  int max_length,  ///< constraints: max length -in
  int min_length, ///< constraints: min length -in
  float max_mass, ///< constraints: max mass -in
  float min_mass, ///< constraints: min mass -in
  PEPTIDE_TYPE_T peptide_type  ///< constraints: peptide type -in
  )
{
  LOOP:
  
  //set kr_idx position
  if(iterator->is_kr){
    if(iterator->seq_marker[iterator->kr_idx-1] < iterator->cur_start){
      iterator->kr_idx = iterator->seq_marker[iterator->kr_idx-1];
    }
    else if(iterator->kr_idx != iterator->first_kr_idx &&
            iterator->kr_idx > iterator->cur_start){
      iterator->kr_idx = iterator->first_kr_idx;
    }
  }

  //check if the smallest mass of a length is larger than max_mass
  if(iterator->cur_length * SMALLEST_MASS > max_mass){
    return FALSE;
  }
  
  //check if out of mass_max idex size
  if(iterator->cur_length > max_length ||
     iterator->cur_length > iterator->protein->length){
    return FALSE;
  }
  
  //check if less than min length
  if(iterator->cur_length < min_length){
    ++iterator->cur_length;
    goto LOOP;
  }
  
  //reached end of length column, check next length
  if((unsigned int)(iterator->cur_start + iterator->cur_length - 1) > iterator->protein->length){
    ++iterator->cur_length;
    iterator->cur_start = 1;
    goto LOOP;
  }
  
  //is mass with in range
  if(iterator->mass_matrix[iterator->cur_length-1][iterator->cur_start-1] < min_mass ||
        iterator->mass_matrix[iterator->cur_length-1][iterator->cur_start-1] > max_mass){
    //does this length have any possibility of having a peptide within mass range?
    if(iterator->mass_matrix[iterator->cur_length-1][iterator->cur_start-1] == 0 ||
       iterator->cur_length * LARGEST_MASS + 19 < min_mass){
      ++iterator->cur_length;
      iterator->cur_start = 1;
    }
    else{
      ++iterator->cur_start;
    }
    goto LOOP;
  }
  
  //examin tryptic type and cleavage
  if(peptide_type != ANY_TRYPTIC){
    if((examine_peptide_type(iterator, 
                             iterator->cur_start, 
                             iterator->cur_length + iterator->cur_start -1) != peptide_type))
      {
        ++iterator->cur_start;
        goto LOOP;
      }
  }

  //examine cleavage
  if(iterator->num_mis_cleavage == 0){
    if(examine_peptide_cleavage(iterator, 
                                iterator->cur_start, 
                                iterator->cur_length + iterator->cur_start -1) != 0)
      {
        ++iterator->cur_start;
        goto LOOP;
      }
  }
    
  return TRUE;
  
}

/**
 * sets the iterator to the next peptide that fits the constraints
 * \returns TRUE if there is a next peptide. FALSE if not.
 */
BOOLEAN_T set_iterator_state(
  PROTEIN_PEPTIDE_ITERATOR_T* iterator  ///< set iterator to next peptide -in
  )
{
  int max_length = get_peptide_constraint_max_length(iterator->peptide_constraint);
  int min_length = get_peptide_constraint_min_length(iterator->peptide_constraint);
  float max_mass = get_peptide_constraint_max_mass(iterator->peptide_constraint);
  float min_mass = get_peptide_constraint_min_mass(iterator->peptide_constraint);
  PEPTIDE_TYPE_T peptide_type = get_peptide_constraint_peptide_type(iterator->peptide_constraint);
  
  return iterator_state_help(iterator, max_length, min_length, max_mass, min_mass, peptide_type);
}


//start_size: total sequence size 
//length_size: max_length from constraint;
//float** mass_matrix = float[length_size][start_size];
/**
 * Dynamically sets the mass of the mass_matrix
 * The mass matrix contains every peptide bellow max length
 * must pass in a heap allocated matrix
 */
void set_mass_matrix(
  float** mass_matrix,  ///< the mass matrix -out
  unsigned int start_size,  ///< the y axis size -in
  unsigned int length_size, ///< the x axis size -in
  PROTEIN_T* protein, ///< the parent protein -in
  MASS_TYPE_T mass_type ///< isotopic mass type (AVERAGE, MONO) -in
  )
{
  unsigned int start_index = 0;
  unsigned int length_index = 1;
  float mass_h2o = MASS_H2O_AVERAGE;

  //set correct H2O mass
  if(mass_type == MONO){
    mass_h2o = MASS_H2O_MONO;
  }
  
  //initialize mass matrix
  for(; start_index < start_size; ++start_index){
    mass_matrix[0][start_index] = 
      get_mass_amino_acid(protein->sequence[start_index], mass_type) + mass_h2o;
  }
  start_index = 0;
  
  //fill in the mass matrix
  for(; start_index < start_size; ++start_index){
    for(length_index = 1; length_index < length_size; ++length_index){
      if(start_index + length_index < protein->length){
        mass_matrix[length_index][start_index] = 
          mass_matrix[length_index - 1][start_index] + mass_matrix[0][start_index + length_index] - mass_h2o; 
      }
    }
  }
}

/**
 * set seq_marker array in protein_peptide_iterator
 * creates an int array the length of the protein sequence.
 * 0 if not 'K | R | P', if 'P' 1, for 'K|R' you store the index of the next incident of 'K|R'  
 * thus, you can look at the index if it is K|R, and if it is you can tell where the next K|R is
 * sets the sequence idx(starts 1) where the first incident of K | R, returns protein_length+1 if no K|R exist
 */
void set_seq_marker(
  PROTEIN_PEPTIDE_ITERATOR_T* protein_peptide_iterator  ///< the iterator to set -out
  )
{
  char* sequence = protein_peptide_iterator->protein->sequence;
  unsigned int protein_length = protein_peptide_iterator->protein->length;
  unsigned int sequence_idx = 0;
  unsigned int previous_kr = 0;
  unsigned int first_kr = protein_length+1;
  BOOLEAN_T first = FALSE; //have we seen a K | R so far

  //create seq_marker
  unsigned int* seq_marker = (unsigned int*)mycalloc(protein_length, sizeof(unsigned int));
  
  //iterate through the sequence
  for(; sequence_idx < protein_length; ++sequence_idx){
    if(sequence[sequence_idx] == 'K' || sequence[sequence_idx] == 'R'){
      if(first){
        seq_marker[previous_kr] = sequence_idx + 1;
      }
      else{
        first = TRUE;
        first_kr = sequence_idx+1;
      }
      previous_kr = sequence_idx;
    }
    else if(sequence[sequence_idx] == 'P'){
      seq_marker[sequence_idx] = 1;
    }
  }
  if(first){
    //last K|R is set to protein length
    seq_marker[previous_kr] = sequence_idx + 1;
  }
  
  //set the complete seq_marker
  protein_peptide_iterator->seq_marker = seq_marker;
  protein_peptide_iterator->first_kr_idx = first_kr;
  protein_peptide_iterator->kr_idx = first_kr;
  protein_peptide_iterator->is_kr = first;
}

/**
 * Instantiates a new peptide_iterator from a peptide.
 * \returns a PROTEIN_PEPTIDE_ITERATOR_T object.
 * assumes that the protein is heavy
 */
PROTEIN_PEPTIDE_ITERATOR_T* new_protein_peptide_iterator(
  PROTEIN_T* protein, ///< the protein's peptide to iterate -in
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraints -in
  )
{
  unsigned int matrix_index = 0;
  unsigned int max_length = get_peptide_constraint_max_length(peptide_constraint);
  MASS_TYPE_T mass_type = get_peptide_constraint_mass_type(peptide_constraint);

  PROTEIN_PEPTIDE_ITERATOR_T* iterator = 
    (PROTEIN_PEPTIDE_ITERATOR_T*)mycalloc(1, sizeof(PROTEIN_PEPTIDE_ITERATOR_T));

  //create mass_matrix
  iterator->mass_matrix = (float**)mycalloc(max_length, sizeof(float*));
  for (; matrix_index < max_length ; ++matrix_index){
    iterator->mass_matrix[matrix_index] = (float*)mycalloc(protein->length, sizeof(float));
  }  
  set_mass_matrix(iterator->mass_matrix, protein->length, max_length, protein, mass_type);
  
  //initialize iterator
  iterator->protein = protein;
  iterator->peptide_idx = 0;
  iterator->peptide_constraint = peptide_constraint;
  iterator->cur_start = 1; // must cur_start-1 for access mass_matrix
  iterator->cur_length = 1;  // must cur_length-1 for access mass_matrix
  iterator->num_mis_cleavage = get_peptide_constraint_num_mis_cleavage(iterator->peptide_constraint);
  set_seq_marker(iterator);
  iterator->has_next = set_iterator_state(iterator);
  return iterator;
}


/**
 * free the heap allocated mass_matrix
 */
void free_mass_matrix(
  float** mass_matrix,  ///< the mass matrix to free -in
  unsigned int length_size ///< the x axis size -in
  )
{
  unsigned int matrix_idx = 0;
  for (; matrix_idx  < length_size; ++matrix_idx){
    free(mass_matrix[matrix_idx]);
  }

  free(mass_matrix);
}

/**
 * Frees an allocated peptide_iterator object.
 */
void free_protein_peptide_iterator(
  PROTEIN_PEPTIDE_ITERATOR_T* protein_peptide_iterator ///< the iterator to free -in
  )
{
  free_mass_matrix(protein_peptide_iterator->mass_matrix, 
                   get_peptide_constraint_max_length(protein_peptide_iterator->peptide_constraint));
  free(protein_peptide_iterator->seq_marker);
  free(protein_peptide_iterator);
}

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptides to iterate over, FALSE if not.
 */
BOOLEAN_T protein_peptide_iterator_has_next(
  PROTEIN_PEPTIDE_ITERATOR_T* protein_peptide_iterator ///< the iterator of interest -in
  )
{
  return protein_peptide_iterator->has_next;
}


/**
 * \returns The next peptide in the protein, in an unspecified order
 * the Peptide is new heap allocated object, user must free it
 */
PEPTIDE_T* protein_peptide_iterator_next(
  PROTEIN_PEPTIDE_ITERATOR_T* protein_peptide_iterator
  )
{
  PEPTIDE_TYPE_T peptide_type;

  if(!protein_peptide_iterator->has_next){
    free_protein_peptide_iterator(protein_peptide_iterator);
    die("ERROR: no more peptides\n");
  }
  
  //set peptide type
  if(get_peptide_constraint_peptide_type(protein_peptide_iterator->peptide_constraint) != ANY_TRYPTIC){
    peptide_type = get_peptide_constraint_peptide_type(protein_peptide_iterator->peptide_constraint);
  }

  //when constraints ANY_TRYPTIC, need to examine peptide
  //possible to skip this step and leave it as ANY_TRYPTIC
  else{
      peptide_type = 
        examine_peptide_type(protein_peptide_iterator,
                             protein_peptide_iterator->cur_start,
                             protein_peptide_iterator->cur_start + protein_peptide_iterator->cur_length -1);
  }
  

  //create new peptide
  PEPTIDE_T* peptide = 
    new_peptide
    (protein_peptide_iterator->cur_length, 
     protein_peptide_iterator->mass_matrix[protein_peptide_iterator->cur_length-1][protein_peptide_iterator->cur_start-1],
     protein_peptide_iterator->protein,
     protein_peptide_iterator->cur_start,
     peptide_type);
  
  //FIXME Not sure what the use delete this field if needed
  ++protein_peptide_iterator->peptide_idx;

  //update poisiton of iterator
  ++protein_peptide_iterator->cur_start;
  protein_peptide_iterator->has_next = set_iterator_state(protein_peptide_iterator);

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

