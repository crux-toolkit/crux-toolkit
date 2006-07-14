/*****************************************************************************
 * \file protein.c
 * $Revision: 1.13 $
 * \brief: Object for representing a single protein.
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include "utils.h"
#include "objects.h"
#include "peptide.h"
#include "protein.h"
#include "protein_peptide_association.h"

/**
 * Constants
 */
#define PROTEIN_ID_LENGTH 100
#define PROTEIN_SEQUENCE_LENGTH 10000
#define PROTEIN_ANNOTATION_LENGTH 100
#define LONGEST_LINE PROTEIN_ID_LENGTH + PROTEIN_ID_LENGTH
#define FASTA_LINE 50

/**
 * \struct protein 
 * \brief A protein sequence.
 */
struct protein{
  char*         id; ///< The protein sequence id.
  char*   sequence; ///< The protein sequence.
  int       length; ///< The length of the protein sequence.
  char* annotation; ///< Optional protein annotation. 
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
  unsigned char cur_length; ///< The length of the current peptide.
  unsigned int peptide_idx; ///< The index of the current peptide.
  PEPTIDE_CONSTRAINT_T* peptide_constraint; ///< The type of peptide to iterate over.
};

//def bellow
static BOOLEAN_T read_title_line
  (FILE* fasta_file,
   char* name,
   char* description);

//def bellow
static BOOLEAN_T read_raw_sequence
  (FILE* fasta_file,   // Input Fasta file.
   char* name,         // Sequence ID (used in error messages).
   int   max_chars,    // Maximum number of allowed characters.
   char* raw_sequence, // Pre-allocated sequence.
   int* sequence_length // the sequence length -chris added
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
 */
PROTEIN_T* new_protein(
  char*         id, ///< The protein sequence id.
  char*   sequence, ///< The protein sequence.
  int       length, ///< The length of the protein sequence.
  char* annotation  ///< Optional protein annotation. 
  )
{
  PROTEIN_T* protein = allocate_protein();
  set_protein_id(id);
  set_protein_sequence(sequence);
  set_protein_length(length);
  set_protein_annotation(annotation);
  return protein;
}         

/**
 * Frees an allocated protein object.
 */
void free_protein(
  PROTEIN_T* protein ///< object to free -in
  )
{
  free(protein->id);
  free(protein->sequence);
  free(protein->annotation);
  free(protein);
}

/**
 * Prints a protein object to file.
 */
void print_protein(
  PROTEIN_T* protein, ///< protein to print -in
  FILE* file ///< output stream -out
  )
{
  int   sequence_index;
  int   sequence_length = get_protein_length(protein);
  char* sequence = get_protein_sequence(protein);
  char* id = get_protein_id(protein);
  char* annotation = get_protein_annotation(protein)

  fprintf(file, ">%s %s\n", id, protein);

  sequence_index = 0;
  while (sequence_length - sequence_index > FASTA_LINE) {
    fprintf(file, "%.*s\n", FASTA_LINE, &(sequence[sequence_index]));
    sequence_index += FASTA_LINE;
  }
  fprintf(outfile, "%s\n\n", &(sequence[sequence_index]));

  free(sequence);
  free(id);
  free(annotation);
}

/**
 * Copies protein object src to dest.
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

  free(id);
  free(sequence);
  free(annotation);
}


//FIXME ID line and annotation might need to be fixed
/**
 * Parses a protein from an open (FASTA) file.
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
  static int sequence_length; //the sequence length

  // Read the title line.
  if (!read_title_line(fasta_file, name, desc)) {
    return(FALSE);
  }
  
  // Read the sequence.
  buffer[0] = '\0';
  if (!read_raw_sequence(fasta_file, name, PROTEIN_SEQUENCE_LENGTH, buffer, &sequence_length)) {
    die("Sequence %s is too long.\n", name);
  }
    
  //update the protein object.
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
   char* description)
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

  // Read the ID and comment line.
  if (fgets(id_line, LONGEST_LINE-1, fasta_file) == NULL) {
    die("Error reading Fasta file.\n");
  }

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
   int   max_chars,    // Maximum number of allowed characters.
   char* raw_sequence, // Pre-allocated sequence.
   int* sequence_length // the sequence length -chris added
   )
{
  // char a_char;
  // tlb; change a_char to integer so it will compile on SGI
  int a_char;
  int i_seq;
  BOOLEAN_T return_value = TRUE;

  // Start at the end of the given sequence.
  i_seq = strlen(raw_sequence);
  assert((int)strlen(raw_sequence) < max_chars);

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
	fprintf(stderr, "Warning: Skipping character %c in sequence %s.\n",
		a_char, name);
      }

    } else {

      // Convert invalid characters to X.
      a_char = toupper((int)a_char);
      if (!char_in_string(get_alphabet(TRUE), a_char)) {
	fprintf(stderr, "Warning: Converting illegal character %c to X ",
		a_char);
	fprintf(stderr, "in sequence %s.\n", name);
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
 */
char* get_protein_id(
  PROTEIN_T* protein ///< the query protein -in 
  )
{
  int id_length = strlen(protein->id) +1; //+\0
  char* copy_id = 
    (char *)mymalloc(sizeof(char)*id_length);
  return strncpy(copy_id, protein->id, id_length); 
}

/**
 * sets the id of the protein
 */
void set_protein_id(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  char* id ///< the sequence to add -in
  )
{
  free(peptide->id);
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
 */
char* get_protein_sequence(
  PROTEIN_T* protein ///< the query protein -in 
  )
{
  int sequence_length = strlen(protein->sequence) +1; //+\0
  char * copy_sequence = 
    (char *)mymalloc(sizeof(char)*sequence_length);
  return strncpy(copy_sequence, protein->sequence, sequence_length);  
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
  int sequence_length = strlen(sequence) +1; //+\0
  char * copy_sequence = 
    (char *)mymalloc(sizeof(char)*sequence_length);

  protein->sequence =
    strncpy(copy_sequence, sequence, sequence_length);  
}

/**
 *\returns the length of the protein
 */
int get_protein_length(
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
  int length ///< the length to add -in
  )
{
  protein->length = length;
}

/**
 *\returns the annotation of the protein
 * returns a heap allocated new copy of the annotation
 * user must free the return annotation
 */
char* get_protein_annotation(
  PROTEIN_T* protein ///< the query protein -in 
  )
{
  int annotation_length = strlen(protein->annotation) +1; //+\0
  char * copy_annotation = 
    (char *)mymalloc(sizeof(char)*annotation_length);
  return strncpy(copy_annotation, protein->annotation, sequence_annotation);  
}

/**
 * sets the annotation of the protein
 */
void set_protein_annotation(
  PROTEIN_T* protein, ///< the protein to set it's fields -out
  char* annotation ///< the sequence to add -in
  )
{
  free(protein->annotation);
  int annotation_length = strlen(annotation) +1; //+\0
  char * copy_annotation = 
    (char *)mymalloc(sizeof(char)*annotation_length);
  protein->annotation =
    strncpy(copy_annotation, annotation, annotation_length);  
}










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
  unsigned char cur_length; ///< The length of the current peptide.
  unsigned int peptide_idx; ///< The index of the current peptide.
  PEPTIDE_CONSTRAINT_T* peptide_constraint; ///< The type of peptide to iterate over.

  unsigned char* KRP_map; ///< shows the location of the residue K|R and P 
  float mass;  ///< the current mass
  PEPTIDE_TYPE_T peptide_type; ///< the current peptide_type
  BOOLEAN_T has_next; ///< is there a next? 

};


BOOLEAN_T set_iterator_state(PROTEIN_PEPTIDE_ITERATOR_T* iterator){
  if(iterator->cur_start == iterator->protein->length-1){
    return FALSE;
  }
  ++iterator->length;
  iterator->mass +=  get_mass_amino_acid_average(iterator->protein->sequence[iterator->cur_start]);

  if(iterator->cur_start != 0){
    if(iterator->KRP_map[iterator->cur_start-1] != 0 &&
       iterator->peptide_type != NOT_TRYPTIC){
      iterator->peptide_type = PARTIALLY_TRYPTIC;
    }
    else if(iterator->peptide_type != NOT_TRYPTIC){
      iterator->peptide_type = TRYPTIC;
    }
  }

  if(iterator->cur_start + iterator->length < iterator->protein->length){
    if(iterator->KRP_map[iterator->cur_start-1] != 0){
      if(iterator->peptide_type == TRYPTIC){
        iterator->peptide_type = PARTIALLY_TRYPTIC;
      }
      else if(iterator->peptide_type == PARTIALLY_TRYPTIC){
        iterator->peptide_type != NOT_TRYPTIC
      }
    }
  } 
  
  //change to new location

  //recursive call to set_iterator_state...
  
}



/**
 * Iterator
 * iterates over the peptides given a partent protein and constraints
 */

/**
 * Instantiates a new peptide_iterator from a peptide.
 * \returns a PROTEIN_PEPTIDE_ITERATOR_T object.
 */
PROTEIN_PEPTIDE_ITERATOR_T* new_protein_peptide_iterator(
  PROTEIN_T* protein, ///< the protein's peptide to iterate -in
  PEPTIDE_CONSTRAINT_T* peptide_constraint ///< the peptide constraints -in
  )
{
  PROTEIN_PEPTIDE_ITERATOR_T* iterator = 
    (PROTEIN_PEPTIDE_ITERATOR_T*)mycalloc(1, sizeof(PROTEIN_PEPTIDE_ITERATOR_T));
  iterator->protein = protein;
  iterator->cur_start = 0;
  iterator->cur_length = 1;
  iterator->peptide_idx = 0;
  iterator->mass = get_mass_amino_acid_average(protein->sequence[iterator->cur_start]);
  iterator->peptide_type = TRYPTIC;
  iterator->peptide_constraint = peptide_constraint;
  iterator->KRP_map = make_KRP_map(protein->sequence, protein->length);
  iterator->has_next = set_iterator_state(iterator);
  return iterator;
}


/**
 * \returns a K|R, P residue map of the protein
 * K|R are 0, P is 1, other residues are 3
 */
unsigned char* make_KRP_map(
  char* sequence,  ///< the sequence to examine
  int length  ///< the length of sequence
  )
{
 int sequence_idx = 0;
 unsigned char* KRP_map = (unsigned char*)mymalloc(sizeof(unsigned char)*length);
  
 for(; sequence_idx < length; ++sequence_idx){
   if(sequence[sequence_idx] == 'K' || sequence[sequence_idx] == 'R'){
     KRP_map[sequence_idx] = 0;
   }
   else if(sequence[sequence_idx] == 'P'){
     KRP_map[sequence_idx] = 1;
   }
   else{
     KRP_map[sequence_idx] = 3;
   }
 }
 return KRP_map;
}

/**
 * Frees an allocated peptide_iterator object.
 */
void free_protein_peptide_iterator(
  PROTEIN_PEPTIDE_ITERATOR_T* protein_peptide_iterator ///< the iterator to free -in
  )
{
  free(protein_peptide_iterator->KRP_map);
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
  BOOLEAN has_next = TRUE;

  do{
    int state = examin_state(iterator);
    if(state == 0){
      ++iterator->cur_start;
      iterator->cur_length = 1;
      iterator->peptide_type = TRYPTIC;
      iterator->mass =  get_mass_amino_acid_average(iterator->protein->sequence[iterator->cur_start]) ;
    }
    else if(state == 1){
      
    }
  }
  while();

  return has_next; 
   

  
  iterator->length = 

}



/**
 * \returns The next peptide in the protein, in an unspecified order
 */
PEPTIDE_T* protein_peptide_iterator_next(
    PROTEIN_PEPTIDE_ITERATOR_T* protein_peptide_iterator);



/*
PEPTIDE_TYPE_T examin_peptide_type(
  char* sequence, ///< the parent protein -in
  int start_idx, ///< the start index of peptide -in
  int end_idx, ///< the end index of peptide -in
  )
{
  int current_idx = start_idx;
  BOOLEAN_T start = TRUE;
  BOOLEAN_T end =TRUE;

  //check start position must be cleaved at K or R residue
  if(start != 0){
    if(sequence[start-1] != 'K' && sequence[start-1] != 'R'){
      start = FALSE;
    }
  }
  for(; current_idx <= end_idx; ++current_idx){
    


  }
}
*/

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

