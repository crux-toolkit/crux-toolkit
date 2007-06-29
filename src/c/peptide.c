/*****************************************************************************
 * \file peptide.c
 * $Revision: 1.52 $
 * \brief: Object for representing a single peptide.
 ****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "crux-utils.h"
#include "objects.h"
#include "mass.h"
#include "peptide.h"
#include "protein.h"
#include "database.h"
#include "carp.h"

/**
 * static global variable
 * determines if the peptide src are created by link lists or array
 * if TRUE, peptides are implented with link list peptide src, else array
 */
static BOOLEAN_T PEPTIDE_SRC_USE_LINK_LIST;

/**
 * \struct peptide
 * \brief A subsequence of a protein.
 */
struct peptide {
  unsigned char length; ///< The length of the peptide
  float peptide_mass;   ///< The peptide's mass.
  PEPTIDE_SRC_T* peptide_src; ///< a linklist of peptide_src   
};

 
/**
 * \struct residue_iterator
 * \brief Object to iterate over the residues in a peptide, starting at the
 * first residue of the peptide, and proceeding in order.
 */
struct residue_iterator {
  PEPTIDE_T*  peptide; ///< The peptide whose residues to iterate over.
  char*   sequence;    ///< The peptide sequence
  int     residue_idx; ///< The index of the current peak
};

/**
 * \struct peptide_src_iterator
 * \brief Object to iterate over the peptide_srcs linklist in a peptide
 */
struct peptide_src_iterator{
  PEPTIDE_T*  peptide; ///< The peptide whose peptide_srcs to iterate over.
  PEPTIDE_SRC_T* current; ///< the current peptide_srcs
};

/**
 * \returns An (empty) peptide object.
 */
PEPTIDE_T* allocate_peptide(void){
  PEPTIDE_T* peptide = (PEPTIDE_T*)mycalloc(1, sizeof(PEPTIDE_T));
  peptide->peptide_src = NULL;
  return peptide;
}

/**
 * \returns The mass of the given peptide.
 */
float calc_peptide_mass(
  PEPTIDE_T* peptide, ///< the query peptide -in
  MASS_TYPE_T mass_type ///< isotopic mass type (AVERAGE, MONO) -in
  )
{
  float peptide_mass = 0;
  RESIDUE_ITERATOR_T * residue_iterator = new_residue_iterator(peptide);
  
  while(residue_iterator_has_next(residue_iterator)){
    peptide_mass += get_mass_amino_acid(residue_iterator_next(residue_iterator), mass_type);
  }
  free_residue_iterator(residue_iterator);

  if(mass_type == AVERAGE){
    return peptide_mass + MASS_H2O_AVERAGE;
  }
  return peptide_mass + MASS_H2O_MONO;
}

//FIXME association part might be need to change
/**
 * \returns A new peptide object, populated with the user specified parameters.
 */
PEPTIDE_T* new_peptide(
  unsigned char length,     ///< The length of the peptide -in
  float peptide_mass,       ///< The neutral mass of the peptide -in
  PROTEIN_T* parent_protein, ///< the parent_protein of this peptide -in
  int start_idx, ///< the start index of this peptide in the protein sequence -in
  PEPTIDE_TYPE_T peptide_type ///<  The type of peptides(TRYPTIC, PARTIALLY_TRYPTIC, NOT_TRYPTIC, ANY_TRYPTIC) -in
  )
{
  PEPTIDE_T* peptide = allocate_peptide();
  set_peptide_length( peptide, length);
  set_peptide_peptide_mass( peptide, peptide_mass);
  peptide->peptide_src =
    new_peptide_src( peptide_type, parent_protein, start_idx );
  
  //increment the database pointer count
  add_database_pointer_count(get_protein_database(parent_protein));
  
  return peptide;
}
  

/** 
 * \returns the neutral mass of the peptide
 */
float get_peptide_neutral_mass(
  PEPTIDE_T* peptide ///< the query peptide -in
  )
{
  return get_peptide_peptide_mass(peptide);
}

/** 
 * \returns the mass of the peptide if it had charge "charge"
 */
float get_peptide_charged_mass(
    PEPTIDE_T* peptide, ///< the query peptide -in
    int charge ///< charge of peptide -in
    )
{
  return get_peptide_mz(peptide, charge) * charge;
}

/** 
 * \returns the m/z of the peptide if it had charge "charge"
 */
float get_peptide_mz(
    PEPTIDE_T* peptide, ///< the query peptide -in
    int charge ///< the charge of peptide -in
    )
{
  return ((get_peptide_peptide_mass(peptide) + MASS_H * charge)/ charge);
}


/**
 * Frees an allocated peptide object.
 * Depending on peptide_src implementation determines how to free srcs
 * This decision is made by global variable PEPTIDE_SRC_USE_LINK_LIST
 */
void free_peptide(
  PEPTIDE_T* peptide ///< peptide to free -in
  )
{
  //decrement the pointer count
  free_database(get_peptide_first_src_database(peptide));

  //check which implementation peptide_src uses
  if(!PEPTIDE_SRC_USE_LINK_LIST){
    //array implementation
    free(peptide->peptide_src);    
  }
  else{
    //link list implementation
    free_peptide_src(peptide->peptide_src);
  }

  free(peptide);
}

/**
 * FIXME, don't need this anymore, may delete
 *
 * Frees an allocated peptide object.
 * This one is used when the peptides is created throuh
 * parsing the index, the peptide_src is not a link list, but
 * an array, thus needs it's own free_peptide version
 */
void free_peptide_for_array(
  PEPTIDE_T* peptide ///< peptide to free -in
  )
{
  //decrement the pointer count
  free_database(get_peptide_first_src_database(peptide));
  
  //array
  free(peptide->peptide_src);
  free(peptide);
}

/**
 * Prints a peptide object to file.
 * prints all information.
 */
void print_peptide(
  PEPTIDE_T* peptide,  ///< the query peptide -in
  FILE* file  ///< the out put stream -out
  )
{
  char* sequence = get_peptide_sequence(peptide);
  
  PEPTIDE_SRC_ITERATOR_T* iterator = 
    new_peptide_src_iterator(peptide);
  fprintf(file,"%s\n","Peptide");
  fprintf(file,"%.2f\t",peptide->peptide_mass);
  fprintf(file,"%d\t",peptide->length);
  fprintf(file,"%s\n",sequence);
  //interate through the linklist of possible parent proteins
  while(peptide_src_iterator_has_next(iterator)){
    print_peptide_src(peptide_src_iterator_next(iterator), file);
  }
  free_peptide_src_iterator(iterator);
  free(sequence);
}

/**
 * Prints a peptide object to file.
 * prints all peptide_src object it's associated 
 * mass \\t protein-id \\t peptide-start \\t peptide-length <\\t peptide-sequence> \n
 *      \\t protein-id \\t peptide-start \\t peptide-length <\\t peptide-sequence> \n
 * prints in correct format for generate_peptide
 */
void print_peptide_in_format(
  PEPTIDE_T* peptide,  ///< the query peptide -in
  BOOLEAN_T flag_out, ///< print peptide sequence? -in
  FILE* file  ///< the out put stream -out
  )
{
  PROTEIN_T* parent = NULL;
  PEPTIDE_SRC_T* next_src = peptide->peptide_src;
  char* id = NULL;
  int start_idx = 0;
  char* sequence = NULL;
  //BOOLEAN_T light = FALSE;

  //print mass of the peptide
  fprintf(file, "%.2f", peptide->peptide_mass);

  //obtain peptide sequence
  if(flag_out){
    parent = get_peptide_src_parent_protein(next_src);
    
    //covnert to heavy protein
    /*
      FIXME, IF use light heavy put back
    if(get_protein_is_light(parent)){
      protein_to_heavy(parent);
      light = TRUE;
    }
    */
    sequence = get_peptide_sequence(peptide);
  }

  //iterate over all peptide src
  while(next_src != NULL){
    //if(!light){
      parent = get_peptide_src_parent_protein(next_src);
      
      //covnert to heavy protein
      /*
        FIXME, IF use light heavy put back
      if(get_protein_is_light(parent)){
        protein_to_heavy(parent);
        light = TRUE;
      }
      */
      //}

    id = get_protein_id_pointer(parent);
    start_idx = get_peptide_src_start_idx(next_src);
        
    fprintf(file, "\t%s\t%d\t%d", id, start_idx, peptide->length);
  
    //print peptide sequence?
    if(flag_out){
      fprintf(file, "\t%s\n", sequence);
    }
    else{
      fprintf(file, "\n");
    }
    next_src = get_peptide_src_next_association(next_src);
    
    /** 
     * uncomment this code if you want to restore a protein to 
     * light after converted to heavy
    //convert back to light
    if(light){
      protein_to_light(parent);
      light = FALSE;
    }
    */
  }

  //free sequence if allocated
  if(flag_out){
    free(sequence);
  }
}

/**
 * Prints a peptide object to file.
 * ONLY prints peptide_src that match the peptide_src
 * mass \\t protein-id \\t peptide-start \\t peptide-length <\\t peptide-sequence> \n
 *      \\t protein-id \\t peptide-start \\t peptide-length <\\t peptide-sequence> \n
 * prints in correct format for generate_peptide
 */
void print_filtered_peptide_in_format(
  PEPTIDE_T* peptide,  ///< the query peptide -in
  BOOLEAN_T flag_out, ///< print peptide sequence? -in
  FILE* file,  ///< the out put stream -out
  PEPTIDE_TYPE_T peptide_type ///< the peptide_type of src to print -in
  )
{
  PROTEIN_T* parent = NULL;
  PEPTIDE_SRC_T* next_src = peptide->peptide_src;
  char* id = NULL;
  int start_idx = 0;
  char* sequence = NULL;
  //BOOLEAN_T light = FALSE;

  //print mass of the peptide
  fprintf(file, "%.2f", peptide->peptide_mass);

  //obtain peptide sequence
  if(flag_out){
    parent = get_peptide_src_parent_protein(next_src);
    
    //covnert to heavy protein
    /*
    FIXME, IF use light heavy put back
    if(get_protein_is_light(parent)){
      protein_to_heavy(parent);
      light = TRUE;
    }
    */
    sequence = get_peptide_sequence(peptide);
  }

  //iterate over all peptide src
  while(next_src != NULL){
    if(peptide_type == ANY_TRYPTIC ||
       peptide_type == get_peptide_src_peptide_type(next_src)){

      //if(!light){
      parent = get_peptide_src_parent_protein(next_src);
        
      //covnert to heavy protein
      /*
      FIXME, IF use light heavy put back
      if(get_protein_is_light(parent)){
        protein_to_heavy(parent);
        light = TRUE;
      }
      */
        //}
      
      id = get_protein_id_pointer(parent);
      start_idx = get_peptide_src_start_idx(next_src);
      
      fprintf(file, "\t%s\t%d\t%d", id, start_idx, peptide->length);
      
      //print peptide sequence?
      if(flag_out){
        fprintf(file, "\t%s\n", sequence);
      }
      else{
        fprintf(file, "\n");
      }
    
      /** 
       * uncomment this code if you want to restore a protein to 
       * light after converted to heavy
      //convert back to light
      if(light){
        protein_to_light(parent);
        light = FALSE;
      }
      */
    }
    next_src = get_peptide_src_next_association(next_src);
  }

  //free sequence if allocated
  if(flag_out){
    free(sequence);
  }
}


//TESTME
/**
 * Copies peptide object src to dest.
 * dest must be a heap allocated peptide
 */
void copy_peptide(
  PEPTIDE_T* src, ///< source peptide -in
  PEPTIDE_T* dest ///< destination peptide -out
  )
{
 
  PEPTIDE_SRC_T* new_association;

  set_peptide_length(dest, get_peptide_length(src));
  set_peptide_peptide_mass(dest, get_peptide_peptide_mass(src));

  //copy all of the peptide_src in the peptide
  new_association = allocate_peptide_src();
  copy_peptide_src(src->peptide_src, new_association);
  set_peptide_peptide_src(dest, new_association);

  //increment the database pointer count
  //since we are creating a new peptide
  add_database_pointer_count(get_peptide_first_src_database(src));
}

//FIXME needs to be rewritten for the new output format -Chris
/**
 * Parses a peptide from file.
 * \returns TRUE if success. FALSE if failure.
 */
/*
BOOLEAN_T parse_peptide_file(
  PEPTIDE_T* peptide,
  FILE* file
  )
{

  char* new_line = NULL;
  int line_length;
  size_t buf_length = 0;
  char* sequence;
  float peptide_mass;
  int peptide_length;

  line_length =  getline(&new_line, &buf_length, file);
  if( (sscanf(new_line,"%f\t%d\t%s\n", // parse peptide line
              &peptide_mass, &peptide_length, sequence) != 3) && 
      (!isalpha((int)sequence[0]))){
    free(new_line);
    return FALSE;
  }
  set_peptide_sequence(sequence);  
  set_peptide_peptide_mass(peptide_mass);
  set_peptide_length(peptide_length);
  free(new_line);
  
  return TRUE;
}
*/

/** 
 * Access routines of the form get_<object>_<field> and set_<object>_<field>. 
 * FIXME Chris, could you create the get and set methods for the object fields?
 */

/**
 * Additional get and set methods
 */

/**
 * get the peptide->first peptide_src->parent protein->database
 */
DATABASE_T* get_peptide_first_src_database(
  PEPTIDE_T* peptide ///< working peptide -in
  )
{
  return get_protein_database(get_peptide_src_parent_protein(peptide->peptide_src));
}

/**
 * \returns the sequence of peptide
 * goes to the first peptide_src to gain sequence, thus must have at least one peptide src
 * returns a char* to a heap allocated copy of the sequence
 * user must free the memory
 */
char* get_peptide_sequence(
 PEPTIDE_T* peptide ///< peptide to query sequence -in
 )
{
  if(peptide->peptide_src == NULL){
    die("ERROR: no peptide_src to retrieve peptide sequence\n");
  }

  char* parent_sequence = 
    get_protein_sequence_pointer(get_peptide_src_parent_protein(peptide->peptide_src));
  int start_idx = get_peptide_src_start_idx(peptide->peptide_src);

  char* copy_sequence = copy_string_part(&parent_sequence[start_idx-1], peptide->length);
 
  return copy_sequence; 
}

/**
 * \returns a pointer to the start of peptide sequence with in it's protein parent sequence
 * goes to the first peptide_src to find the location of start, thus must have at least one peptide src
 * should not print, will result in printing the entire protein sequence
 */
char* get_peptide_sequence_pointer(
  PEPTIDE_T* peptide ///< peptide to query sequence -in
  )
{
  if(peptide->peptide_src == NULL){
    die("ERROR: no peptide_src to retrieve peptide sequence pointer\n");
  }
  char* parent_sequence = 
    get_protein_sequence_pointer(get_peptide_src_parent_protein(peptide->peptide_src));
  int start_idx = get_peptide_src_start_idx(peptide->peptide_src);

  char* pointer_peptide_sequence = &parent_sequence[start_idx-1];
  
  return pointer_peptide_sequence;
}

/**
 * \returns the sequence of peptide with each flanking AA from the specified peptide_src(protein)
 *  template "*.peptide_sequence.*", where "*" are flanking amino acids
 * "*", left empty if no flanking sequence
 * returns a char* to a heap allocated copy of the sequence
 * user must free the memory
 */
char* get_peptide_sequence_from_peptide_src_sqt(
 PEPTIDE_T* peptide, ///< peptide to query sequence -in
 PEPTIDE_SRC_T* peptide_src ///< peptide_src -in 
 )
{
  char* copy_sequence = NULL;
  PROTEIN_T* protein = get_peptide_src_parent_protein(peptide_src);
  //get peptide start idx of protein in prarent protein
  int start_idx = get_peptide_src_start_idx(peptide_src);
  //parent protein length
  int protein_length = get_protein_length(protein);
  //parent protein length
  char* parent_sequence = 
    get_protein_sequence_pointer(protein);

  //allocate peptide memory
  //Template "*.peptide.*", where "*" are flanking amino acids
  copy_sequence = (char*)mycalloc(peptide->length+5, sizeof(char));
  //first copy over the peptide sequences
  strncpy(&copy_sequence[2], &parent_sequence[start_idx-1], peptide->length);
  
  //copy over template
  copy_sequence[1] = '.';
  copy_sequence[peptide->length+2] = '.';
  copy_sequence[peptide->length+4] = '\0';

  //copy flanking AA
  //is there a AA before?
  if(start_idx != 1){
    copy_sequence[0] = parent_sequence[start_idx-2];
  }
  //is there a AA after?
  if((start_idx + peptide->length - 1) < protein_length){
    copy_sequence[peptide->length+3] = parent_sequence[start_idx+peptide->length-1];
  }
  
  //yeah return!!
  return copy_sequence; 
}

/**
 * \returns the sequence of peptide with each flanking AA
 *  template "*.peptide_sequence.*", where "*" are flanking amino acids
 * "*", left empty if no flanking sequence
 * goes to the first peptide_src to gain sequence, thus must have at least one peptide src
 * returns a char* to a heap allocated copy of the sequence
 * user must free the memory
 */
char* get_peptide_sequence_sqt(
 PEPTIDE_T* peptide ///< peptide to query sequence -in
 )
{
  if(peptide->peptide_src == NULL){
    die("ERROR: no peptide_src to retrieve peptide sequence\n");
  }
  
  //yeah return!!
  return get_peptide_sequence_from_peptide_src_sqt(peptide, peptide->peptide_src);
}

/**
 * sets the sequence length of the peptide
 * length maximum of 255
 */
void set_peptide_length( 
  PEPTIDE_T* peptide,  ///< the peptide to set the length -out
  unsigned char length  ///< the length of sequence -in
  )
{
  peptide->length = length;
}

/**
 *\returns the sequence length of the peptide
 */
unsigned char get_peptide_length( 
  PEPTIDE_T* peptide  ///< the peptide to query the length -in
  )
{
  return peptide->length;
}

/**
 * sets the peptide mass
 */
void set_peptide_peptide_mass( 
  PEPTIDE_T* peptide,  ///< the peptide to set -out
  float peptide_mass  ///< the mass of the peptide - in
  )
{
  peptide->peptide_mass = peptide_mass;
}

/**
 * \returns the peptide mass
 */
inline float get_peptide_peptide_mass( 
  PEPTIDE_T* peptide  ///< the peptide to query the mass -in
  )
{
  return peptide->peptide_mass;
}

/**
 * sets the peptide_src field in the peptide
 * this method should be ONLY used when the peptide has no existing list of peptide_src
 * use add_peptide_peptide_src method to add to existing list
 * must pass on a heap allocated peptide_src object
 * does not copy in the object, just the pointer to the object.
 */
void set_peptide_peptide_src(
  PEPTIDE_T* peptide,  ///< the peptide to set -out                                             
  PEPTIDE_SRC_T* new_association ///< new peptide_src -in
  )
{
  peptide->peptide_src = new_association;
}

/**
 * this method adds the new_association to the end of the existing peptide's 
 * linklist of peptide_srcs
 * if no prior existing list, adds it at the front
 * must pass on a heap allocated peptide_src object
 * does not copy in the object, just the pointer to the object.
 */
void add_peptide_peptide_src(
  PEPTIDE_T* peptide,  ///< the peptide to set -out
  PEPTIDE_SRC_T* new_association ///< new peptide_src -in
  )
{
  PEPTIDE_SRC_T* add_association = peptide->peptide_src;
  PEPTIDE_SRC_ITERATOR_T* iterator = NULL;
  
  //is the peptide src list empty?
  if(add_association == NULL){
    peptide->peptide_src = new_association;
    return;
  }

  //create peptide src iterator
  iterator = new_peptide_src_iterator(peptide);

  //find the last peptide_src object in the list
  while(peptide_src_iterator_has_next(iterator)){
    add_association = peptide_src_iterator_next(iterator);
  }
  
  set_peptide_src_next_association(add_association, new_association);
  free_peptide_src_iterator(iterator);
}

/**
 * this method adds the peptide src array to an EMPTY peptide
 * only used in index.c, when the peptide src count for  peptide is known
 * Any existing peptide_src will lose it's reference
 */
void add_peptide_peptide_src_array(
  PEPTIDE_T* peptide,  ///< the peptide to set -out
  PEPTIDE_SRC_T* peptide_src_array ///< new peptide_src -in
  )
{
  //should be empty peptide src list
  peptide->peptide_src = peptide_src_array;

}


/**
 * returns a pointer to the peptide_protein_association field of the peptide
 */
PEPTIDE_SRC_T* get_peptide_peptide_src(
  PEPTIDE_T* peptide  ///< the peptide to query the peptide_peptide_src -in
  )
{
  return peptide->peptide_src;
}

/**
 * returns a pointer to the peptide's first parent protein field of the peptide
 */
PROTEIN_T* get_peptide_parent_protein(
  PEPTIDE_T* peptide  ///< the peptide to query the parent_protein -in
  )
{
  return get_peptide_src_parent_protein(peptide->peptide_src);
}

/**
 *\returns the protein struct size, value of sizeof function
 */
int get_peptide_sizeof(){
  return sizeof(PEPTIDE_T);
}

/**
 * sets the peptide src implementation in the peptide object
 * This should be set only once and not be altered
 */
void set_peptide_src_implementation(
  BOOLEAN_T use_link_list ///< does the peptide use link list peptide src
  )
{  
  PEPTIDE_SRC_USE_LINK_LIST = use_link_list; 
}
                                    
/**
 * Iterator
 */

/**
 * Instantiates a new residue_iterator from a peptide.
 * \returns a RESIDUE_ITERATOR_T object.
 */
RESIDUE_ITERATOR_T* new_residue_iterator(
  PEPTIDE_T* peptide ///< peptide sequence to iterate -in
  )
{
  RESIDUE_ITERATOR_T* residue_iterator =
    (RESIDUE_ITERATOR_T*)mycalloc(1, sizeof(RESIDUE_ITERATOR_T));
  
  residue_iterator->peptide =  peptide;
  residue_iterator->residue_idx = 0;
  residue_iterator->sequence = get_peptide_sequence(peptide);
  return residue_iterator;
}        

/**
 * Frees an allocated residue_iterator object.
 */
void free_residue_iterator(
  RESIDUE_ITERATOR_T* residue_iterator ///< free this object -in
  )
{
  free(residue_iterator->sequence);
  free(residue_iterator);
}

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional residues to iterate over, FALSE if not.
 */
BOOLEAN_T residue_iterator_has_next(
  RESIDUE_ITERATOR_T* residue_iterator ///< the query iterator -in
  )
{
  return (residue_iterator->residue_idx < residue_iterator->peptide->length);
}

/**
 * \returns The next residue (a character) in the peptide.
 */
char residue_iterator_next(
  RESIDUE_ITERATOR_T* residue_iterator  ///< the query iterator -in
  )
{
  ++residue_iterator->residue_idx;
  return residue_iterator->sequence[residue_iterator->residue_idx - 1];
}


/**
 * Protein peptide association Iterator
 */

/**
 * Instantiates a new peptide_src_iterator from a peptide.
 * \returns a PEPTIDE_SRC_T object.
 */
PEPTIDE_SRC_ITERATOR_T* new_peptide_src_iterator(
  PEPTIDE_T* peptide ///< peptide's fields to iterate -in
  )
{
  PEPTIDE_SRC_ITERATOR_T* association_iterator =
    (PEPTIDE_SRC_ITERATOR_T*)mycalloc(1,sizeof(PEPTIDE_SRC_ITERATOR_T));
  association_iterator->peptide = peptide;
  association_iterator->current = peptide->peptide_src;
  return association_iterator;
}

/**
 * Frees an allocated peptide_src_iterator object.
 */
void free_peptide_src_iterator(
  PEPTIDE_SRC_ITERATOR_T* peptide_src_iterator ///< free this object -in
  )
{
  free(peptide_src_iterator);

}

/**
 * The basic iterator functions.
 * \returns TRUE if there are additional peptide_srcs to iterate over, FALSE if not.
 */
BOOLEAN_T peptide_src_iterator_has_next(
  PEPTIDE_SRC_ITERATOR_T* peptide_src_iterator///< the query iterator -in
  )
{
  return !(peptide_src_iterator->current == NULL);
}

/**
 * \returns The next peptide_srcs in the peptide.
 */
PEPTIDE_SRC_T* peptide_src_iterator_next(
  PEPTIDE_SRC_ITERATOR_T* peptide_src_iterator///< the query iterator -in
  )
{
  PEPTIDE_SRC_T* previous = peptide_src_iterator->current;
  if(peptide_src_iterator->current != NULL){
    peptide_src_iterator->current = 
      get_peptide_src_next_association(peptide_src_iterator->current);
  }
  else{
    free(peptide_src_iterator);
    die("ERROR: no more peptide_srcs to iterate\n");
  }
  return previous;
}
 
/////

/*************************************
 * compare peptide method
 ************************************/

/**
 * Compare peptide sequence
 * \returns TRUE if peptide sequence is identical else FALSE
 */
BOOLEAN_T compare_peptide_sequence(
  PEPTIDE_T* peptide_one,  ///< the peptide sequence to compare  -out
  PEPTIDE_T* peptide_two  ///< the peptide sequence to compare  -out
  )
{
  //is Mass and Length identical
  if(compare_float(peptide_one->peptide_mass, peptide_two->peptide_mass) != 0 ||
     peptide_one->length != peptide_two->length){
    return FALSE;
  }
  else{
    int current_idx = 0;
    char* start_one = get_peptide_src_sequence_pointer(peptide_one->peptide_src);
    char* start_two = get_peptide_src_sequence_pointer(peptide_two->peptide_src);
    
    while(current_idx < peptide_one->length){
      if(start_one[current_idx] != start_two[current_idx]){
        return FALSE;
      }
      ++current_idx;
    } 
  }
  return TRUE;
}


/**
 * compares two peptides with the lexical sort type
 * for qsort
 * /returns 1 if peptide_one has lower priority, 0 if equal, -1 if greater priority
 */
int compare_peptide_lexical_qsort(
  PEPTIDE_T** peptide_one, ///< peptide to compare one -in
  PEPTIDE_T** peptide_two ///< peptide to compare two -in
  )
{
  //convert the protein to heavy if needed
  /* uncomment if needed, to use light heavy protein
  protein_to_heavy(get_peptide_parent_protein(peptide_one));
  protein_to_heavy(get_peptide_parent_protein(peptide_two));
  */
  PEPTIDE_T* peptide_1 = *peptide_one;
  PEPTIDE_T* peptide_2 = *peptide_two;
  char* peptide_one_sequence = get_peptide_sequence_pointer(peptide_1);
  char* peptide_two_sequence = get_peptide_sequence_pointer(peptide_2);
  int peptide_one_length = peptide_1->length;
  int peptide_two_length = peptide_2->length;
  int current_idx = 0;
  
  //check if all alphabetically identical
  while(current_idx < peptide_one_length &&
        current_idx < peptide_two_length){
    if(peptide_one_sequence[current_idx] > peptide_two_sequence[current_idx]){        
      return 1;
    }
    else if(peptide_one_sequence[current_idx] < peptide_two_sequence[current_idx]){
      return -1;
    }
    ++current_idx;
  }
  
  //alphabetically identical, check if same length
  if(peptide_one_length > peptide_two_length){
    return 1;
  }
  else if(peptide_one_length < peptide_two_length){
    return -1;
  }
  else{
    return 0;
  }
}

/**
 * compares two peptides with the mass sort type
 * if peptide mass is identical sort by lexicographical order
 * used for qsort function
 * /returns 1 if peptide_one has lower priority, 0 if equal, -1 if greater priority
 */
int compare_peptide_mass_qsort(
  PEPTIDE_T** peptide_one, ///< peptide to compare one -in
  PEPTIDE_T** peptide_two ///< peptide to compare two -in
  )
{
  //determine mass order
  int result = compare_float((*peptide_one)->peptide_mass, 
                             (*peptide_two)->peptide_mass);

  //if identical mass, sort in lexical order
  if(result == 0){
    return compare_peptide_lexical_qsort(peptide_one, peptide_two);
  }
  else{//if not identical
    return result;
  }
}

/**
 * compares two peptides with the length sort type
 * /returns 1 if peptide_one has lower priority, 0 if equal, -1 if greater priority
 */
int compare_peptide_length_qsort(
  PEPTIDE_T** peptide_one, ///< peptide to compare one -in
  PEPTIDE_T** peptide_two ///< peptide to compare two -in
  )
{
  int peptide_one_length = (*peptide_one)->length;
  int peptide_two_length = (*peptide_two)->length;
  
  //alphabetically identical, check if same length
  if(peptide_one_length > peptide_two_length){
    return 1;
  }
  else if(peptide_one_length < peptide_two_length){
    return -1;
  }
  else{
    return 0;
  }
}

/**
 * Compare peptide mass
 * \returns 0 if peptide mass is identical else 1 if peptide_one is larger, -1 if peptide_two is larger
 */
int compare_peptide_mass(
  PEPTIDE_T* peptide_one,  ///< the peptide mass to compare  -out
  PEPTIDE_T* peptide_two  ///< the peptide mass to compare  -out
  )
{
  return compare_float(peptide_one->peptide_mass, peptide_two->peptide_mass);
}

/**
 * Merge to identical peptides, copy all peptide_src into one of the peptide
 * peptide_dest, peptide_bye must have at least one peptide src
 * frees the peptide_bye, once the peptide_src are re-linked to the peptide_dest
 * Assumes that both peptides use linklist implemenation for peptide_src
 * \returns TRUE if merge is successful else FALSE
 */
BOOLEAN_T merge_peptides(
  PEPTIDE_T* peptide_dest, ///< the peptide to merge into  -out
  PEPTIDE_T* peptide_bye ///< the peptide to be merged  -in
  )
{
  PEPTIDE_SRC_T* current_src = peptide_dest->peptide_src;
  PEPTIDE_SRC_T* next_src = get_peptide_src_next_association(current_src);
  
  //does all peptides have at least one peptide_src?
  if(current_src == NULL || peptide_bye->peptide_src == NULL){
    carp(CARP_ERROR, "failed to merge two peptides");
    return FALSE;
  }

  //decrement the database pointer count
  free_database(get_peptide_first_src_database(peptide_bye));

  //find the end of the peptide src link list..
  while(next_src != NULL){
    current_src = next_src;
    next_src =  get_peptide_src_next_association(current_src);
  }
  set_peptide_src_next_association(current_src, peptide_bye->peptide_src);
  free(peptide_bye);
  return TRUE;
}

/*
 * Serialize a peptide to a FILE in binary
 * \returns TRUE if serialization is successful, else FALSE
 *
 * The peptide serialization format looks like this:
 *
 *<PEPTIDE_T: peptide struct><int: number of peptide_src>[<int: protein index><PEPTIDE_TYPE_T: peptide_type><int: peptide start index>]+
 * the bracket peptide src information repeats for the number of peptide src listed before the bracket
 * the protein index is the index of the parent protein in the database DATABASE_T
 *
 */
BOOLEAN_T serialize_peptide(
  PEPTIDE_T* peptide, ///< the peptide to serialize -in
  FILE* file ///< the output file to serlize -out
  )
{

  PEPTIDE_SRC_ITERATOR_T* iterator = 
    new_peptide_src_iterator(peptide);
  PEPTIDE_SRC_T* peptide_src = NULL;
  long num_src_location;
  long original_location;
  int num_src = 0;
  
  //write the peptide struct
  fwrite(peptide, sizeof(PEPTIDE_T), 1, file);
  
  //store number of src location in file
  num_src_location = ftell(file);

  //write dummie peptide src count
  fwrite(&num_src, sizeof(int), 1, file);

  //there must be at least one peptide src
  if(!peptide_src_iterator_has_next(iterator)){
    carp(CARP_WARNING, "no peptide src");
    return FALSE;
  }

  //interate through the peptide src for this peptide, serialize each peptide src
  while(peptide_src_iterator_has_next(iterator)){
    peptide_src = peptide_src_iterator_next(iterator);

    //serialize the peptide src
    serialize_peptide_src(peptide_src, file);
    
    ++num_src;
  }
  free_peptide_src_iterator(iterator);
  
  original_location = ftell(file);
  fseek(file, num_src_location, SEEK_SET);
  
  //over write the dummie peptide_src count with real value
  fwrite(&num_src, sizeof(int), 1, file);
  
  //return to original poistion
  fseek(file, original_location, SEEK_SET);
  return TRUE;
}
/*
BOOLEAN_T serialize_peptide(
  PEPTIDE_T* peptide, ///< the peptide to serialize -in
  FILE* file, ///< the output file to serlize -out
  int num_digits ///< the number of floating point digits to record -in
  )
{
  PEPTIDE_SRC_ITERATOR_T* iterator = 
    new_peptide_src_iterator(peptide);
  long num_src_location;
  long original_location;
  int num_src = 0;
  
  //print mass of the peptide
  if(num_digits == 2){
    fprintf(file, "* %.2f\n", peptide->peptide_mass);
  }
  else{//default if not 2, it's 3
    fprintf(file, "* %.3f\n", peptide->peptide_mass);
  }

  //print length of the peptide
  fprintf(file,"%d\n",peptide->length);
  //print number of src
  num_src_location = ftell(file);
  fprintf(file, "%s\n", "      ");

  //there must be at least one peptide src
  if(!peptide_src_iterator_has_next(iterator)){
    carp(CARP_WARNING, "no peptide src");
    return FALSE;
  }

  //interate through the linklist of possible parent proteins, print each parent protein
  while(peptide_src_iterator_has_next(iterator)){
    PEPTIDE_SRC_T* peptide_src = peptide_src_iterator_next(iterator);
    fprintf(file, "\t%d\n", get_peptide_src_peptide_type(peptide_src));
    fprintf(file, "\t%d\n", get_peptide_src_start_idx(peptide_src));
    fprintf(file, "\t%d\n", get_protein_protein_idx(get_peptide_src_parent_protein(peptide_src)));
    ++num_src;
  }
  free_peptide_src_iterator(iterator);
  
  original_location = ftell(file);
  fseek(file, num_src_location, SEEK_SET);
  fprintf(file, "%d", num_src);
  //update
  fseek(file, original_location, SEEK_SET);
  return TRUE;
}
*/

/*
 * Load a peptide from the FILE
 * \returns TRUE if load is successful, else FALSE
 *
 * See serialize_peptide above for the serialization format
 */
BOOLEAN_T load_peptide(
  PEPTIDE_T* peptide, ///< An allocated peptide
  FILE* file ///< The file pointing to the location of the peptide
  );
 

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

