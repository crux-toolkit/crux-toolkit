#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "check-peak.h"
#include "mass.h"
#include "objects.h"
#include "Spectrum.h"
#include "Peak.h"
#include "Peptide.h"
#include "PeptideSrc.h"
#include "Protein.h"
#include "ProteinPeptideIterator.h"
#include "Database.h"
#include "carp.h"
#include "crux-utils.h"
#include "parameter.h"
using namespace Crux; 
// also in parameter.c
void parse_custom_enzyme(const char* rule_str);

// declare things to setup
static Peptide *peptide1, *peptide2;
static Protein *protein1;
static ProteinPeptideIterator* pp_iterator;
static PeptideConstraint *constraint, *enzyme_constraint;

static PeptideSrc* src;

static Database* db;
//                      0123456789
static char prot1_sequence[256] = "MRVLKFGGTSVANAERFLRVADILESNARQGQVAOOTVLSAPAKITNHLVAMIEKTISGQDALPNISDAERIFAELLTGLAAAQPGFPLAQLKTFWVDQEFAQIKHVLHGISLWLGQC";

void protein_setup(){

  initialize_parameters();
  // a protein must have a database as its source of sequences
  db = new Database("input-data/protein1.fasta", false); //not mem mapped
  db->parse();  // assuming we have already tested database

  protein1 = new Protein("protein1", prot1_sequence, strlen(prot1_sequence),
                         NULL, 0, 0, db);//annotation, file offset, index

  constraint = PeptideConstraint::newFromParameters();
  pp_iterator = new ProteinPeptideIterator(protein1, constraint);

  // For testing enzymes 
  // remove mass constraints, set len 4-100, no missed cleavages
  enzyme_constraint = new PeptideConstraint(TRYPSIN, FULL_DIGEST,
                                             1,70000, 4, 100, 0, AVERAGE);

}

void protein_teardown(){
  Database::freeDatabase(db);
  delete protein1;
  PeptideConstraint::free(constraint);
  delete pp_iterator;
}

START_TEST (test_new){
  fail_unless(protein1 != NULL, "Failed to create protein1.");
  // check getters
}
END_TEST

START_TEST (test_peptide_iterator){
  fail_unless(pp_iterator != NULL, "Failed to create peptide iterator.");
  fail_unless(pp_iterator->hasNext() == true,
              "New iterator should have first peptide.");
  peptide1 = pp_iterator->next();
  fail_unless(peptide1 != NULL, "Iterator failed to return first peptide.");

  // check sequence of all seven peptides
  char* seq = peptide1->getSequence(); 
  fail_unless( strcmp(seq, "FGGTSVANAER") == 0,
               "First peptide should be FGGTSVANAER but is %s.", seq);
  peptide1 = pp_iterator->next();
  seq = peptide1->getSequence(); 
  fail_unless( strcmp(seq, "VADILESNAR") == 0,
               "Second peptide should be VADILESNAR but is %s.", seq);
  peptide1 = pp_iterator->next();
  seq = peptide1->getSequence(); 
  fail_unless( strcmp(seq, "QGQVAOOTVLSAPAK") == 0,
               "Third peptide should be QGQVAOOTVLSAPAK but is %s.", seq);
  peptide1 = pp_iterator->next();
  seq = peptide1->getSequence(); 
  fail_unless( strcmp(seq, "ITNHLVAMIEK") == 0,
               "Fourth peptide should be ITNHLVAMIEK but is %s.", seq);
  peptide1 = pp_iterator->next();
  seq = peptide1->getSequence(); 
  fail_unless( strcmp(seq, "TISGQDALPNISDAER") == 0,
               "Fifth peptide should be TISGQDALPNISDAER but is %s.", seq);
  peptide1 = pp_iterator->next();
  seq = peptide1->getSequence(); 
  fail_unless( strcmp(seq, "IFAELLTGLAAAQPGFPLAQLK") == 0,
               "Sixth peptide should be IFAELLTGLAAAQPGFPLAQLK but is %s.",
               seq);
  peptide1 = pp_iterator->next();
  seq = peptide1->getSequence(); 
  fail_unless( strcmp(seq, "TFWVDQEFAQIK") == 0,
               "Seventh peptide should be TFWVDQEFAQIK but is %s.", seq);

  peptide1 = pp_iterator->next();
  seq = peptide1->getSequence(); 
  fail_unless( strcmp(seq, "HVLHGISLWLGQC") == 0,
               "Eighth peptide should be HVLHGISLWLGQC but is %s.", seq);

  // now there should be no more peptides
  fail_unless( pp_iterator->hasNext() == false,
               "Default iterator should not have more than 7 peptides.");

}
END_TEST

START_TEST (test_elastase){
  // create list of start indexes for protein
  int starts[128];
  int s_idx = 0;
  int idx = 0;
  for(idx=0; idx < (int)strlen(prot1_sequence)-1; idx++){
    if((prot1_sequence[idx] == 'A' || 
        prot1_sequence[idx] == 'L' ||
        prot1_sequence[idx] == 'I' || 
        prot1_sequence[idx] == 'V') 
       && prot1_sequence[idx+1] != 'P'){
      starts[s_idx] = idx+1;
      s_idx++;
    }
  }
  // remove those for which length is too short
  for(idx=1; idx < s_idx; idx++){
    if( starts[idx] - starts[idx-1] < 4){
      starts[idx-1] = -1;
    }
  }

  // create peptide constraint with elastase
  enzyme_constraint->setEnzyme(ELASTASE);

  // create new iterator
  delete pp_iterator;
  pp_iterator = new ProteinPeptideIterator(protein1, enzyme_constraint);

  fail_unless( pp_iterator != NULL, "Failed to create elastase peptide iterator.");

  idx = 0;
  // for each peptide 
  while(pp_iterator->hasNext()){
    // get peptide
    peptide2 = pp_iterator->next();
    fail_unless(peptide2 != NULL, "Failed to get peptide from iterator.");
    char* seq = peptide2->getSequence();

    // get its source
    src = peptide2->getPeptideSrc();
    fail_unless(src != NULL, "Failed to get peptide src from peptide.");

    // get peptide indecies in sequence
    int start_idx = src->getStartIdx();// this is 1-based
    start_idx--;
    int len = peptide2->getLength();
    while( starts[idx] == -1 ){ idx++;}
    fail_unless( start_idx == starts[idx],
                 "Peptide %s should start at %i but starts at %i",
                 seq, starts[idx], start_idx );
    idx++;

    // if at start of protein, skip first end check
    if( start_idx > 0 ){
      fail_unless( (prot1_sequence[start_idx-1] == 'A' ||
                    prot1_sequence[start_idx-1] == 'L' ||
                    prot1_sequence[start_idx-1] == 'I' ||
                    prot1_sequence[start_idx-1] == 'V'),
                   "Peptide %s is flanked by %c but should be ALIV.",
                   seq, prot1_sequence[start_idx-1]); 
      fail_unless( prot1_sequence[start_idx] != 'P',
                   "Peptide %s cleaved before P and should not have.",
                   seq); 
    }
    // if at end of protein, skip second end check
    if( start_idx+len < (int)strlen(prot1_sequence) ){
      fail_unless( (prot1_sequence[start_idx+len-1] == 'A' ||
                    prot1_sequence[start_idx+len-1] == 'L' ||
                    prot1_sequence[start_idx+len-1] == 'I' ||
                    prot1_sequence[start_idx+len-1] == 'V'),
                   "Peptide %s ends in %c but should be ALIV.",
                   seq, prot1_sequence[start_idx+len]); 
      fail_unless( prot1_sequence[start_idx+len] != 'P',
                   "Peptide %s cleaved before P and should not have.",
                   seq); 
    }

    delete peptide2;
    free(seq);
  }
}
END_TEST

START_TEST (test_chymo){
  // create list of start indexes for protein
  int starts[128];
  starts[0] = 0;
  int s_idx = 1;
  int idx = 1;
  for(idx=1; idx < (int)strlen(prot1_sequence)-1; idx++){
    if((prot1_sequence[idx] == 'F' || 
        prot1_sequence[idx] == 'W' ||
        prot1_sequence[idx] == 'Y') 
       && prot1_sequence[idx+1] != 'P'){
      starts[s_idx] = idx+1;
      s_idx++;
    }
  }
  // remove those for which length is too short
  for(idx=1; idx < s_idx; idx++){
    //printf("index %i - %i is length %i\n", starts[idx], starts[idx-1], starts[idx]-starts[idx-1]);
    if( starts[idx] - starts[idx-1] < 4){
      //printf("removing index %i\n", starts[idx-1]);
      starts[idx-1] = -1;
    }
  }

  // create peptide constraint with elastase
  enzyme_constraint->setEnzyme(CHYMOTRYPSIN);

  // create new iterator
  delete pp_iterator;
  pp_iterator = new ProteinPeptideIterator(protein1, enzyme_constraint);

  fail_unless( pp_iterator != NULL, 
               "Failed to create chymotrypsin peptide iterator.");

  idx = 0;
  // for each peptide 
  while(pp_iterator->hasNext()){
    // get peptide
    peptide2 = pp_iterator->next();
    fail_unless(peptide2 != NULL, "Failed to get peptide from iterator.");
    char* seq = peptide2->getSequence();

    // get its source
    src = peptide2->getPeptideSrc();
    fail_unless(src != NULL, "Failed to get peptide src from peptide.");

    // get peptide indecies in sequence
    int start_idx = src->getStartIdx();// this is 1-based
    start_idx--;
    int len = peptide2->getLength();
    while( starts[idx] == -1 ){ idx++;}
    fail_unless( start_idx == starts[idx],
                 "Peptide %s should start at %i but starts at %i",
                 seq, starts[idx], start_idx );
    idx++;

    // if at start of protein, skip first end check
    if( start_idx > 0 ){
      fail_unless( (prot1_sequence[start_idx-1] == 'F' ||
                    prot1_sequence[start_idx-1] == 'W' ||
                    prot1_sequence[start_idx-1] == 'Y'),
                   "Peptide %s is flanked by %c but should be FWY.",
                   seq, prot1_sequence[start_idx-1]); 
      fail_unless( prot1_sequence[start_idx] != 'P',
                   "Peptide %s cleaved before P and should not have.",
                   seq); 
    }
    // if at end of protein, skip second end check
    if( start_idx+len < (int)strlen(prot1_sequence) ){
      fail_unless( (prot1_sequence[start_idx+len-1] == 'F' ||
                    prot1_sequence[start_idx+len-1] == 'W' ||
                    prot1_sequence[start_idx+len-1] == 'Y'),
                   "Peptide %s ends in %c but should be FWY.",
                   seq, prot1_sequence[start_idx+len]); 
      fail_unless( prot1_sequence[start_idx+len] != 'P',
                   "Peptide %s cleaved before P and should not have.",
                   seq); 
    }

    delete peptide2;
    free(seq);
  }
}
END_TEST

START_TEST (test_mod_chymo){
  // create list of start indexes for protein
  int starts[128];
  starts[0] = 0;
  int s_idx = 1;
  int idx = 1;
  for(idx=1; idx < (int)strlen(prot1_sequence)-1; idx++){
    if((prot1_sequence[idx] == 'F' || 
        prot1_sequence[idx] == 'W' ||
        prot1_sequence[idx] == 'Y' ||
        prot1_sequence[idx] == 'L') 
       && prot1_sequence[idx+1] != 'P'){
      starts[s_idx] = idx+1;
      s_idx++;
    }
  }
  // remove those for which length is too short
  for(idx=1; idx < s_idx; idx++){
    //printf("index %i - %i is length %i\n", starts[idx], starts[idx-1], starts[idx]-starts[idx-1]);
    if( starts[idx] - starts[idx-1] < 4){
      //printf("removing index %i\n", starts[idx-1]);
      starts[idx-1] = -1;
    }
  }

  // create peptide constraint with elastase
  enzyme_constraint->setEnzyme(MODIFIED_CHYMOTRYPSIN);

  // create new iterator
  delete pp_iterator;
  pp_iterator = new ProteinPeptideIterator(protein1, enzyme_constraint);

  fail_unless( pp_iterator != NULL, 
               "Failed to create chymotrypsin peptide iterator.");

  idx = 0;
  // for each peptide 
  while(pp_iterator->hasNext()){
    // get peptide
    peptide2 = pp_iterator->next();
    fail_unless(peptide2 != NULL, "Failed to get peptide from iterator.");
    char* seq = peptide2->getSequence();

    // get its source
    src = peptide2->getPeptideSrc();
    fail_unless(src != NULL, "Failed to get peptide src from peptide.");

    // get peptide indecies in sequence
    int start_idx = src->getStartIdx();// this is 1-based
    start_idx--;
    int len = peptide2->getLength();
    while( starts[idx] == -1 ){ idx++;}
    fail_unless( start_idx == starts[idx],
                 "Peptide %s should start at %i but starts at %i",
                 seq, starts[idx], start_idx );
    idx++;

    // if at start of protein, skip first end check
    if( start_idx > 0 ){
      fail_unless( (prot1_sequence[start_idx-1] == 'F' ||
                    prot1_sequence[start_idx-1] == 'W' ||
                    prot1_sequence[start_idx-1] == 'Y' ||
                    prot1_sequence[start_idx-1] == 'L'),
                   "Peptide %s is flanked by %c but should be FWYL.",
                   seq, prot1_sequence[start_idx-1]); 
      fail_unless( prot1_sequence[start_idx] != 'P',
                   "Peptide %s cleaved before P and should not have.",
                   seq); 
    }
    // if at end of protein, skip second end check
    if( start_idx+len < (int)strlen(prot1_sequence) ){
      fail_unless( (prot1_sequence[start_idx+len-1] == 'F' ||
                    prot1_sequence[start_idx+len-1] == 'W' ||
                    prot1_sequence[start_idx+len-1] == 'Y' ||
                    prot1_sequence[start_idx+len-1] == 'L'),
                   "Peptide %s ends in %c but should be FWYL.",
                   seq, prot1_sequence[start_idx+len-1]); 
      fail_unless( prot1_sequence[start_idx+len] != 'P',
                   "Peptide %s cleaved before P and should not have.",
                   seq); 
    }

    delete peptide2;
    free(seq);
  }
}
END_TEST

START_TEST (test_el_tryp_chymo){
  // create list of start indexes for protein
  int starts[128];
  starts[0] = 0;
  int s_idx = 1;
  int idx = 1;
  for(idx=1; idx < (int)strlen(prot1_sequence)-1; idx++){
    if((prot1_sequence[idx] == 'A' || 
        prot1_sequence[idx] == 'L' ||
        prot1_sequence[idx] == 'I' ||
        prot1_sequence[idx] == 'V' ||
        prot1_sequence[idx] == 'K' ||
        prot1_sequence[idx] == 'R' ||
        prot1_sequence[idx] == 'W' ||
        prot1_sequence[idx] == 'F' ||
        prot1_sequence[idx] == 'Y') 
       && prot1_sequence[idx+1] != 'P'){
      starts[s_idx] = idx+1;
      s_idx++;
    }
  }
  // remove those for which length is too short
  for(idx=1; idx < s_idx; idx++){
    //printf("index %i - %i is length %i\n", starts[idx], starts[idx-1], starts[idx]-starts[idx-1]);
    if( starts[idx] - starts[idx-1] < 4){
      //printf("removing index %i\n", starts[idx-1]);
      starts[idx-1] = -1;
    }
  }

  // create peptide constraint with elastase
  enzyme_constraint->setEnzyme(ELASTASE_TRYPSIN_CHYMOTRYPSIN);

  // create new iterator
  delete pp_iterator;
  pp_iterator = new ProteinPeptideIterator(protein1, enzyme_constraint);

  fail_unless( pp_iterator != NULL, 
               "Failed to create chymotrypsin peptide iterator.");

  idx = 0;
  // for each peptide 
  while(pp_iterator->hasNext()){
    // get peptide
    peptide2 = pp_iterator->next();
    fail_unless(peptide2 != NULL, "Failed to get peptide from iterator.");
    char* seq = peptide2->getSequence();

    // get its source
    src = peptide2->getPeptideSrc();
    fail_unless(src != NULL, "Failed to get peptide src from peptide.");

    // get peptide indecies in sequence
    int start_idx = src->getStartIdx();// this is 1-based
    start_idx--;
    int len = peptide2->getLength();
    while( starts[idx] == -1 ){ idx++;}
    fail_unless( start_idx == starts[idx],
                 "Peptide %s should start at %i but starts at %i",
                 seq, starts[idx], start_idx );
    idx++;

    // if at start of protein, skip first end check
    if( start_idx > 0 ){
      fail_unless( (prot1_sequence[start_idx-1] == 'A' || 
                    prot1_sequence[start_idx-1] == 'L' ||
                    prot1_sequence[start_idx-1] == 'I' ||
                    prot1_sequence[start_idx-1] == 'V' ||
                    prot1_sequence[start_idx-1] == 'K' ||
                    prot1_sequence[start_idx-1] == 'R' ||
                    prot1_sequence[start_idx-1] == 'W' ||
                    prot1_sequence[start_idx-1] == 'F' ||
                    prot1_sequence[start_idx-1] == 'Y') ,
                   "Peptide %s is flanked by %c but should be ALIVKRWFY.",
                   seq, prot1_sequence[start_idx-1]); 
      fail_unless( prot1_sequence[start_idx] != 'P',
                   "Peptide %s cleaved before P and should not have.",
                   seq); 
    }
    // if at end of protein, skip second end check
    if( start_idx+len < (int)strlen(prot1_sequence) ){
      fail_unless( (prot1_sequence[start_idx+len-1] == 'A' || 
                    prot1_sequence[start_idx+len-1] == 'L' ||
                    prot1_sequence[start_idx+len-1] == 'I' ||
                    prot1_sequence[start_idx+len-1] == 'V' ||
                    prot1_sequence[start_idx+len-1] == 'K' ||
                    prot1_sequence[start_idx+len-1] == 'R' ||
                    prot1_sequence[start_idx+len-1] == 'W' ||
                    prot1_sequence[start_idx+len-1] == 'F' ||
                    prot1_sequence[start_idx+len-1] == 'Y') ,
                   "Peptide %s ends in %c but should be ALIVKRWFY.",
                   seq, prot1_sequence[start_idx+len-1]); 
      fail_unless( prot1_sequence[start_idx+len] != 'P',
                   "Peptide %s cleaved before P and should not have.",
                   seq); 
    }

    delete peptide2;
    free(seq);
  }
}
END_TEST

START_TEST(test_aspn){
  // create list of start indexes for protein
  int starts[128];
  starts[0] = 0;
  int s_idx = 1;
  int idx = 1;
  for(idx=1; idx < (int)strlen(prot1_sequence)-1; idx++){
    if(prot1_sequence[idx] == 'D' ){
      starts[s_idx] = idx;
      s_idx++;
    }
  }
  // remove those for which length is too short
  for(idx=1; idx < s_idx; idx++){
    //printf("index %i - %i is length %i\n", starts[idx], starts[idx-1], starts[idx]-starts[idx-1]);
    if( starts[idx] - starts[idx-1] < 4){
      //printf("removing index %i\n", starts[idx-1]);
      starts[idx-1] = -1;
    }
  }

  // create peptide constraint with elastase
  enzyme_constraint->setEnzyme(ASPN);

  // create new iterator
  delete pp_iterator;
  pp_iterator = new ProteinPeptideIterator(protein1, enzyme_constraint);

  fail_unless( pp_iterator != NULL, 
               "Failed to create aspn peptide iterator.");

  idx = 0;
  // for each peptide 
  while(pp_iterator->hasNext()){
    // get peptide
    peptide2 = pp_iterator->next();
    fail_unless(peptide2 != NULL, "Failed to get peptide from iterator.");
    char* seq = peptide2->getSequence();

    // get its source
    src = peptide2->getPeptideSrc();
    fail_unless(src != NULL, "Failed to get peptide src from peptide.");

    // get peptide indecies in sequence
    int start_idx = src->getStartIdx();// this is 1-based
    start_idx--;
    int len = peptide2->getLength();
    while( starts[idx] == -1 ){ idx++;}
    fail_unless( start_idx == starts[idx],
                 "Peptide %s should start at %i but starts at %i",
                 seq, starts[idx], start_idx );
    idx++;

    // if at start of protein, skip first end check
    if( start_idx > 0 ){
      fail_unless( (prot1_sequence[start_idx] == 'D' ),
                   "Peptide %s starts with %c but should be D.",
                   seq, prot1_sequence[start_idx]); 
    }
    // if at end of protein, skip second end check
    if( start_idx+len < (int)strlen(prot1_sequence) ){
      fail_unless( prot1_sequence[start_idx+len+1] != 'D',
                   "Peptide %s did not cleave before D and should have.",
                   seq); 
    }

    delete peptide2;
    free(seq);
  }
}
END_TEST

START_TEST(test_other_enzymes){
  int num_enzymes = 5;
  ENZYME_T enzymes[5] = { CLOSTRIPAIN, CYANOGEN_BROMIDE,
                          IODOSOBENZOATE, PROLINE_ENDOPEPTIDASE,
                          STAPH_PROTEASE };
  char residues[5] = { 'R', 'M', 'W', 'P', 'E'};

 int idx = 0;
 for(idx = 0; idx < num_enzymes; idx ++){

   char residue = residues[idx];
   ENZYME_T enzyme = enzymes[idx];
   char* enzyme_name = enzyme_type_to_string(enzyme);

   //printf("enzyme: %s\n", enzyme_name);
   // create list of start indexes for protein
   int starts[128];
   starts[0] = 0;
   int s_idx = 1;
   int idx = 0;
   for(idx=0; idx < (int)strlen(prot1_sequence)-1; idx++){
     if(prot1_sequence[idx] == residue ){
       starts[s_idx] = idx+1;
       s_idx++;
     }
   }
   // remove those for which length is too short
   for(idx=1; idx < s_idx; idx++){
     //printf("index %i - %i is length %i\n", starts[idx], starts[idx-1], starts[idx]-starts[idx-1]);
     if( starts[idx] - starts[idx-1] < 4){
       //printf("removing index %i\n", starts[idx-1]);
       starts[idx-1] = -1;
     }
   }
   
   // create peptide constraint with elastase
  enzyme_constraint->setEnzyme(enzyme);
   
   // create new iterator
   delete pp_iterator;
   pp_iterator = new ProteinPeptideIterator(protein1, enzyme_constraint);
   
   fail_unless( pp_iterator != NULL, 
                "Failed to create %s peptide iterator.", enzyme_name);
   
   idx = 0;
   // for each peptide 
   while(pp_iterator->hasNext()){
     // get peptide
     peptide2 = pp_iterator->next();
     fail_unless(peptide2 != NULL, "Failed to get peptide from iterator.");
     char* seq = peptide2->getSequence();
     
     // get its source
     src = peptide2->getPeptideSrc();
     fail_unless(src != NULL, "Failed to get peptide src from peptide.");
     
     // get peptide indecies in sequence
     int start_idx = src->getStartIdx();// this is 1-based
     start_idx--;
     int len = peptide2->getLength();
     while( starts[idx] == -1 ){ idx++;}
     fail_unless( start_idx == starts[idx],
                  "%s peptide %s should start at %i but starts at %i",
                  enzyme_name, seq, starts[idx], start_idx );
     idx++;
     
     // if at start of protein, skip first end check
     if( start_idx > 0 ){
       fail_unless( (prot1_sequence[start_idx-1] == residue ),
                    "%s peptide %s flanked by %c but should be %c.",
                    enzyme_name, seq, prot1_sequence[start_idx-1], residue); 
     }
     // if at end of protein, skip second end check
     if( start_idx+len < (int)strlen(prot1_sequence) ){
       fail_unless( prot1_sequence[start_idx+len-1] == residue,
                    "%s peptide %s does not end in %c and should.",
                    enzyme_name, seq, residue); 
     }
     
     delete peptide2;
     free(seq);
   }
   free(enzyme_name);
 }// next enzyme
 
}
END_TEST

START_TEST(test_custom_enzyme){
  // define custom rule
  const char* rule = "[MLK]|{P}";
  parse_custom_enzyme(rule);

  // create list of start indexes
  int starts[128];
  starts[0] = 0;
  int s_idx = 1;
  int idx = 0;
  for(idx=0; idx < (int)strlen(prot1_sequence)-1; idx++){
    //    printf("seq[%i] = %c\n", );
   if((prot1_sequence[idx] == 'M' || 
        prot1_sequence[idx] == 'L' ||
        prot1_sequence[idx] == 'K') 
       && prot1_sequence[idx+1] != 'P'){
      starts[s_idx] = idx+1;
      s_idx++;
    }  
  }

  // remove those for which length is too short
  for(idx=1; idx < s_idx; idx++){
    if( starts[idx] - starts[idx-1] < 4){
      starts[idx-1] = -1;
    }
  }

  // create peptide constraint with elastase
  enzyme_constraint->setEnzyme(CUSTOM_ENZYME);

  // create new iterator
  delete pp_iterator;
  pp_iterator = new ProteinPeptideIterator(protein1, enzyme_constraint);

  fail_unless( pp_iterator != NULL, 
               "Failed to create custom enzyme peptide iterator.");

  idx=0;
  // for each peptide 
  while(pp_iterator->hasNext()){
    // get peptide
    peptide2 = pp_iterator->next();
    fail_unless(peptide2 != NULL, "Failed to get peptide from iterator.");
    char* seq = peptide2->getSequence();

    // get its source
    src = peptide2->getPeptideSrc();
    fail_unless(src != NULL, "Failed to get peptide src from peptide.");

    // get peptide indecies in sequence
    int start_idx = src->getStartIdx();// this is 1-based
    start_idx--;
    int len = peptide2->getLength();
    while( starts[idx] == -1 ){ idx++;}
    fail_unless( start_idx == starts[idx],
                 "Peptide %s should start at %i but starts at %i",
                 seq, starts[idx], start_idx );
    idx++;

    // if at start of protein, skip first end check
    if( start_idx > 0 ){
      fail_unless( (prot1_sequence[start_idx-1] == 'M' || 
                    prot1_sequence[start_idx-1] == 'L' ||
                    prot1_sequence[start_idx-1] == 'K') ,
                   "Peptide %s is flanked by %c but should be MLK.",
                   seq, prot1_sequence[start_idx-1]); 
      fail_unless( prot1_sequence[start_idx] != 'P',
                   "Peptide %s cleaved before P and should not have.",
                   seq); 
    }
    // if at end of protein, skip second end check
    if( start_idx+len < (int)strlen(prot1_sequence) ){
      fail_unless( (prot1_sequence[start_idx+len-1] == 'M' || 
                    prot1_sequence[start_idx+len-1] == 'L' ||
                    prot1_sequence[start_idx+len-1] == 'K') ,
                   "Peptide %s ends in %c but should be MLK.",
                   seq, prot1_sequence[start_idx+len-1]); 
      fail_unless( prot1_sequence[start_idx+len] != 'P',
                   "Peptide %s cleaved before P and should not have.",
                   seq); 
    }

    delete peptide2;
    free(seq);
  }

  // define custom rule
  rule = "[X]|[Q]";
  free(pre_cleavage_list);
  free(post_cleavage_list);
  pre_list_size = 0;
  post_list_size = 0;
  parse_custom_enzyme(rule);

  // create list of start indexes
  starts[0] = 0;
  s_idx = 1;
  idx = 0;
  for(idx=0; idx < (int)strlen(prot1_sequence)-1; idx++){
   if( prot1_sequence[idx+1] == 'Q'){
      starts[s_idx] = idx+1;
      s_idx++;
    }  
  }

  // remove those for which length is too short
  for(idx=1; idx < s_idx; idx++){
    //printf("index %i - %i is length %i\n", starts[idx], starts[idx-1], starts[idx]-starts[idx-1]);
    if( starts[idx] - starts[idx-1] < 4){
      //printf("removing index %i\n", starts[idx-1]);
      starts[idx-1] = -1;
    }
  }

  // create peptide constraint with elastase
  enzyme_constraint->setEnzyme(CUSTOM_ENZYME);

  // create new iterator
  delete pp_iterator;
  pp_iterator = new ProteinPeptideIterator(protein1, enzyme_constraint);

  fail_unless( pp_iterator != NULL, 
               "Failed to create custom enzyme peptide iterator.");

  idx=0;
  // for each peptide 
  while(pp_iterator->hasNext()){
    // get peptide
    peptide2 = pp_iterator->next();
    fail_unless(peptide2 != NULL, "Failed to get peptide from iterator.");
    char* seq = peptide2->getSequence();

    // get its source
    src = peptide2->getPeptideSrc();
    fail_unless(src != NULL, "Failed to get peptide src from peptide.");

    // get peptide indecies in sequence
    int start_idx = src->getStartIdx();// this is 1-based
    start_idx--;
    int len = peptide2->getLength();
    while( starts[idx] == -1 ){ idx++;}
    fail_unless( start_idx == starts[idx],
                 "Peptide %s should start at %i but starts at %i",
                 seq, starts[idx], start_idx );
    idx++;

    if( start_idx > 0 ){
    // if at start of protein, skip first end check
    fail_unless( prot1_sequence[start_idx] == 'Q',
                   "Peptide %s should have cleaved before Q and did not.",
                   seq);
    } 

    // if at end of protein, skip second end check
    if( start_idx+len < (int)strlen(prot1_sequence) ){
      fail_unless( prot1_sequence[start_idx+len] == 'Q',
                   "Peptide %s should have cleaved before Q and did not.",
                   seq); 
    }

    delete peptide2;
    free(seq);
  }

}
END_TEST

Suite* protein_suite(void){
  Suite *s = suite_create("Protein");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_new);
  tcase_add_test(tc_core, test_peptide_iterator);
  tcase_add_test(tc_core, test_elastase);
  tcase_add_test(tc_core, test_chymo);
  tcase_add_test(tc_core, test_mod_chymo);
  tcase_add_test(tc_core, test_el_tryp_chymo);
  tcase_add_test(tc_core, test_aspn);
  tcase_add_test(tc_core, test_other_enzymes);
  tcase_add_test(tc_core, test_custom_enzyme);

  //  tcase_add_test(tc_core, test_create);
  tcase_add_checked_fixture(tc_core, protein_setup, protein_teardown);
  return s;
}

/*
START_TEST (test_create){
  char* name = NULL;
  char* name2 = NULL;
  
  //try create a new database
  db = new_database("fasta_file_binary_fasta", false);
  fail_unless(parse_database(db), "failed to parse database");
  fail_unless(strncmp((name = get_database_filename(db)), "fasta_file", 10) == 0, "database filename not set correctly");
  free(name);
  fail_unless(get_database_is_parsed(db), "database parsed field not correctly set");
  fail_unless(get_database_num_proteins(db) == 3, "database number of proteins not set correctly");
  
  //check the fasta file parsing
  FILE* file = fopen("fasta_file", "r");
  protein3 = allocate_protein();
  fail_unless(parse_protein_fasta_file(protein3, file), "failed to parse protein from fasta file");
  set_protein_database(protein3, db);
  print_protein(protein3, stdout);
  
 
  protein1 = new_protein("23 Jordan", "AADAAKAGAAKFFA", 14, "this is a test protein", 45, 3, db);
  print_protein(protein1, stdout);
 
  //try copy protein
  protein2 = allocate_protein();
  copy_protein(protein1, protein2);
  fail_unless(strcmp((name = get_protein_id(protein1)), (name2=get_protein_id(protein2))) == 0, "protein, id not correct");
  free(name); 
  free(name2);
  fail_unless(strcmp((name = get_protein_sequence(protein1)), (name2 = get_protein_sequence(protein2))) == 0, "protein, sequence not correct");
  free(name); 
  free(name2);
  fail_unless(strcmp((name = get_protein_annotation(protein1)), (name2 = get_protein_annotation(protein2))) == 0, "protein, annotation not correct");
  free(name); 
  free(name2);
  fail_unless(get_protein_length(protein1) == get_protein_length(protein2), "protein, protein length not correct");
  fail_unless(get_protein_offset(protein1) == get_protein_offset(protein2), "protein, protein offset not correct");
  fail_unless(get_protein_protein_idx(protein1) == get_protein_protein_idx(protein2), "protein, protein idx not correct");
  fail_unless(get_protein_is_light(protein1) == get_protein_is_light(protein2), "protein, protein is_light not correct");
  
  //peptide constraint
  constraint = new_peptide_constraint(TRYPSIN, FULL_DIGEST, 0, 
                                      1200, 1, 10, 1, AVERAGE);
  
  // ** test, protein_peptide_iterator ** /
  
  //create iterator
  iterator = new_protein_peptide_iterator(protein3, constraint);
  
  //iterate over all possible peptides
  while(protein_peptide_iterator_has_next(iterator)){
    peptide1 = protein_peptide_iterator_next(iterator);
    //print_peptide(peptide1, stdout);
    free_peptide(peptide1);
  }  

  //  print_peptide_in_format(peptide5, true,  stdout);

  //free stuff
  free_protein(protein1);
  free_protein(protein2);
  free_protein(protein3);
  free_protein_peptide_iterator(iterator);
  free_peptide_constraint(constraint);
  free_database(db);
  fclose(file);
}
END_TEST
*/
