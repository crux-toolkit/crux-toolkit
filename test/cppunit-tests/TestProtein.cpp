#include <cppunit/config/SourcePrefix.h>
#include "TestProtein.h"
#include "parameter.h" 
//#include "objects.h"
using namespace Crux;
CPPUNIT_TEST_SUITE_REGISTRATION( TestProtein );

void TestProtein::setUp(){
  initialize_parameters();  // accessing any parameter values requires this

  // initialize variables to use for testing
  tryptic_seq = "MBEFOKFIRPSTKSECONDREND";
  prot1 = new Protein("tryptic_seq",
                      tryptic_seq,
                      strlen(tryptic_seq),
                      "",
                      0, ///< file location in the database
                      0, ///< The index of the protein in it's database.-in  
                      NULL); ///< the database of its origin

}

void TestProtein::tearDown(){
  // delete anything you allocated
  delete prot1;
}

void TestProtein::testPeptideShuffle(){
  prot1->shuffle(PEPTIDE_SHUFFLE_DECOYS);
  const char* new_seq = prot1->getSequence();

  // confirm that the cleavage residues haven't moved
  int no_move[8] = {0, 5, 6, 12, 13, 19, 20, 22};
  for(int i = 0; i < 8; i++){
    CPPUNIT_ASSERT_EQUAL_MESSAGE("Shuffled sequence should have same residue",
                                 tryptic_seq[no_move[i]], new_seq[no_move[i]]);
  }
}










