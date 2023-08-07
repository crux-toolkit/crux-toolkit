#include <cppunit/config/SourcePrefix.h>
#include <stdlib.h>
#include <map>
#include "TestXml.h"
#include "Match.h"
#include "Peptide.h"
#include "modifications.h"
#include "parameter.h"


using namespace std;
bool set_double_parameter(
 const char*     name,  ///< the name of the parameter looking for -in
 double set_value,  ///< the value to be set -in
 double min_value,  ///< the value to be set -in
 double max_value,  ///< the value to be set -in
 const char* usage,      ///< string to print in usage statement
 const char* filenotes,  ///< additional info for param file
 const char* foruser     ///< "true" if should be revealed to user
			       );


CPPUNIT_TEST_SUITE_REGISTRATION( TestXml );

void TestXml::setUp(){
  initialize_parameters();
  set_double_parameter((const char*) "V", 30, 30, 30, "", "", "");
  set_double_parameter((const char*) "P", 40, 40, 40, "", "", "");
  isotopic_type = get_mass_type_parameter("isotopic-mass");
  mass_v = get_mass_amino_acid('V', isotopic_type);
  mass_p = get_mass_amino_acid('P', isotopic_type);
  ord_pep_seq = "VGGAGK"; //ordinary peptide sequence
}

void TestXml::tearDown(){
  var_mods.clear();
  static_mods.clear();
}


/*
 * Checks that the number of internal cleavage should be 
 * zero for ordinary peptide sequence and Trypsin cut
 */
void TestXml::getNumInternalCleavageNone(){
  int num_missed_cleavage = 
    get_num_internal_cleavage((char *)ord_pep_seq.c_str(), TRYPSIN);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Number of internal cleavages is wrong.", 
			       0, num_missed_cleavage);
}

/*
 * Checks that number of internal cleavage or RVGGAGKA 
 * should be two since RV and KA are internal sites
 */
void TestXml::getNumInternalCleavageTwo(){
  string peptide_sequence = "RVGGAGKA"; //RV and KA
  int num_missed_cleavage = 
    get_num_internal_cleavage((char *)peptide_sequence.c_str(), TRYPSIN);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Number of internal cleavages is wrong.", 
			       2, num_missed_cleavage);
}

/*
 * Checks that Dash '-' is counted as terminal cleavage
 */
void TestXml::getNumTerminalCleavageTwoDash(){
  int num = 
    get_num_terminal_cleavage((char*)ord_pep_seq.c_str(), '-', '-', TRYPSIN);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Number of terminal cleavage is wrong.",
			       2, num);
}

/*
 * Checks that it recognizes non-terminal cleavage sites
 * at the beginning of the sequence
 */
void TestXml::getNumTerminalCleavageOnePrev(){
  int num = 
    get_num_terminal_cleavage((char*)ord_pep_seq.c_str(), 'R', 'P', TRYPSIN);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Number of terminal cleavage is wrong.", 
			       1, num);
}

/*
 * Checks that it recognizes non-terminal cleavage sites
 * at the end of the sequence
 */
void TestXml::getNumTerminalCleavageOneNext(){
  string peptide_sequence = "PGGAGK";
  int num = get_num_terminal_cleavage((char*)peptide_sequence.c_str(), 
				      'R', 'A', TRYPSIN);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Number of terminal cleavage is wrong.",
			       1, num);
}

/*
 * Makes sure that the number variable modifications found in
 * sequence is 0
 */
void TestXml::findVariableModificationsNone(){
  find_variable_modifications(var_mods, (char*)ord_pep_seq.c_str());
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Number of variable modifcations is wrong.",
			       0, (int)var_mods.size());
}

/*
 * Makes sure that the correct variable modifications are 
 * identified
 */
void TestXml::findVariableModificationsThree(){
  string mod_seq = "V[100.00]GGA[20.3]K[100.1]";
  find_variable_modifications(var_mods, (char*)mod_seq.c_str());
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Number of variable modifications is wrong.",
			       3, (int)var_mods.size());
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Incorrect variable modification information",
			       100.00, var_mods[1]);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Incorrect variable modification information",
			       20.3, var_mods[4]);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Incorrect variable modification information",
			       100.1, var_mods[5]);
}


/*
 * Makes sure that none of the static modifications are
 * identified because the modifications were already 
 * identified as variable 
 */
void TestXml::findStaticModificationsNoneFromVariable(){
  string mod_seq = "V[100.00]GGA[20.3]K[100.1]";
  find_variable_modifications(var_mods, (char*) mod_seq.c_str());
  find_static_modifications(static_mods, var_mods,  
			    (char*) ord_pep_seq.c_str());
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Number of static modifcations is wrong.", 
			       0, (int)static_mods.size());
}

/*
 * Makes sure that there are no static modifications
 * found in the the sequence RRRRR
 */
void TestXml::findStaticModificationsNone(){
  string mod_seq = "R[100.00]RRR[20.3]R[100.1]";
  string peptide_sequence = "RRRRR";
  find_variable_modifications(var_mods, (char*) mod_seq.c_str());
  find_static_modifications(static_mods, 
			    var_mods, (char*) peptide_sequence.c_str());
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Number of static modifcations is wrong.", 
			       0, (int)static_mods.size());
}


/*
 * Makes sure all static modifications are found
 * and identified correctly
 */
void TestXml::findStaticModificationsThree(){
  string mod_seq = "VGGPA[20.3]GKP";
  string peptide_sequence = "VGGPAGKP";
  find_variable_modifications(var_mods, (char*) mod_seq.c_str());
  find_static_modifications(static_mods, var_mods,  
			    (char*) peptide_sequence.c_str());
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Number of static modifcations is wrong.",
			       3, (int)static_mods.size());
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Incorrect static modification information", 
			       mass_v , static_mods[1]);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Incorrect static modification information", 
			       mass_p , static_mods[4]);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Incorrect static modification information", 
			       mass_p , static_mods[8]);
}

