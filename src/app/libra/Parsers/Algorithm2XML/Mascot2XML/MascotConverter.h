#ifndef MASC_CONV_H
#define MASC_CONV_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <sys/stat.h>

#include "Common/sysdepend.h"
#include "Common/Array.h"
#include "SpectrumPeak.h"
#include "Common/constants.h"
#include "Common/Enzyme/ProteolyticEnzyme/ProteolyticEnzymeFactory/ProteolyticEnzymeFactory.h"
#include "Parsers/mzParser/mzParser.h"
#include "MascotModification.h"
//#include "Common/ResidueMass/ResidueMass.h"
#include "Parsers/Parser/Tag.h"
#include "Parsers/Algorithm2XML/pICalculator.h"
#include "Common/ModificationInfo/ModificationInfo.h"
#include "Common/constants.h"
#include "Parsers/Parser/Parser.h"

#include <map>
#include <utility>
//proteinMap: short id --> protein description
typedef std::map< std::string, std::string, std::less<std::string> > proteinMap;
//idMap: short id --> protein full id
typedef std::map< std::string, std::string, std::less<std::string> > idMap;
typedef std::pair< long, int > cycleExperimentPair;
typedef std::map< cycleExperimentPair, long, std::less<cycleExperimentPair> > scanMap;
typedef std::map< long, double, std::less<long> > scanRTMap;
typedef std::map< std::string, char, std::less<std::string> > writtenDtaMap;


/*

Program       : Mascot2XML converter to pepXML
Author        : Andrew Keller <akeller@systemsbiology.org>                                                       
Date          : 11.27.02 


Copyright (C) 2003 Andrew Keller

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Andrew Keller
Institute for Systems Biology
401 Terry Avenue North 
Seattle, WA  98109  USA
akeller@systemsbiology.org

*/

int comp_specpks(const void* num1, const void* num2);

#include "Parsers/Parser/TagListComparator.h"
#include "Common/tpp_tarball.h" // for constructing tarball without all that disk IO

class MascotConverter {


 public:

  //  MascotConverter(char* filename, Boolean windows);
  //  MascotConverter(char* filename, Boolean windows, char* database);
  MascotConverter(char* filename, Boolean xml, char* database, char* sample_enzyme, Boolean calculatePI, 
       Boolean generateTgz, Boolean generateDta, Boolean generateDescription, Boolean shortProteinId, Boolean verbose, 
       eTagListFilePurpose testmode, std::string &testfilename);
  ~MascotConverter();
  void init(const char* filename, Boolean windows, Boolean generateTgz, Boolean generateDta);
  void databaseSearch();
  void writeHTML(const char* outfile);
  void write_pepXML(const char* outfile);
  bool writeOutfile(int index);
  bool writeDtafile(int index, const char* line, int line_len);
  void tarFiles();
  Boolean isStar(int index);
  Boolean homologous(char* pep1, char* pep2);
  void setFileInfo(const char* file);
  Boolean validate_db();
  Boolean checkResult(int index);
  char decodeURL(char first, char second);
  Array<Tag*>* getModificationTags(char* peptide, char* modstring);
  double getModifiedMass(char res, Boolean staticmod, int index);
  double getModifiedMass(char res, Boolean staticmod, int index, MascotModification*& mod);
  //Boolean enterModification(char* mod, double mass, double error, Boolean variable, int index);
  //Boolean enterModification(char* mod, double mass, double error, Boolean variable, int index, Boolean n_term, Boolean c_term);
  Array<char*> spectra_;
  Array<double> ionscores_;
  Array<double> homologyscores_;
  Array<double> identityscores_;
  Array<double> expectationvalues_;
  Array<char*> peptides_;
  Boolean term_read_;
  Boolean generate_description_;
  Boolean short_protein_id_;
  std::string proteinIDLeftAnchor_;
  std::string proteinIDRightAnchor_;
  proteinMap proteinDescriptionMap_;
  proteinMap missingDescriptionMap_;
  idMap proteinFullIdMap_;
  Boolean verbose_;
  scanMap cycleExptToScanMap_;
  scanRTMap scanRetentionTimeMap_;
  Boolean toLoadScanMap_; //flag to signal loading is required
  Boolean scanMapLoaded_; //the scanMap has been loaded successfully

  Array<char*> proteins_;
  Array<int> pluses_; // record number of additional proteins corresponding to peptide in db
  Array<char> prev_AAs_;
  Array<char> foll_AAs_;
  Array<double> masses_;
  Array<double> massdiffs_;
  Array<int> matchedions_;
  Array<MascotModification*> static_mods_;
  Array<MascotModification*> variable_mods_;

  //rank1_* store other proteins sharing hit#1 position
  Array<char*> rank1_peptides_;  //the peptide is just reference to Array<char*> peptides_
  Array<char*> rank1_proteins_;
  Array<char> rank1_prev_AAs_;
  Array<char> rank1_foll_AAs_;
  Array<int> rank1_inds_; // ref to start index
  
  //Array<char *> written_dtas_; // for tracking DTAs we've written already
  writtenDtaMap written_dtas_;
  tpp_tarball *tarball_; // for constructing tarball without all that disk IO

  char* database_;
  Boolean enzyme_constraint_;
  int mass_type_;
  Boolean calculate_pI_;

  Array<int> alt_inds_; // ref to start index
  Array<char*> alt_peptides_;
  Array<char*> alt_proteins_;
  Array<int> alt_pluses_;
  Array<double> alt_ionscores_;
  Array<char> alt_prev_AAs_;
  Array<char> alt_foll_AAs_;
  Array<double> alt_masses_;
  Array<double> alt_massdiffs_;
  Array<int> alt_matchedions_;
  Array<int> precursor_ioncharges_;
  Array<char*> mods_;
  Array<char*> alt_mods_;

  Array<char*> params_;
  Array<char*> vals_;

  pICalculator* pi_calc_;
  Boolean EXCLUDE_CONST_STATICS_; //True;
  Array<StaticModificationCount> static_consts_;
  char* timestamp_;

  char* outfilename_;
  char* dtafilename_;
  char* dtafilename2_;
  char* tarfilename_;
  char* directoryname_;
  int directoryname_len_; // for speed and storage optimisation in writeDtaFile
  char *outPepXMLfileName_;

  char* fullpath_;
  char* filextn_;

  std::string regressionTestFileName_; // for self test
  std::string regressionTestTarFileName_; // for self test
  eTagListFilePurpose testmode_;
  bool compare_tarred_file(std::string &text,const char *tarredfilename);

  int num_ion_series_;

  char* header_;
  char* enzyme_; // search constraint (if any)
  char* sample_enzyme_; // what the sample was treated with
  ProteolyticEnzyme* sampleEnzymeInstance_; //cached sample enzyme instance
  ProteolyticEnzyme* searchEnzymeInstance_; //cached search enzyme instance

  std::map<char, double> workingAAMasses; // fixed masses, which may or may not be modified
  std::map<char, Boolean> aaMod; // is the corresponding workingAAMasses entry modified?

private:
  char *xmlEscape(const char *string) const;
  const char *getSearchParameter(const char *name) const;
  int getMaxInternalCleavages() const;
  int getMinNumberTermini();
  Boolean isSemiCleavageSearch ();
	char *locatePeptide(const char* aasequence, const char *aasubsequence);
  inline char mapAA(const char aa) {
    return ('*' == aa ? 'X' : aa);
  }
  void updatePeptideTerminus (const char *dbProtein, const char *dbProteinDesc, const char *dbProteinSequence, 
    int &proteinIDType, int &proteinIDGuessedType, int &proteinIDTypeGuessCount, Boolean removeMapEntry, 
    Array<char*> &proteins, Array<char*> &peptides, 
    Array<char> &prevAAs, Array<char> &follAAs, Array<int> &indices);
  void buildProteinPeptideInfoViaLookup (const char *dbProtein, const char *dbProteinDesc, const char *dbProteinSequence, 
    int &proteinIDType, int &proteinIDGuessedType, int &proteinIDTypeGuessCount, 
    Array<int> &indices, Array<int> &rank1_indices, Array<int> &alt_indices);
  void buildProteinPeptideInfo (const char *dbProtein, const char *dbProteinDesc, const char *dbProteinSequence, 
    int &proteinIDType, int &proteinIDGuessedType, int &proteinIDTypeGuessCount, 
    Array<int> &indices, Array<int> &rank1_indices, Array<int> &alt_indices);
  int guessProteinIDType (const char *dbProtein, const char *resultProtein);

  // jmt:
  // helpers for recording mass info
  std::string extractModificationInfo(const std::string & inputString, // input
				      Boolean& nterm, // output
				      Boolean& cterm, // output
				      Boolean& proteinterm // output
				      );

  void recordStaticAAMass(char AA, double mass);
  void recordStaticAAMod(char AA, double fixedDelta, Boolean nterm, Boolean cterm);
  void recordVariableAAMod(char AA, double variableDelta, Boolean nterm, Boolean cterm, int variableIndex);

  void loadCycleExptScanFromMZXML (const char *datFile);
  void parseCyleExperimentForScan (const char *nextspectrum, long &rangeStart, long &rangeEnd);
  long CycleExperimentToScan (long lCycle, long lExperiment) const;
  double scanToRetentionTime (long lScan) const;
  ProteolyticEnzyme* getSampleEnzyme (bool verbose) const;
  ProteolyticEnzyme* getSearchEnzyme (bool verbose) const;
};


#endif
