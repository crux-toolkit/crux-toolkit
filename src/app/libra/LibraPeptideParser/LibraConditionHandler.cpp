#include "LibraConditionHandler.hpp"
#include <algorithm>
#include "Common/TPPVersion.h" // contains version number, name, revision
#include "Parsers/mzParser/mzParser.h"

using namespace mzParser;

LibraConditionHandler::LibraConditionHandler() {
  conditionFileName = NULL;
  m_pseudoCount = 0;
  m_maxPseudoChannels = 99;
}

LibraConditionHandler::~LibraConditionHandler() {}

void LibraConditionHandler::setFileName(const char* infile) {
  conditionFileName = infile;
}

const char* LibraConditionHandler::getFileName() const {
  return conditionFileName;
}

int LibraConditionHandler::getCentroidingPref() const {
  return m_centroidingPref;
}

const vvf &LibraConditionHandler::getMassIsotopes() const {
  return m_massIsotopes;
}

const std::vector<double> &LibraConditionHandler::getReagentMasses() const {
  return m_mass;
}

/**
 * checks and returns the value of m_normalPosition. 
 * If the value is unreasonable it is reset to a default
 * of 0 before return.
 * @return type of normalization to perform 
 * (where [0] = use highest intensity, [-2] = use against TIC,
 * and [> 1] is use this specifed isoptope channel)
 */ 
int LibraConditionHandler::getNormalPosition() const {
  return m_normalPosition;
}

double LibraConditionHandler::getTolerance() const {
  return m_tolerance;
}

int LibraConditionHandler::getIsToNormalize() const {
  return m_isToNormalize;
}

bool LibraConditionHandler::getUseIsotopicCorrection() const {
  return m_useIsotopicCorrection;
}

int LibraConditionHandler::getNumCentroidingIterations() const {
  return m_centroidingIteration;
}

bool LibraConditionHandler::getReporterFromMS3() const {
  return m_reporterFromMS3;
}


/**
 * parse condition xml file
 */
int LibraConditionHandler::readFile() {
  double mass;
  char nextline[10000];
  char* data = NULL;
  Array<Tag*>* tags = NULL;
  Tag* tag = NULL;

  // in  Element tags for elements with children
  int inFragmentMassesElement; // 0 = not read yet; 1= in element; 2 = done reading
  int inIsotopicContributionsElement;
  int inContributingMzElement;
  int aContributingMz=1;
  int anAffectedMz;
  float aCorrection;
  int initIsoMatrix;

  std::ifstream in(conditionFileName);
  int checkread;

  if (dbug) cout << "reading condition file : "<< conditionFileName << endl;

  try  {
    if (!in) throw 100;

    inFragmentMassesElement = 0;
    inIsotopicContributionsElement = 0;
    inContributingMzElement = 0;
    initIsoMatrix = 0;

    while( in.getline(nextline, 10000) ) {
      // check read state:
      checkread = in.rdstate();

      if (checkread & ios::eofbit) int dummyvar = 1;
      else if (checkread & ios::failbit ) throw 101;
      else if (checkread & ios::badbit ) throw 101;

      data = strstr(nextline, "<");

      if (!data) { // old style config file?
	if (inFragmentMassesElement) {
	  mass = atof(nextline);
	  if (mass) {
	    m_mass.push_back( mass );
	  }
	}
      }

      while(data != NULL) {
	tag = new Tag(data);

	// fragmentMasses

	// set flag to in element
	if ( tag->isStart() && ! strcmp(tag->getName(), "fragmentMasses")) {
	  inFragmentMassesElement = 1;
	} 

	// set flag to out of element
	if ( tag->isEnd() && ! strcmp(tag->getName(), "fragmentMasses")) {
	  inFragmentMassesElement = 2;

	  // if didn't find any reagents, throw error
	  if (m_mass.size() == 0 ) throw 111;
	}

	// get reagent mz attributes
	if ( (inFragmentMassesElement == 1) && (tag->isStart())
	     && (! strcmp(tag->getName(), "reagent"))) {
	  mass = atof(tag->getAttributeValue("mz"));

	  m_mass.push_back( mass );
	}

	// initialize isotopic corrections matrix if done reading reagent masses
	if (( inFragmentMassesElement == 2) && (initIsoMatrix == 0) ) {
	  int row, col;

	  // Fill the isotopic contribution matrix with zeros
	  for( row = 0 ; row < (int) m_mass.size() ; row++ ) {
	    m_massIsotopes.push_back( vf() );

	    for( col = 0 ; col < (int) m_mass.size() ; col++ ) {
	      m_massIsotopes[row].push_back( 0 );
	    }
	  }

	  initIsoMatrix = 1;
	}

	// isotopicContributions
	if ( tag->isStart() && ! strcmp(tag->getName(), "isotopicContributions")) {
	  m_useIsotopicCorrection = true;
	  inIsotopicContributionsElement = 1;
	}

	// set flag to out of element
	if ( tag->isEnd() && ! strcmp(tag->getName(), "isotopicContributions")) {
	  inIsotopicContributionsElement = 0;
	}

	if ((inIsotopicContributionsElement == 1) && (initIsoMatrix == 1)) {
	  if ( tag->isStart() && ! strcmp(tag->getName(), "contributingMz")) {
	    aContributingMz = atoi(tag->getAttributeValue("value"));
	  }

	  if ( tag->isStart() && ! strcmp(tag->getName(), "affected")) {
	    anAffectedMz = atoi(tag->getAttributeValue("mz"));
	    aCorrection = atof(tag->getAttributeValue("correction"));
	    m_massIsotopes[aContributingMz-1][anAffectedMz-1] = aCorrection;
	  }
	}

	// remaining elements
	if ( tag->isEnd() && ! strcmp(tag->getName(), "massTolerance")) {
	  m_tolerance = atof(tag->getAttributeValue("value"));

	  // check for smallest mz diff; stop with error if that is less than 2 x tolerance
	  float minMzDiff = 100.0;
	  if(m_mass.size() > 1) {
	    for(int p = 0 ; p < (int) m_mass.size() - 1; p++)
	      for(int q = p+1 ; q < (int) m_mass.size() ; q++)
		if (abs(m_mass[p] - m_mass[q]) < minMzDiff)
		  minMzDiff = abs(m_mass[p] - m_mass[q]);

	    if (2 * m_tolerance > minMzDiff) throw 112;

	  } else {
	    cerr << "Warning: Only one channel specified for Libra quantitation" << endl;
	  }
	}

	if ( tag->isEnd() && ! strcmp(tag->getName(), "centroiding")) {
	  m_centroidingPref = atoi(tag->getAttributeValue("type"));
	  m_centroidingIteration = atoi(tag->getAttributeValue("iterations"));
	}

	if ( tag->isEnd() && ! strcmp(tag->getName(), "normalization")) {
	  m_normalPosition = atoi(tag->getAttributeValue("type"));

	  if( m_normalPosition == -1 )
	    m_isToNormalize = false;
	  else
	    m_isToNormalize = true;

	}

	if ( tag->isEnd() && ! strcmp(tag->getName(), "pseudocount")) {
	  int temp = atoi(tag->getAttributeValue("value"));
	  if(temp < 0)
	    temp = 0;
	  m_pseudoCount = temp;
	  m_maxPseudoChannels = atoi(tag->getAttributeValue("maxChannels"));
	}

	if ( tag->isEnd() && ! strcmp(tag->getName(), "targetMs"))
	  m_requestedMsLevel = atoi(tag->getAttributeValue("level"));

	if ( tag->isEnd() && ! strcmp(tag->getName(), "output"))
	  m_outputPrefs = atoi(tag->getAttributeValue("type"));

	if ( tag->isEnd() && ! strcmp(tag->getName(), "reporterFromMS3"))
	  m_reporterFromMS3 = atoi(tag->getAttributeValue("value"));


	if (tag != NULL) delete tag;
	tag = NULL;

	data = strstr(data+1, "<");
      }
    }
    in.close();

    if(m_pseudoCount > 0 && m_maxPseudoChannels > m_mass.size()) throw 113;

  } catch ( int code) {

    if( code == 100) {
      cerr << "Error opening file '" << conditionFileName << "'" << endl
	   << "Please make sure that the file name and path are correct." << endl;
      return 1; // Failed to open file

    } else if ( code == 101) {
      cerr << "The data in " << conditionFileName << 
	" may not be of type char." << endl;

    } else if ( code == 102) {
      cerr << "Error: Could not find value attribute in contributingMz element\n";

    } else if ( code == 103) {
      cerr << "Error: Could not find contributingMz affected\n";

    } else if ( code == 104) {
      cerr << "Error: Could not find contributingMz correction\n";

    } else if ( code == 105) {
      cerr << "Error: Could not find value attribute in massTolerance element\n";

    } else if ( code == 106) {
      cerr << "Error: Could not find type attribute in centroiding element\n";

    } else if ( code == 107) {
      cerr << "Error: Could not find iterations attribute in centroiding element\n";

    } else if ( code == 108) {
      cerr << "Error: Could not find type attribute in normalization element\n";

    } else if ( code == 109) {
      cerr << "Error: Could not find level attribute in targetMs element\n";

    } else if ( code == 110) {
      cerr << "Error: Could not find type attribute in normalization element\n";

    } else if ( code == 111) {
      cerr << "Error: Could not find reagent masses. condition.xml format incorrect? \n";

    } else if ( code == 112) {
      cerr << "Error: At least two channels have mz values within the specified massTolerance.  Please review your parameters. \n";

    } else if ( code == 113) {
      cerr << "Error: Maximum allowed channels with no signal exceeds total.  Please review your parameters. \n";
    }

    in.close();

    return -2; // failure
  }

  return 0; 
}


/**
 * Fill the remaining positions on a line of a IC matrix
 * with zero values
 * @return 0 for success
 */
int LibraConditionHandler::zeroFill( int n ) {
  for(int p = 0 ; p < (int) m_mass.size() ; p++) {
    m_massIsotopes[n].push_back(0);
  }

  return 0;
}


/**
 * @return the MS level we want to look into
 */
int LibraConditionHandler::getRequestedMsLevel() const {
  return m_requestedMsLevel;
}


/** 
 * give the outfile a name derived from mzML file root, open it as new file
 * @param mzXML file name
 * @return 0
 */
int LibraConditionHandler::openOutFile(const char* inFileName ) {
  const char  *nameStart;
  char *test;

  char  outFileName[512];
  strncpy( outFileName , inFileName , sizeof(outFileName) );

  // Skip the path
  if (( test = findRightmostPathSeperator( outFileName )) != NULL ) {
    nameStart = inFileName + (test-outFileName) + 1;
  } else {
    nameStart = inFileName;
  }

  // replace mzXML/mzData suffix
  strncpy( outFileName , nameStart , sizeof(outFileName) );

  if (!(test = rampValidFileType(outFileName))) { // watch for .mzXML.gz
    test = strrchr(outFileName,'.');
  }

  if (!test) { // no .ext found
    test = outFileName+strlen(outFileName); // attach at end
  }

  strncpy( test , ".aq", sizeof(outFileName)-(test-outFileName) );
  m_fout.open( outFileName );
  return 0;
}


/**
 * close m_fout, previously opened in openOutFile()
 * @return 0
 */
int LibraConditionHandler::closeOutFile() {
  m_fout.close();
  return 0;
}


/**
 * write conditions to output file with name mzXMLFile root .aq
 * @param mzXML file name
 * @return 0
 */
int LibraConditionHandler::writeLibraConditionHandlerToOutFile(const char* mzXMLFile ) {
  openOutFile( mzXMLFile);

  if( !m_fout.good() ) {
    cerr << "Could not open output file\n";
    exit(2);
  }

  vf::iterator massPos;

  // Save condition in the output file
  if( m_useIsotopicCorrection)
    cout << "\n\nUsing the following parameters:" << endl;

  m_fout << ">ScanNum\t";

  for( massPos = m_mass.begin() ; massPos < m_mass.end() ; massPos++ ) {
    if( m_useIsotopicCorrection)
      cout << "\tIC from " << *massPos;

    m_fout << *massPos << "\t";
  }

  cout << endl;
  m_fout << endl;

  if( m_useIsotopicCorrection) {
    massPos = m_mass.begin();

    for( int n = 0 ; n < (int) m_massIsotopes.size() ; n++ ) {
      cout << *massPos << "\t";
      massPos++;

      for( int p = 0 ; p < (int) m_massIsotopes.size() ; p++ ) {
	cout << m_massIsotopes[n][p] << "\t\t";
      }

      cout << endl;
    }

    cout << endl;
  }

  closeOutFile();
  return 0;
}


/**
 * get output preferences:  scanNum (1) or RT (2)
 * @return 1 for scanNum, 2 for RT
 */
int LibraConditionHandler::getOutputPrefs() const {
  return  m_outputPrefs;
}


// format needed for PepXML
Array<Tag*>* LibraConditionHandler::getPepXMLTags() {
  char text[512];
  Array<Tag*>* output = new Array<Tag*>;
  Tag* next = new Tag("libra_summary", True, False);
  snprintf(text,sizeof(text),"LibraPeptideParser (%s)",szTPPVersionInfo);
  next->setAttributeValue("version", text);
  sprintf(text, "%0.3f", m_tolerance);
  next->setAttributeValue("mass_tolerance", text);
  sprintf(text, "%d", m_centroidingPref);
  next->setAttributeValue("centroiding_preference", text);
  sprintf(text, "%d", m_normalPosition);
  next->setAttributeValue("normalization", text);

  if (m_pseudoCount > 0) {
    sprintf(text, "%d", m_pseudoCount);
    next->setAttributeValue("pseudocount", text);
    sprintf(text, "%d", m_maxPseudoChannels);
    next->setAttributeValue("max_pseudo_channels", text);
  }

  sprintf(text, "%d", m_outputPrefs);
  next->setAttributeValue("output_type", text);
  output->insertAtEnd(next);
  int index = output->length()-1; // for adding attribute later
  char code_text[1000];
  code_text[0] = 0;

  for(int k = 0; k < (int) m_mass.size(); k++) {
    next = new Tag("fragment_masses", True, True);
    sprintf(text, "%d", k+1);
    next->setAttributeValue("channel", text);
    sprintf(text, "%0.5f", m_mass[k]);
    strcat(code_text, text);
    next->setAttributeValue("mz", text);
    output->insertAtEnd(next);
  } // next reagent mass

  (*output)[index]->setAttributeValue("channel_code", code_text);

  // isotopic contributions:  
  next = new Tag("isotopic_contributions", True, False);
  output->insertAtEnd(next);
  int row, col;
  for( row = 0 ; row < (int) m_mass.size() ; row++ ) {
    sprintf(text, "%d", row+1);
    next = new Tag("contributing_channel", True, False);

    next->setAttributeValue("channel", text);
    output->insertAtEnd(next);

    for( col = 0 ; col < (int) m_mass.size() ; col++ ) {
      if(col != row) {
	next = new Tag("affected_channel", True, True);
	sprintf(text, "%d", col+1);
	next->setAttributeValue("channel", text);
	sprintf(text, "%0.3f", m_massIsotopes[row][col]);
	next->setAttributeValue("correction", text);
	output->insertAtEnd(next);
      } // if now same as contr
    } // next affected
    next = new Tag("contributing_channel", False, True);
    output->insertAtEnd(next);
  }

  next = new Tag("isotopic_contributions", False, True);
  output->insertAtEnd(next);

  next = new Tag("libra_summary", False, True);
  output->insertAtEnd(next);

  return output;
}


/**
 * print state of class to standard out
 */
int LibraConditionHandler::toString() {
  cout << "fileName: " << getFileName() << endl;
  cout << "centroidingPref: " << getCentroidingPref() << endl;
  cout << "normalPosition: " << getNormalPosition() << endl;
  cout << "numCentroidingIterations: " << getNumCentroidingIterations() << endl;
  cout << "tolerance: " << getTolerance() << endl;

  std::vector<double> tmp2;
  tmp2 = getReagentMasses();

  for (int i = 0; i < (int) tmp2.size(); i++) {
    cout << "reagent mass: " << tmp2[i] << endl;
  }

  return 1;
}
