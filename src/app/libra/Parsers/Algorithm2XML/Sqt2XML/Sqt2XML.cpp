#include "Sqt2XML.h"
#include "Common/TPPVersion.h"
#include "SimpleXMLWriter.h"
#include "Parsers/Algorithm2XML/Sequest2XML/SequestParams.h"
#include "Common/Enzyme/ProteolyticEnzyme/ProteolyticEnzymeFactory/ProteolyticEnzymeFactory.h"
#include "Common/constants.h"
#include <iostream>
#include <iomanip>
#include "gzstream.h"

#ifdef WIN32
#define PATH_DELIMITERS "\\/" // Windows works with both
#else
#define PATH_DELIMITERS "/"
#endif

using namespace mzParser;

static const double protonMass = 1.00727646688;

string trim(string& s, const string& dropChars = " \t\f\r\n") {
  string r = s.erase(s.find_last_not_of(dropChars) + 1);
  return r.erase(0, r.find_first_not_of(dropChars));
}


void writeTagArray(SimpleXMLWriter& writer, Array<Tag*>* tagArray) {
  if( tagArray == NULL )
    return;
  for( int i = 0; i < tagArray->length(); ++i ) {
    Tag* tag = (*tagArray)[i];
    if( tag != NULL ) {
      if( tag->isStart() ) {
	writer.open( tag->getName() );
	for( int j=0; j < tag->getAttributeCount(); ++j )
	  writer.attr( tag->getAttribute(j), tag->getAttributeValue(j) );
      }

      if( tag->isEnd() )
	writer.close();
    }
  }
}

void deleteTagArray( Array<Tag*>* tagArray ) {
  if( tagArray == NULL )
    return;
  for( int i=0; i < tagArray->length(); ++i ) {
    if( (*tagArray)[i] != NULL )
      delete (*tagArray)[i];
  }
  delete tagArray;
}

int main(int argc, char** argv) {
  hooks_tpp handler(argc,argv); // set up install paths etc

  string inputFilename, outputFilename, sourceName;

  int maxResultRank = 10;
  string paramsFilepath = "sequest.params";
  char* sample_enzyme = new char[128];

  int minArgs = 1;
  int flagArgs = 0;

  strcpy(sample_enzyme, "trypsin");
  for( int i=1; i < argc; ++i )
    if( argv[i][0] == '-' ) {
      ++flagArgs;

      if( argv[i][1] == 'r' ) {
	maxResultRank = lexical_cast<int>( argv[i+1] );
	cout << "Setting maximum result rank to " << maxResultRank << endl;
      } else if( argv[i][1] == 'p' ) {
	paramsFilepath = argv[i+1];
	cout << "Setting sequest.params filepath to " << paramsFilepath << endl;
      } else if (argv[i][1] == 'E') {
	strcpy(sample_enzyme, argv[i+1]);
	cout << "Setting sample_enzyme to " << sample_enzyme << endl;
      }

    }

  if( argc < minArgs + flagArgs + 1 ) {
    cerr << "SQT2XML (" << szTPPVersionInfo << ")" << endl
	 << "Usage: " << argv[0] << " [OPTIONS] <filepath to SQT file>" << endl
	 << endl
	 << "   OPTIONS:" << endl
	 << "      -r <max. result rank to keep, default 10>" << endl
	 << "      -p <filepath to sequest.params>" << endl
	 << "      -E <enzyme> " << endl
	 << "       Where <enzyme> is: " << endl
	 << "         trypsin - Cut: KR, No Cut: P, Sense: C-term (default)" << endl
	 << "         stricttrypsin - Cut: KR, No Cut: none, Sense: C-term " << endl
	 << "         argc - Cut: R, No Cut: P, Sense: C-term " << endl
	 << "         aspn - Cut: D, No Cut: none, Sense: N-term " << endl
	 << "         chymotrypsin - Cut: YWFL, No Cut: P, Sense: C-term " << endl
	 << "         cnbr - Cut: M, No Cut: P, Sense: C-term " << endl
	 << "         elastase - Cut: GVLIA, No Cut: P, Sense: C-term " << endl
	 << "         formicacid - Cut: D, No Cut: P, Sense: C-term " << endl
	 << "         gluc - Cut: DE, No Cut: P, Sense: C-term " << endl
	 << "         gluc_bicarb - Cut: E, No Cut: P, Sense: C-term " << endl
	 << "         iodosobenzoate - Cut: W, No Cut: terminal, Sense: C-term " << endl
	 << "         lysc - Cut: K, No Cut: P, Sense: C-term " << endl
	 << "         lysc-p - Cut: K, No Cut: none, Sense: C-term " << endl
	 << "         lysn - Cut: K, No Cut: none, Sense: N-term " << endl
	 << "         lysn_promisc - Cut: KR, No Cut: none, Sense: N-term " << endl
	 << "         nonspecific - Cut: all, No Cut: none, Sense: N/A " << endl
	 << "         pepsina - Cut: FL, No Cut: terminal, Sense: C-term " << endl
	 << "         protein_endopeptidase - Cut: P, No Cut: terminal, Sense: C-term " << endl
	 << "         staph_protease - Cut: E, No Cut: terminal, Sense: C-term " << endl
	 << "         tca - Cut: KR, No Cut: P, Sense: C-term " << endl
	 << "             - Cut: YWFM, No Cut: P, Sense: C-term " << endl
	 << "             - Cut: D, No Cut: none, Sense: N-term " << endl
	 << "         trypsin/cnbr - Cut: KR, No Cut: P, Sense: C-term " << endl
	 << "         trypsin_gluc - Cut: DEKR, No Cut: P, Sense: C-term " << endl
	 << "         trypsin_k - Cut: K, No Cut: P, Sense: C-term " << endl
	 << "         trypsin_r - Cut: R, No Cut: P, Sense: C-term " << endl << endl;

    exit(1);
  }

  inputFilename = argv[argc-1];
  //outputFilename = inputFilename.substr( 0, inputFilename.find_last_of('.')+1 ) + "pep.xml";
  outputFilename = inputFilename.substr( 0, inputFilename.find_last_of('.') ) + get_pepxml_dot_ext();
  sourceName = inputFilename.substr( inputFilename.find_last_of( PATH_DELIMITERS )+1, inputFilename.find_last_of('.')-inputFilename.find_last_of( PATH_DELIMITERS )-1 );

  SequestParams sequestParams( (char*)paramsFilepath.c_str() );

  // write as gzip if filename so indicates
  ogzstream pepXmlFile( outputFilename.c_str(), isDotGZ(outputFilename)?9:0 );
  if( !pepXmlFile ) {
    cerr << "error opening output file " << outputFilename << endl;
    exit(1);
  }

  // Get instrument info from mzXML file
  RAMPFILE *pFI;
  ramp_fileoffset_t *pScanIndex = NULL;
  int iAnalysisLastScan;
  string mzXmlPath_;

  mzXmlPath_ = rampConstructInputFileName(inputFilename.substr( 0, inputFilename.find_last_of('.')));

  InstrumentStruct* ppstruct = NULL;
  ramp_fileoffset_t indexOffset;
  cout <<  "attempting to read mzXML/mzML: "; //add some context to rampOpenFile message "Unknown file type. No file loaded." when mzXML not present
  fflush(stdout);
  if ( (pFI=rampOpenFile(mzXmlPath_.c_str()))==NULL) {
    //suppress warning
    //		cout << " warning: cannot open mzXML file \"" << mzXmlPath_  << "\" for reading MS instrument info and scan times." << endl;
  }
  else {
    cout << "done." << endl;
    ppstruct = getInstrumentStruct(pFI);
    indexOffset = getIndexOffset(pFI);
    pScanIndex = readIndex( pFI , indexOffset, &iAnalysisLastScan );
  }

  SimpleXMLWriter writer;
  Array<Tag*>* tags;
  writer.condenseAttr_ = true;
  writer.setOutputStream( pepXmlFile );
  writer.startDocument();

  writer.open( "msms_pipeline_analysis" );
  writer.attr( "name", sourceName );
  writer.attr( "date", "01/01/2001" );
  writer.attr( "summary_xml", outputFilename );

  writer.open( "msms_run_summary" );
  writer.attr( "base_name", sourceName );

  if (ppstruct!=NULL) {
    writer.attr( "msManufacturer", ppstruct->manufacturer );
    writer.attr( "msModel", ppstruct->model );
    writer.attr( "msIonization", ppstruct->ionisation );
    writer.attr( "msAnalyzer", ppstruct->analyzer );
    writer.attr( "raw_data_type", "raw" );
    writer.attr( "raw_data", strrchr(mzXmlPath_.c_str(),'.'));
  }
  else {
    writer.attr( "raw_data_type", "unknown" );
    writer.attr( "raw_data", "unknown" );
  }

  //TODO: Have to support multiple enzymes
  ProteolyticEnzyme* enzyme = (new ProteolyticEnzymeFactory())->getProteolyticEnzyme(sample_enzyme);
  tags = enzyme->getPepXMLTags();
  writeTagArray( writer, tags );
  deleteTagArray( tags );

  ifstream sqtFile( inputFilename.c_str() );
  if( !sqtFile ) {
    cerr << "error opening input file " << inputFilename << endl;
    exit(1);
  }

  char* engine = strdup("SEQUEST");

  string inputLine;
  while( getline( sqtFile, inputLine ) ) {
    if( inputLine.find("H\tSQTGenerator ") == 0 || inputLine.find("H\tSQTGenerator\t") == 0 )
      engine = strdup( inputLine.substr( 15, inputLine.find_first_of( "\r\n \t", 15 ) ).c_str() );
    else if( inputLine.find("S\t") == 0 )
      break;
  }

  tags = sequestParams.getSearchParamTags( (char*)sourceName.c_str(), engine );
  writeTagArray( writer, tags );
  deleteTagArray( tags );
  free(engine);

  if( !sqtFile ) {
    cerr << "nothing to convert in " << inputFilename << endl;
    exit(1);
  }

  size_t tokenStart, tokenEnd;
  int index = 1;

  while( sqtFile && inputLine[0] == 'S' ) {
    S_entry s;

    tokenStart = 2; tokenEnd = inputLine.find( '\t', tokenStart+1 ); // skip S and \t
    //cout << "\"" << inputLine.substr( tokenStart, tokenEnd-tokenStart ) << "\" ";
    s.firstScan = lexical_cast<int>( inputLine.substr( tokenStart, tokenEnd-tokenStart ) );

    tokenStart = tokenEnd+1; tokenEnd = inputLine.find( '\t', tokenStart+1 );
    //cout << "\"" << inputLine.substr( tokenStart, tokenEnd-tokenStart ) << "\" ";
    s.lastScan = lexical_cast<int>( inputLine.substr( tokenStart, tokenEnd-tokenStart ) );

    tokenStart = tokenEnd+1; tokenEnd = inputLine.find( '\t', tokenStart+1 );
    //cout << "\"" << inputLine.substr( tokenStart, tokenEnd-tokenStart ) << "\" ";
    s.chargeState = lexical_cast<int>( inputLine.substr( tokenStart, tokenEnd-tokenStart ) );

    tokenStart = tokenEnd+1; tokenEnd = inputLine.find( '\t', tokenStart+1 );
    //cout << "\"" << inputLine.substr( tokenStart, tokenEnd-tokenStart ) << "\" ";
    s.processingTime = lexical_cast<float>( inputLine.substr( tokenStart, tokenEnd-tokenStart ) );

    tokenStart = tokenEnd+1; tokenEnd = inputLine.find( '\t', tokenStart+1 );
    //cout << "\"" << inputLine.substr( tokenStart, tokenEnd-tokenStart ) << "\" ";
    //s.processingHostname = inputLine.substr( tokenStart, tokenEnd-tokenStart );

    tokenStart = tokenEnd+1; tokenEnd = inputLine.find( '\t', tokenStart+1 );
    //cout << "\"" << inputLine.substr( tokenStart, tokenEnd-tokenStart ) << "\" ";
    s.observedMassPlus1 = lexical_cast<float>( inputLine.substr( tokenStart, tokenEnd-tokenStart ) );

    tokenStart = tokenEnd+1; tokenEnd = inputLine.find( '\t', tokenStart+1 );
    //cout << "\"" << inputLine.substr( tokenStart, tokenEnd-tokenStart ) << "\" ";
    s.totalIntensity = lexical_cast<float>( inputLine.substr( tokenStart, tokenEnd-tokenStart ) );

    tokenStart = tokenEnd+1; tokenEnd = inputLine.find( '\t', tokenStart+1 );
    //cout << "\"" << inputLine.substr( tokenStart, tokenEnd-tokenStart ) << "\" ";
    s.lowestSp = lexical_cast<float>( inputLine.substr( tokenStart, tokenEnd-tokenStart ) );

    tokenStart = tokenEnd+1; tokenEnd = inputLine.find( '\t', tokenStart+1 );
    //cout << "\"" << inputLine.substr( tokenStart, tokenEnd-tokenStart ) << "\"\n";
    s.numSequenceComparisons = lexical_cast<int>( inputLine.substr( tokenStart, tokenEnd-tokenStart ) );

    getline( sqtFile, inputLine );

    while( sqtFile && inputLine[0] == 'M' ) {
      M_entry m;

      tokenStart = 2; tokenEnd = inputLine.find( '\t', tokenStart+1 ); // skip M and \t
      //cout << "\"" << inputLine.substr( tokenStart, tokenEnd-tokenStart )<< "\" ";
      m.rankByXcorr = lexical_cast<int>( inputLine.substr( tokenStart, tokenEnd-tokenStart ) );
      if( maxResultRank > 0 && m.rankByXcorr > maxResultRank ) {
	do {
	  getline( sqtFile, inputLine );
	} while( sqtFile && inputLine[0] != 'S' );
	break;
      }

      tokenStart = tokenEnd+1; tokenEnd = inputLine.find( '\t', tokenStart+1 );
      //cout << "\"" << inputLine.substr( tokenStart, tokenEnd-tokenStart ) << "\" ";
      m.rankBySp = lexical_cast<size_t>( inputLine.substr( tokenStart, tokenEnd-tokenStart ) );

      tokenStart = tokenEnd+1; tokenEnd = inputLine.find( '\t', tokenStart+1 );
      //cout << "\"" << inputLine.substr( tokenStart, tokenEnd-tokenStart ) << "\" ";
      m.calculatedMassPlus1 = lexical_cast<float>( inputLine.substr( tokenStart, tokenEnd-tokenStart ) );

      tokenStart = tokenEnd+1; tokenEnd = inputLine.find( '\t', tokenStart+1 );
      //cout << "\"" << inputLine.substr( tokenStart, tokenEnd-tokenStart ) << "\" ";
      m.deltacn = lexical_cast<float>( inputLine.substr( tokenStart, tokenEnd-tokenStart ) );

      tokenStart = tokenEnd+1; tokenEnd = inputLine.find( '\t', tokenStart+1 );
      //cout << "\"" << inputLine.substr( tokenStart, tokenEnd-tokenStart ) << "\" ";
      m.xcorr = lexical_cast<float>( inputLine.substr( tokenStart, tokenEnd-tokenStart ) );

      tokenStart = tokenEnd+1; tokenEnd = inputLine.find( '\t', tokenStart+1 );
      //cout << "\"" << inputLine.substr( tokenStart, tokenEnd-tokenStart ) << "\" ";
      m.sp = lexical_cast<float>( inputLine.substr( tokenStart, tokenEnd-tokenStart ) );

      tokenStart = tokenEnd+1; tokenEnd = inputLine.find( '\t', tokenStart+1 );
      //cout << "\"" << inputLine.substr( tokenStart, tokenEnd-tokenStart ) << "\" ";
      m.matchedIons = lexical_cast<int>( inputLine.substr( tokenStart, tokenEnd-tokenStart ) );

      tokenStart = tokenEnd+1; tokenEnd = inputLine.find( '\t', tokenStart+1 );
      //cout << "\"" << inputLine.substr( tokenStart, tokenEnd-tokenStart ) << "\" ";
      m.predictedIons = lexical_cast<int>( inputLine.substr( tokenStart, tokenEnd-tokenStart ) );

      tokenStart = tokenEnd+1; tokenEnd = inputLine.find( '\t', tokenStart+1 );
      //cout << "\"" << inputLine.substr( tokenStart, tokenEnd-tokenStart ) << "\"\n";
      m.sequence = inputLine.substr( tokenStart, tokenEnd-tokenStart );
      // SEQUEST output includes extra whitespace, so must trim here!
      m.sequence = trim( m.sequence ); 
      //m.sequence = m.sequence.substr( 2, m.sequence.length() - 4 ); // trim flanking residue notation
      //m.sequence = ConvertSqtPtmToFreiPtm( r.sequence, &fileResidueMap );
      //cout << m.sequence << endl;

      getline( sqtFile, inputLine );

      while( sqtFile && inputLine[0] == 'L' ) {
	L_entry protein;

	tokenStart = 2; // skip L and \t
	tokenEnd = inputLine.find_first_of( "\r\n\t ", tokenStart+1 );
	//cout << "\"" << inputLine.substr( tokenStart, tokenEnd-tokenStart ) << "\"\n";
	string locus = inputLine.substr( tokenStart, tokenEnd-tokenStart );
	locus = trim( locus );
	m.loci.push_back( L_entry( locus ) ); 
	getline( sqtFile, inputLine );
      }

      if( m.loci.empty() ) {
	cerr << "match " << m.sequence << " in scan " << s.firstScan << "." << s.chargeState << " has no proteins" << endl;
	continue;
      }
      s.matches.push_back(m);
    }

    writer.open( "spectrum_query" );
    stringstream spectrumId;
    spectrumId << sourceName << "." << setfill('0') << setw(5) << s.firstScan << "." << setfill('0') << setw(5) << s.lastScan << "." << s.chargeState;
    writer.attr( "spectrum", spectrumId.str() );
    writer.attr( "start_scan", s.firstScan );
    writer.attr( "end_scan", s.lastScan );
    float precursorNeutralMass = s.observedMassPlus1 - protonMass;
    writer.attr( "precursor_neutral_mass", precursorNeutralMass );
    writer.attr( "assumed_charge", s.chargeState );
    writer.attr( "index", index++ );

    // add retention times if available
    if (pFI != NULL) {
      if(s.firstScan != -1 && s.firstScan <= iAnalysisLastScan) {
	struct ScanHeaderStruct scanHeader;
	readHeader(pFI, pScanIndex[s.firstScan],&scanHeader);
	writer.attr( "retention_time_sec", scanHeader.retentionTime );
      }
    }
    else
      writer.attr( "retention_time_sec", 0 );

    if( !s.matches.empty() ) {
      writer.open( "search_result" );
      writer.attr( "num_comparisons", s.numSequenceComparisons );

      for( size_t matchIndex = 0; matchIndex < s.matches.size(); ++matchIndex ) {
	M_entry& m = s.matches[matchIndex];

	size_t locusIndex = 0;

	string subSequence = m.sequence.substr( 2, m.sequence.length() - 4 );
	string plainSequence;

	Array<Tag*>* modinfo_tags = sequestParams.getModificationInfoTags( (char*)subSequence.c_str() );

	plainSequence.reserve( subSequence.length() );
	for( size_t i = 0; i < subSequence.length(); ++i ) {
	  plainSequence.push_back( subSequence[i] );
	  if( i+1 < subSequence.length()
	      && sequestParams.getMonoisotopicAAMass( subSequence[i+1] ) == 0.0
	      && subSequence[i+1]!='X')
	    {
	      ++i;
	    }
	}

	writer.open( "search_hit" );
	writer.attr( "hit_rank", m.rankByXcorr );
	writer.attr( "peptide", plainSequence );
	writer.attr( "peptide_prev_aa", *m.sequence.begin() );
	writer.attr( "peptide_next_aa", *m.sequence.rbegin() );
	writer.attr( "protein", m.loci[locusIndex] );
	writer.attr( "num_tot_proteins", (int) m.loci.size() );
	writer.attr( "num_matched_ions", m.matchedIons );
	writer.attr( "tot_num_ions", m.predictedIons );
	float sequenceNeutralMass = m.calculatedMassPlus1 - protonMass;
	writer.attr( "calc_neutral_pep_mass", sequenceNeutralMass );
	float precursorToSequenceMassError = s.observedMassPlus1 - m.calculatedMassPlus1;
	writer.attr( "massdiff", precursorToSequenceMassError );
	writer.attr( "num_tol_term", 2 );
	writer.attr( "num_missed_cleavages", enzyme->getNumMissedCleavages( m.sequence.c_str() ) );
	//writer.attr( "is_rejected", 0 );

	writeTagArray( writer, modinfo_tags );
	deleteTagArray( modinfo_tags );

	++locusIndex;
	for( ; locusIndex < m.loci.size(); ++locusIndex ) {
	  writer.open( "alternative_protein" );
	  writer.attr( "protein", m.loci[locusIndex] );
	  writer.close();
	}

	{
	  writer.open( "search_score" );
	  writer.attr( "name", "xcorr" );
	  writer.attr( "value", m.xcorr );
	  writer.close();

	  float nextDeltaCn = 0;   // unless we find better, assume the worse
	  size_t nextRankIndex = matchIndex+1;
	  while( nextRankIndex < s.matches.size() ) {
	    if( s.matches[nextRankIndex].rankByXcorr > m.rankByXcorr ) {
	      nextDeltaCn = s.matches[nextRankIndex].deltacn;
	      break;
	    }
	    ++nextRankIndex;
	  }
	  writer.open( "search_score" );
	  writer.attr( "name", "deltacn" );
	  writer.attr( "value", nextDeltaCn );
	  writer.close();

	  writer.open( "search_score" );
	  writer.attr( "name", "deltacnstar" );
	  writer.attr( "value", 0);
	  writer.close();

	  writer.open( "search_score" );
	  writer.attr( "name", "spscore" );
	  writer.attr( "value", m.sp);
	  writer.close();

	  writer.open( "search_score" );
	  writer.attr( "name", "sprank" );
	  writer.attr( "value", m.rankBySp);
	  writer.close();
	}
	writer.close(); // search_hit
      }
      writer.close(); // search_result
    }
    writer.close(); // spectrum_query
    do {
      getline( sqtFile, inputLine );
    } while( sqtFile && ( inputLine[0] != 'S') );
  }
  writer.close(); // msms_run_summary
  writer.close(); // msms_pipeline_analysis

  if (pFI != NULL)
    rampCloseFile(pFI);
  if (pScanIndex != NULL)
    free(pScanIndex);

  cout << "Done.  File created: " << outputFilename << endl << endl;

  return 0;
}
