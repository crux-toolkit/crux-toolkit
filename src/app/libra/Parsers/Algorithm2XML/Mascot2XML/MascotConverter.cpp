/*
Program       : Mascot2XML converter to pepXML
Author        : Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN Version   : $Id: MascotConverter.cpp 8162 2020-06-01 07:11:22Z real_procopio $

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

#include "MascotConverter.h"
#include "Common/util.h"
#include "Parsers/mzParser/mzParser.h"
#include <vector>
#include "Common/tpp_tarball.h" // tgz handling routines

using std::string;
using std::vector;
using std::pair;

using namespace mzParser;

enum {
  MIN_CONSECUTIVE_CORRECT_GUESS = 10,
  STATE_GIVEUP = -1,
  UNKNOWN = 0,
  ID_TYPE_FULL = 1,
  ID_TYPE_RIGHTMOST_PLUS_LEFT_ANCHOR = 2,
  ID_TYPE_LEFTMOST_PLUS_RIGHT_ANCHOR = 3,
  ID_TYPE_CENTER_PLUS_LEFT_RIGHT_ANCHOR = 4,
  ID_TYPE_VALID = ID_TYPE_FULL
};

/*
  Trim leading spaces
*/
void trimLeft (std::string &str) {
  Boolean changed = False;
  string::iterator iter;
  for (iter = str.begin(); iter != str.end(); iter++) {
    if (!isspace(*iter)) {
      str.erase(str.begin(), iter);
      changed = True;
      break;
    }
  }
  if (!changed && iter == str.end()) {
    str.clear();
  }
}

/*
  Trim trailing spaces
*/
void trimRight (std::string &str) {
  string::iterator iter;
  for (iter = str.end() - 1; ;iter--) {
    if (!isspace(*iter)) {
      str.erase(iter + 1, str.end());
      break;
    }
    if (iter == str.begin()) {
      str.clear();
      break;
    }
  }
}

/*
  Trim spaces from both ends
*/
void trimAll (std::string &str) {
  trimLeft(str);
  trimRight(str);
}

/*
  Keep the cycles and experiments grouped as a single spectra
*/
class ExperimentCycleRecord {
public:
  ExperimentCycleRecord(long lExperiment = 0, long lCycleStart = 0, long lCycleEnd = 0, Boolean bSingleCycle = True, Boolean bRangleCycle = False)
  {
    this->lExperiment=lExperiment;
    this->lCycleStart=lCycleStart;
    this->lCycleEnd=lCycleEnd;
    this->bSingleCycle=bSingleCycle;
    this->bRangleCycle=bRangleCycle;
  }

public:
  long lExperiment;
  long lCycleStart;
  long lCycleEnd;
  Boolean bSingleCycle;
  Boolean bRangleCycle;
};


MascotConverter::MascotConverter(char* filename, Boolean xml, char* database, char* sample_enzyme, Boolean calculate_PI, 
  Boolean generateTgz, Boolean generateDta, Boolean generateDescription, Boolean shortProteinId, Boolean verbose, 
  eTagListFilePurpose testmode, std::string &testfilename) {

  toLoadScanMap_ = True;  //Assume that we need to load scan map
  scanMapLoaded_ = False;

  database_ = NULL;
  sample_enzyme_ = NULL;
  sampleEnzymeInstance_ = NULL;
  searchEnzymeInstance_ = NULL;
  calculate_pI_ = calculate_PI;
  pi_calc_ = NULL;
  if(calculate_pI_)
    pi_calc_ = new pICalculator();
  testmode_ = testmode;
  regressionTestFileName_ = testfilename;
  regressionTestTarFileName_ = testfilename+std::string(".outs_and_dtas.tgz");

  EXCLUDE_CONST_STATICS_ = False; //True;

  timestamp_ = NULL;

  if(database != NULL && database[0]) {
    database_ = strCopy(database);
    validate_db();
    //cout << "database: " << database_ << endl;
  }
  if(sample_enzyme != NULL && strlen(sample_enzyme) > 0) {
    sample_enzyme_ = strCopy(sample_enzyme);
    //cout << "enzyme: " << sample_enzyme_ << endl;
  }

  generate_description_ = generateDescription;
  short_protein_id_ = shortProteinId;
  verbose_ = verbose;

  // initialize the working list AA masses
  for (char c = 'A'; c <= 'Z'; c++) {
    workingAAMasses[c] = -1;
    aaMod[c] = false;
  }

  workingAAMasses['n'] = -1;
  aaMod['n'] = false;

  workingAAMasses['c'] = -1;
  aaMod['c'] = false;


  init(filename, xml, generateTgz, generateDta);
}

#define ARRAYDEL(a) {for(int i=(a).size();i--;){delete (a)[i];}}

const char G_UNKNOWN_AA = '?';

MascotConverter::~MascotConverter() {
  delete[] database_;
  delete[] enzyme_;
  delete sampleEnzymeInstance_;
  delete searchEnzymeInstance_;
  ARRAYDEL(static_mods_);
  ARRAYDEL(variable_mods_);
  ARRAYDEL(vals_);
  delete[] dtafilename2_;
  ARRAYDEL(alt_peptides_);
  ARRAYDEL(alt_proteins_);
  ARRAYDEL(alt_mods_);
  ARRAYDEL(spectra_);
  ARRAYDEL(peptides_);
  ARRAYDEL(proteins_);
  rank1_peptides_.clear();
  ARRAYDEL(rank1_proteins_);
  ARRAYDEL(mods_);
  ARRAYDEL(params_);
  //ARRAYDEL(written_dtas_);
  delete pi_calc_;
  delete[] timestamp_;
  delete[] outfilename_;
  delete[] dtafilename_;
  delete[] tarfilename_;
  delete[] directoryname_;
  delete[] fullpath_;
  delete[] filextn_;
  delete[] header_;
  delete[] outPepXMLfileName_;
}

#define CONTENTTYPE_UNKNOWN -1
#define CONTENTTYPE_HEADER 1
#define CONTENTTYPE_MASSES 2
#define CONTENTTYPE_PARAMETERS 3
#define CONTENTTYPE_PEPTIDES 4
#define CONTENTTYPE_PROTEINS 5
#define CONTENTTYPE_SUMMARY 6
#define CONTENTTYPE_QUERY 7
int mapContentType (const char *name) {
  if (strcmp(name, "header") == 0) {
    return CONTENTTYPE_HEADER;
  } else if  (strcmp(name, "masses") == 0) {
    return CONTENTTYPE_MASSES;
  } else if  (strcmp(name, "parameters") == 0) {
    return CONTENTTYPE_PARAMETERS;
  } else if  (strcmp(name, "peptides") == 0) {
    return CONTENTTYPE_PEPTIDES;
  } else if  (strcmp(name, "decoy_peptides") == 0) {
    return CONTENTTYPE_PEPTIDES;
  } else if  (strcmp(name, "proteins") == 0) {
    return CONTENTTYPE_PROTEINS;
  } else if  (strcmp(name, "summary") == 0) {
    return CONTENTTYPE_SUMMARY;
  } else if  (strcmp(name, "decoy_summary") == 0) {
    return CONTENTTYPE_SUMMARY;
  } else if  (strncmp(name, "query", 5) == 0) {
    return CONTENTTYPE_QUERY;
  } else {
    return CONTENTTYPE_UNKNOWN;
  }
}

void MascotConverter::init(const char* filename, Boolean xml, Boolean generateTgz, Boolean generateDta) {

  Boolean bVerbose = false;
  char delimiter[] = "--gc0p4Jq0M2Yt08jU534c0p";
  char active[500];
  active[0] = 0; // to start with
  int line_width = 1000000;
  char* line = new char[line_width];
  header_ = NULL; // until parsed from parameters
  enzyme_ = NULL;

  term_read_ = False;
  time_t now;
  time(&now);

  char* conversion = ctime(&now);
  timestamp_ = new char[strlen(conversion)];
  strncpy(timestamp_, conversion, strlen(conversion)-1);
  timestamp_[strlen(conversion)-1] = 0;

  // jmt: used for looking for fixed modification section
  // potentially hack-y logic: assumes FixedModResidue<n> always follows FixedMod<n> line
  char fixedModDesc[255];
  bool lookingForFixedModResidues=false;
  double fixedModMass = -1;
  int fixedModNumber = -1;

  fullpath_ = NULL;
  filextn_ = NULL;
  setFileInfo(filename);
  if(fullpath_ == NULL || filextn_ == NULL) {
    cerr << "error extracting file info" << std::endl;
    exit(1);
  }

  outfilename_ = new char[500];
  dtafilename_ = new char[500];
  dtafilename2_ = new char[500];
 
  tarfilename_ = new char[500];
  directoryname_ = new char[500];
  directoryname_len_ = 0;

  strcpy(tarfilename_, fullpath_);
  strcat(tarfilename_, filextn_);
  strcat(tarfilename_, ".tgz");
  strcpy(directoryname_, fullpath_);
  //chee-hong
  //     should we use some random string here as we do not
  //     allow more than 1 instance of Mascot2XML to run.
  //     with "-notgz", this restriction is not there
  strcat(directoryname_, "converter");
  strcat(directoryname_, "/");

  // make tarball?
  if (generateTgz) {
    tarball_=new tpp_tarball(tarfilename_,directoryname_);
  }
  directoryname_len_ = (int)strlen(directoryname_);

  num_ion_series_ = 2; // b and y

  enzyme_constraint_ = False;

  mass_type_ = 0; // average unless proven otherwise

  char query[6];
  int contentType=CONTENTTYPE_UNKNOWN;
  int queryIndex=0; // WCH: TODO: should this be long?

  double mass;
  double massdiff;
  int first;
  char nextpeptide[1024];
  double nextionscore;
  char discard[1024];
  std::string nextprotein;
  char nextc;
  char nextspectrum[5000];
  int nextplus;
  int nextionmatch;
  char nextmods[1024];

  char* nextline = new char[line_width];
  int spectrum_counter = 1;

  char* cleaned = new char[strlen(timestamp_)+1];
  strcpy(cleaned, timestamp_);
  // now replace spaces with underscores, colons with dashes
  int k;
  for(k = 0; cleaned[k]; k++) {
    if(cleaned[k] == ' ')
      cleaned[k] = '_';
    if (strchr("*?[]/\\=+<>:;\", ",cleaned[k])) {
      cleaned[k] = '-';
    }
  }

  RACI fin(filename /* can read gzipped results */
#if defined(_MSC_VER) && !defined(_STLP_WIN32) && (_MSC_VER < 1400)
    , ios::nocreate
#endif
    );
  if(! fin) {
    cerr << "error reading " << filename << endl;
    exit(1);
  }
  int last;
  int line_len=0;
  int active_len=0;
  long rangeStart=-1;
  long rangeEnd=0;
  int z = 0;
  char* ptr;

  while(fin.getline(line, line_width)) {
    //last = strlen(line) - 1;
    last = fin.gcount() - 1;
    while(last > 0 && (line[last] < 33 || line[last] > 126)) 
      line[last--] = 0;

    //line_len = strlen(line);
    line_len = last+1;

    if(strstr(line, delimiter) != NULL) { // have new delimiter
      fin.getline(line, line_width);
      //last = strlen(line) - 1;
      last = fin.gcount() - 1;
      while(last > 0 && (line[last] < 33 || line[last] > 126)) 
        line[last--] = 0;
      if(strstr(line, "Content-Type") != NULL) {
        char* result = strstr(line, "name=\"");
        if(result != NULL && strlen(result) > 6) {
          active[0] = 0;
	  active_len=(int)strlen(result)-6;
          strncat(active, result+6, active_len - 1);
          active[active_len] = 0;
        }
        else {
          cerr << "error parsing name for " << line << ": " << result << endl;
          exit(1);
        }
	if (verbose_) {
	  cerr << "active: " << active << endl;
	}
	contentType = mapContentType(active);

	if (CONTENTTYPE_QUERY == contentType) {
          char *indicator = active;
          queryIndex=atoi(indicator+5);
        }
      }
    } // delimiter
    // peptides
    else if(CONTENTTYPE_PEPTIDES == contentType && line_len > 5) {
      //cerr << line << endl;
      if(line[0] == 'q') {
        query[0] = 0;
        int k = 1;
        while(k < line_len && line[k] != '_') {
          query[k - 1] = line[k];
          k++;
        }
        query[k-1] = 0;

        // avoid reading things like
        // "q29_p1_terms=R,-" 
        // as another hit declaration (although, what SHOULD we do with it?)
        const char *eq=strchr(line+k+1,'=');
        const char *ub=strchr(line+k+1,'_');

        // if top hit or next to top
        if ((eq && ((!ub) || (ub>eq))) &&
            (atoi(query) == peptides_.length() + 1 ||
             (atoi(query) == peptides_.length() && line[k+1] == 'p' && (line[k+2] != '1' || line[k+3] != '=')))) {

          nextline[0] = 0;
          // check here if 10 or greater 
          if(line[k+3] == '=')
            strcat(nextline, line + k + 4);
          else
            strcat(nextline, line + k + 5);
          int strlen_tmp =  (int)strlen(nextline);
          for(int z = 0; z < strlen_tmp; z++) {
            if(nextline[z] == ',')
              nextline[z] = ' ';
            else if (nextline[z] == ';')
              break;
          }

          // initialize
          mass = 0.0;
          massdiff = 0.0;
          nextionmatch = 0;
          nextpeptide[0] = 0;
          nextionscore = 0.0;
          nextprotein.clear();
          nextmods[0] = 0;
          discard[0] = 0;
          nextplus = 0;

          std::vector<std::string> hitProteinIds;
          int result = sscanf(nextline, "%d %lf %lf %d %1024s %d %1024s %lf %1024s %d %d;%c",
            &first, &mass, &massdiff, &nextionmatch, nextpeptide, &first, nextmods, &nextionscore, discard, &first, &first, &nextc);
          if (12==result) {
            // Read first protein match
            const char* left = strchr(nextline, ';');
            const char* right = left;
            while (left != NULL && strlen(left)>1) {
	      //if (left != NULL && strlen(left)>1) {
              left++;
              if (*left == '"') {
		left++;
		right = strchr(left, '"');
              } else {
		right = strchr(left, ':');
              }

              if (right == NULL) {
		right = left + strlen(left);
              }
              nextprotein.clear();
              nextprotein.assign(left, right - left);

              //trim the id
              trimRight(nextprotein);
              trimLeft(nextprotein);

              hitProteinIds.push_back(nextprotein);

              left = strchr(right, ',');
              right = left;
            }
            nextplus = hitProteinIds.size()-1;
            if (nextplus>0) {
              nextprotein = hitProteinIds[0];
            }
          }

          // CONSIDER(brendanx): 1.0 is not the true mass of a proton.
          //                     Should we update to be consistent with other tools.
          mass += 1.007826035; // add proton!

          // the top score
          if(line[k+1] == 'p' && line[k+2] == '1' && line[k+3] == '=') { // the correct one to parse
            char *szSharedPeptide = strCopy(nextpeptide);
            peptides_.insertAtEnd(szSharedPeptide);
            proteins_.insertAtEnd(strCopy(nextprotein.data()));
            ionscores_.insertAtEnd(nextionscore);
            masses_.insertAtEnd(mass);
            massdiffs_.insertAtEnd(massdiff);
            pluses_.insertAtEnd(nextplus);
            alt_inds_.insertAtEnd(alt_proteins_.length()); // start
            matchedions_.insertAtEnd(nextionmatch);
            mods_.insertAtEnd(strCopy(nextmods));

            // if completely failed to match, make sure prev_AAs_ and foll_AAs
            // don't get out of synch.
            if (nextpeptide[0] == 0) {
              prev_AAs_.insertAtEnd(G_UNKNOWN_AA);
              foll_AAs_.insertAtEnd(G_UNKNOWN_AA);
            }

            rank1_inds_.insertAtEnd(rank1_proteins_.length()); // start
            std::vector<std::string>::const_iterator iter = hitProteinIds.begin();
            if (hitProteinIds.end() != iter) {
              // we skip the first id as it is kept in proteins_[]
              iter++;
            }
            for (; iter != hitProteinIds.end(); ++iter) {
              rank1_peptides_.insertAtEnd(szSharedPeptide);
              rank1_proteins_.insertAtEnd(strCopy(iter->data()));
            }

          }
          // alt matches
          else {
            alt_peptides_.insertAtEnd(strCopy(nextpeptide));
            alt_proteins_.insertAtEnd(strCopy(nextprotein.data()));
            alt_ionscores_.insertAtEnd(nextionscore);
            alt_masses_.insertAtEnd(mass);
            alt_massdiffs_.insertAtEnd(massdiff);
            alt_pluses_.insertAtEnd(nextplus);
            alt_matchedions_.insertAtEnd(nextionmatch);
            alt_mods_.insertAtEnd(strCopy(nextmods));
          }
        }
        else if ((eq && ub && ub<eq) &&
                 (atoi(query) == prev_AAs_.length() + 1 ||
		  (atoi(query) == prev_AAs_.length() && line[k+1] == 'p' && (line[k+2] != '1' || line[k+3] != '=')))) {
          // q1_p1_terms=K,N:R,R
          if ('t' == ub[1] && 'e' == ub[2] && 'r' == ub[3] && 'm' == ub[4] && 's' == ub[5] && '=' == ub[6])
	    {
	      term_read_ = True;

	      std::vector<char> hitProteinIdsPrevAA;
	      std::vector<char> hitProteinIdsNextAA;

	      char prev_aa = G_UNKNOWN_AA;
	      char next_aa = G_UNKNOWN_AA;
	      int result = sscanf(eq, "=%c,%c", &prev_aa, &next_aa);
	      if (2==result) {
		hitProteinIdsPrevAA.push_back(prev_aa);
		hitProteinIdsNextAA.push_back(next_aa);
		const char* left = strchr(eq, ':');
		while (NULL!=left) {
		  prev_aa = next_aa = G_UNKNOWN_AA;
		  sscanf(left, ":%c,%c", &prev_aa, &next_aa);
		  hitProteinIdsPrevAA.push_back(prev_aa);
		  hitProteinIdsNextAA.push_back(next_aa);
		  left += 4;
		  left = strchr(left, ':');
		}
		prev_aa = hitProteinIdsPrevAA[0];
		next_aa = hitProteinIdsNextAA[0];
	      }

	      // the top score
	      if(line[k+1] == 'p' && line[k+2] == '1' && line[k+3] == '_') { // the correct one to parse
		prev_AAs_.insertAtEnd(prev_aa);
		foll_AAs_.insertAtEnd(next_aa);

		int numberOfAAs = hitProteinIdsPrevAA.size()-1;
		if (numberOfAAs > 0) {
		  if (pluses_[pluses_.length()-1]<numberOfAAs)
		    numberOfAAs = pluses_[pluses_.length()-1];
		  numberOfAAs++;
		  for (int aaIndex=1; aaIndex<numberOfAAs; ++aaIndex) {
		    rank1_prev_AAs_.insertAtEnd(hitProteinIdsPrevAA[aaIndex]);
		    rank1_foll_AAs_.insertAtEnd(hitProteinIdsNextAA[aaIndex]);
		  }
		}
	      }
	      // alt matches
	      else {
		alt_prev_AAs_.insertAtEnd(prev_aa);
		alt_foll_AAs_.insertAtEnd(next_aa);
	      }
	    }
          // q1_p1_subst=6,X,V 
	}
      } // start with 'q'
    } // if peptides are active
    else if(CONTENTTYPE_QUERY == contentType && active_len > 5) {

      // fill in dummy values for missing result (to be removed later)
      while(spectra_.length() < queryIndex-1) {
        spectra_.insertAtEnd(NULL);
      }


      if (queryIndex == spectra_.length() + 1) {
	
	if(line_len > 6 && line[0] == 'I' && line[1] == 'o' && line[2] == 'n' && line[3] == 's' && line[4] == '1' && line[5] == '=') {
	  if (rangeStart < 0) { // no scan numbers found: set start scan to 0000
	    cout << "Warning: could not find scan numbers of spectrum " << nextspectrum << endl;
	    cout << "Set to 0000" << endl;
	  }
	  sprintf(nextspectrum, "%s.%05ld.%05ld.%d", filextn_, rangeStart, rangeEnd, z); 
	  spectra_.insertAtEnd(strCopy(nextspectrum));
	  
	  nextspectrum[0] = 0;
          rangeStart=-1;
          rangeEnd=0;
	  z=0;

	  if (generateTgz && generateDta) {
	    //if (queryIndex>4600)
	    writeDtafile(spectra_.length() - 1, line + 6, line_len-6);
	  }
	}
	else if ((ptr = strstr(line, "scans=")) != NULL) {
	  if (rangeStart <= 0) {
	    rangeStart = atol(ptr+6);
	  }
	  if (rangeEnd <= 0) {
	    rangeEnd = atol(ptr+6);
	  } 
	  
	}
	else if ((ptr = strstr(line, "charge=")) != NULL) {
	  char* sign = strstr(line,"+");
	  if (sign) *sign = '\0';
	  sign = strstr(line,"-");
	  if (sign) *sign = '\0';
	  
	  if (z <= 0) {
	    z = atol(ptr+7);
	  }
	  
	}
	else if(//line_len > 8 && line[0] == 'n' && line[1] == 'u' && line[2] == 'm' && line[3] == '_' &&
		//line[4] == 'u' && line[5] == 's' && line[6] == 'e' && line[7] == 'd') { // gone too far
		line[0] == '-' && line[1] == '-') { // gone too far
	  
          const char generic_spec[] = "spectrum.00000.00000.";
	  const int next_alloc_len = (int)strlen(generic_spec)+12;
          char* next = new char[next_alloc_len];
          if(precursor_ioncharges_.length() > 0)
            sprintf(next, "%d%s%d", spectra_.length()+1, generic_spec, precursor_ioncharges_[spectra_.length()]);
          else
            sprintf(next, "%d%s", spectra_.length()+1, generic_spec);
          spectra_.insertAtEnd(next);

	  nextspectrum[0] = 0;
          rangeStart=-1;
          rangeEnd=0;
	  z=0;


        }
        else if(line_len > 7 && line[0] == 't' && line[1] == 'i' && line[2] == 't' && line[3] == 'l' && line[4] == 'e' && line[5] == '=') {

          nextspectrum[0] = 0;
          int index = 0;
          char next_decode;
          int strlen_tmp = (int)strlen(line);
	  int strlen_tmp_less_2 = strlen_tmp - 2;
	  bool bNeedsFakeSequestStyleSpectrumName = false;
          int k;
          for(k = 6; k < strlen_tmp; k++) {
            if(line[k] == '%' && k < strlen_tmp_less_2) {
              next_decode = decodeURL(line[k+1], line[k+2]);
              if(next_decode > 0) {
                nextspectrum[index++] = next_decode;
              }
	      k += 2; // whether we convert or not, skip two hex chars (thx Matthias Schlesner)
            }
            else if (line[k] == ',') 
              nextspectrum[index++] = '.';
            else if (line[k] != ' ' && line[k] != ':')
              nextspectrum[index++] = line[k];

          }
          //
          // this section submitted by Matthias Schlesner:
          // get scan numbers from Mascot Distiller info
          //
          nextspectrum[index]=0;
          char tmp_line[500];
          char* rangepos ;
          int mlength; //length of match
          char cset[] = "1234567890";
          //rangeStart=-1;
          //rangeEnd=0;
	  //z=0;

          // find start scan
          mlength = 12;
          rangepos = strstr(nextspectrum, "scansinrange"); //sum of several scans
          if (!rangepos) { 
	    rangepos = strstr(nextspectrum, "scans:");
	    
	  }
	  if (!rangepos) {
            rangepos = strstri(nextspectrum, "scans="); // single scan
            mlength = 6;
          }
	  
	  if (!rangepos) {
            rangepos = strstri(nextspectrum, "Scan"); // single scan
            mlength = 4;
          }
	  
	  if (rangepos && (!strncmp(rangepos, "scans:",6) || !strncmp(rangepos, "scans=",6))) {
	    rangeStart = rangeEnd = atol(rangepos+6); //	   	    
	    while (rangepos = strstr(rangepos, ",")) {
	      long scan = atol(rangepos+1);
	      if (scan && scan < rangeStart) {
		rangeStart = scan;
	      }
	      if (scan && scan > rangeEnd) {
		rangeEnd = scan;
	      }
	      rangepos++;
	    }
	    
	    bNeedsFakeSequestStyleSpectrumName = true;
          }
	  else if (rangepos && (!strncmp(rangepos, "scan:",5) || !strncmp(rangepos, "scan=",5))) {
	    rangeStart = rangeEnd = atol(rangepos+5); //	   	    
	    while (rangepos = strstr(rangepos, ",")) {
	      long scan = atol(rangepos+1);
	      if (scan && scan < rangeStart) {
		rangeStart = scan;
	      }
	      if (scan && scan > rangeEnd) {
		rangeEnd = scan;
	      }
	      rangepos++;
	    }
	    
	    bNeedsFakeSequestStyleSpectrumName = true;
          } 
	  else if (nextspectrum && !strncmp(nextspectrum, "ScanNumber:",11)) {
	    rangeStart = rangeEnd = atol(nextspectrum+11); // single scan
	    bNeedsFakeSequestStyleSpectrumName = true;
          } else if (!strncmp(nextspectrum, "spectrumId=",11) ) {
	    rangepos = NULL;
	    rangeStart = rangeEnd = atol(nextspectrum+11); // single scan
	    bNeedsFakeSequestStyleSpectrumName = true;
	    const char *retentiontime = strstr(nextspectrum,"TimeInSeconds=");
	    if (retentiontime) {
	      scanRetentionTimeMap_.insert(std::make_pair(rangeStart,atof(retentiontime+14)));
	    }

          } else
	    //if we had a match and after the match still some characters and the first of them is a number get start scan
	    if(rangepos != NULL && (int) strlen(rangepos) > mlength && (int)strcspn(rangepos, cset) == mlength) {
	      tmp_line[0] = 0;
	      strncat(tmp_line, rangepos + mlength, strlen(rangepos) - mlength - 1);
	      tmp_line[strlen(rangepos) - mlength] = 0;
	      rangeStart=atol(tmp_line);
	    }
	    else if ((rangepos = strstr(nextspectrum,"FinneganScanNumber:"))!=NULL) {
	      rangeStart = rangeEnd = atol(strchr(rangepos,':')+1);
	    } else {
	      //wch 2007-04-24:
	      //      test out the the scan number mapping
	      loadCycleExptScanFromMZXML (filename);
	      parseCyleExperimentForScan (nextspectrum, rangeStart, rangeEnd);
	      if (-1==rangeStart) {
		// bpratt 1-16-05:
		// look for a title constructed as <anything>.start.end.charge
		// look for a title constructed as or <anything>.start.end.charge.<ext> bpratt 1-20-05
		// look for a title constructed as .start.end.charge bpratt 3-17-06
		//
		char *copy= strdup(nextspectrum);
		char *charge = strrchr(copy,'.');
		if (charge) { // is this really start.end.charge format?
		  int dotcount = 3+(!isdigit(*(charge+1))); // or s.e.c.ext
		  char *dot = charge;
		  for (int i=dotcount;--i;) {
		    while (--dot>copy) {
		      if ('.'==*dot) {
			break;
		      }
		    }
		    if ((dot <= copy) || !isdigit(*(dot+1))) {
		      charge = NULL; // this is not start.end.charge format
		      break;
		    }
		  }
		  if (charge) { // still looks like s.e.c, verify
		    int s,e,c;
		    if (3!=sscanf(dot,".%d.%d.%d",&s,&e,&c)) {
		      charge = NULL; // this is not start.end.charge format
		    }
		  }
		}
		if (charge) {
		  *charge++ = 0;
		  if ((*(charge+1)) && !isdigit(*(charge+1))) { // that was .dta or somesuch
		    charge = strrchr(copy,'.');
		    if (charge) {
		      if (isdigit(*(charge+1))) { // that was charge
			*charge = 0;
		      }
		    }
		  }
		  char *end = strrchr(copy,'.');
		  if (end) {
		    *end++ = 0;
		    char *start = strrchr(copy,'.');
		    if (start++) {
		      rangeStart = atol(start);
		      rangeEnd = atol(end);
		      z = atol(charge);
		    }
		    if ((copy+1)==start) { // of the form .start.end.charge
		      // make form _.start.end.charge, others expect the prefix
		      memmove(nextspectrum+1,nextspectrum,strlen(nextspectrum)+1);
		      nextspectrum[0]='_';
		    }
		  }
		}
		free(copy);
	      }
	    }
   
          if (rangeStart < 0) { // no scan numbers found: set start scan to 0000
	    //cout << "Warning: could not find scan numbers of spectrum " << nextspectrum << endl;
            //cout << "Set to 0000" << endl;


            rangeStart = 0;
            mlength = 0;
          }
          // get end scan if we haven't already
          if ((0==rangeEnd) && (mlength == 12)) {  //was scansinrange: end scan after 'to'
            rangepos = strstr(tmp_line, "to");
            // cout << "rangepos: " << rangepos << endl;
            if (rangepos != NULL && strlen(rangepos) > 6) {
              mlength = 2;
              tmp_line[0] = 0;
              strncat(tmp_line, rangepos + mlength, strlen(rangepos) - mlength - 1);
              // cout << "templ: " << tmp_line << endl;
              rangeEnd = atol(tmp_line);
            }
          }


          //no end scan: either was single scan or "to" not found or no scan numbers at all
          // set end scan to start scan
          if (rangeStart > rangeEnd) {
            rangeEnd = rangeStart;
          }
          //
          // end Matthias' scan number stuff (more below)
          //

          // add on extra suffix if necessary here....
          if(strlen(nextspectrum) > 4 && nextspectrum[index-1] == 'a' && nextspectrum[index-2] == 't' && nextspectrum[index-3] == 'd' && nextspectrum[index-4] == '.')
            nextspectrum[index-4] = 0; // truncate

          if(bNeedsFakeSequestStyleSpectrumName || 
	     ((strlen(nextspectrum) > 12 && nextspectrum[strlen(nextspectrum)-2] != '.') || nextspectrum[strlen(nextspectrum)-1] < '0' || nextspectrum[strlen(nextspectrum)-1] > '9')) {
            int n = (int)strlen(nextspectrum);

            // give a standard name here....


            if(strlen(nextspectrum) > 20 || strstr(nextspectrum, ":") || strstr(nextspectrum, ",") ) {
              //sprintf(nextspectrum, "%s.spectrum%d", cleaned, spectrum_counter++);
	      sprintf(nextspectrum, "%s", filextn_);
	    }

            n = (int)strlen(nextspectrum);
            if (nextspectrum[n-1] != '.') {
              nextspectrum[n++] = '.';
            }
            // cout << "Set Range: " << range << endl;
            sprintf(nextspectrum+n,"%05ld.%05ld.",rangeStart,rangeEnd);

            n = (int)strlen(nextspectrum);
            if(spectra_.length() < precursor_ioncharges_.length())
              nextspectrum[n++] = '0' + precursor_ioncharges_[spectra_.length()];
            else
              nextspectrum[n++] = '5';
            nextspectrum[n] = 0;

            if (verbose_) 
              cout << "\t-->" << nextspectrum << endl;
          }
	  //          spectra_.insertAtEnd(strCopy(nextspectrum));
        }
      }
 
    }
    else if(CONTENTTYPE_SUMMARY == contentType) {
      if(line_len > 7 && line[0] == 'q' && line[1] == 'm' && line[2] == 'a' && line[3] == 't' && line[4] == 'c' && line[5] == 'h') {
        query[0] = 0;
        int k = 6;
        while(k < line_len && line[k] != '=') {
          query[k - 6] = line[k];
          k++;
        }
        query[k-6] = 0;
        if(atoi(query) == identityscores_.length() + 1) {
          char* next = strstr(line, "=");
          if(next != NULL && strlen(next) > 1) 
            identityscores_.insertAtEnd(10 * log(atof(next+1)) / log((double)10));
          else 
            identityscores_.insertAtEnd(100.0);
        }
      } // if qmatch
      else if(line_len > 10 && line[0] == 'q' && line[1] == 'p' && line[2] == 'l' && line[3] == 'u' && line[4] == 'g') {
        query[0] = 0;
        int k = 9;
        while(k < (int) line_len && line[k] != '=') {
          query[k - 9] = line[k];
          k++;
        }
        query[k-9] = 0;
        if(atoi(query) == homologyscores_.length() + 1) {
          char* next = strstr(line, "=");
          if(next != NULL && strlen(next) > 1) 
            homologyscores_.insertAtEnd(atof(next+1));
          else 
            homologyscores_.insertAtEnd(100.0);
        }
      } // homology
      else if(line_len > 10 && line[0] == 'q' && line[1] == 'e' && line[2] == 'x' && line[3] == 'p') {
        query[0] = 0;
        int k = 4;
        while(k < (int) line_len && line[k] != '=') {
          query[k - 4] = line[k];
          k++;
        }
        query[k-4] = 0;
        if(atoi(query) == precursor_ioncharges_.length() + 1) {
          char* next = strstr(line, ",");
          if(next != NULL && strlen(next) > 1) 
            precursor_ioncharges_.insertAtEnd((int)(next[1])-'0');
          else 
            precursor_ioncharges_.insertAtEnd(4); // unknown
        }

      } // charge
    } // summary

    // parameters
    else if(CONTENTTYPE_PARAMETERS == contentType) {
      if(line_len > 5) {
        if(strcmp(line+5, "Monoisotopic") == 0)
          mass_type_ = 1;
      }
      if(line_len > 4 && line[0] == 'C' && line[1] == 'L' && line[2] == 'E' && line[3] == '=') {
        // check that not none
        // WCH: we copy this, but use the flag "enzyme_constraint_" to indicate the constraint presence
        //if(strstr(line+4, "None") == NULL)
	enzyme_ = strCopy(line + 4);
      }
      else if(line_len > 4 && line[0] == 'c' && line[1] == 'l' && line[2] == 'e' && line[3] == '=') {
        // check that not none
        // WCH: we copy this, but use the flag "enzyme_constraint_" to indicate the constraint presence
        //if(strstr(line+4, "None") == NULL)
	enzyme_ = strCopy(line + 4);
      }
      // jmt: commenting this out: variable mods are now picked up in "masses" section below
      /*
	else if(strlen(line) > 8 && line[0] == 'I' && line[1] == 'T' && line[2] == '_' && line[3] == 'M' &&
        line[4] == 'O' && line[5] == 'D' && line[6] == 'S' && line[7] == '=') {
        if(enterModification(line+8, 0.0, MOD_ERROR, True, 0))
	cout << "mod: " << line+8 << endl;
	}
      */
      char* nextcomp = strchr(line, '=');
      if(nextcomp != NULL) {
	int nextcomp_len = (int)strlen(nextcomp);
	int nextparam_len = line_len - nextcomp_len;
        char* nextparam = new char[nextparam_len + 1];
        strncpy(nextparam, line, nextparam_len);
        nextparam[nextparam_len] = 0;
        params_.insertAtEnd(nextparam);
        char* nextval = new char[nextcomp_len];
        strcpy(nextval, nextcomp + 1);
        vals_.insertAtEnd(nextval);
      }

    } // params
    // header
    else if(CONTENTTYPE_HEADER == contentType) {
      if(line_len > 8 && line[0] == 'r' && line[1] == 'e' && line[2] == 'l' && line[3] == 'e' && line[4] == 'a' && line[5] == 's' && line[6] == 'e' && line[7] == '=')
	{
	  if(database_ == NULL)
	    database_ = strCopy(line+8);
	}
      else if(line_len > 5 && line[0] == 'F' && line[1] == 'I' && line[2] == 'L' && line[3] == 'E' && line[4] == '=') {
	header_ = strCopy(line+5);
	// now truncate at last '/' or '\'
	char *p = findRightmostPathSeperator(header_);
	if (p) {
	  *p = 0;
	}

      }

    } // end of header section

    ////
    // masses section ('Content-Type: application/x-Mascot; name="masses"')
    ///////////////////////////////////////////////////////////////////////
    else if(CONTENTTYPE_MASSES == contentType) {

      ////
      // Single AA entry
      //////////////////

      if(line_len > 2 && line[1] == '=' && line[0] >= 'A' && line[0] <= 'Z') {
        //char nextaa[2];
        //nextaa[0] = line[0];
        //nextaa[1] = 0;
        //if(enterModification(nextaa, atof(line+2), MOD_ERROR, False, 0))
        //  cout << "next: " << nextaa << " with " << line+2 << endl;

	// jmt: use these single AA masses instead of ResidueMass.cxx
	// entries although they may be modified values, know that the
	// modifications are listed in the "FixedMod" entries that
	// will follow.  Bookkeeping for fixed mods is dealt with
	// there.
	double mass = atof(line + 2);
	recordStaticAAMass(line[0], mass);
	cout << line[0] << " stored as " << mass << endl;
      }


      ////
      // n or c term BASE (unmodded) mass
      /////////////////////////////////
      else if(line_len > 7 && line[1] == '_' && line[2] == 't' && line[3] == 'e' && line[4] == 'r' &&
	      line[5] == 'm' && line[6] == '=') {
        char term;
        term = (line[0]=='C' ? 'c' : 'n');
	double mass=atof(line+7);
	recordStaticAAMass(term, mass);
	cout << term << " stored as " << mass << endl;
      }


      ////
      // variable modificiations ("delta<n>=...")
      ///////////////////////////////////////////

      /* jmt: Actually store the variable modificiation masses from
	 here, directly.  The old version read IT_MODS from an above
	 section, a comma-separated list of names of variable
	 modifications, and ResidueMass.cxx had to be updated to store
	 the actual values for these mods.
      */
      else if (line_len > 7 && line[0] == 'd' && line[1] == 'e' && line[2] == 'l' && line[3] == 't' && line[4] == 'a') {
	// I assume there are max 9 deltas-- exit with error otherwise
	if (line[6] != '=') {
	  // exit with error
	  cerr << "delta mass mod > 9, exiting." << endl;
	  exit(1);
	}
	// we have a line similar to 
	// "delta4=15.994919,Oxidation (HW)"


	double variableModMass=-1;
	char modName[50];

	sscanf(line+7, "%lf,%s", &variableModMass, modName);
	strcpy(modName, strchr(line+7, ',')+1);

	int modNumber = line[5] - '0';
	cout << "variable modification "<< modNumber <<", delta mass " << variableModMass << " for " << modName << ":" << endl;
	
	// extractModificationInfo
	Boolean nterm;
	Boolean cterm;
	Boolean protterm;
	string residues = extractModificationInfo(modName, nterm, cterm, protterm);
	cout << " residues: " << residues << endl;
	cout << " nterm: " << nterm << endl;
	cout << " cterm: " << cterm << endl;
	cout << " protterm: " << protterm << endl;
	cout << endl;
	
	// go through all the affected residues
	for (int i=0; i<(int)residues.length(); i++) {
	  if (residues[i]=='n' && protterm)  // recordVariableAAMod seems to want first arg 'residue'
	    residues[i]='1';                // to be 1 or 2 for protein terminus modification
	  else if (residues[i]=='c' && protterm)
	    residues[i]='2';
	  recordVariableAAMod(residues[i], variableModMass, nterm, cterm, modNumber);
	  if (verbose_)
	    cout << "added variable modification for " << residues[i] << endl;
	}
      }


      ////
      // fixed modifications
      //////////////////////
      else if (line_len > 7 && 
	       0 == strncmp("FixedMod", line, 8)
	       && line[8] != 'R' // *don't* pick up FixedModResidues here!
	       ) {
	/*
	  This is only here for completeness.  Fixed (static)
	  modifications should have already been directly entered in
	  the single AA list above.  The value for the mod in this
	  location should add to the ResidueMass base AA mass to equal
	  the (hopefully already recorded) static AA mod above.
	 */
	
	// find the equals sign, and grab the number of the fixed mod
	char* eq = strchr(line, '=');
	int digits = eq - (line+8);
	char digitsStr[10];
	strncpy(digitsStr, line+8, digits);
	digitsStr[digits] = 0;
	fixedModNumber = atoi(digitsStr);
	if (verbose_)
	  cout << "got " << digitsStr << "," << fixedModNumber << "th fixed mod" << endl;
	sscanf(eq+1, "%lf", &fixedModMass);
	if (verbose_)
	  cout << "mod value: " << fixedModMass << endl;
	char* comma = strchr(eq, ',');
	strcpy(fixedModDesc,comma+1);
	if (verbose_)
	  cout << "desc: " << fixedModDesc << endl;
	lookingForFixedModResidues = true;
      }
      ////
      // fixed modification residues
      ///////////////////////////////
      else if (line_len > 7 && 
	       0 == strncmp("FixedModResidues", line, 16)
	       ) {

	/*
	  The value for the mod in this location should add to the
	  ResidueMass base AA mass to equal the (hopefully already
	  recorded) static AA mod above.
	*/
	lookingForFixedModResidues = false;
	int fmNum = -1;
	char* eq = strchr(line, '=');
	int digits = eq - (line+16);
	char digitsStr[10];
	strncpy(digitsStr, line+16, digits);
	digitsStr[digits]=0;
	fmNum = atoi(digitsStr);
	if (fmNum != fixedModNumber) {
	  cerr << "error parsing FixedModResidues: " << line << endl;
	  cerr << "unexpected number (expected " << fixedModNumber << ")" << endl;
	  exit(1);
	}

	// extractModificationInfo
	Boolean nterm;
	Boolean cterm;
	Boolean protterm;
	string residues = extractModificationInfo(fixedModDesc, nterm, cterm, protterm);
	if (verbose_) {
	  cout << " residues: " << residues << endl;
	  cout << " nterm: " << nterm << endl;
	  cout << " cterm: " << cterm << endl;
	  cout << " protterm: " << protterm << endl;
	  cout << endl;
	}
	
	// go through all the affected residues
	for (int i=0; i<(int)residues.length(); i++) {
	  recordStaticAAMod(residues[i], fixedModMass, nterm, cterm);
	  if (verbose_)
	    cout << "added fixed modification for " << residues[i] << endl;
	}

      }

    } // end  masses section

    // proteins
    else if(generate_description_ && CONTENTTYPE_PROTEINS == contentType) {
      char *eq=strchr(line,'=');
      if (eq) {
        // let process the line
        // "<id>"=mass,"description"
        char *descriptionStart=strchr(eq+1,'"');
        if (descriptionStart)
	  {
	    char *descriptionEnd=strchr(descriptionStart+1,'"');
	    if (descriptionEnd)
	      {
		*(eq-1) = '\0';
		*descriptionEnd = '\0';
		std::string id(line+1);
		trimAll(id);
		std::string description((descriptionStart+1));
		trimAll(description);
		if (!id.empty())
		  proteinDescriptionMap_.insert(std::make_pair(id,description));
	      }
	  }
      }
    }

  } // next input line


  fin.close();

  assert(proteins_.length() == peptides_.length());

  if(enzyme_ != NULL) {
    char *lcEnzyme = strCopy(enzyme_);
    strlwr(lcEnzyme);
    if (strstr(lcEnzyme, "none") == NULL)
      enzyme_constraint_ = True;
  }

  if(enzyme_constraint_ && verbose_)
    cerr << enzyme_ << "...." << endl;


  if(! validate_db()) {
    cerr << "could not find database ";
    if(database_ != NULL)
      cerr << database_;
    cerr << endl;
    exit(1);
  }
  if (verbose_) {
    for(k = 0; k < static_mods_.length(); k++) {
      cout << "recorded static mods:" << endl;
      static_mods_[k]->print();
      cout << endl;
    }
    for(k = 0; k < variable_mods_.length(); k++) {
      cout << "recorded variable mods:" << endl;
      variable_mods_[k]->print();
      cout << endl;
    }
  }
  // register all static mods here so won't be used in standard modified peptide names
  if(EXCLUDE_CONST_STATICS_) { // must load them with static mod's (not to be displayed)
    StaticModificationCount next_static_count;
    for(k = 0; k < static_mods_.length(); k++) {
      next_static_count.mod = static_mods_[k]->aa;
      next_static_count.mass = static_mods_[k]->mass;
      next_static_count.num = 1;
      static_consts_.insertAtEnd(next_static_count);
    }
  }


  // now remove all the entries without true result
  for(k = 0; k < spectra_.length(); k++) {
    if(spectra_[k] == NULL) {
      peptides_.remove(k);
      proteins_.remove(k);
      ionscores_.remove(k);
      homologyscores_.remove(k);
      identityscores_.remove(k);
      masses_.remove(k);
      massdiffs_.remove(k);
      pluses_.remove(k);
      alt_inds_.remove(k); // start
      matchedions_.remove(k);
      mods_.remove(k);
      spectra_.remove(k);
      precursor_ioncharges_.remove(k);
      rank1_inds_.remove(k);
      k--;
    } // modify spectrum name to be consistent with precursor ion charge
    else if(strlen(spectra_[k]) > 2 && (spectra_[k])[strlen(spectra_[k])-1] != precursor_ioncharges_[k] + '0') {
      // cout << "note: input spectrum name \"" << spectra_[k] << "\" is inconsistent with stated charge " << precursor_ioncharges_[k];
      (spectra_[k])[strlen(spectra_[k])-1] = precursor_ioncharges_[k] + '0';
      // cout << ", changing to \"" << spectra_[k] << "\"" << endl;
    }
  }

  if (identityscores_.length()!=ionscores_.length()) {
    cerr << "Error: " << identityscores_.length() << " identity scores, but " << ionscores_.length() << " ion scores read" << endl;
    exit(1);
  }

  // ASMS Workshop and User Meeting 2005
  // as defined in http://www.matrixscience.com/pdf/2005WKSHP4.pdf
  // E = P{threshold} * (10 ** ((S{threshold}-score)/ 10))
  for(k = 0; k < spectra_.length(); k++) {
    expectationvalues_.insertAtEnd(0.05 * pow(10.0,((identityscores_[k]-ionscores_[k])/10.0)));
  }

  // get the sample enzyme instance as it is used in .out, .pep.xml and .html
  // databaseSearch() could take a long time, fail earlier if the enzymes are not available
  sampleEnzymeInstance_ = getSampleEnzyme (true);
  searchEnzymeInstance_ = getSearchEnzyme (true);

  databaseSearch();

  char announce[400];
  int length=0;

  // skip, if drawing the information from .mzXML file
  if (generateTgz) {
    for(k = 0; k < proteins_.length(); k++) {
      //cout << k << ": " << proteins_[k] << endl;
      if(checkResult(k)) {
        announce[0] = 0;
        length = sprintf(announce, " %d.   opening %s.out    ", k+1, spectra_[k]);
	if (-1==length) length = 0;
        //length = strlen(announce); //17;
        printf("%s", announce);
        writeOutfile(k);
        for(int j = 0; j < length; j++)
          printf("\b");
      }
    }
    printf("\n");

    tarFiles();

  }
  else
    {
      cerr << "Skipping .dta, .out and .tgz generation" << endl;
    }

  char* outfile = new char[strlen(filename) + 16];
  strcpy(outfile, filename);
  char datsuffix[] = ".dat";
  const char* getSuffix = strstr(filename, datsuffix);
  if(getSuffix != NULL && strlen(getSuffix) == strlen(datsuffix)) {
    outfile[strlen(filename) - strlen(datsuffix)] = 0;
  }
  else {
    if (outfile) delete [] outfile;
    outfile = new char[strlen(filename) + 16];
    strcpy(outfile, filename);
  }


  if(xml) {
    strcat(outfile, get_pepxml_dot_ext());
    outPepXMLfileName_ = strCopy(outfile);
    write_pepXML(outfile);
  }
  else {
    strcat(outfile, ".html");
    outPepXMLfileName_ = strCopy(outfile);
    writeHTML(outfile);
  }

  if (bVerbose) {
    cerr << " results written to file " << outfile << endl << endl;
    
    cout << "read in a total of " << spectra_.length() << " spectra, " << ionscores_.length() << " ion scores, and " << peptides_.length() << " peptides and " << pluses_.length() << " pluses, " << homologyscores_.length() << " hom scores, " << identityscores_.length() << " id scores" << endl;
    //exit(1);
    for(k = 0; k < spectra_.length(); k++)
      cout << k+1 << " " << spectra_[k] << " " << ionscores_[k] << " " << identityscores_[k] << " " << homologyscores_[k] << " " << prev_AAs_[k] << "." << peptides_[k] << "." << foll_AAs_[k] << " " << proteins_[k] << " " << pluses_[k] << endl;
  }
  delete [] line;
  delete [] nextline;
  delete [] cleaned;
  delete [] outfile;
}

char MascotConverter::decodeURL(char first, char second) {
  if(first == '2') {
    if(second == 'e')
      return '.';
    if(second == '0')
      return 0; // NULL
    if(second == 'c')
      return ',';
    if(second == 'd')
      return '-';
    if(second == '8')
      return '(';
    if(second == '9')
      return ')';
  }
  if(first == '3') {
    if(second == 'a')
      return ':';
    if(second == 'b')
      return ';';
    if(second == 'c')
      return '<';
    if(second == 'd')
      return '=';

  }
  return 0; // NULL
}



/******************
 extractModificationInfo

 given strings like 
   Me-ester (C-term)

 returns string of residues
*/

string MascotConverter::extractModificationInfo(const string & inputString, // input
						Boolean& nterm, // output
						Boolean& cterm, // output
						Boolean& proteinterm // output
						) {

  proteinterm = inputString.find("Protein") != string::npos;
  nterm= inputString.find("N-term") != string::npos;
  cterm = inputString.find("C-term") != string::npos;

  //Dave@UPENN, Alexandra.panchaud, Fitzgibbon Matthew
  //parenthesis is allowed in modification name
  //E.g.: ICAT-C:13C(9) (C)
  int resStartPos = inputString.find_last_of('(') + 1;
  string affectedRes = inputString.substr(resStartPos, inputString.find_first_of(')', resStartPos) - resStartPos);
  // assert != npos

  string residues; // the string we'll return, containing "ABC" etc

  // if we have exactly "N-term" or "C-term": these translate to AAs "n" or "c"
  if ((affectedRes.length() == 7  || affectedRes.length() == 14) && (proteinterm)) {
    // jmt: following ResidueMass in assuming "Protein" is n term
    // jke: affectedRes could be "Protein N-term"
    residues = "n";
  }
  else if (affectedRes.length() == 6 && (nterm || cterm)) {

    if ((nterm && cterm)) {
      cerr << "error: both c and n term mod detected" << endl;
      exit(1);
    }

    if (nterm) {
      residues = "n";
    }
    else if (cterm) {
      residues = "c";
    }
  }

  // for "N-term Q" we'd just like to have Q,
  else if ( (nterm || cterm) && affectedRes.length() > 6) {
    // so take everything after "N|C-term " text: +7 
    residues = affectedRes.substr(7);
  }
  else {
    residues = affectedRes;
  }
  //cout << inputString << " got following residues: " << residues << endl;

  if(residues.length() == 0) {
    cout << "error: modification " << inputString << " not recognized" << endl;
    exit(1);
  }

  return residues;
}



/*******************
  recordStaticAAMass

  stores the working fixed AA mass-- this may be modified or not (we
  don't know at this point)

  exits with error on faliure
*/
void MascotConverter::recordStaticAAMass(char AA, double mass) {
  if (workingAAMasses[AA] >= 0) {
    cerr << "error: trying to store static AA mass for " << AA  << " when value already stored." << endl;
    exit(1);
  }
  workingAAMasses[AA] = mass;

  switch (AA) {
  case 'B':
  case 'J':
  case 'O':
  case 'U':
  case 'X':
  case 'Z':
    // a value has been recorded for a non-standard AA;
    // record this as a static mod
    // !!!
    // assume it's not nterm or cterm!
    recordStaticAAMod(AA, mass, false, false);
    break;
  default:
    // just return normally
    break;
  }

}


/******************
  recordStaticAAMod

  records a static modification for a specific AA; give the delta from
  the unmodified mass; we will now assume that the working AA mass has
  been modified by this delta, so the base (unmodded) AA mass is
  (working - delta)

  example: Carbamidomethyl (C) will have delta of 57.021469, which
  implies that the working AA for C, 160.030649, is modified, and that
  the base mass for C is (160.030649-57.021469)=103.009180, as we'd
  expect.

  exits with error on faliure.
 */
void MascotConverter::recordStaticAAMod(char AA, double fixedDelta, Boolean nterm, Boolean cterm) {

  // first, make sure we've already recorded *some* mass for the AA
  if (workingAAMasses[AA] < 0) {
    cerr << "error: trying to store static AA modification for " << AA  << " without any prior info for this AA." << endl;
    exit(1);
  }

  // have we tried to modify this AA already?
  if (aaMod[AA] == True) {
    cerr << "error: AA " << AA  << " already modified!" << endl;
    exit(1);
  }

  aaMod[AA] = True;

  double modifiedMass = workingAAMasses[AA];
  // double baseMass = modifiedMass - fixedDelta;
  double delta = fixedDelta;

  // store with index 0, as the index is used to correlate variable mods
  // with the variable mod strings in the peptide ID section and 0 is unused.
  MascotModification* mod = new MascotModification(AA, modifiedMass, delta, False, 0, false, false);
  static_mods_.insertAtEnd(mod);
}



/********************
  recordVariableAAMod

  records a static modification for a specific AA; give the delta from
  the unmodified mass; we will now assume that the working AA mass *can be*
  modified by this delta, so the (variable/optional) modified AA mass would be 
  (working + delta).

  Note: working could be mod or unmod.

  example: Oxidation (M) has a delta of 15.994919, so the variable
  modded mass for M would be base mass 131.040480 + 15.994919 = modded
  mass of 147.035399.

  exits with error on faliure.
 */
void MascotConverter::recordVariableAAMod(char AA, double variableDelta, Boolean nterm, Boolean cterm, int variableIndex) {
  // first, make sure we've already recorded *some* mass for the AA
  if (workingAAMasses[AA] < 0) {
    cerr << "error: trying to store variable AA modification for " << AA  << " without any prior info for this AA." << endl;
    exit(1);
  }

  // note, we record the modified mass as delta+working, where working may be modified or not-- we don't care.
  double baseMass = workingAAMasses[AA]; 
  double modifiedMass = baseMass + variableDelta;

  // jmt: ???
  // jmt: should the delta be based on the working mass if it's a static mod, or the base mass?
  // jmt: ???
  double delta = variableDelta;
  
  // the index is used to correlate variable mods with the variable
  // mod strings in the peptide ID section
  MascotModification* mod = new MascotModification(AA, modifiedMass, delta, True, variableIndex, nterm, cterm);
  variable_mods_.insertAtEnd(mod);  
}


double MascotConverter::getModifiedMass(char res, Boolean staticmod, int index) {
  MascotModification* mod;
  return  getModifiedMass(res, staticmod, index, mod);
}

double MascotConverter::getModifiedMass(char res, Boolean staticmod, int index, MascotModification*& mod) {
  mod = NULL;
  if(staticmod) {
    //cout << "jmt: looking up static mod match to " << res << ": ";
    //cout <<"      len smods: " << static_mods_.length() << endl;
    for(int k = 0; k < static_mods_.length(); k++) {
      //cout << ">>" << static_mods_[k]->aa << "<<" << endl;
      if(static_mods_[k]->aa == res) {
        mod = static_mods_[k]; // was "mod = variable_mods_[k];" 7-27-06 bpratt 
	//cout << "jmt: found mass: " << static_mods_[k]->mass << endl;
        return static_mods_[k]->mass;
      }
    }
    //cout << "jmt: did not find, returning 0" << endl;
    return 0.0; // nothing
  }
  for(int k = 0; k < variable_mods_.length(); k++) {
    if((variable_mods_[k]->aa == res || 
	(variable_mods_[k]->n_terminal && res == 'n') || 
	(variable_mods_[k]->c_terminal && res == 'c' )) && 
       (index < 0 || variable_mods_[k]->index == index)) {
      mod = variable_mods_[k];
      return variable_mods_[k]->mass;
    }
  }
  return 0.0; // nothing
}

// Section modified by Matthias Schlesner to perform a single Mascot search in quantitative proteomics experiments
// by defining light mod as fixed and heavy mod as variable and having Xpress and ASAPRatio working
//
Array<Tag*>* MascotConverter::getModificationTags(char* peptide, char* modstring) {
  Boolean verbose = false; 
  if(verbose)
    cout << "pep: " << peptide << " mod: " << modstring << " ";

  if(! strlen(peptide))
    return NULL;

  if(strlen(peptide) != strlen(modstring) - 2) {
    cout << "error with mods: " << peptide << " and " << modstring << endl;
    return NULL;
    //exit(1);
  }
  if(static_mods_.length() == 0 && variable_mods_.length() == 0)
    return NULL; // no mods possible


  // first check to see whether could be null
  Boolean modified = getModifiedMass('n', True, -1) || getModifiedMass('c', True, -1);
  int k;
  for(k = 0; k < static_mods_.length(); k++) {
    if(static_mods_[k]->aa == 'n' || static_mods_[k]->aa == 'c' || 
       strchr(peptide, static_mods_[k]->aa) != NULL) {
      modified = True;
      k = static_mods_.length();
    }
  }
  if(! modified) { // now check variable
    for(int k = 0; k < (int) strlen(modstring); k++)
      if(modstring[k] != '0') {
        modified = True;
        k = strlen(modstring);
      }
  }
  if(! modified)
    return NULL; // nothing to do
  char text[100];
  Array<Tag*>* output = new Array<Tag*>;
  Tag* next = NULL; 
  Tag* modInfo = new Tag("modification_info", True, False);
  output->insertAtEnd(modInfo);
  //next = new Tag("modification_info", True, False);
  MascotModification* mod;
  double nextmod;  
  if(modstring[0] != '0') { // get it
    nextmod = getModifiedMass('n', False, modstring[0]-'0', mod); // was "modstring[0]+48", thx Matthias Schlesner
    if(nextmod > 0.0) { //variable
      sprintf(text, "%f", nextmod);
      if (mod && mod->n_terminal && mod->aa != 'n') {
	next = new Tag("mod_aminoacid_mass", True, True);
	sprintf(text, "%d", 1);
	next->setAttributeValue("position", text);
	sprintf(text, "%f", nextmod);
	next->setAttributeValue("mass", text);
      }
      else {
	//next = new Tag("modification_info", True, True);
	sprintf(text, "%f", nextmod);
	modInfo->setAttributeValue("mod_nterm_mass", text);
      }
    }
    //    if(nextmod > 0.0) {
    //  sprintf(text, "%f", nextmod);
    //  next->setAttributeValue("mod_nterm_mass", text);
    //}
  } else {
    nextmod = getModifiedMass('n', True, 0, mod);
    next = NULL;
    if(nextmod > 0.0) { //static
      sprintf(text, "%f", nextmod);
      if (mod && mod->n_terminal && mod->aa != 'n') {
        next = new Tag("mod_aminoacid_mass", True, True);
        sprintf(text, "%d", 1);
        next->setAttributeValue("position", text);
        sprintf(text, "%f", nextmod);
        next->setAttributeValue("mass", text);
      }
      else {
        //next = new Tag("modification_info", True, True);
        sprintf(text, "%f", nextmod);
        modInfo->setAttributeValue("mod_nterm_mass", text);
      }
    }
  }

  if(modstring[strlen(modstring)-1] != '0') { // get it
    nextmod = getModifiedMass('c', False, modstring[strlen(modstring)-1]-'0', mod);
    if(nextmod > 0.0) { //variable
      sprintf(text, "%f", nextmod);
      if (mod->c_terminal && mod->aa != 'c') {
	next = new Tag("mod_aminoacid_mass", True, True);
	sprintf(text, "%d", (int)(strlen(modstring)-2));
	next->setAttributeValue("position", text);
	sprintf(text, "%f", nextmod);
	next->setAttributeValue("mass", text);
      }
      else {
	//next = new Tag("modification_info", True, True);
	sprintf(text, "%f", nextmod);
	modInfo->setAttributeValue("mod_cterm_mass", text);
      }
    }
    //    if(nextmod > 0.0) {
    //  sprintf(text, "%f", nextmod);
    //  next->setAttributeValue("mod_cterm_mass", text);
    //}
  } else {	
    //  nextmod = getModifiedMass('c', True, 0, mod);
    nextmod = getModifiedMass('c', True, 0, mod);
    if(nextmod > 0.0) { //static
      sprintf(text, "%f", nextmod);
      if (mod->c_terminal && mod->aa != 'c') {
        next = new Tag("mod_aminoacid_mass", True, True);
        sprintf(text, "%d", (int)(strlen(modstring)-2));
        next->setAttributeValue("position", text);
        sprintf(text, "%f", nextmod);
        next->setAttributeValue("mass", text);
      }
      else {
        //next = new Tag("modification_info", True, True);
        sprintf(text, "%f", nextmod);
        modInfo->setAttributeValue("mod_cterm_mass", text);
      }
    }
  }

  if (next != NULL) 
    output->insertAtEnd(next);
  
  // now look for mod aas
  for(k = 1; k < (int) strlen(modstring)-1; k++) {
    if(modstring[k] != '0') { // variable
      nextmod = getModifiedMass(peptide[k-1], False, modstring[k]-'0');
      if(nextmod > 0.0) {
	if(verbose_)
	  cout << "mod at position: " << k << " ";
	next = new Tag("mod_aminoacid_mass", True, True);
	sprintf(text, "%d", k);
	next->setAttributeValue("position", text);
	sprintf(text, "%f", nextmod);
	next->setAttributeValue("mass", text);
	output->insertAtEnd(next);
      }
    } else {
      nextmod = getModifiedMass(peptide[k-1], True, 0);
      if(nextmod > 0.0) { //static
        if(verbose_)
          cout << "mod at position: " << k << " ";
        next = new Tag("mod_aminoacid_mass", True, True);
        sprintf(text, "%d", k);
        next->setAttributeValue("position", text);
        sprintf(text, "%f", nextmod);
        next->setAttributeValue("mass", text);
        output->insertAtEnd(next);
      }
    }

  } //  next position in peptide
  next = new Tag("modification_info", False, True);
  output->insertAtEnd(next);
  return output;
}

// End of modified section by Matthias Schlesner
//


void MascotConverter::setFileInfo(const char* file) {
  //cerr << "file: " << file << endl;
  char suff[] = ".dat";
  
  // check if already have full path
  if(isAbsolutePath(file)) {
    fullpath_ = new char[strlen(file)+1]; // max
    strcpy(fullpath_, file);
    // strip off filename
    char *slash = findRightmostPathSeperator(fullpath_);
    int stop = slash?(slash-fullpath_):-1;
    fullpath_[stop+1] = 0; // set it here
    if(strlen(file) > 4) {
      const char* match = strstr(file, suff);
      if(match != NULL && strlen(match) == strlen(suff)) {
        filextn_ = new char[strlen(file) - stop - strlen(suff)];
        strncpy(filextn_, file + stop + 1, strlen(file) - stop - strlen(suff) - 1);
        filextn_[strlen(file) - stop - strlen(suff) - 1] = 0;
      }
      else { // use the whole thing
        filextn_ = new char[strlen(file) - stop];
        strcpy(filextn_, file + stop + 1);
      }
    }
    else { // use the whole thing
      filextn_ = new char[strlen(file) - stop];
      strcpy(filextn_, file + stop + 1);
    }
  }
  else {
    char curr_dir[5000];
    safepath_getcwd(curr_dir, sizeof(curr_dir));
    fullpath_ = new char[strlen(curr_dir) + 2];
    strcpy(fullpath_, curr_dir);
    strcat(fullpath_, "/");

    if(strlen(file) > 4) {
      const char* match = strstr(file, suff);
      if(match != NULL && strlen(match) == strlen(suff)) {
        filextn_ = new char[strlen(file) - strlen(suff) + 1];
        strncpy(filextn_, file, strlen(file) - strlen(suff));
        filextn_[strlen(file) - strlen(suff)] = 0;
      }
      else { // use the whole thing
        filextn_ = new char[strlen(file) + 1];
        strcpy(filextn_, file);
      }
    }
    else { // use the whole thing
      filextn_ = new char[strlen(file) + 1];
      strcpy(filextn_, file);
    }

  }
  unCygwinify(fullpath_); // no effect in Cygwin builds
  cerr << "filepath: " << fullpath_ << ", extn: " << filextn_ << endl;
  //exit(1);
}

void MascotConverter::tarFiles() {
  tarball_->commit(tarfilename_,directoryname_);

  if (LEARN_TEST==testmode_) { // copy tar file for later regression test
    copy_file(tarfilename_,regressionTestTarFileName_.c_str());
  }

}

bool MascotConverter::writeDtafile(int index, const char* line, int line_len) {
  bool result = false;
  char suffix[] = ".dta";
  if(index < 0 || index >= spectra_.length()) {
    cerr << "error with index " << index << endl;
    exit(1);
  }

  if(! checkResult(index))
    return result;

  dtafilename_[0] = 0;
  dtafilename2_[0] = 0;
  strcpy(dtafilename_, directoryname_);
  strcat(dtafilename_, spectra_[index]);
  int dtafilename_len = directoryname_len_ + (int)strlen(spectra_[index]);
  
  //DDS: make dta for 2+ and 3+
  if (dtafilename_[dtafilename_len-2] == '.' &&
      dtafilename_[dtafilename_len-1] == '2') {
    strcpy(dtafilename2_, dtafilename_);
    dtafilename2_[dtafilename_len-1] = '3';
    strcat(dtafilename2_, suffix);
  }
  strcat(dtafilename_, suffix);
  // check if exists already...
  if (written_dtas_.find(dtafilename_+directoryname_len_)!=written_dtas_.end()) {
    return result; // nothing to do
  }
  
  if (dtafilename2_[0]) {
    // does a legit +3 file already exist?
    if (written_dtas_.find(dtafilename2_+directoryname_len_)!=written_dtas_.end()) {
      dtafilename2_[0] = 0;  // don't overwrite legit +3 file
    }
  }

  // track what we've written
  written_dtas_.insert(std::make_pair(dtafilename_+directoryname_len_,' '));
  char buf[2048];

  snprintf(buf, sizeof(buf), "%0.2f %c\n", masses_[index] + massdiffs_[index], (spectra_[index])[strlen(spectra_[index])-1]);
  std::string header1(buf);
  if (dtafilename2_[0]) {
    snprintf(buf, sizeof(buf), "%0.2f %c\n", masses_[index] + massdiffs_[index], (spectra_[index])[strlen(spectra_[index])-1]+1);
  }
  std::string header2(buf);
  std::string output;
  
  int start = 0;
  Array<double>* masses = new Array<double>;
  Array<double>* intens = new Array<double>;
  //int line_len=strlen(line);
  int k;
  for(k = 0; line[k]; k++) {
    if(line[k] == ':') {
      masses->insertAtEnd(atof(line+start));
      assert(masses->length() == intens->length()+1);
      start = k+1;
    }
    else if(line[k] == ',' || line[k] == '\n') {
      intens->insertAtEnd(atof(line+start));
      start = k+1;
      assert(masses->length() == intens->length());
    }
  } // next position
  if(start < line_len) { // one last intens
    intens->insertAtEnd(atof(line+start));
  }
  assert(masses->length() == intens->length());
  SpectrumPeak** spec_pks = new SpectrumPeak*[masses->length()];
  for(k = 0; k < masses->length(); k++) {
    spec_pks[k] = new SpectrumPeak((*masses)[k], (*intens)[k]);
  }
  qsort(spec_pks, masses->length(), sizeof(SpectrumPeak*), (int(*)(const void*, const void*)) comp_specpks);

  for(k = 0; k < masses->length(); k++) {
    // chee hong: for Analyst centroid, the intensity needs 4 decimal places
    snprintf(buf, sizeof(buf), "%0.4f %0.4f\n", spec_pks[k]->mass_, spec_pks[k]->intens_);
    output.append(buf);
  }

  header1+=output;
  tarball_->create_memberfile(tarfilename_,dtafilename_,header1);
  if (dtafilename2_[0]) {
    header2+=output;
    // this may be overwritten by an actual +3 later...
    tarball_->create_tentative_memberfile(tarfilename_,dtafilename2_,header2);
  }

#define DTA_TAR_CHECK_FREQ 500
#define DTA_TAR_CHECK_FIRST 50
  if (((RUN_TEST==testmode_) || (FORCE_TEST==testmode_)) && 
      ((index < DTA_TAR_CHECK_FIRST) || !(index%DTA_TAR_CHECK_FREQ))) // doing each one is too slow...
    { // compare with tar file saved earlier for regression test
      compare_tarred_file(header1,dtafilename_);
    }

  for(k = 0; k < masses->length(); k++)
    if(spec_pks[k] != NULL)
      delete spec_pks[k];
  if(spec_pks != NULL)
    delete [] spec_pks;

  if(masses != NULL)
    delete masses;
  if(intens != NULL)
    delete intens;

  return true;
}

bool MascotConverter::compare_tarred_file(std::string &text,const char *tarredfilename) {
  bool result = false;
#define BBUFLEN 1024
  char tarfile[BBUFLEN];
  strncpy(tarfile,regressionTestTarFileName_.c_str(),BBUFLEN);
  char tarredfile[BBUFLEN];
  strncpy(tarredfile,tarredfilename,BBUFLEN);
  FILE *ppIn = 
    read_dta_or_out_from_tgz_file(tarfile,tarredfile, 
				  tarredfile,BBUFLEN);
  if (!ppIn) {
    printf(" Error - can't read file %s from %s\n\n", tarredfilename,tarredfile);
    exit(EXIT_FAILURE);
  }
#define INBUFSIZE 30000
  char *buf = (char *)malloc(INBUFSIZE);
  int ntotal=0,nlines=0,lastpos=0;
  while (fgets(buf,INBUFSIZE,ppIn) && !result) {
    nlines++;
    int nread = (int)strlen(buf);
    std::string subtext = text.substr(lastpos);
    int step = (int)subtext.find('\n')+1;
    subtext.erase(step);
    lastpos+=step;
    // remove double whitespace
    while (char *cp=strstr(buf,"  ")) {
      memmove(cp,cp+1,strlen(cp));
    }
    while (char *cp=strstr(buf,"--")) {
      memmove(cp,cp+1,strlen(cp));
    }
    while ((step=(int)subtext.find("  "))>=0) {
      subtext.erase(step,1);
    }
    while ((step=(int)subtext.find("--"))>=0) {
      subtext.erase(step,1);
    }

    if (strcasecmp(buf,subtext.c_str())) {
      // they differ - just a path style difference, perhaps?
      char *txt = strdup(subtext.c_str());
      char * diff = buf;
      char *cmp = txt;
      Boolean hope;
      do {
	// advance to next difference
	while (*diff && *cmp && (tolower(*diff) == tolower(*cmp))) {
	  diff++;
	  cmp++;
	}
	// file path difference?
	force_unCygwinify(diff); // this works even if this is a cygwin build
	hope = false;
	if (*diff && (diff[1]==':')&&(diff[2]=='/')) {
	  memmove(diff,diff+2,strlen(diff+1));
	  hope = true;
	}
	force_unCygwinify(cmp); // this works if this is a cygwin build
	if (*cmp && (cmp[1]==':')&&(cmp[2]=='/')) {
	  memmove(cmp,cmp+2,strlen(cmp+1));
	  hope = true;
	}
	if (!hope) { // just an install path difference?
	  char *diff_cgi = strstr(diff,".cgi");
	  char *cmp_cgi = strstr(cmp,".cgi");
	  if (diff_cgi && cmp_cgi) {
	    memmove(diff,diff_cgi,strlen(diff_cgi)+1);
	    memmove(cmp,cmp_cgi,strlen(cmp_cgi)+1);
	    hope = 1;
	  }
	}
      } while (hope && *diff && *cmp && strcasecmp(diff,cmp));
      bool bad = (0!=strcasecmp(buf,txt));
      if (bad) {
	bad = (0!=TagListComparator::check_difference(tarredfilename,buf,txt));
      }
      if (bad) { // look for numeric difference
	char *b=buf;
	char *t=txt;
	char *lastnondigit=b-1;
	while (*b && *t && (*b==*t)) {
	  if (!strchr("+-.",*b) && !isdigit(*b)) {
	    lastnondigit = b;
	  }
	  b++;
	  t++;
	}
	double bd=atof(lastnondigit+1);
	double bt=atof(txt+(lastnondigit-buf)+1);
	if (bd && bt) {
	  bad = ((fabs(bd-bt)/Max(bd,bt)) > .01);
	}
      }
      if (bad) {
	cerr << "regression test failure in " << tarredfilename << " at line #"
	     << nlines << endl << " expected: " << endl << buf << " got:" << endl
	     << txt <<  "difference begins at:" << endl << diff << endl;
	result = true;
      }
      free(txt);
    }
    ntotal+=nread;
  }
  free(buf);
  close_dta_or_out_from_tgz_file(ppIn);
  return result;
}

int comp_specpks(const void* num1, const void* num2) {
  SpectrumPeak** dd1 = (SpectrumPeak**)num1;
  SpectrumPeak** dd2 = (SpectrumPeak**)num2;
  if((*dd1)->mass_ < (*dd2)->mass_)
    return -1;
  if((*dd1)->mass_ == (*dd2)->mass_)
    return 0;
  return 1;
}

bool MascotConverter::writeOutfile(int index) {
  bool result = false;
  char suffix[] = ".out";
  if(index < 0 || index >= spectra_.length()) {
    cerr << "error with index " << index << endl;
    exit(1);
  }
  outfilename_[0] = 0;
  strcat(outfilename_, directoryname_);
  strcat(outfilename_, spectra_[index]);
  strcat(outfilename_, suffix);

  ModificationInfo* modinfo = NULL;
  Array<Tag*>* modtags = NULL;
  char* stdpep = NULL;

  int peplen = (int)strlen(peptides_[index]);
  if(! peplen)
    return result;
  if (!sampleEnzymeInstance_) {
    sampleEnzymeInstance_ = getSampleEnzyme (true);
  }

  int stop = alt_proteins_.length();
  if(index < alt_inds_.length() - 1)
    stop = alt_inds_[index+1];

  int max_prot_len = (int)strlen(proteins_[index]);

  double calc_pi=-1;

  int k;
  int alt_protein_len;
  for(k = alt_inds_[index]; k < stop; k++) {
    alt_protein_len = (int)strlen(alt_proteins_[k]);
    if(alt_protein_len  > max_prot_len)
      max_prot_len = alt_protein_len;
  }

  char *spaces = (char *)malloc(max_prot_len+1);
  memset(spaces,' ',max_prot_len);
  spaces[max_prot_len] = 0;
  char *dashes = (char *)malloc(max_prot_len+1);
  memset(dashes,'-',max_prot_len);
  dashes[max_prot_len] = 0;

  char buf[2048];
  
  std::string output("<HTML><BODY BGCOLOR=\"#FFFFFF\"><PRE>\n");

  snprintf(buf,sizeof(buf), "<font color=\"green\">%s.out</font>\n\n", spectra_[index]);
  output.append(buf);


  snprintf(buf,sizeof(buf), "%s, ", database_);
  output.append(buf);
  if(mass_type_ == 1)
    snprintf(buf,sizeof(buf), "MONO/MONO");
  else
    snprintf(buf,sizeof(buf), "AVG/AVG");
  output.append(buf);
  snprintf(buf,sizeof(buf), "\n\n");
  output.append(buf);

  // now list modifications
  if(static_mods_.length() > 0) {
    output.append("Static Modifications\n");
    for(k = 0; k < static_mods_.length(); k++) {
      snprintf(buf,sizeof(buf), "   %d: %c=%0f\n", k+1, static_mods_[k]->aa, static_mods_[k]->mass);
      output.append(buf);
    }
    output.append("\n");
  }
  if(variable_mods_.length() > 0) {
    output.append("Variable Modifications\n");
    for(k = 0; k < variable_mods_.length(); k++) {
      snprintf(buf,sizeof(buf), "   %d: %c=%0f\n", k+1, variable_mods_[k]->aa, variable_mods_[k]->mass);
      output.append(buf);
    }
    output.append("\n");
  }
  output.append("\n\n");
  snprintf(buf,sizeof(buf), "Identity: %0.2f\tHomology: %0.2f\n", identityscores_[index], homologyscores_[index]);
  output.append(buf);
  output.append("\n\n");

  output.append("<font color=\"green\">");
  output.append(" #    ");
  output.append(" MH+          IonSc   Ions    Ref");
#define WRITESPACES(n) {if ((n)>0) {spaces[n]=0;output.append(spaces);spaces[n]=' ';}}
  WRITESPACES(max_prot_len);
  if(calculate_pI_)
    output.append("  CalcPI ");


  output.append("   Sequence\n");

  output.append("---  ");
  output.append("-------------  -----  -------  ---");
#define WRITEDASHES(n) {if ((n)>0) {dashes[n]=0;output.append(dashes);dashes[n]='-';}}
  WRITEDASHES(max_prot_len - 3);
  output.append("    ");
  if(calculate_pI_)
    output.append("------   ");

  output.append(" ");

  output.append("--------------------------");

  output.append("</font>\n");

  output.append("1    ");
  if(masses_[index] < 1000)
    output.append(" ");
  snprintf(buf,sizeof(buf), "%0.1f (", masses_[index]);
  output.append(buf);
  if(massdiffs_[index] >= 0)
    output.append("+");
  snprintf(buf,sizeof(buf), "%0.1f)  ", massdiffs_[index]);
  output.append(buf);
  if(ionscores_[index] < 10)
    output.append(" ");
  if(ionscores_[index] >= identityscores_[index])
    snprintf(buf,sizeof(buf), "<font color=\"#DD00DD\">%0.2f</font>  ", ionscores_[index]);
  else
    snprintf(buf,sizeof(buf), "%0.2f  ", ionscores_[index]);
  output.append(buf);

  int totnumions = (int)(strlen(peptides_[index]) - 1) * num_ion_series_;
  if((spectra_[index])[strlen(spectra_[index])-1] > '2')
    totnumions *= 2;

  modtags = getModificationTags(peptides_[index], mods_[index]);

  if (modtags) {
    modinfo = new ModificationInfo(modtags);
  }

  snprintf(buf,sizeof(buf),"<A TARGET=\"Win1\" HREF=\"%splot-msms-js.cgi?Dta=%s%s/%s.dta&amp;MassType=%d&amp;NumAxis=1&amp;ISeries=010000011&amp;Pep=%s", getCgiUrl(), fullpath_, filextn_, spectra_[index], mass_type_, peptides_[index]);
 
  output.append(buf);
  if (modinfo != NULL) {
    for (int pepidx=1; pepidx <= peplen; pepidx++) {
      if (modinfo->getModifiedResidueMass(pepidx-1) > 0.0 || modinfo->getModifiedResidueMass(pepidx-1) < 0.0) {
	snprintf(buf,sizeof(buf), "&amp;Mod%d=%f",pepidx, modinfo->getModifiedResidueMass(pepidx-1));
        output.append(buf);
      }
    }
  }
  output.append("\">");

  if(matchedions_[index] < 100)
    output.append(" ");
  if(matchedions_[index] < 10)
    output.append(" ");
  snprintf(buf,sizeof(buf), "%d/", matchedions_[index]);
  output.append(buf);
  if(totnumions < 100)
    output.append(" ");
  if(totnumions < 10)
    output.append(" ");
  snprintf(buf,sizeof(buf), "%d", totnumions);
  output.append(buf);
  output.append("</A>  ");

  snprintf(buf,sizeof(buf), "<A TARGET=\"Win1\" HREF=\"%scomet-fastadb.cgi?Ref=%s&amp;Db=%s&amp;Pep=%s&amp;MassType=%d\">%s", getCgiUrl(), proteins_[index], database_, peptides_[index], mass_type_, proteins_[index]);
  output.append(buf);
  output.append("</A>");
  int plen = (int)strlen(proteins_[index]);
  WRITESPACES(max_prot_len - plen);

  if(pluses_[index]) {
    if(pluses_[index] < 100)
      output.append(" ");
    if(pluses_[index] < 10)
      output.append(" ");

    int min_ntt = getMinNumberTermini();
    snprintf(buf,sizeof(buf), "<A TARGET=\"Win1\" HREF=\"%scomet-fastadb.cgi?Db=%s&amp;Pep=%s&amp;MassType=%d&amp;sample_enzyme=%s&amp;min_ntt=%d\">+%d</A>", getCgiUrl(), database_, peptides_[index], mass_type_, sampleEnzymeInstance_->getName(), min_ntt, pluses_[index]);
    output.append(buf);
  }
  else
    output.append("    ");
  output.append(" ");

  if(modtags != NULL) {
    if(modinfo != NULL) {
      stdpep = modinfo->getStandardModifiedPeptide(peptides_[index], EXCLUDE_CONST_STATICS_?&static_consts_:NULL, MOD_ERROR, "<font size=\"-2\">", "</font>");
      // convert to html friendly format
      // in the future, use modinfo to output codes for spectrum view CGI

      if(calculate_pI_ && pi_calc_ != NULL)
	calc_pi = pi_calc_->Peptide_pI(peptides_[index], modinfo);
    } // if modinfo
  } // if modtags
  else if(calculate_pI_)
    calc_pi = pi_calc_->Peptide_pI(peptides_[index], NULL);

  if(calculate_pI_) {
    if(calc_pi < 10)
      output.append(" ");
    snprintf(buf,sizeof(buf), "%0.2f", calc_pi);
    output.append(buf);
    output.append("    ");
  }

  if(stdpep != NULL) {
    if(isSemiCleavageSearch())
      snprintf(buf,sizeof(buf), "<A TARGET=\"Win1\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&amp;LAYOUT=TwoWindows&amp;AUTO_FORMAT=Semiauto&amp;ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;CDD_SEARCH=on&amp;CLIENT=web&amp;COMPOSITION_BASED_STATISTICS=on&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=(none)&amp;EXPECT=1000&amp;FILTER=L&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;I_THRESH=0.005&amp;MATRIX_NAME=BLOSUM62&amp;NCBI_GI=on&amp;PAGE=Proteins&amp;PROGRAM=blastp&amp;SERVICE=plain&amp;SET_DEFAULTS.x=41&amp;SET_DEFAULTS.y=5&amp;SHOW_OVERVIEW=on&amp;END_OF_HTTPGET=Yes&amp;SHOW_LINKOUT=yes&amp;QUERY=%s\">%s</A>", peptides_[index], stdpep);
    else
      snprintf(buf,sizeof(buf), "%c.<A TARGET=\"Win1\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&amp;LAYOUT=TwoWindows&amp;AUTO_FORMAT=Semiauto&amp;ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;CDD_SEARCH=on&amp;CLIENT=web&amp;COMPOSITION_BASED_STATISTICS=on&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=(none)&amp;EXPECT=1000&amp;FILTER=L&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;I_THRESH=0.005&amp;MATRIX_NAME=BLOSUM62&amp;NCBI_GI=on&amp;PAGE=Proteins&amp;PROGRAM=blastp&amp;SERVICE=plain&amp;SET_DEFAULTS.x=41&amp;SET_DEFAULTS.y=5&amp;SHOW_OVERVIEW=on&amp;END_OF_HTTPGET=Yes&amp;SHOW_LINKOUT=yes&amp;QUERY=%s\">%s</A>.%c", prev_AAs_[index], peptides_[index], stdpep, foll_AAs_[index]);
    output.append(buf);
    delete[] stdpep;
  }
  else {

    if(isSemiCleavageSearch())
      snprintf(buf,sizeof(buf), "<A TARGET=\"Win1\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&amp;LAYOUT=TwoWindows&amp;AUTO_FORMAT=Semiauto&amp;ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;CDD_SEARCH=on&amp;CLIENT=web&amp;COMPOSITION_BASED_STATISTICS=on&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=(none)&amp;EXPECT=1000&amp;FILTER=L&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;I_THRESH=0.005&amp;MATRIX_NAME=BLOSUM62&amp;NCBI_GI=on&amp;PAGE=Proteins&amp;PROGRAM=blastp&amp;SERVICE=plain&amp;SET_DEFAULTS.x=41&amp;SET_DEFAULTS.y=5&amp;SHOW_OVERVIEW=on&amp;END_OF_HTTPGET=Yes&amp;SHOW_LINKOUT=yes&amp;QUERY=%s\">%s</A>", peptides_[index], peptides_[index]);
    else
      snprintf(buf,sizeof(buf), "%c.<A TARGET=\"Win1\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&amp;LAYOUT=TwoWindows&amp;AUTO_FORMAT=Semiauto&amp;ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;CDD_SEARCH=on&amp;CLIENT=web&amp;COMPOSITION_BASED_STATISTICS=on&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=(none)&amp;EXPECT=1000&amp;FILTER=L&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;I_THRESH=0.005&amp;MATRIX_NAME=BLOSUM62&amp;NCBI_GI=on&amp;PAGE=Proteins&amp;PROGRAM=blastp&amp;SERVICE=plain&amp;SET_DEFAULTS.x=41&amp;SET_DEFAULTS.y=5&amp;SHOW_OVERVIEW=on&amp;END_OF_HTTPGET=Yes&amp;SHOW_LINKOUT=yes&amp;QUERY=%s\">%s</A>.%c", prev_AAs_[index], peptides_[index], peptides_[index], foll_AAs_[index]);
    output.append(buf);
  }

  output.append("\n");
  
  if (modinfo != NULL) {
    delete modinfo;
    modinfo = NULL;
  }
  if (modtags != NULL) {
    ARRAYDEL(*modtags);
    delete modtags;
    modtags = NULL;
  }

  for(k = alt_inds_[index]; k < stop; k++) {
    peplen = (int)strlen(alt_peptides_[k]);
    modtags = getModificationTags(alt_peptides_[k], alt_mods_[k]);
    if (modtags != NULL) 
      modinfo = new ModificationInfo(modtags);
    
    snprintf(buf,sizeof(buf), "%d ", k - alt_inds_[index] + 2);
    output.append(buf);
    if(k - alt_inds_[index] + 2 < 1000)
      output.append(" ");
    if(k - alt_inds_[index] + 2 < 100)
      output.append(" ");
    if(k - alt_inds_[index] + 2 < 10)
      output.append(" ");
    if(alt_masses_[k] < 1000)
      output.append(" ");
    snprintf(buf,sizeof(buf), "%0.1f (", alt_masses_[k]);
    output.append(buf);
    if(alt_massdiffs_[k] >= 0)
      output.append("+");
    snprintf(buf,sizeof(buf), "%0.1f)  ", alt_massdiffs_[k]);
    output.append(buf);
    if(alt_ionscores_[k] < 10)
      output.append(" ");
    if(alt_ionscores_[k] >= identityscores_[index])
      snprintf(buf,sizeof(buf), "<font color=\"#DD00DD\">%0.2f</font>  ", alt_ionscores_[k]);
    else
      snprintf(buf,sizeof(buf), "%0.2f  ", alt_ionscores_[k]);
    output.append(buf);

 
    totnumions = (peplen - 1) * num_ion_series_;
    if((spectra_[index])[strlen(spectra_[index])-1] > '2')
      totnumions *= 2;

    snprintf(buf,sizeof(buf),"<A TARGET=\"Win1\" HREF=\"%splot-msms-js.cgi?Dta=%s%s/%s.dta&amp;MassType=%d&amp;NumAxis=1&amp;ISeries=010000011&amp;Pep=%s", getCgiUrl(), fullpath_, filextn_, spectra_[index], mass_type_, alt_peptides_[k]);
    output.append(buf);
  
    if (modinfo != NULL) {
      for (int pepidx=1; pepidx <= peplen; pepidx++) {
	if (modinfo->getModifiedResidueMass(pepidx-1) > 0.0 || modinfo->getModifiedResidueMass(pepidx-1) < 0.0) {
	  snprintf(buf,sizeof(buf),"&amp;Mod%d=%f",pepidx, modinfo->getModifiedResidueMass(pepidx-1));
	  output.append(buf);
	}
      }
    }
    output.append("\">");
    
    if(alt_matchedions_[k] < 100)
      output.append(" ");
    if(alt_matchedions_[k] < 10)
      output.append(" ");
    snprintf(buf,sizeof(buf), "%d/", alt_matchedions_[k]);
    output.append(buf);
    if(totnumions < 100)
      output.append(" ");
    if(totnumions < 10)
      output.append(" ");
    snprintf(buf,sizeof(buf), "%d", totnumions);
    output.append(buf);
    output.append("</A>  ");

    snprintf(buf,sizeof(buf), "<A TARGET=\"Win1\" HREF=\"%scomet-fastadb.cgi?Ref=%s&amp;Db=%s&amp;Pep=%s&amp;MassType=%d\">%s", getCgiUrl(), alt_proteins_[k], database_, alt_peptides_[k], mass_type_, alt_proteins_[k]);
    output.append(buf);
    output.append("</A>");
    int alen = (int)strlen(alt_proteins_[k]);
    WRITESPACES(max_prot_len - alen);

    if(alt_pluses_[k]) {
      if(alt_pluses_[k] < 100)
	output.append(" ");
      if(alt_pluses_[k] < 10)
	output.append(" ");
      
      int min_ntt = getMinNumberTermini();
      snprintf(buf,sizeof(buf), "<A TARGET=\"Win1\" HREF=\"%scomet-fastadb.cgi?Db=%s&amp;Pep=%s&amp;MassType=%d&amp;sample_enzyme=%s&amp;min_ntt=%d\">+%d</A>", getCgiUrl(), database_, alt_peptides_[k], mass_type_, sampleEnzymeInstance_->getName(), min_ntt, alt_pluses_[k]);
      output.append(buf);
    }
    else
      output.append("    ");
    output.append(" ");

    if(calculate_pI_) {
      if(calc_pi < 10)
	output.append(" ");
      snprintf(buf,sizeof(buf), "%0.2f", calc_pi);
      output.append(buf);
      output.append("    ");
    }

    stdpep = NULL;
    
    if(modtags != NULL) {
      
      if(modinfo != NULL) {
	stdpep = modinfo->getStandardModifiedPeptide(alt_peptides_[k], EXCLUDE_CONST_STATICS_?&static_consts_:NULL, MOD_ERROR, "<font size=\"-2\">", "</font>");
	//stdpep = modinfo->getStandardModifiedPeptide(alt_peptides_[k], static_consts_, MOD_ERROR);
  
        if(calculate_pI_ && pi_calc_ != NULL)
          calc_pi = pi_calc_->Peptide_pI(alt_peptides_[k], modinfo);
      } // if modinfo
    } // if modtags

    if(stdpep != NULL) {
      if(isSemiCleavageSearch())
        snprintf(buf,sizeof(buf), "<A TARGET=\"Win1\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&amp;LAYOUT=TwoWindows&amp;AUTO_FORMAT=Semiauto&amp;ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;CDD_SEARCH=on&amp;CLIENT=web&amp;COMPOSITION_BASED_STATISTICS=on&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=(none)&amp;EXPECT=1000&amp;FILTER=L&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;I_THRESH=0.005&amp;MATRIX_NAME=BLOSUM62&amp;NCBI_GI=on&amp;PAGE=Proteins&amp;PROGRAM=blastp&amp;SERVICE=plain&amp;SET_DEFAULTS.x=41&amp;SET_DEFAULTS.y=5&amp;SHOW_OVERVIEW=on&amp;END_OF_HTTPGET=Yes&amp;SHOW_LINKOUT=yes&amp;QUERY=%s\">%s</A>", alt_peptides_[k], stdpep);
      else
	snprintf(buf,sizeof(buf), "%c.<A TARGET=\"Win1\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&amp;LAYOUT=TwoWindows&amp;AUTO_FORMAT=Semiauto&amp;ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;CDD_SEARCH=on&amp;CLIENT=web&amp;COMPOSITION_BASED_STATISTICS=on&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=(none)&amp;EXPECT=1000&amp;FILTER=L&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;I_THRESH=0.005&amp;MATRIX_NAME=BLOSUM62&amp;NCBI_GI=on&amp;PAGE=Proteins&amp;PROGRAM=blastp&amp;SERVICE=plain&amp;SET_DEFAULTS.x=41&amp;SET_DEFAULTS.y=5&amp;SHOW_OVERVIEW=on&amp;END_OF_HTTPGET=Yes&amp;SHOW_LINKOUT=yes&amp;QUERY=%s\">%s</A>.%c", alt_prev_AAs_[k], alt_peptides_[k], stdpep, alt_foll_AAs_[k]);
      output.append(buf);
      delete[] stdpep;
    }
    else {
      if(isSemiCleavageSearch())
        snprintf(buf,sizeof(buf),  "<A TARGET=\"Win1\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&amp;LAYOUT=TwoWindows&amp;AUTO_FORMAT=Semiauto&amp;ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;CDD_SEARCH=on&amp;CLIENT=web&amp;COMPOSITION_BASED_STATISTICS=on&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=(none)&amp;EXPECT=1000&amp;FILTER=L&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;I_THRESH=0.005&amp;MATRIX_NAME=BLOSUM62&amp;NCBI_GI=on&amp;PAGE=Proteins&amp;PROGRAM=blastp&amp;SERVICE=plain&amp;SET_DEFAULTS.x=41&amp;SET_DEFAULTS.y=5&amp;SHOW_OVERVIEW=on&amp;END_OF_HTTPGET=Yes&amp;SHOW_LINKOUT=yes&amp;QUERY=%s\">%s</A>", alt_peptides_[k], alt_peptides_[k]);
      else
	snprintf(buf,sizeof(buf), "%c.<A TARGET=\"Win1\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&amp;LAYOUT=TwoWindows&amp;AUTO_FORMAT=Semiauto&amp;ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;CDD_SEARCH=on&amp;CLIENT=web&amp;COMPOSITION_BASED_STATISTICS=on&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=(none)&amp;EXPECT=1000&amp;FILTER=L&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;I_THRESH=0.005&amp;MATRIX_NAME=BLOSUM62&amp;NCBI_GI=on&amp;PAGE=Proteins&amp;PROGRAM=blastp&amp;SERVICE=plain&amp;SET_DEFAULTS.x=41&amp;SET_DEFAULTS.y=5&amp;SHOW_OVERVIEW=on&amp;END_OF_HTTPGET=Yes&amp;SHOW_LINKOUT=yes&amp;QUERY=%s\">%s</A>.%c", alt_prev_AAs_[k], alt_peptides_[k], alt_peptides_[k], alt_foll_AAs_[k]);
      output.append(buf);
    } // no stdpep

    output.append("\n");
    if (modinfo != NULL) { 
      delete modinfo;
      modinfo = NULL;
    }
    if (modtags != NULL) {
      ARRAYDEL(*modtags);
      delete modtags;
      modtags = NULL;
    }
  } // next

  tarball_->create_memberfile(tarfilename_,outfilename_,output);
 
#define OUT_TAR_CHECK_FREQ 500
#define OUT_TAR_CHECK_FIRST 50
  if (((RUN_TEST==testmode_)||(FORCE_TEST==testmode_)) && 
      ((index < OUT_TAR_CHECK_FIRST) || !(index%OUT_TAR_CHECK_FREQ))) // doing each one is too slow...
    { // compare with tar file saved earlier for regression test
      result |= compare_tarred_file(output,outfilename_);
    }
  free(spaces);
  free(dashes);
  return result;
}

// go into database and retrieve full protein name
// and in the case of unconstrained searches, the preceding and following aa
void MascotConverter::databaseSearch() {

  if (generate_description_) {
    // Let's build up the missing protein description map's lookup
    // we will fill them in against the sequence database
    int k;
    for (k=0; k<proteins_.length(); k++) {
      std::string id(proteins_[k]);
      if (!id.empty()) {
        if (proteinDescriptionMap_.find (id) == proteinDescriptionMap_.end()) {
          if (missingDescriptionMap_.find(id) == missingDescriptionMap_.end()) {
            std::string description;
            missingDescriptionMap_.insert(std::make_pair(id,description));
          }
        }
      }
    }
    for (k=0; k<rank1_proteins_.length(); k++) {
      std::string id(rank1_proteins_[k]);
      if (!id.empty()) {
        if (proteinDescriptionMap_.find (id) == proteinDescriptionMap_.end()) {
          if (missingDescriptionMap_.find(id) == missingDescriptionMap_.end()) {
            std::string description;
            missingDescriptionMap_.insert(std::make_pair(id,description));
          }
        }
      }
    }
    for (k=0; k<alt_proteins_.length(); k++) {
      std::string id(alt_proteins_[k]);
      if (!id.empty()) {
        if (proteinDescriptionMap_.find (id) == proteinDescriptionMap_.end()) {
          if (missingDescriptionMap_.find(id) == missingDescriptionMap_.end()) {
            std::string description;
            missingDescriptionMap_.insert(std::make_pair(id,description));
          }
        }
      }
    }
  }

  if (term_read_) {
    //if terminus have been read, let's check if they tally in term of total count
    if (!(prev_AAs_.length() == spectra_.length() && rank1_prev_AAs_.length() == rank1_proteins_.length() && alt_prev_AAs_.length() == alt_proteins_.length())) {
      term_read_ = False; // trigger database lookup later
    }
  }

  if (term_read_ && generate_description_ && missingDescriptionMap_.empty()) {
    // make sure all the prev_aa and next_aa terms were read correctly.
    if (prev_AAs_.length() == spectra_.length() && rank1_prev_AAs_.length() == rank1_proteins_.length() && alt_prev_AAs_.length() == alt_proteins_.length()) {
      if (short_protein_id_) {
        cerr << "prev_aa and next_aa and protein descriptions have already been read from .dat file, database search skipped." << endl;
        return;
      }
    }
    else {
      // trigger database lookup later
      term_read_ = False;
    }
  }
  else if (term_read_ && !generate_description_) {
    if (prev_AAs_.length() == spectra_.length() && rank1_prev_AAs_.length() == rank1_proteins_.length() && alt_prev_AAs_.length() == alt_proteins_.length()) {
      if (short_protein_id_) {
        cerr << "prev_aa and next_aa have already been read from .dat file, database search skipped." << endl;
        return;
      }
    }
    else {
      // trigger database lookup later
      term_read_ = False;
    }
  }

  if (!generate_description_ || !short_protein_id_) {
    // if we are not using short id, we have to re-run thru'
    // all the proteinDescription(Map_) that we have gathered
    // so that the full protein id can be retrieved from the database

    // we use this as lookup to speed up checking 
    // whether we need to fill in the terminus
    int k;
    for (k=0; k<proteins_.length(); k++) {
      std::string id(proteins_[k]);
      if (!id.empty()) {
        proteinMap::iterator iter = proteinDescriptionMap_.find (id);
        std::string description;
        if (iter!=proteinDescriptionMap_.end()) {
          description = iter->second;
          proteinDescriptionMap_.erase(iter);
        }
        missingDescriptionMap_.insert(std::make_pair(id,description));
      }
    }
    for (k=0; k<rank1_proteins_.length(); k++) {
      std::string id(rank1_proteins_[k]);
      if (!id.empty()) {
        proteinMap::iterator iter = proteinDescriptionMap_.find (id);
        std::string description;
        if (iter!=proteinDescriptionMap_.end()) {
          description = iter->second;
          proteinDescriptionMap_.erase(iter);
        }
        missingDescriptionMap_.insert(std::make_pair(id,description));
      }
    }
    for (k=0; k<alt_proteins_.length(); k++) {
      std::string id(alt_proteins_[k]);
      if (!id.empty()) {
        proteinMap::iterator iter = proteinDescriptionMap_.find (id);
        std::string description;
        if (iter!=proteinDescriptionMap_.end()) {
          description = iter->second;
          proteinDescriptionMap_.erase(iter);
        }
        missingDescriptionMap_.insert(std::make_pair(id,description));
      }
    }
  }

  if (!term_read_) {
    cerr << "Completing " << (spectra_.length()+rank1_proteins_.length()+alt_proteins_.length()) << " prev_aa and next_aa attributes by parsing through the database..." << endl;
  }
  else {
    cerr << "prev_aa and next_aa have already been read from .dat file." << endl;
  }
  if (generate_description_) {
    if (!missingDescriptionMap_.empty()) {
      cerr << "Completing " << missingDescriptionMap_.size() << " protein descriptions by parsing through the database..." << endl;
    }
    else {
      cerr << "Protein descriptions have already been read from .dat file." << endl;
    }
  }
  if (!short_protein_id_) {
    if (!missingDescriptionMap_.empty()) {
      cerr << "Replacing " << missingDescriptionMap_.size() << " protein id with full ids by parsing through the database..." << endl;
    }
  }

  int k;
  Array<int> inds;
  inds.reserve(spectra_.length());
  Array<int> rank1_inds;
  rank1_inds.reserve(rank1_proteins_.length());
  Array<int> alt_inds;
  alt_inds.reserve(alt_proteins_.length());

  if (!term_read_) {
    prev_AAs_.reserve(spectra_.length());
    foll_AAs_.reserve(spectra_.length());
    for(k = 0; k < spectra_.length(); k++) {
      prev_AAs_.insertAtEnd(G_UNKNOWN_AA);
      foll_AAs_.insertAtEnd(G_UNKNOWN_AA);
      if(checkResult(k))
        inds.insertAtEnd(k);
    }

    rank1_prev_AAs_.reserve(rank1_proteins_.length());
    rank1_foll_AAs_.reserve(rank1_proteins_.length());
    for(k = 0; k < rank1_proteins_.length(); k++) {
      rank1_prev_AAs_.insertAtEnd(G_UNKNOWN_AA);
      rank1_foll_AAs_.insertAtEnd(G_UNKNOWN_AA);
      rank1_inds.insertAtEnd(k);
    }

    alt_prev_AAs_.reserve(alt_proteins_.length());
    alt_foll_AAs_.reserve(alt_proteins_.length());
    for(k = 0; k < alt_proteins_.length(); k++) {
      alt_prev_AAs_.insertAtEnd(G_UNKNOWN_AA);
      alt_foll_AAs_.insertAtEnd(G_UNKNOWN_AA);
      alt_inds.insertAtEnd(k);
    }
  }
  else {
    for(k = 0; k < spectra_.length(); k++) {
      if(checkResult(k))
        inds.insertAtEnd(k);
    }

    for(k = 0; k < rank1_proteins_.length(); k++) {
      rank1_inds.insertAtEnd(k);
    }

    for(k = 0; k < alt_proteins_.length(); k++) {
      alt_inds.insertAtEnd(k);
    }
  }

  RACI fdata(database_); // can read gzipped results

  if(! fdata) {
    cerr << "cannot open database " << database_ << endl;
    exit(1);
  }
  // get length of file:
  cerr << "searching " << database_;
  fdata.seekg (0, ios::end);
  double dbase_length = (double)fdata.tellg();
  fdata.seekg (0, ios::beg);  
  int last_pct = 0;
  int last_pct_chunk = 0;

  int prot_size = 50000;
  char* seq =  new char[prot_size];
  int prot_name_size = 20000; //200;
  char* prot = new char[prot_name_size];
  int prot_desc_size = 20000; //200;
  char* prot_desc = new char[prot_desc_size];
  int line_size = prot_name_size * 2;
  char* nextline = new char[line_size];
  seq[0] = 0;
  prot[0] = 0;
  prot_desc[0] = 0;
  Boolean first = True;

  //wch: 2006-09-23
  int proteinIDTypeGuessCount = 0;
  int proteinIDGuessedType = UNKNOWN;
  int proteinIDType = UNKNOWN;
  //END-wch: 2006-09-23

  Boolean fCheckAllIndices = (!term_read_ || !short_protein_id_);

  while(!fdata.eof() && 
	(
	 (fCheckAllIndices && 
	  (inds.length() > 0 || rank1_inds.length() > 0 || alt_inds.length() > 0))
	 || !missingDescriptionMap_.empty())
	) {
    fdata.getline(nextline, line_size);
    if (fdata.fail()) {
      if ((int)strlen(nextline) == (line_size - 1)) {
        // line was too large for the buffer.
        
        // if it was the ID line (most likely case), only the part before
        // the first space or tab is required, so just read to the end of
        // the line discarding the rest.
        if (nextline[0] == '>') {
          char trash[4096];
          while (fdata.fail()) {
            fdata.clear();
            fdata.getline(trash, sizeof(trash) / sizeof(char));
          }
        } else {
          cerr << "Sequence line exceeds " << line_size << "characters" << endl;
          exit(1);
        }
      } else if (fdata.eof()) {
	break;
      } else {
        cerr << "Fail to get line from sequence database" << endl;
        exit(1);
      }
    }

#if 1
    // WCH: there seems to be a special case in Windows
    //      where the first aa in a sequence is missed
    //      when ifstream is in text mode due to the conversion
    //      of CRLF to just LF
    int nextlinelen = (int)strlen(nextline)-1;
    while (nextlinelen >= 0) {
      if ('\r' == nextline[nextlinelen] ||
	  '\n' == nextline[nextlinelen]) {
        nextline[nextlinelen] = '\0';
        nextlinelen--;
      } else {
        break;
      } 
    }
#endif

    if(strlen(nextline) > 0 && nextline[0] == '>') { // protein
      if(strlen(seq) > 0 && strlen(prot) > 0) { // process

        // let's process the read entry before this new one
        // but we have to trim to get the id of the protein

        // start processing now
        double dbase_pos = (double)fdata.tellg();
        int pct = (int)(100.0*(dbase_pos/dbase_length));
        int pctInterval=5;
        if (pct > (last_pct_chunk+pctInterval)) { 
          // show percent complete once in a while
          last_pct_chunk = pct;
          cerr << pctInterval*(pct/pctInterval) << "%";
        } else if (pct > last_pct) {
          cerr << ".";
        }
        last_pct = pct;

        buildProteinPeptideInfo (prot, prot_desc, seq, proteinIDType, proteinIDGuessedType, proteinIDTypeGuessCount, inds, rank1_inds, alt_inds);

      } // process

      // reset for this new incoming entry
      seq[0] = 0;
      prot[0] = 0; // reset
      prot_desc[0] = 0; // reset
      int last = 0;
      first = True;
      int nextline_len = (int)strlen(nextline);
      while(last < nextline_len && nextline[last] != ' ' && nextline[last] != '\t')
        last++;

      if(last-1 > prot_name_size) {
        cerr << "prot name for " << nextline << " exceeds " << prot_name_size << endl;
        exit(1);
      }
      //fasta uses Ctrl+A character to concatenate comment lines
      //we exclude anything after that to get a proper comment line
      int desclast = last;
      /*int*/ nextline_len = (int)strlen(nextline);
      while(desclast < nextline_len  && nextline[desclast] != '\x01')
        desclast++;
      if ((last+1) < nextline_len) {
        // check that we do not copy beyond our buffer
        int charToCopy = desclast-last-1;
        if (charToCopy > prot_desc_size) {
          charToCopy = (prot_desc_size-1);
        }

        strncat(prot_desc, nextline+last+1, charToCopy);
        prot_desc[charToCopy] = 0;
      }
      nextline[last] = 0; // set it to end
      strncat(prot, nextline+1, last-1);
      prot[last-1] = 0;

    }
    else {
      //we do not want a case where line_size failed,
      //but a string of residue of much smaller size can get pass
      //on future check before the next protein entry
      if(first && (int) (strlen(seq) + strlen(nextline)) <= prot_size) {
        strcat(seq, nextline);
      }
      else {
        if(first) {
          cerr << "warning, size of protein " << prot << " exceeds " << prot_size << endl;
          first = False;
        }
      }
    }
  } // next line of db

  if(strlen(seq) > 0 && strlen(prot) > 0) { // process last db entry
    buildProteinPeptideInfo (prot, prot_desc, seq, proteinIDType, proteinIDGuessedType, proteinIDTypeGuessCount, inds, rank1_inds, alt_inds);
  } // process

  cerr << "done" << endl;

  if (generate_description_ && !missingDescriptionMap_.empty()) {
    cerr << "Protein ID type: " << proteinIDType << " (guessed)" << endl;
    cerr << "WARNING: Failed to fill remaining " << missingDescriptionMap_.size() << " protein descriptions" << endl;
    cerr << "         Your database might not be the exact version used or the IDs are not before the first space in the description" << endl;
    int k=1;
    for (proteinMap::const_iterator it = missingDescriptionMap_.begin();
	 it != missingDescriptionMap_.end(); it++, k++) {
      cerr << "WARNING: #" << k << ", ->" << it->first.c_str() << "<-" << endl;
    }
  }
  if (!term_read_) {
    if (inds.length()>0) {
      cerr << "WARNING: Failed to locate terminus for the remaining " << inds.length() << " peptides " << endl;
      cerr << "         Your database might not be the exact version used or the sequence has changed" << endl;
      int k;
      for(k = 0; k < inds.length(); k++) {
        cerr << "WARNING: #" << (k+1) << ", " << inds[k] << " | " << proteins_[inds[k]] << " | " << peptides_[inds[k]] << endl;
      }
    }
    if (rank1_inds.length()>0) {
      cerr << "WARNING: Failed to locate terminus for the remaining " << rank1_inds.length() << " peptides " << endl;
      cerr << "         Your database might not be the exact version used or the sequence has changed" << endl;
      int k;
      for(k = 0; k < rank1_inds.length(); k++) {
        cerr << "WARNING: #" << (k+1) << ", " << rank1_inds[k] << " | " << rank1_proteins_[rank1_inds[k]] << endl;
      }
    }
    if (alt_inds.length()>0) {
      cerr << "WARNING: Failed to locate terminus for the remaining " << alt_inds.length() << " peptides " << endl;
      cerr << "         Your database might not be the exact version used or the sequence has changed" << endl;
      int k;
      for(k = 0; k < alt_inds.length(); k++) {
        cerr << "WARNING: #" << (k+1) << ", " << alt_inds[k] << " | " << alt_proteins_[alt_inds[k]] << " | " << alt_peptides_[alt_inds[k]] << endl;
      }
    }
  }
  
  fdata.close();
  delete [] nextline;
  delete [] seq;
  delete [] prot;
  delete [] prot_desc;
}

void MascotConverter::writeHTML(const char* outfile) {

  //char prefix[] = "/users/akeller/MASCOT/";
  char suffix[] = ".out";
  char outfilecgi[] = "/mascotout.pl";

  FILE* fout;

  if( (fout = fopen(outfile, "w")) == NULL) {
    cerr << "cannot open outfile " << outfile << endl;
    exit(1);
  }


  int max_prot_len = 0;
  int max_spec_len = 0;
  int j,k;
  for(k = 0; k < spectra_.length(); k++) {
    if((int) strlen(proteins_[k]) > max_prot_len)
      max_prot_len = (int)strlen(proteins_[k]);
    if((int) strlen(spectra_[k]) > max_spec_len)
      max_spec_len = (int)strlen(spectra_[k]);
  }


  time_t now;
  time(&now);
  struct tm* tmstruct = localtime(&now);

  fprintf(fout, "<HTML><BODY BGCOLOR=\"#FFFFFF\"><PRE>\n");
  fprintf(fout, "\n<font color=\"green\">MASCOT-SUMMARY v.1   Copyright 2003</font>\n");
  fprintf(fout, "Institute for Systems Biology, A.Keller\n");
  
  char* mins = new char[3];
  if(tmstruct->tm_min < 10)
    sprintf(mins, "0%d", tmstruct->tm_min);
  else
    sprintf(mins, "%d", tmstruct->tm_min);

  if(tmstruct->tm_hour > 12) 
    fprintf(fout, "%d/%d/%d, %d:%s PM, ", tmstruct->tm_mon + 1, tmstruct->tm_mday, tmstruct->tm_year + 1900, tmstruct->tm_hour - 12, mins);
  else
    fprintf(fout, "%d/%d/%d, %d:%s AM, ", tmstruct->tm_mon + 1, tmstruct->tm_mday, tmstruct->tm_year + 1900, tmstruct->tm_hour, mins);
  fprintf(fout, "%s, ", database_);

  delete [] mins;
  if(mass_type_ == 1)
    fprintf(fout, "MONO/MONO");
  else
    fprintf(fout, "AVG/AVG");
  fprintf(fout, "\n\n\n");

  fprintf(fout, "<font color=\"green\">");
  fprintf(fout, " #    File");
  for(j = 0; j < max_spec_len - 2; j++)
    fprintf(fout, " ");
  fprintf(fout, " MH+          IonSc    Id    Homol   Ions    Ref");
  for(j = 0; j < max_prot_len - 1; j++)
    fprintf(fout, " ");
  fprintf(fout, "    Sequence\n");
  fprintf(fout, "---  ----");
  for(j = 0; j < max_spec_len - 4; j++)
    fprintf(fout, "-");
  fprintf(fout, "   -----         -----   -----  -----  -------  ---");
  for(j = 0; j < max_prot_len - 3; j++)
    fprintf(fout, "-");
  fprintf(fout, "      --------------------");
  fprintf(fout, "</font>\n");

  char slash = '/';
  int counter = 1;

  for(k = 0; k < spectra_.length(); k++) {

    if(checkResult(k)) {
      fprintf(fout, "%d ", counter);
      if(counter < 1000)
        fprintf(fout, " ");
      if(counter < 100)
        fprintf(fout, " ");
      if(counter < 10)
        fprintf(fout, " ");
      counter++;
      fprintf(fout, "<A TARGET=\"Win1\" HREF=\"%s%s?OutFile=", getCgiUrl(), outfilecgi);
      fprintf(fout, "%s%s%c%s%s", fullpath_, filextn_, slash, spectra_[k], suffix);
      fprintf(fout, "\">%s</A>", spectra_[k]);
      int stopIndex=max_spec_len - (int) strlen(spectra_[k]);
      for(j = 0; j < stopIndex; j++)
        fprintf(fout, " ");

      fprintf(fout, "  "); // space
      if(masses_[k] < 1000)
        fprintf(fout, " ");
      fprintf(fout, "%0.1f (", masses_[k]);
      if(massdiffs_[k] >= 0)
        fprintf(fout, "+");
      fprintf(fout, "%0.1f)  ", massdiffs_[k]);
      if(ionscores_[k] < 10)
        fprintf(fout, " ");
      if(ionscores_[k] >= identityscores_[k])
        fprintf(fout, "<font color=\"#DD00DD\">%0.2f</font>", ionscores_[k]);
      else
        fprintf(fout, "%0.2f", ionscores_[k]);
      if(isStar(k))
        fprintf(fout, "*");
      else
        fprintf(fout, " ");
      fprintf(fout, "  ");

      if(identityscores_[k] < 10)
        fprintf(fout, " ");
      fprintf(fout, "%0.2f  ", identityscores_[k]);
      if(homologyscores_[k] < 10)
        fprintf(fout, " ");
      fprintf(fout, "%0.2f  ", homologyscores_[k]);

      int totnumions = ((int)strlen(peptides_[k]) - 1) * num_ion_series_;
      if((spectra_[k])[strlen(spectra_[k])-1] > '2')
        totnumions *= 2;

      fprintf(fout, "<A TARGET=\"Win1\" HREF=\"%ssequest-tgz-out.cgi?Dta=%s%s%c%s.dta&amp;MassType=%d&amp;NumAxis=1&amp;ISeries=010000011&amp;Pep=%s\">", getCgiUrl(), fullpath_, filextn_, slash, spectra_[k], mass_type_, peptides_[k]);
      if(matchedions_[k] < 100)
        fprintf(fout, " ");
      if(matchedions_[k] < 10)
        fprintf(fout, " ");
      fprintf(fout, "%d/", matchedions_[k]);
      if(totnumions < 100)
        fprintf(fout, " ");
      if(totnumions < 10)
        fprintf(fout, " ");
      fprintf(fout, "%d", totnumions);
      fprintf(fout, "</A>  ");

      fprintf(fout, "<A TARGET=\"Win1\" HREF=\"%scomet-fastadb.cgi?Ref=%s&amp;Db=%s&amp;Pep=%s&amp;MassType=%d\">%s", getCgiUrl(), proteins_[k], database_, peptides_[k], mass_type_, proteins_[k]);
      fprintf(fout, "</A>");
      /*int*/ stopIndex=max_prot_len - (int) strlen(proteins_[k]);
      for(j = 0; j < stopIndex; j++)
        fprintf(fout, " ");

      if(pluses_[k]) {
        if(pluses_[k] < 100)
          fprintf(fout, " ");
        if(pluses_[k] < 10)
          fprintf(fout, " ");

        fprintf(fout, "<A TARGET=\"Win1\" HREF=\"%scomet-fastadb.cgi?Db=%s&amp;Pep=%s&amp;MassType=%d\">+%d</A>", getCgiUrl(), database_, peptides_[k], mass_type_, pluses_[k]);
      }
      else
        fprintf(fout, "    ");
      fprintf(fout, "  ");
      if(isSemiCleavageSearch())
        fprintf(fout, "<A TARGET=\"Win1\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&amp;LAYOUT=TwoWindows&amp;AUTO_FORMAT=Semiauto&amp;ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;CDD_SEARCH=on&amp;CLIENT=web&amp;COMPOSITION_BASED_STATISTICS=on&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=(none)&amp;EXPECT=1000&amp;FILTER=L&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;I_THRESH=0.005&amp;MATRIX_NAME=BLOSUM62&amp;NCBI_GI=on&amp;PAGE=Proteins&amp;PROGRAM=blastp&amp;SERVICE=plain&amp;SET_DEFAULTS.x=41&amp;SET_DEFAULTS.y=5&amp;SHOW_OVERVIEW=on&amp;END_OF_HTTPGET=Yes&amp;SHOW_LINKOUT=yes&amp;QUERY=%s\">%s</A>\n", peptides_[k], peptides_[k]);
      else
        fprintf(fout, "%c.<A TARGET=\"Win1\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&amp;LAYOUT=TwoWindows&amp;AUTO_FORMAT=Semiauto&amp;ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;CDD_SEARCH=on&amp;CLIENT=web&amp;COMPOSITION_BASED_STATISTICS=on&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=(none)&amp;EXPECT=1000&amp;FILTER=L&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;I_THRESH=0.005&amp;MATRIX_NAME=BLOSUM62&amp;NCBI_GI=on&amp;PAGE=Proteins&amp;PROGRAM=blastp&amp;SERVICE=plain&amp;SET_DEFAULTS.x=41&amp;SET_DEFAULTS.y=5&amp;SHOW_OVERVIEW=on&amp;END_OF_HTTPGET=Yes&amp;SHOW_LINKOUT=yes&amp;QUERY=%s\">%s</A>.%c\n", prev_AAs_[k], peptides_[k], peptides_[k], foll_AAs_[k]);
    } // if ok result
  } // next

  fclose(fout);

}

const char *MascotConverter::getSearchParameter(const char *name) const {
  for(int k = 0; k < params_.length(); k++) {
    if (0==strcmp(name, params_[k])) {
      return vals_[k];
    }
  }
  return NULL;
}

int MascotConverter::getMaxInternalCleavages() const {
  const char *value = getSearchParameter("PFA");
  if (NULL==value) value = getSearchParameter("pfa");

  if (NULL==value) {
    // SHOULD NOT happen
    return -1;
  } else {
    return atoi(value);
  }
}

int MascotConverter::getMinNumberTermini() {
  if (NULL==enzyme_) {
    // could not get from .dat file
    // should not happen, and we are not assuming anything
    //return -1;
    // but want to keep to the old behaviour, so return 0 instead of -1
    return 0;
  } else {

    // there is no search constraint, equivalent to "nonspecific" and "none"
    if (!enzyme_constraint_)
      return 0;

    if (!searchEnzymeInstance_) {
      searchEnzymeInstance_ = getSearchEnzyme(true);
    }

    return searchEnzymeInstance_->getMinNumberTermini();
  }
}

Boolean MascotConverter::isSemiCleavageSearch () {
  if (!enzyme_constraint_) return False;

  if (!searchEnzymeInstance_) {
    searchEnzymeInstance_ = getSearchEnzyme (true);
  }
  return searchEnzymeInstance_->isSemiSpecific();
}

void MascotConverter::write_pepXML(const char* outfile) {

  // mzXML file name, for reading of instrument info purposes only
  int len;
  char* mzXMLfile = new char[len=(int)strlen(outfile)+11];
  const char *rawDatFileExt="?";

  const char* match = hasValidPepXMLFilenameExt(outfile);
  if(match != NULL) {
    strncpy(mzXMLfile, outfile, strlen(outfile) - strlen(match));
    mzXMLfile[strlen(outfile) - strlen(match)] = 0;
    rampConstructInputFileName(mzXMLfile,len,mzXMLfile);
    struct stat statbuf;
    if (stat(mzXMLfile,&statbuf)) { // no exist - try using a scan name
      if (spectra_.length()) {
	char *copy = strdup(spectra_[0]);
	char *dot = strchr(copy,'.');
	if (dot) {
	  *dot = 0;
	}
	int len;
	char *buf = (char *)malloc(len=((int)strlen(copy)+30));
	strcpy(buf,copy);
	rampConstructInputFileName(buf,len,buf);
	if (stat(buf,&statbuf)) { // no exist - path thing?
	  char *slash = findRightmostPathSeperator(mzXMLfile);
	  if (slash) {
	    char *buf2 = (char *)malloc(len=((int)strlen(buf)+(int)strlen(mzXMLfile)+30));
	    strcpy(buf2,mzXMLfile);
	    *(findRightmostPathSeperator(buf2)+1) = 0;
	    strcat(buf2,copy);
	    free(buf);
	    buf = buf2; // no dangling pointers!
	    rampConstructInputFileName(buf,len,buf); // try again with path
	  }
	}
	free(copy);

	if (!stat(buf,&statbuf)) { // exists!
	  delete[] mzXMLfile;
	  mzXMLfile = strCopy(buf);
	}
	free(buf);
      }
    }
    rawDatFileExt = strrchr(mzXMLfile,'.');
  }
  else {
    cerr << "error, cannot parse output xmlfile: " << outfile << endl;
    //exit(1);
  }

  InstrumentStruct* ppstruct = NULL;

  RAMPFILE* pFI;
  if ( (pFI=rampOpenFile(mzXMLfile))==NULL) {
    cout << " warning: cannot open \"" << mzXMLfile << "\" for reading MS instrument info." << endl;
    //exit(1);
  }
  else {
    ppstruct = getInstrumentStruct(pFI);
    rampCloseFile(pFI);
  }

  FILE* fout;

  if( (fout = fopen(outfile, "w")) == NULL) {
    cerr << "cannot open outfile " << outfile << endl;
    exit(1);
  }

  Boolean nucleic_acid_db = False;
  if (!sampleEnzymeInstance_) {
    sampleEnzymeInstance_ = getSampleEnzyme (true);
  }
  ModificationInfo* modinfo = NULL;

  fprintf(fout, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
 
  fprintf(fout, "<?xml-stylesheet type=\"text/xsl\" href=\"%s%s\"?>\n", getPepXML_std_xsl_web_path(), "pepXML_std.xsl");

#ifndef USING_RELATIVE_WEBSERVER_PATH 
  char *dt;
  fprintf(fout, "<msms_pipeline_analysis date=\"%s\" xmlns=\"%s\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"%s %s%s\" summary_xml=\"%s\">\n", dt=Parser::getDateTime(), PEPXML_NAMESPACE, PEPXML_NAMESPACE, getPepXML_std_xsl_web_path(), PEPXML_SCHEMA, outfile);
  delete[] dt;
#else
  char* szWebserverRoot = new char[256];
  char* szCommand = new char[256];
  char* szBuf = new char[SIZE_BUF];
  char *pStr=getenv("WEBSERVER_ROOT");
  if (pStr==NULL)
    {
      printf("<PRE> Environment variable WEBSERVER_ROOT does not exist.\n\n");
      printf(" For Windows users, you can set this environment variable\n");
      printf(" through the Advanced tab under System Properties when you\n");
      printf(" right-mouse-click on your My Computer icon.\n\n");
      
      printf(" Set this environment variable to your webserver's document\n");
      printf(" root directory such as c:\\inetpub\\wwwroot for IIS or\n");
      printf(" c:\\website\\htdocs or WebSite Pro.\n\n");
      printf(" Exiting.\n");
      exit(0);
    }
  else {
    strcpy(szWebserverRoot, pStr);
  }
  char *datetime;
  fprintf(fout, "<msms_pipeline_analysis date=\"%s\" xmlns=\"%s\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"%s %s%s\" summary_xml=\"%s\">\n", datetime=Parser::getDateTime(), PEPXML_NAMESPACE, PEPXML_NAMESPACE, szWebserverRoot, PEPXML_SCHEMA, outfile);
  delete [] datetime;
  delete [] szWebserverRoot;
  delete [] szBuf;
  delete [] szCommand;
#endif

  int tlen;
  char* full_file_path = new char[tlen=((int)strlen(filextn_)+(int)strlen(fullpath_)+50)];
  char* full_rawfile_path = new char[tlen];;
  if(strlen(filextn_) > 0 && filextn_[0] == '/')
    strcpy(full_file_path, filextn_);
  else
    sprintf(full_file_path, "%s%s", fullpath_, filextn_);
  strcpy(full_rawfile_path,full_file_path);

  // is this a suitable basename, or is it a rename by Mascot?
  char *testpath = new char[tlen+1];
  rampConstructInputFileName(testpath,tlen,full_rawfile_path);
  struct stat s;
  if (stat(testpath,&s)) {
    // try to construct from query names
    if (spectra_.size()) {
      char *common = strdup(spectra_[spectra_.size()-1]);
      for (int i=spectra_.size();i--;) {
	char *s = spectra_[i];
	for (char *c=common;*c && *s;c++) {
	  if (*c != *s++) {
	    *c = 0;
	    break;
	  }
	}
      }
      if (*common) { // there was some commonality in spectra names
	if ('.'==common[strlen(common)-1]) {
	  common[strlen(common)-1] = 0;
	}
	delete[] testpath;
	tlen=(int)strlen(fullpath_)+(int)strlen(common)+12;
	testpath = new char[tlen];
	sprintf(testpath, "%s%s", fullpath_, common);
	rampConstructInputFileName(testpath,tlen,testpath);
	if (!stat(testpath,&s)) {
	  delete[] full_rawfile_path;
	  full_rawfile_path = new char[tlen+1];
	  sprintf(full_rawfile_path, "%s%s", fullpath_, common);
	}
      }
      free(common);
    }
  }
  delete[] testpath;

  if(ppstruct != NULL) 
    fprintf(fout, "  <msms_run_summary base_name=\"%s\" msManufacturer=\"%s\" msModel=\"%s\" msIonization=\"%s\" msMassAnalyzer=\"%s\" msDetector=\"%s\" raw_data_type=\"%s\" raw_data=\"%s\">\n", full_file_path, ppstruct->manufacturer, ppstruct->model, ppstruct->ionisation, ppstruct->analyzer, ppstruct->detector, "raw", rawDatFileExt);
  else // leave out mass spec info
    fprintf(fout, "  <msms_run_summary base_name=\"%s\" raw_data_type=\"%s\" raw_data=\"%s\">\n", full_file_path, "raw", rawDatFileExt);

  ////////////////////////////////////////////////////////////////////////////////////
  // THIS MUST BE FILLED IN
  sampleEnzymeInstance_->writeTraditionalPepXMLTags(fout);
  if(mass_type_ == 0)
    fprintf(fout, "<search_summary base_name=\"%s\" search_engine=\"MASCOT\" precursor_mass_type=\"average\" fragment_mass_type=\"average\" out_data_type=\"out\" out_data=\".tgz\" search_id=\"1\">\n", full_file_path);
  else
    fprintf(fout, "<search_summary base_name=\"%s\" search_engine=\"MASCOT\" precursor_mass_type=\"monoisotopic\" fragment_mass_type=\"monoisotopic\" out_data_type=\"out\" out_data=\".tgz\" search_id=\"1\">\n", full_file_path);

  fprintf(fout,"<search_database local_path=\"%s\" type=\"%s\"/>\n", database_, nucleic_acid_db ? "NA" : "AA");
  // THIS MUST BE FILLED IN

  if(enzyme_constraint_) {
    if (!searchEnzymeInstance_) {
      searchEnzymeInstance_ = getSearchEnzyme (true);
    }
    if(searchEnzymeInstance_ == NULL) {
      cerr << "error: enzyme in search constraint \"" << sample_enzyme_ << "\" is not recognized" << endl;
      exit(1);
    }
    fprintf(fout, "<enzymatic_search_constraint enzyme=\"%s\" max_num_internal_cleavages=\"%d\" min_number_termini=\"%d\"/>\n", searchEnzymeInstance_->getName(), getMaxInternalCleavages(), getMinNumberTermini());
  }

  // MUST PUT MODIFICATION STUFF HERE]
  int k;
  for(k = 0; k < static_mods_.length(); k++) {
    if(static_mods_[k]->aa == 'n' || static_mods_[k]->aa == 'c')
      fprintf(fout, "<terminal_modification terminus=\"%c\" mass=\"%f\" massdiff=\"%f\" variable=\"N\" protein_terminus=\"%c\"/>\n", static_mods_[k]->aa, static_mods_[k]->mass, static_mods_[k]->massdiff, static_mods_[k]->prot_term ? 'Y' : 'N');
    else {
      if (static_mods_[k]->n_terminal ||static_mods_[k]->c_terminal ) {
        if (static_mods_[k]->n_terminal && static_mods_[k]->c_terminal ) {
          fprintf(fout, "<aminoacid_modification aminoacid=\"%c\" mass=\"%f\" massdiff=\"%f\" peptide_terminus=\"nc\" variable=\"N\"/>\n", static_mods_[k]->aa, static_mods_[k]->mass, static_mods_[k]->massdiff);
        }
        else if (static_mods_[k]->n_terminal) {
          fprintf(fout, "<aminoacid_modification aminoacid=\"%c\" mass=\"%f\" massdiff=\"%f\" peptide_terminus=\"n\" variable=\"N\"/>\n", static_mods_[k]->aa, static_mods_[k]->mass, static_mods_[k]->massdiff);
        }
        else {
          fprintf(fout, "<aminoacid_modification aminoacid=\"%c\" mass=\"%f\" massdiff=\"%f\" peptide_terminus=\"c\" variable=\"N\"/>\n", static_mods_[k]->aa, static_mods_[k]->mass, static_mods_[k]->massdiff);
        }
      }
      else {
        fprintf(fout, "<aminoacid_modification aminoacid=\"%c\" mass=\"%f\" massdiff=\"%f\" variable=\"N\"/>\n", static_mods_[k]->aa, static_mods_[k]->mass, static_mods_[k]->massdiff);
      }
    }
  }
  for(k = 0; k < variable_mods_.length(); k++) {
    if(variable_mods_[k]->aa == 'n' || variable_mods_[k]->aa == 'c')
      fprintf(fout, "<terminal_modification terminus=\"%c\" mass=\"%f\" massdiff=\"%f\" variable=\"Y\" protein_terminus=\"%c\"/>\n", variable_mods_[k]->aa, variable_mods_[k]->mass, variable_mods_[k]->massdiff, variable_mods_[k]->prot_term ? 'Y' : 'N');
    else {
      if (variable_mods_[k]->n_terminal ||variable_mods_[k]->c_terminal ) {
        if (variable_mods_[k]->n_terminal && variable_mods_[k]->c_terminal ) {
          fprintf(fout, "<aminoacid_modification aminoacid=\"%c\" mass=\"%f\" massdiff=\"%f\" peptide_terminus=\"nc\" variable=\"Y\"/>\n", variable_mods_[k]->aa, variable_mods_[k]->mass, variable_mods_[k]->massdiff);
        }
        else if (variable_mods_[k]->n_terminal) {
          fprintf(fout, "<aminoacid_modification aminoacid=\"%c\" mass=\"%f\" massdiff=\"%f\" peptide_terminus=\"n\" variable=\"Y\"/>\n", variable_mods_[k]->aa, variable_mods_[k]->mass, variable_mods_[k]->massdiff);
        }
        else {
          fprintf(fout, "<aminoacid_modification aminoacid=\"%c\" mass=\"%f\" massdiff=\"%f\" peptide_terminus=\"c\" variable=\"Y\"/>\n", variable_mods_[k]->aa, variable_mods_[k]->mass, variable_mods_[k]->massdiff);
        }
      }
      else {
        fprintf(fout, "<aminoacid_modification aminoacid=\"%c\" mass=\"%f\" massdiff=\"%f\" variable=\"Y\"/>\n", variable_mods_[k]->aa, variable_mods_[k]->mass, variable_mods_[k]->massdiff);
      }
    }

  }
  // put the Mascot specific parameters here....

  for(k = 0; k < params_.length(); k++) {
    char *escapedVals = xmlEscape(vals_[k]);
    fprintf(fout, "<parameter name=\"%s\" value=\"%s\"/>\n", params_[k], escapedVals);
    delete [] escapedVals;
  }
  fprintf(fout, "</search_summary>\n");

  int counter = 1;
  int spec_divisors[] = {-1, -1, -1};

  Array<Tag*>* modtags = NULL;

  for(k = 0; k < spectra_.length(); k++) {

    if(checkResult(k)) {

      if((spectra_[k])[strlen(spectra_[k])-1] != precursor_ioncharges_[k] + '0') {
        cout << "error: spectrum name " << spectra_[k] << " inconsistent with charge " << precursor_ioncharges_[k] << endl;
        exit(1);
      }

      fprintf(fout, "    <spectrum_query spectrum=\"%s\"", spectra_[k]);
    
      int index = 2;
      int j = (int)strlen(spectra_[k])-1;
      while((index+1) && j) {
        if((spectra_[k])[j] == '.')
          spec_divisors[index--] = j;
        j--;
      }
      if(index+1) {
        cerr << " error: only found " << (2 - index) << " periods in " << spectra_[k] << endl;
        exit(1);
      }
      
      fprintf(fout, " start_scan=\"");
      for(j = spec_divisors[0]+1; j < spec_divisors[1]; j++)
        fprintf(fout, "%c", (spectra_[k])[j]);
      fprintf(fout, "\" end_scan=\"");
      for(j = spec_divisors[1]+1; j < spec_divisors[2]; j++)
        fprintf(fout, "%c", (spectra_[k])[j]);
      fprintf(fout, "\" precursor_neutral_mass=\"%0.4f\" assumed_charge=\"%d\" index=\"%d\"", masses_[k] - 1.0 + massdiffs_[k], precursor_ioncharges_[k], counter++);

      {
        double rtSec = -1.0;
        // let's get the start and end scan numbers from the spectra label
        long lStartScan=0;
        const char *szLabel=spectra_[k];
        szLabel+=(spec_divisors[0]+1);
        sscanf(szLabel,"%ld", &lStartScan);
        rtSec = scanToRetentionTime(lStartScan);
        if (-1!=rtSec) {
          fprintf(fout, " retention_time_sec=\"%f\"", rtSec);
        }
      }
      fprintf(fout, ">\n");
      
      int num_ion_series_ = 2;
      int totnumions = ((int)strlen(peptides_[k]) - 1) * num_ion_series_;
      if((spectra_[k])[strlen(spectra_[k])-1] > '2')
        totnumions *= 2;
#if 0 //WCH: this is likely translated from writeOutFile
      // now format matched ions to 3 digits
      char matched[] = "   ";
      char totions[] = "   ";

      if(matchedions_[k] < 10)
        sprintf(matched+2, "%d", matchedions_[k]);
      else if(matchedions_[k] < 100)
        sprintf(matched+1, "%d", matchedions_[k]);
      else 
        sprintf(matched, "%d", matchedions_[k]);
      if(totnumions < 10)
        sprintf(totions+2, "%d", totnumions);
      else if(totnumions < 100)
        sprintf(totions+1, "%d", totnumions);
      else 
        sprintf(totions, "%d", totnumions);
#endif
      fprintf(fout, "    <search_result>\n");
     
      fprintf(fout, "      <search_hit hit_rank=\"%d\" peptide=\"%s\"", 1, peptides_[k]);

      if(prev_AAs_[k] != G_UNKNOWN_AA) 
        fprintf(fout, " peptide_prev_aa=\"%c\"", prev_AAs_[k]);
      if(foll_AAs_[k] != G_UNKNOWN_AA)
        fprintf(fout, " peptide_next_aa=\"%c\"", foll_AAs_[k]);

      // int num_compat_term = isSemiCleavageSearch() ? 2 : enzyme->getNumTolTerm(prev_AAs_[k], peptides_[k], foll_AAs_[k]);
      int num_compat_term = sampleEnzymeInstance_ == NULL ? 2 : sampleEnzymeInstance_->getNumTolTerm(prev_AAs_[k], peptides_[k], foll_AAs_[k]);

      const char *szProtein = proteins_[k];
      if (!short_protein_id_) {
        idMap::const_iterator it = proteinFullIdMap_.find (proteins_[k]);
        if (it != proteinFullIdMap_.end()) {
          szProtein = it->second.c_str();
        }
      }
      // watch out for FASTA with reserved XML characters
      std::string escapedProtein(XMLEscape(std::string(szProtein)));
      szProtein = escapedProtein.c_str();

      if(massdiffs_[k] < 0.0) 
        //fprintf(fout, " protein=\"%s\" num_tot_proteins=\"%d\" num_matched_ions=\"%s\" tot_num_ions=\"%s\" calc_neutral_pep_mass=\"%0.4f\" massdiff=\"%0.1f\" num_tol_term=\"%d\" num_missed_cleavages=\"%d\"", szProtein, pluses_[k]+1, matched, totions, masses_[k] - 1.0 /*+massdiffs_[k]*/, massdiffs_[k], num_compat_term, enzyme->getNumMissedCleavages(peptides_[k]));
        fprintf(fout, " protein=\"%s\" num_tot_proteins=\"%d\" num_matched_ions=\"%d\" tot_num_ions=\"%d\" calc_neutral_pep_mass=\"%0.4f\" massdiff=\"%0.4f\" num_tol_term=\"%d\" num_missed_cleavages=\"%d\"", szProtein, pluses_[k]+1, matchedions_[k], totnumions, masses_[k] - 1.0 /*+massdiffs_[k]*/, massdiffs_[k], num_compat_term, sampleEnzymeInstance_->getNumMissedCleavages(peptides_[k]));

      else 
        //fprintf(fout, " protein=\"%s\" num_tot_proteins=\"%d\" num_matched_ions=\"%s\" tot_num_ions=\"%s\" calc_neutral_pep_mass=\"%0.4f\" massdiff=\"+%0.1f\" num_tol_term=\"%d\" num_missed_cleavages=\"%d\"", szProtein, pluses_[k]+1, matched, totions, masses_[k] - 1.0 /*+massdiffs_[k]*/, massdiffs_[k], num_compat_term, enzyme->getNumMissedCleavages(peptides_[k]));
        fprintf(fout, " protein=\"%s\" num_tot_proteins=\"%d\" num_matched_ions=\"%d\" tot_num_ions=\"%d\" calc_neutral_pep_mass=\"%0.4f\" massdiff=\"+%0.4f\" num_tol_term=\"%d\" num_missed_cleavages=\"%d\"", szProtein, pluses_[k]+1, matchedions_[k], totnumions, masses_[k] - 1.0 /*+massdiffs_[k]*/, massdiffs_[k], num_compat_term, sampleEnzymeInstance_->getNumMissedCleavages(peptides_[k]));

      modtags = getModificationTags(peptides_[k], mods_[k]);

      if(calculate_pI_ && modtags != NULL && pi_calc_ != NULL) {
        modinfo = new ModificationInfo(modtags);
        fprintf(fout, " calc_pI=\"%0.2f\"", pi_calc_->Peptide_pI(peptides_[k], modinfo));
        delete modinfo;
        modinfo = NULL;
      }
      else if(calculate_pI_)
        fprintf(fout, " calc_pI=\"%0.2f\"", pi_calc_->Peptide_pI(peptides_[k], NULL));

      fprintf(fout, " is_rejected=\"%d\"", 0);

      if (generate_description_) {
        proteinMap::const_iterator it = proteinDescriptionMap_.find (proteins_[k]);
        if (it != proteinDescriptionMap_.end()) {
          fprintf(fout, " protein_descr=\"%s\"", XMLEscape(it->second).c_str());
        }
        else {
          fprintf(fout, " protein_descr=\"NON_EXISTENT PROTEIN DESCRIPTION\"");
        }
      }

      fprintf(fout, ">\n");

      // print out alternative protein sharing hit/rank#1
      int rank1Index=rank1_inds_[k];
      for (int plusesIndex=0; plusesIndex<pluses_[k]; ++plusesIndex,++rank1Index) {

        int num_compat_term = sampleEnzymeInstance_ == NULL ? 2 : sampleEnzymeInstance_->getNumTolTerm(rank1_prev_AAs_[rank1Index], peptides_[k], rank1_foll_AAs_[rank1Index]);
        const char *szProtein = rank1_proteins_[rank1Index];
        if (!short_protein_id_) {
          idMap::const_iterator it = proteinFullIdMap_.find (rank1_proteins_[rank1Index]);
          if (it != proteinFullIdMap_.end()) {
            szProtein = it->second.c_str();
          }
        }
	// watch out for FASTA with reserved XML characters
	std::string escapedProtein(XMLEscape(std::string(szProtein)));
	szProtein = escapedProtein.c_str();

        fprintf(fout, "      <alternative_protein");
        fprintf(fout, " protein=\"%s\"", szProtein);
        if (generate_description_) {
          proteinMap::const_iterator it = proteinDescriptionMap_.find (rank1_proteins_[rank1Index]);
          if (it != proteinDescriptionMap_.end()) {
	    fprintf(fout, " protein_descr=\"%s\"", XMLEscape(it->second).c_str());
          }
          else {
            fprintf(fout, " protein_descr=\"NON_EXISTENT PROTEIN DESCRIPTION\"");
          }
        }
        fprintf(fout, " num_tol_term=\"%d\"", num_compat_term);
        fprintf(fout, " peptide_prev_aa=\"%c\"", rank1_prev_AAs_[rank1Index]);
        fprintf(fout, " peptide_next_aa=\"%c\"", rank1_foll_AAs_[rank1Index]);
        fprintf(fout, "/>\n");
      }

      if(modtags != NULL) {
        for(int k = 0; k < modtags->length(); k++)
          if((*modtags)[k] != NULL) {
            (*modtags)[k]->writeTraditional(fout);
            delete (*modtags)[k];
          }

        delete modtags;

      }

      fprintf(fout, "        <search_score name=\"ionscore\" value=\"%0.2f\"/>\n", ionscores_[k]);
      fprintf(fout, "        <search_score name=\"identityscore\" value=\"%0.2f\"/>\n", identityscores_[k]);

      int star = isStar(k) ? 1 : 0;
      fprintf(fout, "        <search_score name=\"star\" value=\"%d\"/>\n", star);
      fprintf(fout, "        <search_score name=\"homologyscore\" value=\"%0.2f\"/>\n", homologyscores_[k]);
      fprintf(fout, "        <search_score name=\"expect\" value=\"%0.4f\"/>\n", expectationvalues_[k]);

      fprintf(fout, "      </search_hit>\n");
      fprintf(fout, "    </search_result>\n");
      fprintf(fout, "    </spectrum_query>\n");

    } // if ok result
  } // next
  fprintf(fout, "  </msms_run_summary>\n");
  fprintf(fout, "</msms_pipeline_analysis>\n");

  fclose(fout);

  delete[] mzXMLfile;
  delete[] full_file_path;

}

Boolean MascotConverter::isStar(int index) {
  double min_ionscore_pct = 0.0; //0.9; // 0.96;
  if(index < spectra_.length() - 1 && alt_inds_[index] == alt_inds_[index+1])
    return False;
  if(index == spectra_.length() - 1 && alt_inds_[index] >= alt_proteins_.length())
    return False;
  int stop = alt_proteins_.length();
  if(index < alt_inds_.length() - 1)
    stop = alt_inds_[index+1];
  int k = alt_inds_[index];
  while(k < stop && alt_ionscores_[k] >= ionscores_[index] * min_ionscore_pct) 
    if(homologous(peptides_[index], alt_peptides_[k++])) 
      return True;

  return False;
}

// J.Eng
Boolean MascotConverter::homologous(char* pep1, char* pep2) {
  if(pep1 == NULL || pep2 == NULL)
    return False;

  int peplen1 = (int)strlen(pep1);
  int iDiffInSequence = abs(peplen1-(int)strlen(pep2));

  /*
   * These sequences contain preceeding and trailing amino acids
   */
 
  double MIN_PERCENTAGE = 0.75;

  for (int i=0; pep1[i]&&pep2[i]; i++)
    {
      /*
       * K/Q and I/L don't count as differences
       */
      if (pep1[i] != pep2[i])
	{
	  if (!((pep1[i]=='K'|| pep1[i]=='Q') &&
		(pep2[i]=='K'|| pep2[i]=='Q')) &&
	      !((pep1[i]=='I'|| pep1[i]=='L')&&
		(pep2[i]=='I'|| pep2[i]=='L')))
            {
	      iDiffInSequence++;
            }
	}
    }

  /*
   * Calculate deltCn only if sequences are less than PERCENTAGE similar
   * I'm not sure why there's a -3 here (vs. -4 to take care of the 2
   * preceeding and 2 trailing characters)
   *
   * PERCENTAGE is defined as 0.75
   */

  return peplen1 > iDiffInSequence && (double)(peplen1 - iDiffInSequence)/peplen1 >= MIN_PERCENTAGE;

}

Boolean MascotConverter::validate_db() {
  Boolean result = False;
  if(!(database_ == NULL || strcmp(database_, "") == 0)) {
    unCygwinify(database_); // no effect in cygwin builds
    RACI ftest(database_); // can read gzipped files
    if( ftest.is_open()) {
      result = True;
      ftest.close();
    }
  }
  if(! result ) {
    cerr << "could not find database ";
    if(database_ != NULL)
      cerr << database_;
    cerr << endl;
    exit(1);
  } else {
    if (!findRightmostPathSeperator(database_)) {
      // things downstream work better with path info
      const int nlen=5000;
      char *newdb = new char[nlen];
      safepath_getcwd(newdb, nlen);
      strcat(newdb, "/");
      strcat(newdb,database_);
      delete[] database_;
      database_ = newdb;
    }
  }
  return result;
}

Boolean MascotConverter::checkResult(int index) {
  return index >= 0 && index < peptides_.length() &&
    peptides_[index] != NULL && strlen(peptides_[index]) > 0;
}

char *MascotConverter::xmlEscape(const char *string) const {
  // first pass to get the final size
  int finalSize=0;
  int originalLen=(int)strlen(string);
  for (int i=0; i<originalLen; i++) {
    switch(string[i]) {
    case '&': finalSize+=5; break;
    case '<': finalSize+=4; break;
    case '>': finalSize+=4; break;
    case '"': finalSize+=6; break;
    case '\'': finalSize+=6; break;
    default: finalSize++; break;
    }
  }
  char *escapedString = new char[finalSize+1];
  escapedString[0] = 0;
  if (originalLen==finalSize) {
    strcpy (escapedString, string);
  } else {
    finalSize=0;
    for (int i=0; i<originalLen; i++) {
      switch(string[i]) {
      case '&':
	escapedString[finalSize++]='&';
	escapedString[finalSize++]='a';
	escapedString[finalSize++]='m';
	escapedString[finalSize++]='p';
	escapedString[finalSize++]=';';
	break;
      case '<':
	escapedString[finalSize++]='&';
	escapedString[finalSize++]='l';
	escapedString[finalSize++]='t';
	escapedString[finalSize++]=';';
	break;
      case '>':
	escapedString[finalSize++]='&';
	escapedString[finalSize++]='g';
	escapedString[finalSize++]='t';
	escapedString[finalSize++]=';';
	break;
      case '"':
	escapedString[finalSize++]='&';
	escapedString[finalSize++]='q';
	escapedString[finalSize++]='u';
	escapedString[finalSize++]='o';
	escapedString[finalSize++]='t';
	escapedString[finalSize++]=';';
	break;
      case '\'':
	escapedString[finalSize++]='&';
	escapedString[finalSize++]='a';
	escapedString[finalSize++]='p';
	escapedString[finalSize++]='o';
	escapedString[finalSize++]='s';
	escapedString[finalSize++]=';';
	break;
      default:
	escapedString[finalSize++]=string[i];
	break;
      }
    }
    escapedString[finalSize]=0;
  }

  return escapedString;
}

// compare the sequences in Mascot style
// * treated as X
// X treated as any of the 20 aa
// B matches N and D
// Z matches E and Q
char *MascotConverter::locatePeptide(const char* aasequence, const char *aasubsequence) {
  char *aa = (char *) aasequence;
  char *aa1, *aa2;

  if ( !*aasubsequence )
    return((char *)aasequence);

  while (*aa)
    {
      aa1 = aa;
      aa2 = (char *) aasubsequence;

      while ( *aa1 && *aa2 && 
	      (((*aa1)==(*aa2)) || 
	       ('*' == *aa1 || '*' == *aa2 || 'X' == *aa1 || 'X' == *aa2 ) || 
	       ('B' == *aa1 && ('N' == *aa2 || 'D' == *aa2)) || 
	       ('B' == *aa2 && ('N' == *aa1 || 'D' == *aa1)) || 
	       ('Z' == *aa1 && ('E' == *aa2 || 'Q' == *aa2)) || 
	       ('Z' == *aa2 && ('E' == *aa1 || 'Q' == *aa1))
	       )
	      )
	aa1++, aa2++;

      if (!*aa2)
	return(aa);

      aa++;
    }

  return(NULL);
}


// update the peptide terminus
// optionally, we update the protein description (when we are still guessing)
inline void MascotConverter::updatePeptideTerminus (
  const char *dbProtein, const char *dbProteinDesc, const char *dbProteinSequence, 
  int &proteinIDType, int &proteinIDGuessedType, int &proteinIDTypeGuessCount, 
  Boolean removeMapEntry, 
  Array<char*> &proteins, Array<char*> &peptides, 
  Array<char> &prevAAs, Array<char> &follAAs, Array<int> &indices) {

  int k;
  Boolean firstGuess = True;
  for(k = 0; k < indices.length(); k++) {
    if(strstr(dbProtein, proteins[indices[k]]) != NULL) { // match
      // let's guess and check our guess
      if (firstGuess && UNKNOWN == proteinIDType) {
        int guessType = guessProteinIDType (dbProtein, proteins[indices[k]]);
        if (UNKNOWN == guessType) {
          // we do not know, and let's forget about future check then
          proteinIDType = STATE_GIVEUP;
        }
        else {
          if (UNKNOWN == proteinIDGuessedType) {
            proteinIDGuessedType = guessType;
          }
          else if (proteinIDGuessedType == guessType) {
            proteinIDTypeGuessCount++;
            if (proteinIDTypeGuessCount >= MIN_CONSECUTIVE_CORRECT_GUESS)
	      {
		proteinIDType = proteinIDGuessedType;
	      }
          }
          else {
            // new guess is not the same as previous conclusion!
            // we should give up guessing
            proteinIDType = STATE_GIVEUP;
          }
        }

        // we do want to count multiple peptide from the same protein
        // as multiple protein ID match
        firstGuess = False;
      }

      if (removeMapEntry && !missingDescriptionMap_.empty()) {
        // let's check if we need this description?
        proteinMap::iterator it = missingDescriptionMap_.find (proteins[indices[k]]);
        if (it != missingDescriptionMap_.end()) {
          // no description available, let's use the supplied one
          std::string id(proteins[indices[k]]);
          std::string description(dbProteinDesc);
          proteinDescriptionMap_.insert(std::make_pair(id,description));
          missingDescriptionMap_.erase(it);
        }
      }

      //match = strstr(dbProteinSequence, peptides[indices[k]]);
      char *match = locatePeptide(dbProteinSequence, peptides[indices[k]]);
      if(match != NULL) {

        if (!short_protein_id_) {
          idMap::iterator it = proteinFullIdMap_.find (proteins[indices[k]]);
          if (it == proteinFullIdMap_.end()) {
            // no full ids version available, let's cache it
            std::string shortid(proteins[indices[k]]);
            std::string fullids(dbProtein);
            proteinFullIdMap_.insert(std::make_pair(shortid,fullids));
          }
        }
#if 0
        // cheehong: original code here, we use look up approach to conserve memory
        if(strcmp(dbProtein, proteins[indices[k]]) != 0) { // correct...
          if(proteins[indices[k]] != NULL)
            delete[] proteins[indices[k]];
          char* newprot = new char[strlen(dbProtein)+1];
          strcpy(newprot, dbProtein);
          proteins.replace(indices[k], newprot);
        }
#endif

        if (!term_read_) {
          if(strlen(match) < strlen(dbProteinSequence)) 
            prevAAs.replace(indices[k], mapAA(dbProteinSequence[strlen(dbProteinSequence) - strlen(match) - 1]));
          else
            prevAAs.replace(indices[k], '-');
          if(strlen(match) > strlen(peptides[indices[k]]))
            follAAs.replace(indices[k], mapAA(match[strlen(peptides[indices[k]])]));
          else
            follAAs.replace(indices[k], '-');
        }

        indices.remove(k);
        k--;
      } // match for peptide
      else if (term_read_) {
        indices.remove(k);
        k--;
      }
    }
  }
}


// Extract the protein id from fasta description
//   perform lookup approach for protein description and peptide terminus 
// Otherwise, original system is used
inline void MascotConverter::buildProteinPeptideInfoViaLookup (
  const char *dbProtein, const char *dbProteinDesc, const char *dbProteinSequence, 
  int &proteinIDType, int &proteinIDGuessedType, int &proteinIDTypeGuessCount, 
  Array<int> &indices, Array<int> &rank1_indices, Array<int> &alt_indices) {

  if (proteinIDType < ID_TYPE_VALID)
    return;

  // we have established the protein id format
  // we can now use a faster processing approach
  std::string proteinId(dbProtein);
  if (ID_TYPE_FULL == proteinIDType) {
    // do nothing, we use the whole string
  }
  else if (ID_TYPE_RIGHTMOST_PLUS_LEFT_ANCHOR == proteinIDType) {
    const char *szStart = strrchr (dbProtein, proteinIDLeftAnchor_[0]);
    if (NULL==szStart) {
      //TODO: what do we want to do? take the whole?
    } else {
      proteinId.erase (0, szStart-dbProtein+1);
    }
  }
  else if (ID_TYPE_LEFTMOST_PLUS_RIGHT_ANCHOR == proteinIDType) {
    const char *szEnd = strchr (dbProtein, proteinIDRightAnchor_[0]);
    if (NULL==szEnd) {
      //TODO: what do we want to do? take the whole?
    } else {
      proteinId.erase (szEnd-dbProtein);
    }
  }
  else if (ID_TYPE_CENTER_PLUS_LEFT_RIGHT_ANCHOR == proteinIDType) {
    const char *szStart = strchr (dbProtein, proteinIDLeftAnchor_[0]);
    if (NULL==szStart) {
      //TODO: what do we want to do? take the whole?
    } else {
      proteinId.erase (0, szStart-dbProtein+1);
    }
    const char *szNewStart = proteinId.c_str ();
    const char *szEnd = strchr (szNewStart, proteinIDRightAnchor_[0]);
    if (NULL==szEnd) {
      //TODO: what do we want to do? take the whole?
    } else {
      proteinId.erase (szEnd-szNewStart);
    }
  }
  else {
    cerr << "ERROR: unknown protein id type " << proteinIDType << endl;
    exit(1);
  }

  Boolean matched = False;
  if (!missingDescriptionMap_.empty()) {
    // let's check if we need this description?
    proteinMap::iterator it = missingDescriptionMap_.find (proteinId);
    if (it != missingDescriptionMap_.end()) {
      // no description available, let's build one
      std::string description;
      description.append(dbProteinDesc);
      proteinDescriptionMap_.insert(std::make_pair(proteinId,description));
      missingDescriptionMap_.erase(it);
      matched = True;
    }
  }

  if (!matched) {
    proteinMap::iterator it = proteinDescriptionMap_.find (proteinId);
    matched = (it != proteinDescriptionMap_.end());
  }

  if (matched) {
    if (!term_read_) {
      updatePeptideTerminus (dbProtein, dbProteinDesc, dbProteinSequence, 
			     proteinIDType, proteinIDGuessedType, proteinIDTypeGuessCount, False, 
			     proteins_, peptides_, prev_AAs_, foll_AAs_, indices);
      updatePeptideTerminus (dbProtein, dbProteinDesc, dbProteinSequence, 
			     proteinIDType, proteinIDGuessedType, proteinIDTypeGuessCount, False, 
			     rank1_proteins_, rank1_peptides_, rank1_prev_AAs_, rank1_foll_AAs_, rank1_indices);
      updatePeptideTerminus (dbProtein, dbProteinDesc, dbProteinSequence, 
			     proteinIDType, proteinIDGuessedType, proteinIDTypeGuessCount, False, 
			     alt_proteins_, alt_peptides_, alt_prev_AAs_, alt_foll_AAs_, alt_indices);
    }

    if (!short_protein_id_) {
      idMap::iterator it = proteinFullIdMap_.find (proteinId);
      if (it == proteinFullIdMap_.end()) {
        // no full ids version available, let's cache it
        std::string shortid(proteinId);
        std::string fullids(dbProtein);
        proteinFullIdMap_.insert(std::make_pair(shortid,fullids));
      }
    }
  }
}


// if we are certain of the protein id system, 
//   perform lookup approach for protein description and peptide terminus 
// Otherwise, original system is used
inline void MascotConverter::buildProteinPeptideInfo (
  const char *dbProtein, const char *dbProteinDesc, const char *dbProteinSequence, 
  int &proteinIDType, int &proteinIDGuessedType, int &proteinIDTypeGuessCount, 
  Array<int> &indices, Array<int> &rank1_indices, Array<int> &alt_indices) {

  if (proteinIDType >= ID_TYPE_VALID) {
    buildProteinPeptideInfoViaLookup (dbProtein, dbProteinDesc, dbProteinSequence, 
				      proteinIDType, proteinIDGuessedType, proteinIDTypeGuessCount, 
				      indices, rank1_indices, alt_indices);
  }
  else {
    updatePeptideTerminus (dbProtein, dbProteinDesc, dbProteinSequence, 
			   proteinIDType, proteinIDGuessedType, proteinIDTypeGuessCount, True, 
			   proteins_, peptides_, prev_AAs_, foll_AAs_, indices);
    updatePeptideTerminus (dbProtein, dbProteinDesc, dbProteinSequence, 
			   proteinIDType, proteinIDGuessedType, proteinIDTypeGuessCount, True, 
			   rank1_proteins_, rank1_peptides_, rank1_prev_AAs_, rank1_foll_AAs_, rank1_indices);
    updatePeptideTerminus (dbProtein, dbProteinDesc, dbProteinSequence, 
			   proteinIDType, proteinIDGuessedType, proteinIDTypeGuessCount, True, 
			   alt_proteins_, alt_peptides_, alt_prev_AAs_, alt_foll_AAs_, alt_indices);
  }
}


// try to guess the type of proteinID system in the fasta description
inline int MascotConverter::guessProteinIDType (
  const char *dbProtein, 
  const char *resultProtein) {

  // let's try to guess it
  const char *idStart = strstr(dbProtein, resultProtein);
  if (!idStart) {
    return UNKNOWN;
  }

  size_t resultLen = strlen(resultProtein);
  if (idStart == dbProtein) {
    // started from left most
    if (strlen(idStart) == resultLen) {
      // ended right most, ==> full id
      return ID_TYPE_FULL;
    }

    // have to terminated before the right most
    // but how is it terminated?
    char anchor = idStart[resultLen];
    if ('|' == anchor || ';' == anchor) {
      // <resultProtein>|
      // <resultProtein>;
      // if we are using this to search, it cannot appear in <resultProtein>
      if (strchr(resultProtein, anchor)) {
        return UNKNOWN;
      }
      proteinIDRightAnchor_.erase ();
      proteinIDRightAnchor_.insert (proteinIDRightAnchor_.begin (), anchor);
      return ID_TYPE_LEFTMOST_PLUS_RIGHT_ANCHOR;
    }
    else {
      // some other anchor that we did not expect
      return UNKNOWN;
    }
  }

  // we did not start from the left-most
  if (strlen(idStart) == resultLen) {
    // ended right most, ==> can use right boundary as the anchor
    char anchor = idStart[-1];
    if ('|' == anchor || ';' == anchor) {
      // <resultProtein>|
      // <resultProtein>;
      // if we are using this to search, it cannot appear in <resultProtein>
      if (strchr(resultProtein, anchor)) {
        return UNKNOWN;
      }
      proteinIDLeftAnchor_.erase ();
      proteinIDLeftAnchor_.insert (proteinIDLeftAnchor_.begin (), anchor);
      return ID_TYPE_RIGHTMOST_PLUS_LEFT_ANCHOR;
    }
    else {
      // some other anchor that we did not expect
      return UNKNOWN;
    }
  }

  // we did not start from the left-most, and did not end right-most either
  // i.e. id is in the middle of the dbProtein
  char anchor = idStart[-1];
  if ('|' == anchor || ';' == anchor || ':' == anchor) {
    const char *anchorPos = strchr(dbProtein, anchor);
    if ((idStart-1) != anchorPos) {
      return UNKNOWN;
    }
    proteinIDLeftAnchor_.erase ();
    proteinIDLeftAnchor_.insert (proteinIDLeftAnchor_.begin (), anchor);

    char anchor = idStart[resultLen];
    if ('|' == anchor || ';' == anchor || '.' == anchor) {
      const char *anchorPos = strchr(idStart, anchor);
      if ((idStart+resultLen) != anchorPos) {
        return UNKNOWN;
      }
      proteinIDRightAnchor_.erase ();
      proteinIDRightAnchor_.insert (proteinIDRightAnchor_.begin (), anchor);

      return ID_TYPE_CENTER_PLUS_LEFT_RIGHT_ANCHOR;
    }
  }

  // FOR FUTURE:
  // the above handling are meant for very simple id system in fasta description
  // we should improve the handling to detect more id type so that we can get 
  // the speed improvement from looking up a map for the id rather than
  // looping the the entry in the map and check that id is a substr of entry

  return UNKNOWN;
}

/*
  Parse mzXML file for the associated scan number for each <cycle, experiment> pair
  cycleExptToScanMap_ is updated to get scan number (via <cycle, experiment>)
  scanRetentionTimeMap_ is updated to get retention time ( via <scan>)
*/
void MascotConverter::loadCycleExptScanFromMZXML (const char *datFile)
{
  if (!toLoadScanMap_) {
    return;
  }

  if (scanMapLoaded_) {
    toLoadScanMap_ = False;
    return;
  }

  // let's try to load the mapping from pair<cycle,experiment> --> scan
  cycleExptToScanMap_.clear();
  scanRetentionTimeMap_.clear();

  int len;
  char* mzXMLfile = new char[len=(int)strlen(datFile)+11]; /*.mzData, .mzXML*/
  char datsuffix[] = ".dat";
  const char* getSuffix = strstr(datFile, datsuffix);
  if(getSuffix != NULL) {
    int numCharToCopy = (int)( strlen(datFile) - strlen(getSuffix) );
    strncpy(mzXMLfile, datFile, numCharToCopy);
    mzXMLfile[numCharToCopy] = 0;
    rampConstructInputFileName(mzXMLfile,len,mzXMLfile);
  }
  else {
    cerr << "error, cannot parse .dat file: " << datFile << endl;
    return;
  }

  RAMPFILE* pFI;
  if ( (pFI=rampOpenFile(mzXMLfile))==NULL) {
    cout << " warning: cannot open \"" << mzXMLfile << "\" for reading scan numbers." << endl;
    toLoadScanMap_ = False;
  }
  else {
    //read the scan number and generate cycle,experiment association

    //Read the offset of the index
    ramp_fileoffset_t indexOffset = getIndexOffset(pFI);
    if (!indexOffset)
      {
	cerr << "Could not find index offset in " << mzXMLfile << ". File may be corrupted." << endl;
	rampCloseFile(pFI);
	toLoadScanMap_ = False;
	return;
      }

    //Read the scan index into a vector, get LastScan
    int iFirstScan = 1;
    int iLastScan = 0;
    long lCycle = 0;
    int lExperiment = 0;
    ramp_fileoffset_t* pScanIndices = readIndex(pFI , indexOffset, &iLastScan);
    for (long lScan = iFirstScan; lScan <= iLastScan; lScan++)
      {
	//read scan header
	struct ScanHeaderStruct scanHeader;
	readHeader(pFI, pScanIndices[lScan], &scanHeader);
	if (1==scanHeader.msLevel)
	  {
	    lCycle++;
	    lExperiment=1;
	  } else if (2==scanHeader.msLevel) {
	  lExperiment++;
	  cycleExperimentPair key(lCycle,lExperiment);
	  cycleExptToScanMap_.insert(std::make_pair(key, lScan));
	  scanRetentionTimeMap_.insert(std::make_pair(lScan,scanHeader.retentionTime));
	  /*
	    cout << lCycle << "\t" << lExperiment << "\t" << lScan 
	    << "\t" << scanHeader.retentionTime << endl;
          */
	}
      }

    rampCloseFile(pFI);

    scanMapLoaded_ = True; // the map has been loaded
  }
}

/*
  Parse the title value from the query<X> block which is generated
  by Analyst(QS) mascot.dll script
*/
void MascotConverter::parseCyleExperimentForScan (const char *nextspectrum, long &rangeStart, long &rangeEnd)
{
  if (!scanMapLoaded_) return;

  const char cycleTag[] = ",Cycle(s):";
  const char *cycleInfo = strstr(nextspectrum, cycleTag);

  if (NULL==cycleInfo) return;
  cycleInfo+=strlen(cycleTag);

  std::vector<ExperimentCycleRecord> experimentCycles;

  const char experimentTag[] = "(Experiment";
  const char experimentDelimeter[] = ",";

  if (verbose_)
    cout << nextspectrum << endl;

  while (cycleInfo!=NULL) {

    const char *experimentInfo = strstr(cycleInfo, experimentTag);

    if (NULL==experimentInfo) return;
    experimentInfo+=strlen(experimentTag);

    ExperimentCycleRecord record;
    record.lExperiment = atol(experimentInfo);
    char szDelimetered='\0';
    int iNum = sscanf(cycleInfo, "%ld-%ld%c", &(record.lCycleStart), &record.lCycleEnd, &szDelimetered);
    if (1==iNum) {
      record.lCycleEnd = record.lCycleStart;
      record.bRangleCycle = False;
    } else {
      record.bRangleCycle = True;
    }
    record.bSingleCycle=(','==szDelimetered);

    if (iNum>0) {
      experimentCycles.insert(experimentCycles.end(), record);
      if (verbose_)
        cout << "\tExpt:" << record.lExperiment << 
          "\tCycle:" << record.lCycleStart << "," << record.lCycleEnd << 
          "\tSingle:" << record.bSingleCycle << "\tRange:" << record.bRangleCycle << 
          endl;
    }

    const char *nextCycleInfo = strstr (experimentInfo, experimentDelimeter);
    if (NULL!=nextCycleInfo) {
      cycleInfo=nextCycleInfo+strlen(experimentDelimeter);
    } else {
      cycleInfo=NULL;
    }
  }

  // next, we need to select the representative cycle and experiment for this query
  /*
  the selection criteria are:
   NOTE: n(i)<n(i+1), c(i)<c(i+1)
   1. if only one single experiment is present
      a. (n1[,n2[,...]],c) --> (n1,c),(n1,c)
      b. (n1-n2[,n3[,...]],c) --> (n1,c),(n2,c)
   2. if two experiments are present
     [1,1]
      a. (n1,c1) and (n2,c2) --> (n1,c1),(n2,c2)
      b. (n2,c1) and (n1,c2) --> (n1,c1),(n2,c2)
      c. (n1,c1) and (n1,c2) --> (n1,c1),(n1,c1)
     [1,x]
      d. (n1,c1) and (n2,n3,...,c2) --> (n1,c1),(n2,c2)
      e. (n2,c1) and (n1,n3,...,c2) --> (n1,c1),(n2,c2)
      f. (n1,c1) and (n1,n2,...,c2) --> (n1,c1),(n1,c1) [NOT SEEN CASES]
     [x,1],[x,y]
      g. (n1,n3,...,c1) and (n2,c2) --> (n1,c1),(n1,c1)
      h. (n1-n2,n4,...,c1) and (n3,c2) --> (n1,c1),(n2,c2)
   2. if three experiments
      a. (n1, c1), (n2, c2), (n3,c3) --> ???
      b. (n1, c1), (n2, c2) --> (n1,c1),(n2,c2) or (n1,c1),(n2,c3) ???
      c. (n1, c1), (n3,c3) --> (n1,c1),(n3,c3)???
  */

  long minCycle=0;
  long minCycleExperiment=0;
  long maxCycle=0;
  long maxCycleExperiment=0;
  if (0==experimentCycles.size()) {
    // there is no cycle for us select, do nothing
  } else if (1==experimentCycles.size()) {
    // single experiment with possibly a few cycles
    minCycle=experimentCycles[0].lCycleStart;
    minCycleExperiment=experimentCycles[0].lExperiment;
    maxCycle=experimentCycles[0].lCycleEnd;
    maxCycleExperiment=experimentCycles[0].lExperiment;
  } else {
    //2<=experimentCycleStart.size()
    //two or more experiments
    ExperimentCycleRecord experiment1 = experimentCycles[0];
    ExperimentCycleRecord experiment2 = experimentCycles[1];

    if (experiment1.bSingleCycle) {
      // (1,1), (1,n)
      if (experiment1.bRangleCycle) {
        minCycle = experiment1.lCycleStart;
        maxCycle = experiment1.lCycleEnd;
        minCycleExperiment = experiment1.lExperiment;
        maxCycleExperiment = experiment2.lExperiment;
      } else {
        if (experiment2.lCycleStart<experiment1.lCycleStart) {
          minCycle = experiment2.lCycleStart;
          maxCycle = experiment1.lCycleStart;
        } else {
          minCycle = experiment1.lCycleStart;
          maxCycle = experiment2.lCycleEnd;
        }
        minCycleExperiment = experiment1.lExperiment;
        if (experiment1.lCycleStart==experiment2.lCycleStart) {
          maxCycleExperiment = experiment1.lExperiment;
        } else {
          maxCycleExperiment = experiment2.lExperiment;
        }
      }
    } else {
      // all others cases: (n,1), (n,m)
      minCycle = experiment1.lCycleStart;
      minCycleExperiment = experiment1.lExperiment;
      maxCycle = experiment1.lCycleEnd;
      if (experiment1.bRangleCycle) {
        //TODO: WCH: we have not seen this case, so NEED TO PROVE this!
        maxCycleExperiment = experiment2.lExperiment;
      } else {
        maxCycleExperiment = experiment1.lExperiment;
      }
    }
  }

  // let's look up for the scan number
  if (0!=minCycle) {
    rangeStart = CycleExperimentToScan (minCycle, minCycleExperiment);
    if (-1==rangeStart)
      cerr << "Warning: Fail to locate scan number of cycle " << minCycle << " experiment " << minCycleExperiment << endl;
    rangeEnd = CycleExperimentToScan (maxCycle, maxCycleExperiment);
    if (-1==rangeEnd)
      cerr << "Warning: Fail to locate scan number of cycle " << maxCycle << " experiment " << maxCycleExperiment << endl;
  }
  if (verbose_)
    cout << "\tStart(" << minCycle << "," << minCycleExperiment << ")->Scan(" << rangeStart << ")"
	 << "\tEnd(" << maxCycle << "," << maxCycleExperiment << ")->Scan(" << rangeEnd << ")"
	 << endl;
}

/*
  return the associated scan number of <cycle, experiment> pair
*/
long MascotConverter::CycleExperimentToScan (long lCycle, long lExperiment) const
{
  scanMap::const_iterator iter = cycleExptToScanMap_.find(std::make_pair(lCycle, lExperiment));
  if (cycleExptToScanMap_.end()!=iter) {
    return (*iter).second;
  } else {
    if (2==lExperiment) {
      // hmm... how can this be possible?!
      // fall thru'
    } else {
      scanMap::const_iterator iter = cycleExptToScanMap_.find(std::make_pair(lCycle, lExperiment));
      if (cycleExptToScanMap_.end()!=iter) {
        return ((*iter).second+(lExperiment-2));
      } else {
        // hmm... how can this be possible?!
        // fall thru'
      }
    }
  }
  return -1;
}

/*
  return the associated retention time of <scan>
*/
double MascotConverter::scanToRetentionTime (long lScan) const
{
  scanRTMap::const_iterator iter = scanRetentionTimeMap_.find(lScan);
  if (scanRetentionTimeMap_.end()!=iter) {
    return (*iter).second;
  } else {
    return -1.0;
  }
}

/*
  return the sample enzyme instance
*/
ProteolyticEnzyme* MascotConverter::getSampleEnzyme (bool verbose) const
{
  char sample_enzyme[] = "Trypsin";
  ProteolyticEnzymeFactory factory;
  //ProteolyticEnzyme* enzyme = factory.getProteolyticEnzyme(sample_enzyme);
  //WCH: default to "Trypsin" if no "sample enzyme" has been provided
  //     May be we should default to the enzyme used in the search?
  //ProteolyticEnzyme* enzyme = factory.getProteolyticEnzyme((NULL==sample_enzyme_) ? sample_enzyme : sample_enzyme_);
  const char *finalEnzyme;
  verbose |= (verbose_!=False);
  if (NULL==sample_enzyme_) {
    // user did not provide the sample enzyme, we try to use CLE in .dat file
    if (enzyme_constraint_ && NULL!=enzyme_) {
      // whatever user used in the Mascot search 
      finalEnzyme = enzyme_;
      if (verbose) {
	cerr << "Sample enzyme has not been specified, \"" << enzyme_ << "\" from .dat file assumed" << endl;
      }
    } else {
      //there is no specification or "None" in CLE, we default to "Trypsin"
      finalEnzyme = sample_enzyme;
      if (verbose) {
	cerr << "Sample enzyme has not been specified and could not guess from .dat file, \"" << finalEnzyme << "\" assumed" << endl;
      }
    }
  } else {
    // we use what user provided as the sample enzyme
    finalEnzyme = sample_enzyme_;
  }
  
  ProteolyticEnzyme* enzyme = factory.getProteolyticEnzyme(finalEnzyme);
  if(enzyme == NULL) {
    if (sample_enzyme_) {
      cerr << "error: enzyme " << sample_enzyme_ << " is not recognized" << endl;
    } else {
      cerr << "error: enzyme " << enzyme_ << " in .dat is not recognized, and none is provided by user" << endl;
    }
    exit(1);
  } else {
    enzyme->write(cerr);
  }

  return enzyme;
}

/*
  return the sample enzyme instance
*/
ProteolyticEnzyme* MascotConverter::getSearchEnzyme (bool verbose) const
{
  if(!enzyme_constraint_)
    return NULL;

  ProteolyticEnzymeFactory factory;
  ProteolyticEnzyme* searchEnzyme = factory.getProteolyticEnzyme(enzyme_);
  if(searchEnzyme == NULL) {
    cerr << "error: enzyme in search constraint \"" << enzyme_ << "\" is not recognized" << endl;
    exit(1);
  }
  return searchEnzyme;
}
