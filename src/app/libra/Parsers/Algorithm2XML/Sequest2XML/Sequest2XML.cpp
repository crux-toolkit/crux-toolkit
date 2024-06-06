
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <errno.h>

#include "Parsers/Algorithm2XML/SearchResult/SequestResult.h"
#include "Common/sysdepend.h"
#include "Common/Array.h"
#include "SequestParams.h"
#include "Parsers/mzParser/mzParser.h"
#include "Common/ModificationInfo/ModificationInfo.h"
#include "Common/ModifiedResidueMass/ModifiedResidueMass.h"
#include "Parsers/Algorithm2XML/pICalculator.h"
#include "Common/Enzyme/ProteolyticEnzyme/ProteolyticEnzymeFactory/ProteolyticEnzymeFactory.h"
#include "Parsers/Parser/Parser.h"
#include "Common/util.h"

double isoElectricPt(char* pep);


/*

Program       : DiscriminantCalculator for discr_calc of PeptideProphet 
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


#include "Common/TPPVersion.h" // contains version number, name, revision
#include <ostream>

#include "gzstream.h"

using namespace mzParser;

int main(int argc, char** argv) {
  hooks_tpp handler(argc,argv); // set up install paths etc

  if(argc < 2) {
     cerr << " " << argv[0] << "(" << szTPPVersionInfo << ")" << std::endl;
    cerr << " usage: Sequest2XML htmlfile" << endl;
    cerr << "    or: Sequest2XML htmlfile (prob) (-Eenzyme) (-Pparamsfile) (-Sspectrometer) (-Vschemafile) (-pI) (-all)" << endl;
    exit(1);
  }
  Boolean pre_exist_probs = False;
  char default_enzyme[] = "Trypsin";
  char sample_enzyme[100];
  sample_enzyme[0] = 0; // init
  SequestParams* sequest_params = NULL;
  char massspec[100];
  strcpy(massspec, "LCQ");
  char schema[500];
  schema[0] = 0;
  Boolean write_isoelectric = False;
  Boolean write_all = False;

  Boolean WRITE_MODIFICATION_INFO = True; // explicit

  Boolean MALDI = False;
  int mass_type = 0; // 0 for default, 1 for average, 2 for monoisotopic

  for(int k = 2; k < argc; k++) {
    if(! strcmp(argv[k], "prob"))
      pre_exist_probs = True;
    else if(strlen(argv[k]) > 2 && argv[k][0] == '-' && argv[k][1] == 'E') // enzyme
      strcpy(sample_enzyme, argv[k]+2);
    else if(strlen(argv[k]) > 2 && argv[k][0] == '-' && argv[k][1] == 'P') // parameters
      sequest_params = new SequestParams(argv[k]+2);
    else if(strlen(argv[k]) > 2 && argv[k][0] == '-' && argv[k][1] == 'S') // parameters
      strcpy(massspec, argv[k]+2);
    else if(strlen(argv[k]) > 2 && argv[k][0] == '-' && argv[k][1] == 'V') // schema
      strcpy(schema, argv[k]+2);
    else if(! strcmp(argv[k], "-pI")) // schema
       write_isoelectric = True;
    else if(! strcmp(argv[k], "-M")) // schema
       MALDI = True;
    else if(! strcmp(argv[k], "-m")) // schema
       mass_type = 2;
    else if(! strcmp(argv[k], "-a")) // schema
       mass_type = 1;
    else if(! strcmp(argv[k], "-all")) // schema
       write_all = True;

  }

  if(sequest_params == NULL) { // use default local copy of sequest.params
    sequest_params = new SequestParams("sequest.params");
  }

  if(mass_type) 
    sequest_params->setMassType(mass_type == 2); // override what's in the sequest params file


  // override params masstype if user-specified


  // now check for non default enzyme in parameters
  if(strlen(sample_enzyme) == 0 && sequest_params != NULL 
     && strcmp(sequest_params->getParamsEnzyme(), "Nonspecific"))
    strcpy(sample_enzyme, sequest_params->getParamsEnzyme());

  if(strlen(sample_enzyme) == 0) // default 
     strcpy(sample_enzyme, default_enzyme);

  //cout << "enzyme: " << sample_enzyme << endl; exit(1);

  char* htmlfile = new char[strlen(argv[1])+1];
  strcpy(htmlfile, argv[1]);

  char* outfile = new char[strlen(argv[1])+1];
  // now use ramp's function to extract mass spec info from mzXML file
  int mzXMLfilenamelen;
  char* mzXMLfile = new char[mzXMLfilenamelen=(strlen(argv[1])+12)];
  const char* rawFileExt = "?";


  char htm_suf[] = ".htm";
  char* match = strstr(htmlfile, htm_suf);
  if(match != NULL && 
     (strlen(match) == strlen(htm_suf) || (strlen(match) == strlen(htm_suf)+1 && match[strlen(htm_suf)] == 'l'))) {
    strncpy(outfile, htmlfile, strlen(htmlfile) - strlen(match));
    outfile[strlen(htmlfile) - strlen(match)] = 0;
    strcat(outfile, get_pepxml_dot_ext());
    strncpy(mzXMLfile, htmlfile, strlen(htmlfile) - strlen(match));
    mzXMLfile[strlen(htmlfile) - strlen(match)] = 0;
    rampConstructInputFileName(mzXMLfile,mzXMLfilenamelen,mzXMLfile);
	if (!(rawFileExt=rampValidFileType(mzXMLfile))) { // find .mzxml, .mzxml.gz etc
		rawFileExt = strrchr(mzXMLfile,'.');
	}
  }
  else {
    cerr << "error, cannot parse input htmlfile: " << htmlfile << endl;
    exit(1);
  }

  InstrumentStruct* ppstruct = NULL;
  RAMPFILE* pFI;
  ramp_fileoffset_t *pScanIndex = NULL;
  ramp_fileoffset_t indexOffset;
  int iAnalysisLastScan;

  if ( (pFI=rampOpenFile(mzXMLfile))==NULL) {
    cout << " warning: cannot open \"" << mzXMLfile << "\" for reading MS instrument info and scan times." << endl;
  }
  else {
    ppstruct = getInstrumentStruct(pFI);
    indexOffset = getIndexOffset(pFI);
    pScanIndex = readIndex( pFI , indexOffset, &iAnalysisLastScan );
  }


  char xslfile[] = ""; 


  const int dirsize = 1000;
  char curr_dir[dirsize];



  ProteolyticEnzyme* enzyme = (new ProteolyticEnzymeFactory())->getProteolyticEnzyme(sample_enzyme);
  if(enzyme == NULL) {
    cerr << "error: enzyme " << sample_enzyme << " not recognized" << endl;
    exit(1);
  }

  SequestResult* result = NULL;
  const int datasize = 10000;

  char* data = new char[datasize];
  char* spec = NULL;
  char* database = NULL;

  char peptide[100];
  char prev_aa;
  char next_aa;

  int max_prot_len = 0;
  int max_spec_len = 0;
  int max_pep_len = 0;

  char basename[5000];
  basename[0] = 0;


  Boolean maldi_ = False;
  Boolean exclusion_list = False;

  Array<char*>* excluded_prots = NULL;
  RACI fhtml(htmlfile); // can read gzipped html

  if(! fhtml) {
    cerr << "could not open html file " << htmlfile << endl;
    exit(1);
  }

  time_t now;
  time(&now);

  char* conversion = ctime(&now);
  char* cleaned = new char[strlen(conversion)];
  strncpy(cleaned, conversion, strlen(conversion)-1);
  cleaned[strlen(conversion)-1] = 0;

  ModificationInfo* modinfo = NULL;
  
  pICalculator* pi_calc = NULL;
  if(write_isoelectric)
    pi_calc = new pICalculator();



  FILE* fout;

  unlink(outfile); // delete any previous copy (sometimes needed before fopen can work)
  if( (fout = fopen(outfile, "w")) == NULL) {
    cerr << "cannot open outfile " << outfile << ": " << strerror(errno) << endl;
    exit(1);
  }

  int format = -1;
  int result_index = 1;



  while(fhtml.getline(data, datasize)) {
  // then go through each line, if valid data line, then create SequestResult object and pass on down the line....
    //cout << data << endl;

    result = new SequestResult(data, pre_exist_probs, format);
    if(result == NULL) {
      cerr << "error parsing " << data << endl;
      exit(1);
    }
    if(result->isProcessed()) { // && result->xcorr_ > 0.0) {

      format = result->format_;

      if(maldi_)
	result->charge_ = 1;

	int spec_divisors[] = {-1, -1, -1};
	int index = 2;
	int k = strlen(result->spectrum_)-1;
	  
	while((index+1) && k) {
	  if(result->spectrum_[k] == '.')
	    spec_divisors[index--] = k;
	  k--;
	}
	int maldi_spot = 0;
	while(MALDI && (index < 0) && k) {
	  if(result->spectrum_[k] == '.')
	    maldi_spot = k;
	  k--;
	}

	if(index+1) {
	  cerr << " error: only found " << (2 - index) << " periods in " << result->spectrum_ << endl;
	  exit(1);
	}

	if(database == NULL) {
	  database = result->extractDatabase(data);
	  if(database == NULL) {
	    cerr << "cannot extract database from " << data << endl;
	    exit(1);
	  }
	  

	  if(strlen(basename) == 0) {
	    // first put full path on here
	    safepath_getcwd(curr_dir, dirsize);
	    strcpy(basename, curr_dir);
	    strcat(basename, "/");
	    if(MALDI) {
	      strncat(basename, result->spectrum_, maldi_spot);
	      basename[strlen(curr_dir)+1+maldi_spot] = 0;
	    }
	    else {

	      strncat(basename, result->spectrum_, spec_divisors[0]);
	      basename[strlen(curr_dir)+1+spec_divisors[0]] = 0;
	    }
	    // must make it full path
	  }

	  fprintf(fout, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");


	  
	  fprintf(fout, "<?xml-stylesheet type=\"text/xsl\" href=\"%s%s\"?>\n", PEPXML_NAMESPACE, "_std.xsl");

	  // for now, make summary.xml same as xml
#ifndef USING_RELATIVE_WEBSERVER_PATH 
	  fprintf(fout, "<msms_pipeline_analysis date=\"%s\" xmlns=\"%s\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"%s %s%s\" summary_xml=\"%s\">\n", Parser::getDateTime(), PEPXML_NAMESPACE, PEPXML_NAMESPACE, getPepXML_std_xsl_web_path(), PEPXML_SCHEMA, outfile);
#else
	  char* szWebserverRoot = new char[256];
	  char* szCommand = new char[256];
	  char* szBuf = new char[SIZE_BUF];
	  const char *pStr=getWebserverRoot();
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
	 cygwinify(szWebserverRoot,256); // no effect in non-cygwin builds
	  } 
	  fprintf(fout, "<msms_pipeline_analysis date=\"%s\" xmlns=\"%s\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"%s %s%s\" summary_xml=\"%s\">\n", Parser::getDateTime(), PEPXML_NAMESPACE, PEPXML_NAMESPACE, szWebserverRoot, PEPXML_SCHEMA, outfile);
	  delete [] szWebserverRoot;
	  delete [] szBuf;  
	  delete [] szCommand;
#endif
	  if(ppstruct != NULL) 
	    fprintf(fout, "<msms_run_summary base_name=\"%s\" msManufacturer=\"%s\" msModel=\"%s\" msIonization=\"%s\" msMassAnalyzer=\"%s\" msDetector=\"%s\" raw_data_type=\"%s\" raw_data=\"%s\">\n", basename, ppstruct->manufacturer, ppstruct->model, ppstruct->ionisation, ppstruct->analyzer, ppstruct->detector, "raw", rawFileExt);
	  else // leave out mass spec info
	    fprintf(fout, "<msms_run_summary base_name=\"%s\" raw_data_type=\"%s\" raw_data=\"%s\">\n", basename, "raw", rawFileExt);

	  enzyme->writeTraditionalPepXMLTags(fout);

	  if(sequest_params != NULL) {
	    Array<Tag*>* tags = sequest_params->getSearchParamTags(basename, "SEQUEST", database);
	    for(int k = 0; k < tags->length(); k++) 
	      if((*tags)[k] != NULL) {
		(*tags)[k]->writeTraditional(fout);
		delete (*tags)[k];
	      }
	    delete tags;
	  }


	}

	peptide[0] = 0;
	strncat(peptide, result->peptide_+2, strlen(result->peptide_)-4);
	peptide[strlen(result->peptide_)-3] = 0;
	prev_aa = result->peptide_[0];

	// make any substitutions to symbols here...
	sequest_params->modifyAminoacidModifications(peptide);
	sequest_params->modifyTerminalModifications(peptide);

	next_aa = result->peptide_[strlen(result->peptide_)-1];


	char* next_str = enzyme->strip(result->peptide_, 1); //result->stripped();
	
	if (!write_all && strchr(next_str, 'X') != NULL) {
	  cerr << "WARNING: Peptide " << next_str << " contains an X, skipping (use -all option to override) ..." << endl;
	  if (next_str != NULL) {
	    delete next_str;
	  }
	}
	else {


	  if((int) strlen(result->spectrum_) > max_spec_len)
	    max_spec_len = strlen(result->spectrum_);
	  
	  if((int) strlen(result->protein_) > max_prot_len)
	    max_prot_len = strlen(result->protein_);
	  
	  if((int) strlen(result->peptide_) > max_pep_len)
	    max_pep_len = strlen(result->peptide_);
	  
	  fprintf(fout, "    <spectrum_query spectrum=\"%s\"", result->spectrum_);
	  
	  fprintf(fout, " start_scan=\"");
	  Boolean scan_start = False;
	  for(k = spec_divisors[0]+1; k < spec_divisors[1]; k++) {
	    if(! scan_start && result->spectrum_[k] != '0')
	      scan_start = True;
	    if(scan_start)
	      fprintf(fout, "%c", result->spectrum_[k]);
	  }
	  scan_start = False;
	  fprintf(fout, "\" end_scan=\"");
	  for(k = spec_divisors[1]+1; k < spec_divisors[2]; k++) {
	    if(! scan_start && result->spectrum_[k] != '0')
	      scan_start = True;
	    if(scan_start)
	      fprintf(fout, "%c", result->spectrum_[k]);
	  }
	  
	  fprintf(fout, "\" precursor_neutral_mass=\"%0.4f\" assumed_charge=\"%d\" index=\"%d\"",
		  result->neutral_mass_+result->massdiff_, result->charge_, result_index++);
	  
	  /* add RT from mzXML file */
	  if (pFI != NULL)
	    {
	      int iStartScan = -1;
	      char *pStr;
	      char szBaseName[256];
	      
	      strcpy(szBaseName, result->spectrum_);
	      
	      /* read scan number encoded in .dta name */
	      if ((pStr=strrchr(szBaseName, '.'))!=NULL)
		{
		  *pStr = '\0';
		  if ((pStr=strrchr(szBaseName, '.'))!=NULL)
		    if ((pStr=strrchr(szBaseName, '.'))!=NULL)
		      sscanf(pStr+1, "%d", &iStartScan);
		}
	      
	      if(iStartScan != -1 && iStartScan <= iAnalysisLastScan)
		{
		  struct ScanHeaderStruct scanHeader;
		  readHeader(pFI, pScanIndex[iStartScan],&scanHeader);
		  fprintf(fout, " retention_time_sec=\"%0.1f\"", scanHeader.retentionTime);
		}
	    }
	  fprintf(fout, ">\n");
	  
	  char matched[4];
	  char totions[4];
	  sprintf(matched, "%d", result->num_matched_ions_);
	  sprintf(totions, "%d", result->tot_num_ions_);
	  
	  fprintf(fout, "    <search_result>\n");
	  
	  if(WRITE_MODIFICATION_INFO) {
	    
	    double pepmass = result->neutral_mass_;
	    
	    // modification info here....
	    Array<Tag*>* modinfo_tags = sequest_params->getModificationInfoTags(peptide);
	    if(modinfo_tags != NULL) {
	      modinfo = new ModificationInfo(modinfo_tags);
	    }
	    
	    if(mass_type > 0) 
	      pepmass = ModifiedResidueMass::getModifiedPeptideMass(modinfo, next_str, mass_type == 2); // monoisotopic
	    
	    if(result->massdiff_ < 0.0)
			fprintf(fout, "      <search_hit hit_rank=\"%d\" peptide=\"%s\" peptide_prev_aa=\"%c\" peptide_next_aa=\"%c\" protein=\"%s\" num_tot_proteins=\"%d\" num_matched_ions=\"%s\" tot_num_ions=\"%s\" calc_neutral_pep_mass=\"%0.4f\" massdiff=\"%f\" num_tol_term=\"%d\" num_missed_cleavages=\"%d\"",  1, next_str, prev_aa, next_aa, XMLEscape(std::string(result->protein_)).c_str(), result->degen_+1, matched, totions, pepmass, result->massdiff_, enzyme->getNumTolTerm(prev_aa, next_str, next_aa), enzyme->getNumMissedCleavages(result->peptide_));
	    else 
	      fprintf(fout, "      <search_hit hit_rank=\"%d\" peptide=\"%s\" peptide_prev_aa=\"%c\" peptide_next_aa=\"%c\" protein=\"%s\" num_tot_proteins=\"%d\" num_matched_ions=\"%s\" tot_num_ions=\"%s\" calc_neutral_pep_mass=\"%0.4f\" massdiff=\"+%f\" num_tol_term=\"%d\" num_missed_cleavages=\"%d\"",  1, next_str, prev_aa, next_aa, XMLEscape(std::string(result->protein_)).c_str(), result->degen_+1, matched, totions, pepmass, result->massdiff_, enzyme->getNumTolTerm(prev_aa, next_str, next_aa), enzyme->getNumMissedCleavages(result->peptide_));
	    
	    
	    
	    if(write_isoelectric){
	      if(WRITE_MODIFICATION_INFO && pi_calc != NULL){
		fprintf(fout, " calc_pI=\"%0.2f\"", pi_calc->Peptide_pI(peptide, modinfo));
	      }else{
		fprintf(fout, " calc_pI=\"%0.2f\"", isoElectricPt(peptide));
              }
	    }
	    fprintf(fout, " is_rejected=\"%d\">\n", 0);
	    
	    if(modinfo_tags != NULL) {
	      for(int k = 0; k < modinfo_tags->length(); k++) 
		if((*modinfo_tags)[k] != NULL) {
		  fprintf(fout, "        ");
		  (*modinfo_tags)[k]->writeTraditional(fout);
		  delete (*modinfo_tags)[k];
		}
	      delete modinfo_tags;
	    }
	    if(modinfo != NULL)
	      delete modinfo;
	    
	  } // if write mod info
	  else {
	    if(result->massdiff_ < 0.0)
	      fprintf(fout, "      <search_hit hit_rank=\"%d\" peptide=\"%s\" peptide_prev_aa=\"%c\" peptide_next_aa=\"%c\" protein=\"%s\" num_tot_proteins=\"%d\" num_matched_ions=\"%s\" tot_num_ions=\"%s\" calc_neutral_pep_mass=\"%0.1f\" massdiff=\"%f\" num_tol_term=\"%d\" num_missed_cleavages=\"%d\" stripped_peptide=\"%s\"",  1, peptide, prev_aa, next_aa, XMLEscape(std::string(result->protein_)).c_str(), result->degen_+1, matched, totions, result->neutral_mass_, result->massdiff_, enzyme->getNumTolTerm(prev_aa, peptide, next_aa), enzyme->getNumMissedCleavages(result->peptide_), next_str);
	    else 
	      fprintf(fout, "      <search_hit hit_rank=\"%d\" peptide=\"%s\" peptide_prev_aa=\"%c\" peptide_next_aa=\"%c\" protein=\"%s\" num_tot_proteins=\"%d\" num_matched_ions=\"%s\" tot_num_ions=\"%s\" calc_neutral_pep_mass=\"%0.1f\" massdiff=\"+%f\" num_tol_term=\"%d\" num_missed_cleavages=\"%d\" stripped_peptide=\"%s\"",  1, peptide, prev_aa, next_aa, XMLEscape(std::string(result->protein_)).c_str(), result->degen_+1, matched, totions, result->neutral_mass_, result->massdiff_, enzyme->getNumTolTerm(prev_aa, peptide, next_aa), enzyme->getNumMissedCleavages(result->peptide_), next_str);
	    if(write_isoelectric)
	      fprintf(fout, " isoelectric_point=\"%0.1f\"", isoElectricPt(peptide));
	    fprintf(fout, " is_rejected=\"%d\">\n", 0);
	  }
	  
	  if(next_str != NULL)
	    delete next_str;
	  
	  fprintf(fout, "        <search_score name=\"xcorr\" value=\"%0.3f\"/>\n", result->xcorr_);
	  fprintf(fout, "        <search_score name=\"deltacn\" value=\"%0.3f\"/>\n", result->delta_);
	  int delta_star = result->deltastar_ ? 1 : 0;
	  fprintf(fout, "        <search_score name=\"deltacnstar\" value=\"%d\"/>\n", delta_star);
	  fprintf(fout, "        <search_score name=\"spscore\" value=\"%0.1f\"/>\n", result->sp_score_);
	  fprintf(fout, "        <search_score name=\"sprank\" value=\"%d\"/>\n", result->rank_);
	  
	  
	  if(pre_exist_probs) {
	    fprintf(fout, "        <peptideprophet_result probability=\"%0.4f\"/>\n", result->probability_);
	  }
	  
	  
	  fprintf(fout, "      </search_hit>\n");
	  fprintf(fout, "    </search_result>\n");
	  fprintf(fout, "    </spectrum_query>\n");
	}

	if(! maldi_ && spec != NULL)
	  delete spec;

    } // if real data line

    if(result != NULL)
      delete result;

  } // next html file line

  if(exclusion_list) {
    for(int k = 0; k < excluded_prots->length(); k++)
      if((*excluded_prots)[k] != NULL)
	delete (*excluded_prots)[k];
    delete excluded_prots;
  }

  if (database != NULL) {
    fprintf(fout, "  </msms_run_summary>\n");
    fprintf(fout, "</msms_pipeline_analysis>\n");
  }

  fhtml.close();
  fclose(fout);

  if (pFI != NULL)
     rampCloseFile(pFI);
     free(pScanIndex);

  return 0;
}

// to be implemented
double isoElectricPt(char* pep) {
  return 0.0;
}
