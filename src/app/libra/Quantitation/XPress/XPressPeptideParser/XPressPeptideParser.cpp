/*
Program       : XPressPeptideParser
Author        : J.Eng and Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN Info      : $Id: XPressPeptideParser.cpp 8870 2023-02-25 06:40:57Z real_procopio $

Additional work for mzData handling Copyright (C) Brian Pratt Insilicos LLC 2005

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

#include "XPressPeptideParser.h"
#include "Common/TPPVersion.h"         // contains version number, name, revision
#include "Parsers/Parser/TagListComparator.h"  // regression test stuff - bpratt Insilicos LLC, Nov 2005
#include "XPressPeptideParserTagListComparator.h" // regression test stuff - bpratt Insilicos LLC, Nov 2005

#include <string>
#define CACHED_RAMP_HOME
#include "Parsers/mzParser/cached_ramp.h"

#define MAX_MASS_FUDGE 0.2
#define SCAN_BUFFER 300        // buffer for # of scans to read before/after CID scan for filtering & finding endpoints

using namespace mzParser;

XPressPeptideParser::XPressPeptideParser(const char *xmlfile, const InputStruct & options, const char *testMode):Parser("xpress")
{
  pInput_ = options;
  m_XMLfile_state = -1; // no mzxml file yet
  index_ = NULL;
  dTmpFilter_ = NULL;
  dLightMS_ = NULL;
  dHeavyMS_ = NULL;
  dLightFilteredMS_ = NULL;
  dHeavyFilteredMS_ = NULL;
  for (int n=MAX_ISOTOPES+1; n--;) {
     dIsotopeLight[n] = NULL;
     dIsotopeHeavy[n] = NULL;
  }
  fp_ = NULL;
  
#ifdef USE_STD_MODS
  modinfo_ = NULL;
#endif
  
  testMode_ = testMode ? strdup(testMode) : NULL;      // regression test stuff - bpratt Insilicos LLC, Nov 2005
  
  init(xmlfile);   // init from Parsers/Parser.cxx will call parse()
}


XPressPeptideParser::~XPressPeptideParser()
{
  free(testMode_);
}


void XPressPeptideParser::parse(const char *xmlfile)
{
    readAllModMasses(xmlfile);
    char  *engine = NULL;
    char  *enzyme = NULL;
    char  *massspec = NULL;
    Array < Tag * >*tags = NULL;
    Tag   *tag = NULL;

    char  *data = NULL;

    double MIN_PROB = modelOpts_.minprob_;        //0.0; //0.05; // for now

    char text3[100];
    char szAnalysisName[100];

    if (pInput_.bLabelFreeMode)
       strcpy(szAnalysisName, "xpresslabelfree");
    else
       strcpy(szAnalysisName, getName());

    //
    // regression test stuff - bpratt Insilicos LLC, Nov 2005
    //
    Array < Tag * >test_tags;
    eTagListFilePurpose testType;
    char  *testFileName = NULL;
    checkRegressionTestArgs(testMode_, testType);
    if (testType!=NO_TEST)
    {
       std::string options;
       testFileName = constructTagListFilename(xmlfile,  // input file
             testMode_,        // program args
             "XPressPeptideParser",    // program name
             testType);        // user info output
    }

#define RECORD(tag) {(tag)->write(fout);if (testType!=NO_TEST) {test_tags.insertAtEnd(new Tag(*(tag)));}}

    Tag   *timestamp_start = new Tag("analysis_timestamp", True, False);

    timestamp_start->setAttributeValue("analysis", szAnalysisName);

    timestamp_start->setAttributeValue("time", time_);
    timestamp_start->setAttributeValue("id", "1");
    Tag   *timestamp_stop = new Tag("analysis_timestamp", False, True);
    Tag   *timestamp;

    Tag   *result_start = new Tag("analysis_result", True, False);

    if (pInput_.bLabelFreeMode)
       timestamp = new Tag("xpresslabelfree_timestamp", True, True);
    else
    {
       timestamp = new Tag("xpressratio_timestamp", True, True);
       sprintf(text3, "%d", pInput_.bXpressLight1);
       timestamp->setAttributeValue("xpress_light", text3);
    }
    result_start->setAttributeValue("analysis", szAnalysisName);

    Tag   *analysis_start = new Tag("analysis_summary", True, False);
    analysis_start->setAttributeValue("analysis", szAnalysisName);
    analysis_start->setAttributeValue("time", time_);
    Tag   *analysis_stop = new Tag("analysis_summary", False, True);
    Tag   *summary = getSummaryTag(pInput_);

    Tag   *result_stop = new Tag("analysis_result", False, True);
    TagFilter *xpress_filter = new TagFilter("analysis_result");
    xpress_filter->enterRequiredAttributeVal("analysis", szAnalysisName);
    TagFilter *xpress_summ_filter = new TagFilter("analysis_timestamp");
    xpress_summ_filter->enterRequiredAttributeVal("analysis", szAnalysisName);
    TagFilter *xpress_summ = new TagFilter("analysis_summary");
    xpress_summ->enterRequiredAttributeVal("analysis", szAnalysisName);

    // construct a tmpfile name based on xmlfile
    std::string outfile = make_tmpfile_name(xmlfile);
    ofstream fout(outfile.c_str());
    if (!fout)
    {
       cerr << "cannot write output to file " << outfile << endl;
       exit(1);
    }

    Boolean first = False;
    Tag   *xpress_ratio = NULL;
    int ratio_tags_written = 0;

    Boolean collected = False;

#ifdef USE_STD_MODS
    monoisotopic_ = False;       // unless proven otherwise
    Boolean mod_on = False;
    Array < Tag * >*modifications = NULL;
#endif
    Boolean top_hit = False;
    double dPeptideProb=1.0;
    
    RACI fin(xmlfile); // can read gzipped xml
    if (!fin)
    {
       cerr << "XPressPeptideParser: error opening " << xmlfile << endl;
       exit(1);
    }

    char  *nextline = new char[line_width_];
    while (fin.getline(nextline, line_width_))
    {
      //cout << "next: " << nextline << endl;

       data = strstr(nextline, "<");
       while (data != NULL)
       {
          tag = new Tag(data);

          collected = False;

          setFilter(tag);

         if ((!xpress_filter->filter(tag) && !xpress_summ_filter->filter(tag) && !xpress_summ->filter(tag)))
         {
            if (tag->isStart() && !strcmp(tag->getName(), "msms_pipeline_analysis"))
            {
               RECORD(tag);
               RECORD(analysis_start);
               RECORD(summary);
               RECORD(analysis_stop);
               delete analysis_start;
               analysis_start = NULL;
               delete summary;
               summary = NULL;
               delete analysis_stop;
               analysis_stop = NULL;
            }
            else if (tag->isStart() && !strcmp(tag->getName(), "msms_run_summary"))
            {
               rampConstructInputPath(mzXMLfile_, sizeof(mzXMLfile_),
                     pInput_.szMzXMLDir, tag->getAttributeValue("base_name"));

               if ((fp_ = cached_ramp_rampOpenFile(mzXMLfile_)) == NULL)
               {
                  printf("XPRESS error - cannot open file from basename %s, will try to derive from scan names\n",
                        tag->getAttributeValue("base_name"));
               } else {
                  // set index here.....
                  index_ = cached_ramp_readIndex(fp_,  cached_ramp_getIndexOffset(fp_), &(pInput_.iAnalysisLastScan));
               }

               const struct ScanHeaderStruct* pHeader = cached_ramp_readHeader(fp_, index_[pInput_.iAnalysisLastScan]);

               // iAnalysisLastScan can either be the real last scan or some crazy native ID number.
               // In order to deal with this, I'm going to map everything back to a sequential scans.
               // piSequentialScan[]:  index is each sequential scan number; value in array is native ID
               // piReverseSequentialScan[]:  index is nativeID; value in array is sequential scan number
               pInput_.iScanCount=1;
               for (int x=0; x<=pInput_.iAnalysisLastScan; x++)
                  if (index_[x] != -1)
                     pInput_.iScanCount++;

               if ( (piSequentialScan = (int *)malloc(pInput_.iScanCount * sizeof(int)))==NULL)
               {
                  printf(" Error cannot malloc piSequentialScan[%d]\n\n", pInput_.iScanCount);
                  exit(1);
               }

               int iLastScanPlus1 = pInput_.iAnalysisLastScan+1;

               if ( (piReverseSequentialScan = (int *)malloc(iLastScanPlus1 * sizeof(int)))==NULL)
               {
                  printf(" Error cannot malloc piReverseSequentialScan[%d]\n\n", iLastScanPlus1);
                  exit(1);
               }

	       piSequentialScan[0] = 0; // fallback, when no signal is found in scan range
               int iTmp=1;
               for (int x=0; x<=pInput_.iAnalysisLastScan; x++)
               {
                  piReverseSequentialScan[x] = -1;
                  if (index_[x] != -1)
                  {
                     piReverseSequentialScan[x] = iTmp;
                     piSequentialScan[iTmp++] = x;
                  }
               }

               pInput_.iAnalysisFirstScan = 1;
               pInput_.iAnalysisLastScan = pInput_.iScanCount - 1;

//             for (int x=1; x<=iScanCount; x++)
//                printf("sequential scan#=%d, nativId=%d\n", x, piSequentialScan[x]);

               if (dTmpFilter_ != NULL)
                  delete  dTmpFilter_;
               if ( (dTmpFilter_ = new Array<double> (pInput_.iScanCount) )==NULL)
               {
                  printf(" Error - cannot malloc dTmpFilter[%d]\n", pInput_.iScanCount);
                  exit(1);
               }

               if (dLightMS_ != NULL)
                  delete [] dLightMS_;
               if ( (dLightMS_ = new double [pInput_.iScanCount] )==NULL)
               {
                  printf(" Error - cannot malloc dLightMS[%d]\n", pInput_.iScanCount);
                  exit(1);
               }

               if (dHeavyMS_ != NULL)
                 delete [] dHeavyMS_;
               if ( (dHeavyMS_ = new double [pInput_.iScanCount] )==NULL)
               {
                  printf(" Error - cannot malloc dHeavyMS[%d]\n", pInput_.iScanCount);
                  exit(1);
               }

               if (dLightFilteredMS_ != NULL)
                  delete dLightFilteredMS_;
               if ( (dLightFilteredMS_ = new Array<double> (pInput_.iScanCount) )==NULL)
               {
                  printf(" Error - cannot malloc dHeavyMS[%d]\n", pInput_.iScanCount);
                  exit(1);
               }

               if (dHeavyFilteredMS_ != NULL)
                  delete dHeavyFilteredMS_;
               if ( (dHeavyFilteredMS_ = new Array<double> (pInput_.iScanCount) )==NULL)
               {
                  printf(" Error - cannot malloc dHeavyMS[%d]\n", pInput_.iScanCount);
                  exit(1);
               } 

               dLightFilteredMS_->nullify();
               dHeavyFilteredMS_->nullify();
               dTmpFilter_->nullify();

               for (int n=0 ; n<=MAX_ISOTOPES; n++)
               {
                  if (dIsotopeLight[n] != NULL)
                     delete [] dIsotopeLight[n];
                  if ((dIsotopeLight[n] = new double [pInput_.iScanCount])==NULL)
                  {
                     printf("Error malloc dIsotopeLight[%d][%d]\n", n , pInput_.iScanCount);
                     exit(1);
                  }

                  if (dIsotopeHeavy[n] != NULL)
                     delete [] dIsotopeHeavy[n];
                  if ((dIsotopeHeavy[n] = new double [pInput_.iScanCount])==NULL)
                  {
                     printf("Error malloc dIsotopeHeavy[%d][%d]\n", n , pInput_.iScanCount);
                     exit(1);
                  }
               }
               
               first = True;
               RECORD(tag);
            }                   // msms_summ
            else if (tag->isEnd() && !strcmp(tag->getName(), "search_summary"))
            {
               RECORD(tag);
               RECORD(timestamp_start);
               RECORD(timestamp);
               RECORD(timestamp_stop);
            }
            else if (tag->isStart() && !strcmp(tag->getName(), "search_summary"))
            {
               monoisotopic_ = !strcmp(tag->getAttributeValue("precursor_mass_type"), "monoisotopic");
               RECORD(tag);
            }
            // have a modification worth recording here
            else if (tag->isStart()
                  && !strcmp(tag->getName(), "terminal_modification")
                  && strchr(pInput_.szXpressResidues, tag->getAttributeValue("terminus")[0]) != NULL)
            {
               RECORD(tag);
            }
            else if (tag->isStart()
                  && !strcmp(tag->getName(), "aminoacid_modification")
                  && strchr(pInput_.szXpressResidues, tag->getAttributeValue("aminoacid")[0]) != NULL)
            {
                RECORD(tag);
            }
            else if (filter_)
            {
               if (tag->isStart() && !strcmp("spectrum_query", tag->getName()))
               {
                  pInput_.iFirstScan = atoi(tag->getAttributeValue("start_scan"));
                  pInput_.iLastScan = atoi(tag->getAttributeValue("end_scan"));
                  pInput_.iChargeState = atoi(tag->getAttributeValue("assumed_charge"));
                  pInput_.dPeptideMass = (double) (atof (tag->getAttributeValue("precursor_neutral_mass")) + PROTON_MASS);
               }
               else if (tag->isStart() && !strcmp("search_hit", tag->getName()) && !strcmp("1", tag->getAttributeValue("hit_rank")))
               {
                  strcpy(pInput_.szPeptide, tag->getAttributeValue("peptide"));
                  pInput_.dPeptideMass = (double) (atof (tag-> getAttributeValue("calc_neutral_pep_mass")) + PROTON_MASS);
                  top_hit = True;
                  pInput_.iModCount = 0;
               }
               else if (tag->isStart() && !strcmp("search_hit", tag->getName()) && strcmp("1", tag->getAttributeValue("hit_rank")))
               {
                  top_hit = False;
               }
               else if (tag->isStart() && !strcmp("peptideprophet_result", tag->getName()))
               {
                  dPeptideProb = (double) (atof (tag-> getAttributeValue("probability")));
               }
               else if (pInput_.iMetabolicLabeling && top_hit && tag->isStart() && !strcmp("mod_aminoacid_mass", tag->getName()))
               {
                  // This just counts # of mods; for metabolic quant, should equal pep length
                  // as every residue is statically modified.  So mod count == sequence length
                  // is a quick/dirty way of determing a heavy peptide.
                  pInput_.iModCount += 1;
               }
#ifdef USE_STD_MODS
               else if (top_hit && tag->isStart() && !strcmp("modification_info", tag->getName()))
               {
                  if (modifications == NULL)
                     modifications = new Array < Tag * >;
                  modifications->insertAtEnd(tag);
                  mod_on = !tag->isEnd();
                  if (!mod_on) { // tag already closed, process it now
                     modinfo_ = new ModificationInfo(modifications);
                  }
               }
               else if (top_hit && mod_on && tag->isEnd() && !strcmp("modification_info", tag->getName()))
               {
                  modifications->insertAtEnd(tag);
                  modinfo_ = new ModificationInfo(modifications);
                  mod_on = False;
               }
               else if (top_hit && mod_on)
               {
                  modifications->insertAtEnd(tag);
               }
#endif

               if (tags == NULL)
                  tags = new Array < Tag * >;
               tags->insertAtEnd(tag);
               collected = True;
            }
            else {
	      if (tag->isEnd() && !strcmp(tag->getName(), "msms_run_summary")) {
		if (fp_ != NULL)
		  cached_ramp_rampCloseFile(fp_);

		if (index_ != NULL) {
		  free(index_);
		  index_ = NULL;
		}

		if (piSequentialScan != NULL) {
		  free(piSequentialScan);
		  piSequentialScan = NULL;
		}
		if (piReverseSequentialScan != NULL) {
		  free(piReverseSequentialScan);
		  piReverseSequentialScan = NULL;
		}
	      }

	      if (tag != NULL)
		RECORD(tag);
            }

            if (filter_memory_) {                   // process
	      if (tags != NULL) {
		if (pInput_.iChargeState >= 0) {
		  xpress_ratio = NULL;

		  if (dPeptideProb >= pInput_.dMinPprob)
		    xpress_ratio = getRatio();
		}

		for (int k = 0; k < tags->length(); k++) {
		  if ((*tags)[k] != NULL) {
		    if (!xpress_filter->filter((*tags)[k])) {
		      // here check for correct time to write xpress tag
		      if ((*tags)[k]->isEnd() && !strcmp((*tags)[k]->getName(), "search_hit") && xpress_ratio != NULL) {
			ratio_tags_written++;
			if (ratio_tags_written % 1000 == 0)
			  cout << " " << ratio_tags_written/1000 << "k" << endl;
			else if (ratio_tags_written % 100 == 0)
			  cout << ":";
			else if (ratio_tags_written % 10 == 0) {
			  cout << ".";
			  fflush(stdout);
			}

			RECORD(result_start);
			RECORD(xpress_ratio);
			RECORD(result_stop);
			delete xpress_ratio;
			xpress_ratio = NULL;
		      }
		      RECORD((*tags)[k]);
		    }
		    delete(*tags)[k];
		  }
		}             // next tag
		delete tags;
		tags = NULL;
#ifdef USE_STD_MODS
		if (modifications != NULL)
                  {
		    delete modifications;
		    modifications = NULL;
                  }
#endif
	      }

	      pInput_.iChargeState = -1;       // reset

#ifdef USE_STD_MODS
	      if (modinfo_ != NULL)
		delete modinfo_;
	      modinfo_ = NULL;
#endif

	      //xpress_ratio = NULL;
            }
         }                      // if not filtered
         else
         {

         }                      // filtered
         if (!collected && tag != NULL)
	   delete tag;
         data = strstr(data + 1, "<");
       }                         // next tag

    }                            // next line
    fin.close();
    fout.close();

    if (dTmpFilter_ != NULL)
      delete dTmpFilter_;
    if (dLightMS_ != NULL)
      delete [] dLightMS_;
    if (dHeavyMS_ != NULL)
      delete [] dHeavyMS_;
    if (dLightFilteredMS_ != NULL)
      delete dLightFilteredMS_;
    if (dHeavyFilteredMS_ != NULL)
      delete dHeavyFilteredMS_;

    dTmpFilter_ = NULL;
    dLightMS_ = NULL;
    dHeavyMS_ = NULL;
    dLightFilteredMS_ = NULL;
    dHeavyFilteredMS_ = NULL;

    for (int n=0 ; n<=MAX_ISOTOPES; n++) {
      if (dIsotopeLight[n] != NULL)
	delete [] dIsotopeLight[n] ;
      if (dIsotopeHeavy[n] != NULL)
	delete [] dIsotopeHeavy[n] ;
    }

    if (!overwrite(xmlfile, outfile.c_str(), "</msms_pipeline_analysis>"))
      cerr << "error: no xpress data written to file " << xmlfile << endl;

    if (testType!=NO_TEST)
      {
	//
	// regression test stuff - bpratt Insilicos LLC, Nov 2005
	//
	TagListComparator("XPressPeptideParser", testType, test_tags, testFileName);
	delete [] testFileName;
	for (int k = test_tags.length(); k--;)
	  {
	    delete test_tags[k];
	  }
      }

    delete [] nextline;

    delete timestamp_start;
    delete timestamp_stop;
    delete timestamp;

    if (analysis_start != NULL)
      delete analysis_start;

    if (analysis_stop != NULL)
      delete analysis_stop;

    if (summary != NULL)
      delete summary;

    delete result_start;
    delete result_stop;
    delete xpress_filter;
    delete xpress_summ_filter;
    delete xpress_summ;

    cout << " Total: " << ratio_tags_written << endl;
}

void XPressPeptideParser::readAllModMasses(const char *xmlfile)
{
   char  *engine = NULL;
   char  *enzyme = NULL;
   char  *massspec = NULL;
   Array < Tag * >*tags = NULL;
   Tag   *tag = NULL;

   char  *data = NULL;

   double MIN_PROB = modelOpts_.minprob_;        //0.0; //0.05; // for now

   char text3[100];
   char szAnalysisName[100];

   if (pInput_.bLabelFreeMode)
      strcpy(szAnalysisName, "xpresslabelfree");
   else
      strcpy(szAnalysisName, getName());

   Tag   *timestamp_start = new Tag("analysis_timestamp", True, False);
   timestamp_start->setAttributeValue("analysis", szAnalysisName);
   timestamp_start->setAttributeValue("time", time_);
   timestamp_start->setAttributeValue("id", "1");
   Tag   *timestamp_stop = new Tag("analysis_timestamp", False, True);
   Tag   *timestamp;

   if (pInput_.bLabelFreeMode)
   {
      timestamp = new Tag("xpresslabelfree_timestamp", True, True);
   }
   else
   {
      timestamp = new Tag("xpressratio_timestamp", True, True);
      sprintf(text3, "%d", pInput_.bXpressLight1);
      timestamp->setAttributeValue("xpress_light", text3);
   }

   Tag   *analysis_start = new Tag("analysis_summary", True, False);
   analysis_start->setAttributeValue("analysis", szAnalysisName);
   analysis_start->setAttributeValue("time", time_);
   Tag   *analysis_stop = new Tag("analysis_summary", False, True);

   Tag   *summary = getSummaryTag(pInput_);

   Tag   *result_start = new Tag("analysis_result", True, False);
   result_start->setAttributeValue("analysis", szAnalysisName);

   Tag   *result_stop = new Tag("analysis_result", False, True);

   TagFilter *xpress_filter = new TagFilter("analysis_result");
   xpress_filter->enterRequiredAttributeVal("analysis", szAnalysisName);
   TagFilter *xpress_summ_filter = new TagFilter("analysis_timestamp");
   xpress_summ_filter->enterRequiredAttributeVal("analysis", szAnalysisName);
   TagFilter *xpress_summ = new TagFilter("analysis_summary");
   xpress_summ->enterRequiredAttributeVal("analysis", szAnalysisName);

   Boolean first = False;
   Tag   *xpress_ratio = NULL;
   int ratio_tags_written = 0;

   Boolean collected = False;

#ifdef USE_STD_MODS
   monoisotopic_ = False;       // unless proven otherwise
   Boolean mod_on = False;
   Array < Tag * >*modifications = NULL;
#endif
   Boolean top_hit = False;
    
   RACI fin(xmlfile); // can read gzipped xml
   if (!fin)
   {
      cerr << "error opening " << xmlfile << endl;
      exit(1);
   }
   char  *nextline = new char[line_width_];
   while (fin.getline(nextline, line_width_))
   {
      //cout << "next: " << nextline << endl;

      data = strstr(nextline, "<");
      while (data != NULL)
      {
         tag = new Tag(data);

         collected = False; 
         
         setFilter(tag);
       
         if ((!xpress_filter->filter(tag) && !xpress_summ_filter->filter(tag) && !xpress_summ->filter(tag)))
         {
            if (tag->isStart() && !strcmp(tag->getName(), "msms_pipeline_analysis"))
            {
               delete analysis_start;
               analysis_start = NULL;
               delete summary;
               summary = NULL;
               delete analysis_stop;
               analysis_stop = NULL;

#ifdef USE_STD_MODS
               for (int k = 0; k < 26; k++)
               {
                  light_label_masses_[k] = 0.0;
                  heavy_label_masses_[k] = 0.0;
               }
               heavy_nterm_mass_ = 0.0;
               light_nterm_mass_ = 0.0;
               heavy_cterm_mass_ = 0.0;
               light_cterm_mass_ = 0.0;
#endif
            }
            else if (tag->isStart() && !strcmp(tag->getName(), "msms_run_summary"))
            {
               pInput_.iAnalysisFirstScan = 1;
               first = True;
            }                   // msms_summ
            else if (tag->isEnd() && !strcmp(tag->getName(), "search_summary"))
            {
#ifdef USE_STD_MODS
               // check label tags here
               for (int k = 0; pInput_.szXpressResidues[k]; k++)
               {
                  if (pInput_.szXpressResidues[k] == 'n')
                  {
                     if (heavy_nterm_mass_ > 0.0 && light_nterm_mass_ == 0.0)
                     {
                        if (pInput_.dXpressMassDiff1 != 0.0 && strchr(pInput_.szXpressResidues1, 'n') != 0)
                           light_nterm_mass_ = heavy_nterm_mass_ - pInput_.dXpressMassDiff1;
                        else if (pInput_.dXpressMassDiff2 != 0.0 && strchr(pInput_.szXpressResidues2, 'n') != 0)
                           light_nterm_mass_ = heavy_nterm_mass_ - pInput_.dXpressMassDiff2;
                        else if (pInput_.dXpressMassDiff3 != 0.0 && strchr(pInput_.szXpressResidues3, 'n') != 0)
                           light_nterm_mass_ = heavy_nterm_mass_ - pInput_.dXpressMassDiff3;
                        else
                           light_nterm_mass_ = ResidueMass::getMass('n', monoisotopic_);
                     }
                     else if (heavy_nterm_mass_ == 0.0 && light_nterm_mass_ > 0.0)
                     {          
                        //must switch them around
                        if (pInput_.dXpressMassDiff1 > 0.0 && strchr(pInput_.szXpressResidues1, 'n') != 0)
                           heavy_nterm_mass_ = light_nterm_mass_ + pInput_.dXpressMassDiff1;
                        else if (pInput_.dXpressMassDiff2 > 0.0 && strchr(pInput_.szXpressResidues2, 'n') != 0)
                           heavy_nterm_mass_ = light_nterm_mass_ + pInput_.dXpressMassDiff2;
                        else if (pInput_.dXpressMassDiff3 > 0.0 && strchr(pInput_.szXpressResidues3, 'n') != 0)
                           heavy_nterm_mass_ = light_nterm_mass_ + pInput_.dXpressMassDiff3;
                        else if (pInput_.dXpressMassDiff1 < 0.0 && strchr(pInput_.szXpressResidues1, 'n') != 0) {
                           heavy_nterm_mass_ = light_nterm_mass_;
                           light_nterm_mass_ = heavy_nterm_mass_ + pInput_.dXpressMassDiff1;
                        }
                        else if (pInput_.dXpressMassDiff2 < 0.0 && strchr(pInput_.szXpressResidues2, 'n') != 0) {
                           heavy_nterm_mass_ = light_nterm_mass_;
                           light_nterm_mass_ = heavy_nterm_mass_ + pInput_.dXpressMassDiff2;
                        }
                        else if (pInput_.dXpressMassDiff3 < 0.0 && strchr(pInput_.szXpressResidues3, 'n') != 0) {
                           heavy_nterm_mass_ = light_nterm_mass_;
                           light_nterm_mass_ = heavy_nterm_mass_ + pInput_.dXpressMassDiff3;
                        }
                        else {
                           heavy_nterm_mass_ = light_nterm_mass_;
                           light_nterm_mass_ = ResidueMass::getMass('n', monoisotopic_);
                        }
                     }
                     else if (heavy_nterm_mass_ == 0.0 && light_nterm_mass_ == 0.0)
                     {
                        if (pInput_.dXpressMassDiff1 != 0.0 && strchr(pInput_.szXpressResidues1, 'n') != 0)
                        {
                           light_nterm_mass_ = ResidueMass::getMass('n', monoisotopic_);
                           heavy_nterm_mass_ = light_nterm_mass_ + pInput_.dXpressMassDiff1;
                        }
                        else if (pInput_.dXpressMassDiff2 != 0.0 && strchr(pInput_.szXpressResidues2, 'n') != 0)
                        {
                           light_nterm_mass_ = ResidueMass::getMass('n', monoisotopic_);
                           heavy_nterm_mass_ = light_nterm_mass_ + pInput_.dXpressMassDiff2;
                        }
                        else if (pInput_.dXpressMassDiff3 != 0.0 && strchr(pInput_.szXpressResidues3, 'n') != 0)
                        {
                           light_nterm_mass_ = ResidueMass::getMass('n', monoisotopic_);
                           heavy_nterm_mass_ = light_nterm_mass_ + pInput_.dXpressMassDiff3;
                        }
                        else
                        {
                           // error
                           cout << "nterm: " << light_nterm_mass_ << " vs " << heavy_nterm_mass_ << endl;
                           exit(1);
                        }
                     }
                     // must override if user defined
                     else
                     {
                        if (pInput_.dXpressMassDiff1 != 0.0 && strchr(pInput_.szXpressResidues1, 'n') != 0)
                           heavy_nterm_mass_ = light_nterm_mass_ + pInput_.dXpressMassDiff1;
                        else if (pInput_.dXpressMassDiff2 != 0.0 && strchr(pInput_.szXpressResidues2, 'n') != 0)
                           heavy_nterm_mass_ = light_nterm_mass_ + pInput_.dXpressMassDiff2;
                        else if (pInput_.dXpressMassDiff3 != 0.0 && strchr(pInput_.szXpressResidues3, 'n') != 0)
                           heavy_nterm_mass_ = light_nterm_mass_ + pInput_.dXpressMassDiff3;
                     }
                  }             // n case
                  else if (pInput_.szXpressResidues[k] == 'c')
                  {
                     if (heavy_cterm_mass_ > 0.0 && light_cterm_mass_ == 0.0)
                     {
                        if (pInput_.dXpressMassDiff1 != 0.0 && strchr(pInput_.szXpressResidues1, 'c') != 0)
                           light_cterm_mass_ = heavy_cterm_mass_ - pInput_.dXpressMassDiff1;
                        else if (pInput_.dXpressMassDiff2 != 0.0 && strchr(pInput_.szXpressResidues2, 'c') != 0)
                           light_cterm_mass_ = heavy_cterm_mass_ - pInput_.dXpressMassDiff2;
                        else if (pInput_.dXpressMassDiff3 != 0.0 && strchr(pInput_.szXpressResidues3, 'c') != 0)
                           light_cterm_mass_ = heavy_cterm_mass_ - pInput_.dXpressMassDiff3;
                        else
                           light_cterm_mass_ = ResidueMass::getMass('c', monoisotopic_); 
                     }
                     else if (heavy_cterm_mass_ == 0.0 && light_cterm_mass_ > 0.0)
                     {          //must switch them around
                        if (pInput_.dXpressMassDiff1 > 0.0 && strchr(pInput_.szXpressResidues1, 'c') != 0)
                           heavy_cterm_mass_ = light_cterm_mass_ + pInput_.dXpressMassDiff1;
                        else if (pInput_.dXpressMassDiff2 > 0.0 && strchr(pInput_.szXpressResidues2, 'c') != 0)
                           heavy_cterm_mass_ = light_cterm_mass_ + pInput_.dXpressMassDiff2;
                        else if (pInput_.dXpressMassDiff3 > 0.0 && strchr(pInput_.szXpressResidues3, 'c') != 0)
                           heavy_cterm_mass_ = light_cterm_mass_ + pInput_.dXpressMassDiff3;
                        else if (pInput_.dXpressMassDiff1 < 0.0 && strchr(pInput_.szXpressResidues1, 'c') != 0) {
                           heavy_cterm_mass_ = light_cterm_mass_;
                           light_cterm_mass_ = heavy_cterm_mass_ + pInput_.dXpressMassDiff1;
                        }
                        else if (pInput_.dXpressMassDiff2 < 0.0 && strchr(pInput_.szXpressResidues2, 'c') != 0) {
                           heavy_cterm_mass_ = light_cterm_mass_;
                           light_cterm_mass_ = heavy_cterm_mass_ + pInput_.dXpressMassDiff2;
                        }
                        else if (pInput_.dXpressMassDiff3 < 0.0 && strchr(pInput_.szXpressResidues3, 'c') != 0) {
                           heavy_cterm_mass_ = light_cterm_mass_;
                           light_cterm_mass_ = heavy_cterm_mass_ + pInput_.dXpressMassDiff3;
                        }
                        else {
                           heavy_cterm_mass_ = light_cterm_mass_;
                           light_cterm_mass_ = ResidueMass::getMass('c', monoisotopic_);
                        }
                     }
                     else if (heavy_cterm_mass_ == 0.0 && light_cterm_mass_ == 0.0)
                     {
                        if (pInput_.dXpressMassDiff1 != 0.0 && strchr(pInput_.szXpressResidues1, 'c') != 0)
                        {
                           light_cterm_mass_ = ResidueMass::getMass('c', monoisotopic_);
                           heavy_cterm_mass_ = light_cterm_mass_ + pInput_.dXpressMassDiff1;
                        }
                        else if (pInput_.dXpressMassDiff2 != 0.0 && strchr(pInput_.szXpressResidues2, 'c') != 0)
                        {
                           light_cterm_mass_ = ResidueMass::getMass('c', monoisotopic_);
                           heavy_cterm_mass_ = light_cterm_mass_ + pInput_.dXpressMassDiff2;
                        }
                        else if (pInput_.dXpressMassDiff3 != 0.0 && strchr(pInput_.szXpressResidues3, 'c') != 0)
                        {
                           light_cterm_mass_ = ResidueMass::getMass('c', monoisotopic_);
                           heavy_cterm_mass_ = light_cterm_mass_ + pInput_.dXpressMassDiff3;
                        }
                        else
                        {
                           // error
                           cout << "cterm: " << light_cterm_mass_ << " vs " << heavy_cterm_mass_ << endl;
                           exit(1);
                        }
                     }
                     // must override if user defined
                     else
                     {
                        if (pInput_.dXpressMassDiff1 != 0.0 && strchr(pInput_.szXpressResidues1, 'c') != 0)
                           heavy_cterm_mass_ = light_nterm_mass_ + pInput_.dXpressMassDiff1;
                        else if (pInput_.dXpressMassDiff2 != 0.0 && strchr(pInput_.szXpressResidues2, 'c') != 0)
                           heavy_cterm_mass_ = light_nterm_mass_ + pInput_.dXpressMassDiff2;
                        else if (pInput_.dXpressMassDiff3 != 0.0 && strchr(pInput_.szXpressResidues3, 'c') != 0)
                           heavy_cterm_mass_ = light_nterm_mass_ + pInput_.dXpressMassDiff3;
                     
                     }
                  }             // c case
                  else
                  {
                     char next_res = pInput_.szXpressResidues[k];
                     if (heavy_label_masses_[next_res - 'A'] > 0.0 && light_label_masses_[next_res - 'A'] == 0.0)
                     {
                        if (pInput_.dXpressMassDiff1 != 0.0 && strchr(pInput_.szXpressResidues1, next_res) != 0)
                           light_label_masses_[next_res - 'A'] = heavy_label_masses_[next_res - 'A'] - pInput_.dXpressMassDiff1;
                        else if (pInput_.dXpressMassDiff2 != 0.0 && strchr(pInput_.szXpressResidues2, next_res) != 0)
                           light_label_masses_[next_res - 'A'] = heavy_label_masses_[next_res - 'A'] - pInput_.dXpressMassDiff2;
                        else if (pInput_.dXpressMassDiff3 != 0.0 && strchr(pInput_.szXpressResidues3, next_res) != 0)
                           light_label_masses_[next_res - 'A'] = heavy_label_masses_[next_res - 'A'] - pInput_.dXpressMassDiff3;
                        else
                           light_label_masses_[next_res - 'A'] = ResidueMass::getMass(next_res, monoisotopic_);
                     }
                     else if (heavy_label_masses_[next_res - 'A'] == 0.0 && light_label_masses_[next_res - 'A'] > 0)
                     {
                        if (pInput_.dXpressMassDiff1 > 0.0 && strchr(pInput_.szXpressResidues1, next_res) != 0)
                           heavy_label_masses_[next_res - 'A'] = light_label_masses_[next_res - 'A'] + pInput_.dXpressMassDiff1;
                        else if (pInput_.dXpressMassDiff2 > 0.0 && strchr(pInput_.szXpressResidues2, next_res) != 0)
                           heavy_label_masses_[next_res - 'A'] = light_label_masses_[next_res - 'A'] + pInput_.dXpressMassDiff2;
                        else if (pInput_.dXpressMassDiff3 > 0.0 && strchr(pInput_.szXpressResidues3, next_res) != 0)
                           heavy_label_masses_[next_res - 'A'] = light_label_masses_[next_res - 'A'] + pInput_.dXpressMassDiff3;
                        else if (pInput_.dXpressMassDiff1 < 0.0 && strchr(pInput_.szXpressResidues1, next_res) != 0) {
                           heavy_label_masses_[next_res - 'A'] = light_label_masses_[next_res - 'A'];
                           light_label_masses_[next_res - 'A'] = heavy_label_masses_[next_res - 'A'] + pInput_.dXpressMassDiff1;
                        }
                        else if (pInput_.dXpressMassDiff2 < 0.0 && strchr(pInput_.szXpressResidues2, next_res) != 0) {
                           heavy_label_masses_[next_res - 'A'] = light_label_masses_[next_res - 'A'];
                           light_label_masses_[next_res - 'A'] = heavy_label_masses_[next_res - 'A'] + pInput_.dXpressMassDiff2;
                        }
                        else if (pInput_.dXpressMassDiff3 < 0.0 && strchr(pInput_.szXpressResidues3, next_res) != 0) {
                           heavy_label_masses_[next_res - 'A'] = light_label_masses_[next_res - 'A'];
                           light_label_masses_[next_res - 'A'] = heavy_label_masses_[next_res - 'A'] + pInput_.dXpressMassDiff3;
                        }
                        else {
                           heavy_label_masses_[next_res - 'A'] = light_label_masses_[next_res - 'A'];
                           light_label_masses_[next_res - 'A'] = ResidueMass::getMass(next_res, monoisotopic_);
                        }
                     }
                     else if (heavy_label_masses_[next_res - 'A'] == 0.0 && light_label_masses_[next_res - 'A'] == 0)
                     {
                        if (pInput_.dXpressMassDiff1 != 0.0 && strchr(pInput_.szXpressResidues1, next_res) != 0)
                        {
                           light_label_masses_[next_res - 'A'] = ResidueMass::getMass(next_res, monoisotopic_);
                           heavy_label_masses_[next_res - 'A'] = light_label_masses_[next_res - 'A'] + pInput_.dXpressMassDiff1;
                        }
                        else if (pInput_.dXpressMassDiff2 != 0.0 && strchr(pInput_.szXpressResidues2, next_res) != 0)
                        {
                           light_label_masses_[next_res - 'A'] = ResidueMass::getMass(next_res, monoisotopic_);
                           heavy_label_masses_[next_res - 'A'] = light_label_masses_[next_res - 'A'] + pInput_.dXpressMassDiff2;
                        }
                        else if (pInput_.dXpressMassDiff3 != 0.0 && strchr(pInput_.szXpressResidues3, next_res) != 0)
                        {
                           light_label_masses_[next_res - 'A'] = ResidueMass::getMass(next_res, monoisotopic_);
                           heavy_label_masses_[next_res - 'A'] = light_label_masses_[next_res - 'A'] + pInput_.dXpressMassDiff3;
                        }
                        else
                        {
                           // error
                           cout << "label " << next_res << ": " <<
                              light_label_masses_[next_res - 'A'] << " vs " <<
                              heavy_label_masses_[next_res - 'A'] << endl;
                           exit(1);
                        }
                     }
                     // must override if user defined
                     else if (heavy_label_masses_[next_res - 'A'] <= 0.0)
                     {
                        if (pInput_.dXpressMassDiff1 != 0.0 && strchr(pInput_.szXpressResidues1, next_res) != 0)
                           heavy_label_masses_[next_res - 'A'] = light_label_masses_[next_res - 'A'] + pInput_.dXpressMassDiff1;
                        else if (pInput_.dXpressMassDiff2 != 0.0 && strchr(pInput_.szXpressResidues2, next_res) != 0)
                           heavy_label_masses_[next_res - 'A'] = light_label_masses_[next_res - 'A'] + pInput_.dXpressMassDiff2;
                        else if (pInput_.dXpressMassDiff3 != 0.0 && strchr(pInput_.szXpressResidues3, next_res) != 0)
                           heavy_label_masses_[next_res - 'A'] = light_label_masses_[next_res - 'A'] + pInput_.dXpressMassDiff3;
                     }
 //                  cout << "setting light to " <<  light_label_masses_[next_res-'A'] << " and heavy to " << heavy_label_masses_[next_res-'A'] << " for " << next_res << endl;
                  }             // aa case
               }                // next label
#endif
            }
#ifdef USE_STD_MODS
            else if (tag->isStart() && !strcmp(tag->getName(), "search_summary"))
            {
               monoisotopic_ = !strcmp(tag->getAttributeValue("precursor_mass_type"), "monoistotopic");
            }
            // have a modification worth recording here
            else if (tag->isStart()
                  && !strcmp(tag->getName(), "terminal_modification")
                  && strchr(pInput_.szXpressResidues, tolower(tag->getAttributeValue("terminus")[0])) != NULL)
            {
               char next_res = tolower(tag->getAttributeValue("terminus")[0]);
               double nextmass = atof(tag->getAttributeValue("mass"));
               double massdiff = atof(tag->getAttributeValue("massdiff"));
               
               bool bTandem = false;  // check if special Tandem mods by looking for "symbol" attribute

               // X!Tandem special mods use attribute symbol="^"; may be too lax but if
               // "symbol" attribute is present, assume this is the special Tandem mods
               if (tag->getAttributeValue("symbol"))
                  bTandem = true;

               if (next_res == 'n')
               {
                  if (!strcmp(tag->getAttributeValue("variable"), "Y"))
                  {
                     // static mods define the "light" mass & variable mods define "heavy"
                     // must be heavy

                     // rule out Tandem's default +42 terminal mod
                     if (!bTandem || fabs(massdiff - 42.0106) > MAX_MASS_FUDGE)
                     {
                        if (heavy_nterm_mass_ > 0.0 && fabs(nextmass - heavy_nterm_mass_)>MAX_MASS_FUDGE)
                           cout << "WARNING: Found more than one variable mod on N-terminus." << endl; 

                        heavy_nterm_mass_ = nextmass;
                     }

                  }             // variable
                  else
                  {
                     // place in light for now
                     light_nterm_mass_ = nextmass;
                  }
               }                // n terminal
               else if (next_res == 'c')
               {
                  if (!strcmp(tag->getAttributeValue("variable"), "Y"))
                  {             // must be heavy
                     if (heavy_cterm_mass_ > 0.0 && fabs(nextmass - heavy_cterm_mass_)>MAX_MASS_FUDGE)
                        cout << "WARNING: Found more than one variable mod on C-terminus." << endl; 

                     heavy_cterm_mass_ = nextmass;
                  }             // variable
                  else
                  {
                     // place in light for now
                     light_cterm_mass_ = nextmass;
                  }
               }                // c
            }                   // term modification
            else if (tag->isStart()
                  && !strcmp(tag->getName(), "aminoacid_modification")
                  && strchr(pInput_.szXpressResidues, tag->getAttributeValue("aminoacid")[0]) != NULL)
            {
               char next_res = tag->getAttributeValue("aminoacid")[0];
               double nextmass = atof(tag->getAttributeValue("mass"));

               if (!strcmp(tag->getAttributeValue("variable"), "Y"))
               {                // must be heavy
                  if (heavy_label_masses_[next_res - 'A'] > 0.0
                        && fabs(nextmass - heavy_label_masses_[next_res - 'A'])>MAX_MASS_FUDGE)
                  {
                     cout << "WARNING: Found more than one variable mod on \'" << next_res << "\'." << endl; 
                  }

                  heavy_label_masses_[next_res - 'A'] = nextmass;
               }                // variable
               else
               {
                  // place in light for now
                  light_label_masses_[next_res - 'A'] = nextmass;
               }
            }                   // aa modification

            if (modifications != NULL)
            {
               delete modifications;
               modifications = NULL;
            }

            if (modinfo_ != NULL)
               delete modinfo_;
            modinfo_ = NULL;
#endif

            pInput_.iChargeState = -1;       // reset

         }
         // if not filtered
         if (tag != NULL)
            delete tag;
         data = strstr(data + 1, "<");
      }                     // next tag
   }                            // next line
   fin.close();

   delete [] nextline;

   delete timestamp_start;
   delete timestamp_stop;
   delete timestamp;
   if (analysis_start != NULL)
      delete analysis_start;
   if (analysis_stop != NULL)
      delete analysis_stop;
   if (summary != NULL)
      delete summary;
   delete result_start;
   delete result_stop;
   delete xpress_filter;
   delete xpress_summ_filter;
   delete xpress_summ;
}


void XPressPeptideParser::setFilter( Tag * tag) {
   if (tag == NULL)
      return;

   if (filter_memory_) {
      filter_memory_ = False;
      filter_ = False;
   }

   if (!strcmp(tag->getName(), "spectrum_query")) {
      if (tag->isStart())
	filter_ = True;
      else
	filter_memory_ = True;
   }
}


Tag *XPressPeptideParser::getSummaryTag( const InputStruct &opts) {
  char text[1024];
  Tag *output;
  
  if (pInput_.bLabelFreeMode)
    output = new Tag("xpresslabelfree_summary", True, True);
  else
    output = new Tag("xpressratio_summary", True, True);

  char version[200];
  snprintf(version, sizeof(version), "%s (%s)", PROGRAM_VERSION, szTPPVersionInfo);
  output->setAttributeValue("version", version);
  output->setAttributeValue("author", PROGRAM_AUTHOR);

  if (!pInput_.bLabelFreeMode) {
    if (opts.bUseSameScanRange)
      output->setAttributeValue("same_scan_range", "Y");
    else
      output->setAttributeValue("same_scan_range", "N");
    output->setAttributeValue("labeled_residues", opts.szXpressResidues);

    sprintf(text, "%d", opts.bXpressLight1);
    output->setAttributeValue("xpress_light", text);
   
    if (opts.dXpressMassDiff3 > 0.0) {
      sprintf(text, "%s,%0.6f %s,%0.6f %s,%0.6f",
	      opts.szXpressResidues1, opts.dXpressMassDiff1,
	      opts.szXpressResidues2, opts.dXpressMassDiff2,
	      opts.szXpressResidues3, opts.dXpressMassDiff3);
    }
    else if (opts.dXpressMassDiff2 > 0.0) {
      sprintf(text, "%s,%0.6f %s,%0.6f",
	      opts.szXpressResidues1, opts.dXpressMassDiff1,
	      opts.szXpressResidues2, opts.dXpressMassDiff2);
    }
    else
      sprintf(text, "%s,%0.6f", opts.szXpressResidues1, opts.dXpressMassDiff1);

    output->setAttributeValue("massdiff", text);
  }

  sprintf(text, "%0.6f", opts.dMassTol);
  output->setAttributeValue("masstol", text);

  sprintf(text, "%d", opts.bPpmMassTol);
  output->setAttributeValue("ppmtol", text);

  sprintf(text, "%d", opts.iMinNumChromatogramPoints);
  output->setAttributeValue("min_num_chromatogram_points", text);

  sprintf(text, "%d", opts.iMinNumIsotopePeaks);
  output->setAttributeValue("min_num_isotope_peaks", text);

  return output;
}

#ifdef USE_STD_MODS

Tag *XPressPeptideParser::getRatio()
{
   char szBuf[SIZE_BUF];
   // check modifications only to see if heavy or light or none
   Boolean light = False;
   Boolean heavy = False;
   Boolean unmod = False;       // could be counted as light

   int negMassDiff = 1;
   // calculate if correct, and mass diff between heavy and light
   double massdiff = 0.0;

   if (pInput_.bLabelFreeMode)
   {
      pXpress_.dLightPeptideMass = pInput_.dPeptideMass;
      this->XPRESS_ANALYSIS_LABELFREE();
      Tag *output = new Tag("xpresslabelfree_result", True, True);

      sprintf(szBuf, "%d", pInput_.iChargeState);
      output->setAttributeValue("charge", szBuf);

      sprintf(szBuf, "%d", piSequentialScan[pXpress_.iLightFirstScan]);
      output->setAttributeValue("first_scan", szBuf);
      sprintf(szBuf, "%d", piSequentialScan[pXpress_.iLightLastScan]);
      output->setAttributeValue("last_scan", szBuf);

      sprintf(szBuf, "%0.0lf", pXpress_.dLightFirstScanRT);
      output->setAttributeValue("first_scan_RT_seconds", szBuf);
      sprintf(szBuf, "%0.0lf", pXpress_.dLightLastScanRT);
      output->setAttributeValue("last_scan_RT_seconds", szBuf);

      sprintf(szBuf, "%0.4f", (PROTON_MASS * (pInput_.iChargeState - 1) + pXpress_.dLightPeptideMass) / pInput_.iChargeState);
      output->setAttributeValue("precursor_mz", szBuf);

      sprintf(szBuf, "%0.2e", pXpress_.dLightArea);
      output->setAttributeValue("peak_area", szBuf);
      sprintf(szBuf, "%0.2e", pXpress_.dLightIntensity);
      output->setAttributeValue("peak_intensity", szBuf);
      sprintf(szBuf, "%0.0lf", pXpress_.dLightIntensityRT);
      output->setAttributeValue("peak_intensity_RT_seconds", szBuf);
      sprintf(szBuf, "%d", piSequentialScan[pXpress_.iLightIntensityScan]);
      output->setAttributeValue("peak_intensity_scan", szBuf);

      return output;
   }

   // MUST CHECK FOR N AND C TERMINAL MODS HERE FIRST............
   if (modinfo_ != NULL)
   {
      if (strchr(pInput_.szXpressResidues, 'n') != NULL)
      {
         if (modinfo_->getNtermModMass() > 0.0)
         {
            double nextmass = modinfo_->getNtermModMass();
            if (nextmass - light_nterm_mass_ <= MOD_ERROR && light_nterm_mass_ - nextmass <= MOD_ERROR)
               light = True;
            else if (nextmass - heavy_nterm_mass_ <= MOD_ERROR && heavy_nterm_mass_ - nextmass <= MOD_ERROR)
               heavy = True;
         }                      // have modified terminus
         else
         {
            if (fabs(light_nterm_mass_ - ResidueMass::getMass('n', monoisotopic_)) < 0.001)
               light = True;
            else
               unmod = True;    // illegal
         }

         if (pInput_.dXpressMassDiff1 != 0.0 && strchr(pInput_.szXpressResidues1, 'n') != 0)
            massdiff += pInput_.dXpressMassDiff1;
         else if (pInput_.dXpressMassDiff2 != 0.0 && strchr(pInput_.szXpressResidues2, 'n') != 0)
            massdiff += pInput_.dXpressMassDiff2;
         else if (pInput_.dXpressMassDiff3 != 0.0 && strchr(pInput_.szXpressResidues3, 'n') != 0)
            massdiff += pInput_.dXpressMassDiff3;
         else
            massdiff += heavy_nterm_mass_ - light_nterm_mass_;

      }                         // n
      if (strchr(pInput_.szXpressResidues, 'c') != NULL)
      {
         if (modinfo_->getCtermModMass() > 0.0)
         {
            double nextmass = modinfo_->getCtermModMass();
            if (nextmass - light_cterm_mass_ <= MOD_ERROR && light_cterm_mass_ - nextmass <= MOD_ERROR)
               light = True;
            else if (nextmass - heavy_cterm_mass_ <= MOD_ERROR && heavy_cterm_mass_ - nextmass <= MOD_ERROR)
               heavy = True;
         }                      // have modified terminus
         else
         {
            if (fabs(light_cterm_mass_ - ResidueMass::getMass('c', monoisotopic_)) < 0.001)
               light = True;
            else
               unmod = True;    // illegal
         }
         if (pInput_.dXpressMassDiff1 != 0.0 && strchr(pInput_.szXpressResidues1, 'c') != 0)
            massdiff += pInput_.dXpressMassDiff1;
         else if (pInput_.dXpressMassDiff2 != 0.0 && strchr(pInput_.szXpressResidues2, 'c') != 0)
            massdiff += pInput_.dXpressMassDiff2;
         else if (pInput_.dXpressMassDiff3 != 0.0 && strchr(pInput_.szXpressResidues3, 'c') != 0)
            massdiff += pInput_.dXpressMassDiff3;
         else
            massdiff += heavy_cterm_mass_ - light_cterm_mass_;

      }                         // c
   }                            // some mods to look at

   for (int k = 0; pInput_.szPeptide[k]; k++)
   {
      if (strchr(pInput_.szXpressResidues, pInput_.szPeptide[k]) != NULL)
      {                         // have a labeled aa
         double nextmass = modinfo_ == NULL ? 0.0 : modinfo_->getModifiedResidueMass(k);
         if (nextmass > 0.0)
         {                      // really is modified, find out if light or heavy

            if (nextmass - light_label_masses_[pInput_.szPeptide[k] - 'A'] <=
                MOD_ERROR && light_label_masses_[pInput_.szPeptide[k] - 'A'] - nextmass <= MOD_ERROR)
               light = True;
            else if (nextmass - heavy_label_masses_[pInput_.szPeptide[k] - 'A'] <= MOD_ERROR
                     && heavy_label_masses_[pInput_.szPeptide[k] - 'A'] - nextmass <= MOD_ERROR)
               heavy = True;
         }                      // have modified aa
         else
         {
            light = True;
         }
         if (pInput_.dXpressMassDiff1 != 0.0 && strchr(pInput_.szXpressResidues1, pInput_.szPeptide[k]) != 0)
            massdiff += pInput_.dXpressMassDiff1;
         else if (pInput_.dXpressMassDiff2 != 0.0 && strchr(pInput_.szXpressResidues2, pInput_.szPeptide[k]) != 0)
            massdiff += pInput_.dXpressMassDiff2;
         else if (pInput_.dXpressMassDiff3 != 0.0 && strchr(pInput_.szXpressResidues3, pInput_.szPeptide[k]) != 0)
            massdiff += pInput_.dXpressMassDiff3;
         else
            massdiff += heavy_label_masses_[pInput_.szPeptide[k] - 'A'] - light_label_masses_[pInput_.szPeptide[k] - 'A'];

      }                         // modified res
   }

   //DDS: Check if the light mod is on the terminals
   if (!heavy && !light) {
      massdiff=0;
      if (strchr(pInput_.szXpressResidues, 'c') != NULL) {
         if (pInput_.dXpressMassDiff1 != 0.0 && strchr(pInput_.szXpressResidues1, 'c') != 0)
            massdiff += pInput_.dXpressMassDiff1;
         else if (pInput_.dXpressMassDiff2 != 0.0 && strchr(pInput_.szXpressResidues2, 'c') != 0)
            massdiff += pInput_.dXpressMassDiff2;
         else if (pInput_.dXpressMassDiff3 != 0.0 && strchr(pInput_.szXpressResidues3, 'c') != 0)
            massdiff += pInput_.dXpressMassDiff3;
         else
            massdiff += heavy_cterm_mass_ - light_cterm_mass_;
         light = True;
         unmod = False;
      }
      if (strchr(pInput_.szXpressResidues, 'n') != NULL) {
         if (pInput_.dXpressMassDiff1 != 0.0 && strchr(pInput_.szXpressResidues1, 'n') != 0)
            massdiff += pInput_.dXpressMassDiff1;
         else if (pInput_.dXpressMassDiff2 != 0.0 && strchr(pInput_.szXpressResidues2, 'n') != 0)
            massdiff += pInput_.dXpressMassDiff2;
         else if (pInput_.dXpressMassDiff3 != 0.0 && strchr(pInput_.szXpressResidues3, 'n') != 0)
            massdiff += pInput_.dXpressMassDiff3;
         else
            massdiff += heavy_nterm_mass_ - light_nterm_mass_;
         light = True;
         unmod = False;
      } 
   }
   if (heavy && massdiff<0) {
      massdiff = -1.*massdiff;
   }
   if (pInput_.iMetabolicLabeling || pInput_.bLabelFreeMode
       || (!unmod && (light || heavy) && (!light || !heavy) && massdiff != 0.0))
   {                            // valid quant data to process
      if (pInput_.iMetabolicLabeling == 1)
      {
         if (strlen(pInput_.szPeptide) == pInput_.iModCount)  // N15 heavy since every residue is modified
         {
            pXpress_.dHeavyPeptideMass = pInput_.dPeptideMass;
            pXpress_.dLightPeptideMass = pInput_.dPeptideMass - ELEMENT_COUNT(pInput_.szPeptide, 'N');
         }
         else
         {
            pXpress_.dLightPeptideMass = pInput_.dPeptideMass;
            pXpress_.dHeavyPeptideMass = pInput_.dPeptideMass + ELEMENT_COUNT(pInput_.szPeptide, 'N');
         }
      }
      else if (pInput_.iMetabolicLabeling == 2)
      {
         if (strlen(pInput_.szPeptide) == pInput_.iModCount)  // C13 heavy since every residue is modified
         {
            pXpress_.dHeavyPeptideMass = pInput_.dPeptideMass;
            pXpress_.dLightPeptideMass = pInput_.dPeptideMass - ELEMENT_COUNT(pInput_.szPeptide, 'C');
         }
         else
         {
            pXpress_.dLightPeptideMass = pInput_.dPeptideMass;
            pXpress_.dHeavyPeptideMass = pInput_.dPeptideMass + ELEMENT_COUNT(pInput_.szPeptide, 'C');
         }
      }
      else
      {
         pXpress_.dLightPeptideMass = light ? pInput_.dPeptideMass : pInput_.dPeptideMass - massdiff;
         pXpress_.dHeavyPeptideMass = heavy ? pInput_.dPeptideMass : pInput_.dPeptideMass + massdiff;
      }

      // continue to process now....
      this->XPRESS_ANALYSIS((int)light);

      Tag *output = new Tag("xpressratio_result", True, True);

      pXpress_.iChargeState = pInput_.iChargeState;
      pXpress_.bXpressLight1 = pInput_.bXpressLight1;
      pXpress_.iMetabolicLabeling = pInput_.iMetabolicLabeling;
      pXpress_.dMassTol = pInput_.dMassTol;

      sprintf(szBuf, "%d", piSequentialScan[pXpress_.iLightFirstScan]);
      output->setAttributeValue("light_firstscan", szBuf);
      sprintf(szBuf, "%d", piSequentialScan[pXpress_.iLightLastScan]);
      output->setAttributeValue("light_lastscan", szBuf);
      sprintf(szBuf, "%0.4f", pXpress_.dLightPeptideMass);
      output->setAttributeValue("light_mass", szBuf);
      sprintf(szBuf, "%d", piSequentialScan[pXpress_.iHeavyFirstScan]);
      output->setAttributeValue("heavy_firstscan", szBuf);
      sprintf(szBuf, "%d", piSequentialScan[pXpress_.iHeavyLastScan]);
      output->setAttributeValue("heavy_lastscan", szBuf);
      sprintf(szBuf, "%0.4f", pXpress_.dHeavyPeptideMass);
      output->setAttributeValue("heavy_mass", szBuf);

      sprintf(szBuf, "%0.4f", pXpress_.dMassTol);
      output->setAttributeValue("mass_tol", szBuf);
      sprintf(szBuf, "%s", pXpress_.szQuan);
      output->setAttributeValue("ratio", szBuf);

      // now the flipped guy....
      flipRatio(szBuf, szBuf);
      output->setAttributeValue("heavy2light_ratio", szBuf);

      sprintf(szBuf, "%0.2e", pXpress_.dLightArea);
      output->setAttributeValue("light_area", szBuf);
      sprintf(szBuf, "%0.2e", pXpress_.dHeavyArea);
      output->setAttributeValue("heavy_area", szBuf);

      if (pXpress_.dLightArea == 0.0)
      {
         if (pXpress_.dHeavyArea > 0.0)
            sprintf(szBuf, "%0.3f", 0.001);
         else
            sprintf(szBuf, "%0.1f", -1.0);
      }
      else
      {
         if (pXpress_.dHeavyArea != 0.0)
            sprintf(szBuf, "%0.2f", pXpress_.dLightArea / pXpress_.dHeavyArea);
         else
            sprintf(szBuf, "%0.1f", 999.0);
      }
      output->setAttributeValue("decimal_ratio", szBuf);

      if (pInput_.bIntensityRatio)
      {
         sprintf(szBuf, "%0.2e", pXpress_.dLightIntensity);
         output->setAttributeValue("light_intensity", szBuf);
         sprintf(szBuf, "%0.0lf", pXpress_.dLightIntensityRT);
         output->setAttributeValue("light_intensity_RT", szBuf);
         sprintf(szBuf, "%d", piSequentialScan[pXpress_.iLightIntensityScan]);
         output->setAttributeValue("light_intensity_scan", szBuf);

         sprintf(szBuf, "%0.2e", pXpress_.dHeavyIntensity);
         output->setAttributeValue("heavy_intensity", szBuf);
         sprintf(szBuf, "%0.0lf", pXpress_.dHeavyIntensityRT);
         output->setAttributeValue("heavy_intensity_RT", szBuf);
         sprintf(szBuf, "%d", pXpress_.iHeavyIntensityScan);
         output->setAttributeValue("heavy_intensity_scan", szBuf);

         sprintf(szBuf, "%0.4f", pXpress_.dIntensityRatio);
         output->setAttributeValue("intensity_ratio", szBuf);
      }

      return output;
   }

   return NULL;
}


#endif
#ifndef USE_STD_MODS

Tag *XPressPeptideParser::getRatio()
{
   int bHasXpress;
   int iLen;
   int ii;
   char szBuf[SIZE_BUF];
   int negMassDiff = 1;

   if (pInput_.bLabelFreeMode)
   {
      pXpress_.dLightPeptideMass = dLightPeptideMass;
      this->XPRESS_ANALYSIS_LABELFREE();
      Tag *output = new Tag("xpresslabelfree_result", True, True);

      sprintf(szBuf, "%d", pXpress_.iChargeState);
      output->setAttributeValue("charge", szBuf);
      sprintf(szBuf, "%d", piSequentialScan[pXpress_.iLightFirstScan]);
      output->setAttributeValue("firstscan", szBuf);
      sprintf(szBuf, "%d", piSequentialScan[pXpress_.iLightLastScan]);
      output->setAttributeValue("lastscan", szBuf);

      output->setAttributeValue("firstscan_time_sec", szBuf);
      output->setAttributeValue("lastscan_time_sec", szBuf);

      sprintf(szBuf, "%0.4f", pXpress_.dLightPeptideMass);
      output->setAttributeValue("precursor_mz", szBuf);

      sprintf(szBuf, "%0.2e", pXpress_.dLightIntensity);
      output->setAttributeValue("peak_intensity", szBuf);
      sprintf(szBuf, "%0.0lf", pXpress_.dLightIntensityRT);
      output->setAttributeValue("peak_intensity_RT", szBuf);
      sprintf(szBuf, "%d", piSequentialScan[pXpress_.iLightIntensityScan]);
      output->setAttributeValue("peak_intensity_scan", szBuf);

      return output;
   }

   // check to see if it has a labeled residue
   bHasXpress = FALSE;
   iLen = (int) strlen(pInput_.szPeptide);
   for (ii = 0; ii < iLen; ii++)
   {
      if (strchr(pInput_.szXpressResidues1, pInput_.szPeptide[ii])
          || strchr(pInput_.szXpressResidues2, pInput_.szPeptide[ii])
          || strchr(pInput_.szXpressResidues3, pInput_.szPeptide[ii]))
      {
         bHasXpress = TRUE;
         break;
      }
   }

   // check if peptide contains Cys, is light or heavy, and valid (not both light & heavy)
   if (bHasXpress)
   {
      int iNumModRes1, iNumModRes2, iNumModRes3, bLightPeptide = TRUE, bError;

      bError = FALSE;

      if ((ii + 1) == iLen || isalpha(pInput_.szPeptide[ii + 1]))
         bLightPeptide = TRUE;
      else
         bLightPeptide = FALSE;

      if (pInput_.dXpressMassDiff1 + pInput_.dXpressMassDiff2 + pInput_.dXpressMassDiff3 < 0) {
         negMassDiff = -1;
      }

      if (bLightPeptide && negMassDiff < 0) {
         bLightPeptide = FALSE;
      }


      iNumModRes1 = 0;
      iNumModRes2 = 0;
      iNumModRes3 = 0;

      for (ii = 0; ii < iLen; ii++)
      {
         if (strchr(pInput_.szXpressResidues1, pInput_.szPeptide[ii]))
         {
            iNumModRes1++;

            // Make sure not both heavy/light residue in peptide
            if (bLightPeptide)
            {
               if (ii + 1 < iLen && !isalpha(pInput_.szPeptide[ii + 1]))
               {
                  bError = TRUE;
                  break;
               }
            }
            else
            {
               if (ii + 1 == iLen || isalpha(pInput_.szPeptide[ii + 1]))
               {
                  bError = TRUE;
                  break;
               }
            }
         }
         else if (strchr(pInput_.szXpressResidues2, pInput_.szPeptide[ii]))
         {
            iNumModRes2++;

            // Make sure not both heavy/light residue in peptide
            if (bLightPeptide)
            {
               if (ii + 1 < iLen && !isalpha(pInput_.szPeptide[ii + 1]))
               {
                  bError = TRUE;
                  break;
               }
            }
            else
            {
               if (ii + 1 == iLen || isalpha(pInput_.szPeptide[ii + 1]))
               {
                  bError = TRUE;
                  break;
               }
            }
         }
         else if (strchr(pInput_.szXpressResidues3, pInput_.szPeptide[ii]))
         {
            iNumModRes3++;

            // Make sure not both heavy/light residue in peptide
            if (bLightPeptide)
            {
               if (ii + 1 < iLen && !isalpha(pInput_.szPeptide[ii + 1]))
               {
                  bError = TRUE;
                  break;
               }
            }
            else
            {
               if (ii + 1 == iLen || isalpha(pInput_.szPeptide[ii + 1]))
               {
                  bError = TRUE;
                  break;
               }
            }
         }
      }

      // Peptide is a valid, isotopically labeled peptide with only heavy or light residues
      if (!bError)
      {
         char szTmp[SIZE_BUF];
         int iLightFirstScan,
             iLightLastScan,
             iHeavyFirstScan,
             iHeavyLastScan;
         double dLightPeptideMass,
                dHeavyPeptideMass;

         // Calculate the corresponding heavy/light MH+ peptide mass
         if (bLightPeptide)
         {
            dLightPeptideMass = pInput_.dPeptideMass;
            dHeavyPeptideMass =
               pInput_.dPeptideMass + pInput_.dXpressMassDiff1 * iNumModRes1 +
               pInput_.dXpressMassDiff2 * iNumModRes2 +
               pInput_.dXpressMassDiff3 * iNumModRes3;
         }
         else
         {
            dHeavyPeptideMass = pInput_.dPeptideMass;
            dLightPeptideMass =
               pInput_.dPeptideMass - negMassDiff * pInput_.dXpressMassDiff1 * iNumModRes1 -
               negMassDiff * pInput_.dXpressMassDiff2 * iNumModRes2 -
               negMassDiff * pInput_.dXpressMassDiff3 * iNumModRes3;
         }


         pXpress_.dLightPeptideMass = dLightPeptideMass;
         pXpress_.dHeavyPeptideMass = dHeavyPeptideMass;

         XPRESS_ANALYSIS(bLightPeptide);

         pXpress_.iChargeState = pInput_.iChargeState;
         pXpress_.bXpressLight1 = pInput_.bXpressLight1;
         pXpress_.dMassTol = pInput_.dMassTol;

         // print the XPRESS link

         Tag   *output = new Tag("xpressratio_result", True, True);
         sprintf(szBuf, "%d", piSequentialScan[pXpress_.iLightFirstScan]);
         output->setAttributeValue("light_firstscan", szBuf);
         sprintf(szBuf, "%d", piSequentialScan[pXpress_.iLightLastScan]);
         output->setAttributeValue("light_lastscan", szBuf);
         sprintf(szBuf, "%0.3f", pXpress_.dLightPeptideMass);
         output->setAttributeValue("light_mass", szBuf);
         sprintf(szBuf, "%d", piSequentialScan[pXpress_.iHeavyFirstScan]);
         output->setAttributeValue("heavy_firstscan", szBuf);
         sprintf(szBuf, "%d", piSequentialScan[pXpress_.iHeavyLastScan]);
         output->setAttributeValue("heavy_lastscan", szBuf);
         sprintf(szBuf, "%0.3f", pXpress_.dHeavyPeptideMass);
         output->setAttributeValue("heavy_mass", szBuf);

         sprintf(szBuf, "%0.3f", pXpress_.dMassTol);
         output->setAttributeValue("mass_tol", szBuf);
         sprintf(szBuf, "%s", pXpress_.szQuan);
         output->setAttributeValue("ratio", szBuf);

         // now the flipped guy....
         flipRatio(szBuf, szBuf);
         output->setAttributeValue("heavy2light_ratio", szBuf);

         sprintf(szBuf, "%0.2e", pXpress_.dLightArea);
         output->setAttributeValue("light_area", szBuf);
         sprintf(szBuf, "%0.2e", pXpress_.dHeavyArea);
         output->setAttributeValue("heavy_area", szBuf);

         if (pXpress_.dLightArea == 0.0) {
	   if (pXpress_.dHeavyArea > 0.0)
	     sprintf(szBuf, "%0.3f", 0.001);
	   else
	     sprintf(szBuf, "%0.1f", -1.0);
         }
         else {
	   if (pXpress_.dHeavyArea != 0.0)
	     sprintf(szBuf, "%0.2f", pXpress_.dLightArea / pXpress_.dHeavyArea);
	   else
	     sprintf(szBuf, "%0.1f", 999.0);
         }
         output->setAttributeValue("decimal_ratio", szBuf);

         if (verbose_)
	   output->write(cout);

         return output;
      }
   }
   return NULL;
}

#endif


// Return mass difference based on number of nitrogens/carbons in peptide
// sequence for metabolic labeling
double XPressPeptideParser::ELEMENT_COUNT(char *szPeptide,
      char cElement)
{
   int i;
   int iMassDiff;
   int iLen;
   double dDiff;
   if (cElement=='N')
      dDiff = 0.997034893;
   else if (cElement=='C')
      dDiff = 1.0033548378;
   else
   {
      printf(" Error - unknown element %c for ELEMENT_COUNT\n", cElement);
      exit(1);
   }

   // http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=N&ascii=html&isotype=some
   // 14N = 14.003074005, 15N = 15.000108898, diff = 0.997034893

   // http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=C&ascii=html&isotype=some
   // 12C = 12.0000000, 13C = 13.0033548378, diff = 1.0033548378

   // retired:
   // 14N = 14.003074020, 15N = 15.00010897, diff = 0.99703477
   // http://www.webelements.com/webelements/elements/text/N/isot.html

   iMassDiff = 0;
   iLen = (int) strlen(szPeptide);

   if (cElement=='N')
   {
      for (i = 0; i < iLen; i++)
      {
         if (szPeptide[i] == 'A')
            iMassDiff += 1;
         else if (szPeptide[i] == 'R')
            iMassDiff += 4;
         else if (szPeptide[i] == 'N')
            iMassDiff += 2;
         else if (szPeptide[i] == 'D')
            iMassDiff += 1;
         else if (szPeptide[i] == 'C')
            iMassDiff += 1;
         else if (szPeptide[i] == 'E')
            iMassDiff += 1;
         else if (szPeptide[i] == 'Q')
            iMassDiff += 2;
         else if (szPeptide[i] == 'G')
            iMassDiff += 1;
         else if (szPeptide[i] == 'H')
            iMassDiff += 3;
         else if (szPeptide[i] == 'I')
            iMassDiff += 1;
         else if (szPeptide[i] == 'L')
            iMassDiff += 1;
         else if (szPeptide[i] == 'K')
            iMassDiff += 2;
         else if (szPeptide[i] == 'M')
            iMassDiff += 1;
         else if (szPeptide[i] == 'F')
            iMassDiff += 1;
         else if (szPeptide[i] == 'P')
            iMassDiff += 1;
         else if (szPeptide[i] == 'S')
            iMassDiff += 1;
         else if (szPeptide[i] == 'T')
            iMassDiff += 1;
         else if (szPeptide[i] == 'W')
            iMassDiff += 2;
         else if (szPeptide[i] == 'Y')
            iMassDiff += 1;
         else if (szPeptide[i] == 'V')
            iMassDiff += 1;
      }
   }
   else if (cElement=='C')
   {
      for (i = 0; i < iLen; i++)
      {
         if (szPeptide[i] == 'A')
            iMassDiff += 3;
         else if (szPeptide[i] == 'R')
            iMassDiff += 6;
         else if (szPeptide[i] == 'N')
            iMassDiff += 4;
         else if (szPeptide[i] == 'D')
            iMassDiff += 4;
         else if (szPeptide[i] == 'C')
            iMassDiff += 3;
         else if (szPeptide[i] == 'E')
            iMassDiff += 5;
         else if (szPeptide[i] == 'Q')
            iMassDiff += 5;
         else if (szPeptide[i] == 'G')
            iMassDiff += 2;
         else if (szPeptide[i] == 'H')
            iMassDiff += 6;
         else if (szPeptide[i] == 'I')
            iMassDiff += 6;
         else if (szPeptide[i] == 'L')
            iMassDiff += 6;
         else if (szPeptide[i] == 'K')
            iMassDiff += 6;
         else if (szPeptide[i] == 'M')
            iMassDiff += 5;
         else if (szPeptide[i] == 'F')
            iMassDiff += 9;
         else if (szPeptide[i] == 'P')
            iMassDiff += 5;
         else if (szPeptide[i] == 'S')
            iMassDiff += 3;
         else if (szPeptide[i] == 'T')
            iMassDiff += 4;
         else if (szPeptide[i] == 'W')
            iMassDiff += 11;
         else if (szPeptide[i] == 'Y')
            iMassDiff += 9;
         else if (szPeptide[i] == 'V')
            iMassDiff += 5;
      }
   }

   return (iMassDiff * dDiff);
}

// verify mzXMLfile goodness
void XPressPeptideParser::validate_mzXMLfile() {
   if (-1 == m_XMLfile_state) {
      // first pass
      if (fp_) { // basename-derived file opened ok
         m_XMLfile_state = 1; // normal operation
      } else {
         m_XMLfile_state = 0; // custom operation
      }
   }
   if (!m_XMLfile_state) { // deriving filenames from spectrum names
      // are we still looking at the same raw data as last time through?
      char *specname = (char *)malloc(strlen(pInput_.szSpectrumName)+strlen(mzXMLfile_)+2);
      strcpy(specname,mzXMLfile_);
      char *slash = findRightmostPathSeperator(specname);
      strcpy(slash?slash+1:specname,pInput_.szSpectrumName);
      char *dot = strchr(specname,'.');
      if (dot) {
         *dot = 0;
      }

      char *oldname = strdup(mzXMLfile_);
      rampConstructInputPath(mzXMLfile_, sizeof(mzXMLfile_), pInput_.szMzXMLDir, specname);
      if (strcmp(oldname,mzXMLfile_))
      {
         // changing files
         cached_ramp_rampCloseFile(fp_);
         free(index_);
         index_ = NULL;

         if ((fp_ = cached_ramp_rampOpenFile(mzXMLfile_)) == NULL)
         {
            printf("scan-derived scan file %s (from %s) not found, cannot proceed...\n",
                  mzXMLfile_,pInput_.szSpectrumName);
            exit(1);
         }
      }
      free(specname);
      free(oldname);
   } // end deriving filenames from spectrum names
   if (!fp_) {
      printf("scan file %s not found, cannot proceed...\n", mzXMLfile_);
      exit(1);
   }
}


// Reads mzXML files and get quantitation numbers
void XPressPeptideParser::XPRESS_ANALYSIS(Boolean bLightPeptide, int iScanLocalBuffer)
{

   int ii=0,
       ctScan=0,
       iLightStartScan=0,
       iLightEndScan=0,
       iHeavyStartScan=0,
       iHeavyEndScan=0,
       iStart=0,
       iEnd=0;
   int iNumChromatogramPoints;
   double dLightMass,
          dHeavyMass,
          dLightQuanValue,
          dHeavyQuanValue;
   validate_mzXMLfile(); // do we have a valid scan data file?

   if (pInput_.iChargeState < 1) 
   {
      printf(" Error, charge state = %d\n\n", pInput_.iChargeState);
      exit(EXIT_FAILURE);
   }
   else
   {
      dLightMass = (PROTON_MASS * (pInput_.iChargeState - 1) + pXpress_.dLightPeptideMass) / pInput_.iChargeState;
      dHeavyMass = (PROTON_MASS * (pInput_.iChargeState - 1) + pXpress_.dHeavyPeptideMass) / pInput_.iChargeState;
   }

   iStart = piReverseSequentialScan[pInput_.iFirstScan] - (SCAN_BUFFER + iScanLocalBuffer);
   iEnd = piReverseSequentialScan[pInput_.iLastScan] + (SCAN_BUFFER + iScanLocalBuffer);

   if (iStart < pInput_.iAnalysisFirstScan)
      iStart = pInput_.iAnalysisFirstScan;
   if (iEnd > pInput_.iAnalysisLastScan)
      iEnd = pInput_.iAnalysisLastScan;

   // Clear all values
   memset(dLightMS_ + iStart, 0, (iEnd - iStart + 1) * sizeof(double));
   memset(dHeavyMS_ + iStart, 0, (iEnd - iStart + 1) * sizeof(double));

   // Read all MS scan values

   // Expensive loop, cache values.
   double dLightMassLower;
   double dLightMassUpper;
   double dHeavyMassLower;
   double dHeavyMassUpper;

   if( pInput_.bPpmMassTol )
     {
       dLightMassLower = dLightMass - dLightMass * pInput_.dMassTol / 1000000;
       dLightMassUpper = dLightMass + dLightMass * pInput_.dMassTol / 1000000;
       dHeavyMassLower = dHeavyMass - dHeavyMass * pInput_.dMassTol / 1000000;
       dHeavyMassUpper = dHeavyMass + dHeavyMass * pInput_.dMassTol / 1000000;
     }
   else
     {
       dLightMassLower = dLightMass - pInput_.dMassTol;
       dLightMassUpper = dLightMass + pInput_.dMassTol;
       dHeavyMassLower = dHeavyMass - pInput_.dMassTol;
       dHeavyMassUpper = dHeavyMass + pInput_.dMassTol;
     }

   const RAMPREAL *pPeaks;
   RAMPREAL fMass;
   RAMPREAL fInten;
   int n;

   // only track up to first 2 isotope peaks (in addition to monoisotopic peak)
   double dC13diff = 1.00335483;

   for (n=0 ; n<=pInput_.iMinNumIsotopePeaks; n++)
   {
      memset(dIsotopeLight[n]+iStart, 0, (iEnd - iStart +1)*sizeof(double));
      memset(dIsotopeHeavy[n]+iStart, 0, (iEnd - iStart +1)*sizeof(double));
   }

   // pInput_.iMinNumIsotopePeaks only valid if pInput_.dMassTol<=0.5 otherwise
   // makes no sense to try to reconstruct individual isotope chromatograms
   if (pInput_.iMinNumIsotopePeaks<0 || (pInput_.iMinNumIsotopePeaks>=1 && pInput_.dMassTol>0.5))
      pInput_.iMinNumIsotopePeaks=0;
   if (pInput_.iMinNumIsotopePeaks>MAX_ISOTOPES)
      pInput_.iMinNumIsotopePeaks=MAX_ISOTOPES;

   // dLightMS_/dHeavyMS_ contains the main chromatogram which in the case
   // of highres data (i.e. iMinNumIsotopePeaks>=1) means the summed
   // intensities of the individual isotopic peaks.
   //
   // To extend better highres MS1 support, dIsotopeLight and dIsotopeHeavy
   // track individual isotope chromatograms.

   double dSmallestMass = dLightMassLower;
   if (dHeavyMassLower < dSmallestMass)
      dSmallestMass = dHeavyMassLower;

   double dLargestMass = dHeavyMassUpper + 2*dC13diff;
   if (dLightMassUpper+2*dC13diff > dLargestMass)
      dLargestMass = dLightMassUpper+2*dC13diff;

   for (ctScan = iStart; ctScan <= iEnd; ctScan++)
   {
      const struct ScanHeaderStruct* pHeader = cached_ramp_readHeader(fp_, index_[piSequentialScan[ctScan]]);
      if (pHeader->msLevel != 1)
          continue;

      // Open a scan
      pPeaks = cached_ramp_readPeaks_const(fp_, index_[piSequentialScan[ctScan]]);
      if (pPeaks == NULL)
          continue;

      // Binary search to the heaviest peak <= dSmallestMass.
      int l = 0, r = pHeader->peaksCount;
      while (l < r)
      {
          int m = l + (r - l) / 2;
          fMass = pPeaks[m*2];
          if (fMass == dSmallestMass)
          {
              l = r = m;
              break;
          }
          else if (fMass > dSmallestMass)
               r = m - 1;
          else
               l = m + 1;
      }
      n = l * 2;


      // Step forward to lightest peak >= dSmallestMass, if necessary.
      fMass = pPeaks[n];
      if (fMass != -1 && fMass < dSmallestMass)
          n += 2;

      // Scan peak intensities between dLightMassLower and dLightMassUpper.
      while (pPeaks[n] != -1)
      {
         fMass = pPeaks[n];
         n++;

         if (fMass>dLargestMass)
            break;

         fInten = pPeaks[n];
         n++;

         for (ii=0; ii<=pInput_.iMinNumIsotopePeaks; ii++)
         {
            if (dLightMassLower+ii*dC13diff <= fMass && fMass <= dLightMassUpper+ii*dC13diff)
            {
               if (fInten > dIsotopeLight[ii][ctScan])
                  dIsotopeLight[ii][ctScan] = (double) fInten;
            }

            if (dHeavyMassLower+ii*dC13diff <= fMass && fMass <= dHeavyMassUpper+ii*dC13diff)
            {
               if (fInten > dIsotopeHeavy[ii][ctScan])
                  dIsotopeHeavy[ii][ctScan] = (double) fInten;
            }
         }
      }

      for (ii=0; ii<=pInput_.iMinNumIsotopePeaks; ii++)
      {
         dLightMS_[ctScan] += dIsotopeLight[ii][ctScan];
         dHeavyMS_[ctScan] += dIsotopeHeavy[ii][ctScan];
      }

   }

   // Now that we have an MS profile of each of the two
   // light/heavy masses, we can calculate their intensities.
   // First, need to determine the start & end scan of the
   // entire peptide peak by smoothing MS profile.
   FILTER_MS(dLightMS_, *dLightFilteredMS_, iStart, iEnd);
   FILTER_MS(dHeavyMS_, *dHeavyFilteredMS_, iStart, iEnd);

   // Starting from the start and end scans read from .out
   // files, need to see the real start/end scan of eluting
   // peptide by looking at smoothed/filtered MS profile.

   // Get peptide start & end scans
   iLightStartScan = piReverseSequentialScan[pInput_.iFirstScan];
   iLightEndScan = piReverseSequentialScan[pInput_.iLastScan];

   // Backtrack to last MS scan for start
   for (ctScan = iLightStartScan; ctScan >= pInput_.iAnalysisFirstScan; ctScan--)
   {
      const struct ScanHeaderStruct* scanHeader = cached_ramp_readHeader(fp_, index_[piSequentialScan[ctScan]]);

      if (scanHeader->msLevel == 1)
         break;
   }

   iLightStartScan = ctScan;

   // Backtrack to last MS scan for end
   for (ctScan = iLightEndScan; ctScan >= pInput_.iAnalysisFirstScan; ctScan--)
   {
     const struct ScanHeaderStruct* scanHeader = cached_ramp_readHeader(fp_, index_[piSequentialScan[ctScan]]);;

      if (scanHeader->msLevel == 1)
         break;
   }

   iLightEndScan = ctScan;

   if (iLightStartScan < pInput_.iAnalysisFirstScan)
      iLightStartScan = pInput_.iAnalysisFirstScan;
   if (iLightEndScan > pInput_.iAnalysisLastScan)
      iLightEndScan = pInput_.iAnalysisLastScan;

   iHeavyStartScan = iLightStartScan;
   iHeavyEndScan = iLightEndScan;

   if (!(pInput_.bUseSameScanRange))    // don't use same scan range
   {
      if (bLightPeptide)   // try and force offset H before L
      {
         iHeavyStartScan -= 10;
         iHeavyEndScan -= 10;
      }
      else
      {
         iLightStartScan += 10;
         iLightEndScan += 10;
      }

      FIND_ENDPOINTS(*dLightFilteredMS_, iLightStartScan, iLightEndScan,
                     pInput_.iAnalysisFirstScan, pInput_.iAnalysisLastScan);
      FIND_ENDPOINTS(*dHeavyFilteredMS_, iHeavyStartScan, iHeavyEndScan,
                     pInput_.iAnalysisFirstScan, pInput_.iAnalysisLastScan);					 

      // Make sure not at zero point since using filtered data to find endpoints
      if (dLightMS_[iLightStartScan] != 0.0)
         iLightStartScan++;
      if (dLightMS_[iLightEndScan] == 0.0)
         iLightEndScan--;
      if (dHeavyMS_[iHeavyStartScan] != 0.0)
         iHeavyStartScan++;
      if (dHeavyMS_[iHeavyEndScan] == 0.0)
         iHeavyEndScan--;
   }
   else if (pInput_.bUseFixedScanRange == 1)  // take +- .iFixedScanRange scans from ID scan
   {
      if (bLightPeptide)
      {
         iLightStartScan -= pInput_.iFixedScanRange;
         if (iLightStartScan < pInput_.iAnalysisFirstScan)
            iLightStartScan = pInput_.iAnalysisFirstScan;

         iLightEndScan += pInput_.iFixedScanRange;
         if (iLightEndScan > pInput_.iAnalysisLastScan)
            iLightEndScan = pInput_.iAnalysisLastScan;

         iHeavyStartScan = iLightStartScan;
         iHeavyEndScan = iLightEndScan;
      }
      else
      {
         iHeavyStartScan -= pInput_.iFixedScanRange;
         if (iHeavyStartScan < pInput_.iAnalysisFirstScan)
            iHeavyStartScan = pInput_.iAnalysisFirstScan;

         iHeavyEndScan += pInput_.iFixedScanRange;
         if (iHeavyEndScan > pInput_.iAnalysisLastScan)
            iHeavyEndScan = pInput_.iAnalysisLastScan;

         iLightStartScan = iHeavyStartScan;
         iLightEndScan = iHeavyEndScan;
      }
   }
   else if (pInput_.bUseFixedScanRange == 2)  // take +- .iFixedScanRange scans from apex
   {
      if (bLightPeptide)
      {
         FIND_ENDPOINTS_FIX(*dLightFilteredMS_, iLightStartScan,
                            iLightEndScan, pInput_.iAnalysisFirstScan,
                            pInput_.iAnalysisLastScan);

         iHeavyStartScan = iLightStartScan;
         iHeavyEndScan = iLightEndScan;
      }
      else
      {
         FIND_ENDPOINTS_FIX(*dHeavyFilteredMS_, iHeavyStartScan,
                            iHeavyEndScan, pInput_.iAnalysisFirstScan,
                            pInput_.iAnalysisLastScan);

         iLightStartScan = iHeavyStartScan;
         iLightEndScan = iHeavyEndScan;
      }
   }
   else if (bLightPeptide)      // light identified, use same scan range
   {
      FIND_ENDPOINTS(*dLightFilteredMS_, iLightStartScan, iLightEndScan,
                     pInput_.iAnalysisFirstScan, pInput_.iAnalysisLastScan);

      if (dLightMS_[iLightStartScan] != 0.0)
         iLightStartScan++;
      if (dLightMS_[iLightEndScan] == 0.0)
         iLightEndScan--;

      iHeavyStartScan = iLightStartScan;
      iHeavyEndScan = iLightEndScan;
   }
   else                         // heavy identified, use same scan range
   {
      FIND_ENDPOINTS(*dHeavyFilteredMS_, iHeavyStartScan, iHeavyEndScan,
                     pInput_.iAnalysisFirstScan, pInput_.iAnalysisLastScan);

      if (dHeavyMS_[iHeavyStartScan] != 0.0)
         iHeavyStartScan++;

      if (dHeavyMS_[iHeavyEndScan] == 0.0)
         iHeavyEndScan--;

      iLightStartScan = iHeavyStartScan;
      iLightEndScan = iHeavyEndScan;
   }

   // Make sure the end point did not overrun the buffered scan range.
   // If they did increase buffer range for this scan and repeat analysis.
   if ((iLightStartScan < iStart) || (iHeavyStartScan < iStart) || (iLightEndScan > iEnd) || (iHeavyEndScan > iEnd))
   {
      iScanLocalBuffer += 500;
      this->XPRESS_ANALYSIS(bLightPeptide, iScanLocalBuffer);
      //cout << "Upping buffer for scan: " << pInput_.iFirstScan << endl;
      return;
   }

   // Not needed for XPRESS but this is option is being added for external user
   // Force start and end scans to point to MS1 scan
   if (pInput_.bForceMS1Endpoints)
   {
      for (ctScan = iLightEndScan; ctScan >= pInput_.iAnalysisFirstScan; ctScan--)
      {
        const struct ScanHeaderStruct* scanHeader = cached_ramp_readHeader(fp_, index_[piSequentialScan[ctScan]]);;
   
         if (scanHeader->msLevel == 1)
            break;
      }
      iLightEndScan = ctScan;
   
      for (ctScan = iHeavyEndScan; ctScan >= pInput_.iAnalysisFirstScan; ctScan--)
      {
        const struct ScanHeaderStruct* scanHeader = cached_ramp_readHeader(fp_, index_[piSequentialScan[ctScan]]);;
   
         if (scanHeader->msLevel == 1)
            break;
      }
      iHeavyEndScan = ctScan;
   
      for (ctScan = iLightStartScan; ctScan <= pInput_.iAnalysisLastScan; ctScan++)
      {
        const struct ScanHeaderStruct* scanHeader = cached_ramp_readHeader(fp_, index_[piSequentialScan[ctScan]]);;

         if (scanHeader->msLevel == 1)
            break;
      }
      iLightStartScan = ctScan;
   
      for (ctScan = iHeavyStartScan; ctScan <= pInput_.iAnalysisLastScan; ctScan++)
      {
        const struct ScanHeaderStruct* scanHeader = cached_ramp_readHeader(fp_, index_[piSequentialScan[ctScan]]);;

         if (scanHeader->msLevel == 1)
            break;
      }
      iHeavyStartScan = ctScan;
   }


   dLightQuanValue = 0.0;
   dHeavyQuanValue = 0.0;

   pXpress_.dLightIntensity = 0.0;     // intensity of most intense peak
   pXpress_.dLightIntensityRT = 0.0;   // retentionTime of most intense peak
   pXpress_.iLightIntensityScan = 0;   // scan number of most intense peak
   pXpress_.dHeavyIntensity = 0.0;
   pXpress_.dHeavyIntensityRT = 0.0;
   pXpress_.iHeavyIntensityScan = 0;

   iNumChromatogramPoints=0;
   for (ii = iLightStartScan; ii <= iLightEndScan; ii++) {
     dLightQuanValue += dLightMS_[ii];
     if (bLightPeptide && dLightMS_[ii]>0.0)
       iNumChromatogramPoints++;

     if (dLightMS_[ii] > pXpress_.dLightIntensity) {
       pXpress_.dLightIntensity = dLightMS_[ii];
       pXpress_.iLightIntensityScan = ii;
     }
   }
   if (dLightQuanValue>0) {
     const struct ScanHeaderStruct* pHeader = cached_ramp_readHeader(fp_, index_[piSequentialScan[pXpress_.iLightIntensityScan]]);;
     pXpress_.dLightIntensityRT = pHeader->retentionTime;
   }

   for (ii = iHeavyStartScan; ii <= iHeavyEndScan; ii++) {
     dHeavyQuanValue += dHeavyMS_[ii];
     if (!bLightPeptide && dHeavyMS_[ii]>0.0)
       iNumChromatogramPoints++;

     if (dHeavyMS_[ii] > pXpress_.dHeavyIntensity) {
       pXpress_.dHeavyIntensity = dHeavyMS_[ii];
       pXpress_.iHeavyIntensityScan = ii;
     }
   }
   if (dHeavyQuanValue>0) {
       const struct ScanHeaderStruct* pHeader = cached_ramp_readHeader(fp_, index_[piSequentialScan[pXpress_.iHeavyIntensityScan]]);;
       pXpress_.dHeavyIntensityRT = pHeader->retentionTime;
   }

   if (iNumChromatogramPoints < pInput_.iMinNumChromatogramPoints) {
     pXpress_.dLightArea = 0.0;
     pXpress_.dHeavyArea = 0.0;
     sprintf(pXpress_.szQuan, "?:?");

     pXpress_.dHeavyIntensity = 0.0;
     pXpress_.dLightIntensity = 0.0;
     pXpress_.dIntensityRatio = -1.0;
   }
   else {
     // todo: calculate KL here by summing dIsotopeLight and dIsotopeHeavy over scan range

     if (pInput_.bXpressLight1 == 1) {
       if (dLightQuanValue == 0.0)
	 sprintf(pXpress_.szQuan, "1:INF");
       else
	 sprintf(pXpress_.szQuan, "1:%0.2f", dHeavyQuanValue / dLightQuanValue);

       if (pXpress_.dLightIntensity == 0.0)
	 pXpress_.dIntensityRatio = 999.0;
       else if (pXpress_.dHeavyIntensity == 0.0)
	 pXpress_.dIntensityRatio = 0.001;
       else
	 pXpress_.dIntensityRatio = pXpress_.dHeavyIntensity / pXpress_.dLightIntensity;

     }
     else if (pInput_.bXpressLight1 == 2) {
       if (dHeavyQuanValue == 0.0)
	 sprintf(pXpress_.szQuan, "INF:1");
       else
	 sprintf(pXpress_.szQuan, "%0.2f:1", dLightQuanValue / dHeavyQuanValue);

       if (pXpress_.dHeavyIntensity == 0.0)
	 pXpress_.dIntensityRatio = 999.0;
       else if (pXpress_.dLightIntensity == 0.0)
	 pXpress_.dIntensityRatio = 0.001;
       else
	 pXpress_.dIntensityRatio = pXpress_.dLightIntensity / pXpress_.dHeavyIntensity;
     }
     else {
       if (dLightQuanValue == 0.0 && dHeavyQuanValue == 0.0)
	 sprintf(pXpress_.szQuan, "?:?");
       else if (dLightQuanValue > dHeavyQuanValue)
	 sprintf(pXpress_.szQuan, "1:%0.2f", dHeavyQuanValue / dLightQuanValue);
       else
	 sprintf(pXpress_.szQuan, "%0.2f:1", dLightQuanValue / dHeavyQuanValue);

       if (pXpress_.dLightIntensity == 0.0 && pXpress_.dHeavyIntensity == 0.0)
	 pXpress_.dIntensityRatio = -1.0;
       else if (pXpress_.dLightIntensity > pXpress_.dHeavyIntensity)
	 pXpress_.dIntensityRatio = pXpress_.dHeavyIntensity / pXpress_.dLightIntensity;
       else
	 pXpress_.dIntensityRatio = pXpress_.dLightIntensity / pXpress_.dHeavyIntensity;
     }

     pXpress_.dLightArea = dLightQuanValue;
     pXpress_.dHeavyArea = dHeavyQuanValue;
   }

   pXpress_.iLightFirstScan = iLightStartScan;
   pXpress_.iLightLastScan = iLightEndScan;
   pXpress_.iHeavyFirstScan = iHeavyStartScan;
   pXpress_.iHeavyLastScan = iHeavyEndScan;

}  // XPRESS_ANALYSIS 


// Reads mzXML files and get precursor ion information
void XPressPeptideParser::XPRESS_ANALYSIS_LABELFREE(int iScanLocalBuffer)
{
   validate_mzXMLfile(); // do we have a valid scan data file?

   int ii,
       ctScan,
       iLightStartScan,
       iLightEndScan,
       iStart,
       iEnd;
   int iNumChromatogramPoints;
   double dLightMass;
   const struct ScanHeaderStruct* scanHeader;

   if (pInput_.iChargeState < 1) 
   {
      printf(" Error, charge state = %d\n\n", pInput_.iChargeState);
      exit(EXIT_FAILURE);
   }
   else
   {
      dLightMass = (PROTON_MASS * (pInput_.iChargeState - 1) + pXpress_.dLightPeptideMass) / pInput_.iChargeState;
   }

   iStart = piReverseSequentialScan[pInput_.iFirstScan] - (SCAN_BUFFER + iScanLocalBuffer);
   iEnd = piReverseSequentialScan[pInput_.iLastScan] + (SCAN_BUFFER + iScanLocalBuffer);

   if (iStart < pInput_.iAnalysisFirstScan)
      iStart = pInput_.iAnalysisFirstScan;
   if (iEnd > pInput_.iAnalysisLastScan)
      iEnd = pInput_.iAnalysisLastScan;

   // Clear all values
   memset(dLightMS_ + iStart, 0, (iEnd - iStart + 1) * sizeof(double));

   // Read all MS scan values

   // Expensive loop, cache values.
   double dLightMassLower;
   double dLightMassUpper;

   if( pInput_.bPpmMassTol )
   {
      dLightMassLower = dLightMass - dLightMass * pInput_.dMassTol / 1000000;
      dLightMassUpper = dLightMass + dLightMass * pInput_.dMassTol / 1000000;
   }
   else
   {
      dLightMassLower = dLightMass - pInput_.dMassTol;
      dLightMassUpper = dLightMass + pInput_.dMassTol;
   }
   
   const RAMPREAL *pPeaks;
   RAMPREAL fMass;
   RAMPREAL fInten;
   int n;

   // only track up to first 2 isotope peaks (in addition to monoisotopic peak)
   double dC13diff = 1.00335483;

   for (n=0 ; n<=pInput_.iMinNumIsotopePeaks; n++)
      memset(dIsotopeLight[n]+iStart, 0, (iEnd - iStart +1)*sizeof(double));

   // pInput_.iMinNumIsotopePeaks only valid if pInput_.dMassTol<=0.5 otherwise
   // makes no sense to try to reconstruct individual isotope chromatograms
   // This tolerance cutoff should really scale with charge state
   if (pInput_.iMinNumIsotopePeaks<0 || (pInput_.iMinNumIsotopePeaks>=1 && pInput_.dMassTol>0.5))
      pInput_.iMinNumIsotopePeaks=0;
   if (pInput_.iMinNumIsotopePeaks>MAX_ISOTOPES)
      pInput_.iMinNumIsotopePeaks=MAX_ISOTOPES;

   // dLightMS_ contains the main chromatogram which in the case
   // of highres data (i.e. iMinNumIsotopePeaks>=1) means the summed
   // intensities of the individual isotopic peaks.
   //
   // To extend better highres MS1 support, dIsotopeLight 
   // tracks individual isotope chromatograms.

   double dSmallestMass = dLightMassLower;
   double dLargestMass = dLightMassUpper+2*dC13diff;

   for (ctScan = iStart; ctScan <= iEnd; ctScan++) {
     scanHeader = cached_ramp_readHeader(fp_, index_[piSequentialScan[ctScan]]);
     if (scanHeader->msLevel != 1)
       continue;

      // Open a scan
      pPeaks =  cached_ramp_readPeaks_const(fp_, index_[piSequentialScan[ctScan]]);
      if (pPeaks == NULL)
          continue;

      // Binary search to the heaviest peak <= dSmallestMass.
      int l = 0, r = scanHeader->peaksCount;
      while (l < r) {
	int m = l + (r - l) / 2;
	fMass = pPeaks[m*2];
	if (fMass == dSmallestMass) {
	  l = r = m;
	  break;
	}
	else if (fMass > dSmallestMass)
	  r = m - 1;
	else
	  l = m + 1;
      }
      n = l * 2;

      // Step forward to lightest peak >= dSmallestMass, if necessary.
      fMass = pPeaks[n];
      if (fMass != -1 && fMass < dSmallestMass)
	n += 2;

      // Scan peak intensities between dLightMassLower and dLightMassUpper.
      while (pPeaks[n] != -1) {
	fMass = pPeaks[n];
	n++;

	if (fMass>dLargestMass)
	  break;

	fInten = pPeaks[n];
	n++;

	for (ii=0; ii<=pInput_.iMinNumIsotopePeaks; ii++) {
	  if (dLightMassLower+ii*dC13diff <= fMass && fMass <= dLightMassUpper+ii*dC13diff) {
	    if (fInten > dIsotopeLight[ii][ctScan])
	      dIsotopeLight[ii][ctScan] = (double) fInten;
	  }
	}
      }

      for (ii=0; ii<=pInput_.iMinNumIsotopePeaks; ii++) {
	dLightMS_[ctScan] += dIsotopeLight[ii][ctScan];
      }
   }

   // Now that we have an MS profile of each of the 
   // precursor mass, we can calculate its intensities.
   // First, need to determine the start & end scan of the
   // entire peptide peak by smoothing MS profile.
   FILTER_MS(dLightMS_, *dLightFilteredMS_, iStart, iEnd);

   // Starting from the start and end scans read from .out
   // files, need to see the real start/end scan of eluting
   // peptide by looking at smoothed/filtered MS profile.

   // Get peptide start & end scans
   iLightStartScan = piReverseSequentialScan[pInput_.iFirstScan];
   iLightEndScan = piReverseSequentialScan[pInput_.iLastScan];

   // Backtrack to last MS scan for start
   for (ctScan = iLightStartScan; ctScan >= pInput_.iAnalysisFirstScan; ctScan--) {
     scanHeader = cached_ramp_readHeader(fp_, index_[piSequentialScan[ctScan]]);
     if (scanHeader->msLevel == 1)
       break;
   }

   iLightStartScan = ctScan;

   // Backtrack to last MS scan for end
   for (ctScan = iLightEndScan; ctScan >= pInput_.iAnalysisFirstScan; ctScan--) {
     scanHeader = cached_ramp_readHeader( fp_, index_[piSequentialScan[ctScan]]);
     if (scanHeader->msLevel == 1)
       break;
   }

   iLightEndScan = ctScan;

   if (iLightStartScan < pInput_.iAnalysisFirstScan)
     iLightStartScan = pInput_.iAnalysisFirstScan;
   if (iLightEndScan > pInput_.iAnalysisLastScan)
     iLightEndScan = pInput_.iAnalysisLastScan;


   FIND_ENDPOINTS(*dLightFilteredMS_, iLightStartScan, iLightEndScan,
                  pInput_.iAnalysisFirstScan, pInput_.iAnalysisLastScan);

   if (dLightMS_[iLightStartScan] != 0.0)
     iLightStartScan++;
   if (dLightMS_[iLightEndScan] == 0.0)
     iLightEndScan--;

   // Make sure the end point did not overrun the buffered scan range.
   // If they did increase buffer range for this scan and repeat analysis.
   if ((iLightStartScan < iStart) || (iLightEndScan > iEnd)) {
     iScanLocalBuffer += 100;
     this->XPRESS_ANALYSIS_LABELFREE(iScanLocalBuffer);
     //cout << "Upping buffer for scan: " << pInput_.iFirstScan << endl;
     return;
   }

   double dLightQuanValue = 0.0;

   pXpress_.dLightIntensity = 0.0;     // intensity of most intense peak
   pXpress_.dLightIntensityRT = 0.0;   // retentionTime of most intense peak
   pXpress_.iLightIntensityScan = 0;   // scan number of most intense peak

   iNumChromatogramPoints=0;
   for (ii = iLightStartScan; ii <= iLightEndScan; ii++) {
     dLightQuanValue += dLightMS_[ii];
     if (dLightMS_[ii]>0.0)
       iNumChromatogramPoints++;

     if (dLightMS_[ii] > pXpress_.dLightIntensity) {
       pXpress_.dLightIntensity = dLightMS_[ii];
       pXpress_.iLightIntensityScan = ii;
     }
   }

   if (dLightQuanValue>0) {
     scanHeader = cached_ramp_readHeader(fp_, index_[piSequentialScan[pXpress_.iLightIntensityScan]]);
     pXpress_.dLightIntensityRT = scanHeader->retentionTime;
   }
   scanHeader = cached_ramp_readHeader(fp_, index_[piSequentialScan[iLightStartScan]]);
   pXpress_.dLightFirstScanRT = scanHeader->retentionTime;
   scanHeader = cached_ramp_readHeader(fp_, index_[piSequentialScan[iLightEndScan]]);
   pXpress_.dLightLastScanRT = scanHeader->retentionTime;

   pXpress_.dLightArea = dLightQuanValue;
   pXpress_.iLightFirstScan = iLightStartScan;
   pXpress_.iLightLastScan = iLightEndScan;

} // XPRESS_ANALYSIS_LABELFREE


void XPressPeptideParser::flipRatio(char *ratio, char *flipped)
{
   // find the : and 1
   double left, right;
   if (strstr(ratio, "?:?") != NULL)
      return;
   sscanf(ratio, "%lf:%lf", &left, &right);
   if (left == 1.0)
      ;
   else if (left == 0.0)
      left = 999;
   else if (left >= 999.)
      left = 0.0;
   else
      left = 1.0 / left;
   if (right == 1.0)
      ;
   else if (right >= 999.)
      right = 0.0;
   else if (right == 0.0)
      right = 999;
   else
      right = 1.0 / right;

   if (left == 1.0)
      sprintf(flipped, "1:%0.2f", right);
   else
      sprintf(flipped, "%0.2f:1", left);
}


void XPressPeptideParser::FIND_ENDPOINTS(
   Array<double> & pdFiltered,
   int &iPepStartScan,
   int &iPepEndScan,
   int iAnalysisFirstScan,
   int iAnalysisLastScan)
{
   int i, bLeftFalling, bRightFalling;

   // Find start and end scan of peak in the MS profile

   // check to see if pdFiltered is rising or falling as iPepStartScan is decreased.
   i = 1;
   while (iPepStartScan >= 0
	  && iPepStartScan < iAnalysisLastScan
	  && iPepStartScan - i < iAnalysisLastScan
	  && iPepStartScan - i >= 0
	  && pdFiltered[iPepStartScan] == pdFiltered[iPepStartScan - i] )
   {
      i++;
   }

   if (iPepStartScan - i >= 0 && pdFiltered[iPepStartScan] > pdFiltered[iPepStartScan - i])
      bLeftFalling = TRUE;
   else
      bLeftFalling = FALSE;

   i = 1;
   while (iPepEndScan + i < iAnalysisLastScan && pdFiltered[iPepEndScan] == pdFiltered[iPepEndScan + i] )
   {
      i++;
   }
   if (pdFiltered[iPepEndScan] > pdFiltered[iPepEndScan + i])
      bRightFalling = TRUE;
   else
      bRightFalling = FALSE;

   if (bLeftFalling == FALSE && bRightFalling == FALSE)
   {
      // error in valley of 2 peaks - just sum up from start to end scan
   }
   else if (bLeftFalling == TRUE && bRightFalling == TRUE)
   {
      // at a peak ... continue down left side to find end
      while (iPepStartScan - 1 > iAnalysisFirstScan
             && pdFiltered[iPepStartScan] >= pdFiltered[iPepStartScan - 1])
      {
         (iPepStartScan)--;
      }

      // Backtrack if the end point is at a flat plateau
      while (pdFiltered[iPepStartScan] == pdFiltered[iPepStartScan + 1])
         (iPepStartScan)++;

      // at a peak ... continue down right side to find end
      while (iPepEndScan + 1 < iAnalysisLastScan
             && pdFiltered[iPepEndScan] >= pdFiltered[iPepEndScan + 1])
      {
         (iPepEndScan)++;
      }

      // Backtrack if the end point is at a flat plateau
      while (pdFiltered[iPepEndScan] == pdFiltered[iPepEndScan - 1])
         (iPepEndScan)--;
   }
   else if (bLeftFalling == TRUE && bRightFalling == FALSE)
   {
      // Walk down to the left
      while (iPepStartScan > 0
             && pdFiltered[iPepStartScan] >= pdFiltered[iPepStartScan - 1])
      {
         (iPepStartScan)--;
      }

      // Backtrack if the end point is at a flat plateau
      while (pdFiltered[iPepStartScan] == pdFiltered[iPepStartScan + 1])
         (iPepStartScan)++;

      // Need to walk over right hump and down other side of peak
      while (iPepEndScan + 1 < iAnalysisLastScan
             && pdFiltered[iPepEndScan] <= pdFiltered[iPepEndScan + 1])
      {
         (iPepEndScan)++;
      }
      while (iPepEndScan + 1 < iAnalysisLastScan
             && pdFiltered[iPepEndScan] >= pdFiltered[iPepEndScan + 1])
      {
         (iPepEndScan)++;
      }

      // Backtrack if the end point is at a flat plateau
      while (pdFiltered[iPepEndScan] == pdFiltered[iPepEndScan - 1])
         (iPepEndScan)--;
   }
   else if (bLeftFalling == FALSE && bRightFalling == TRUE)
   {
      // Need to walk over left hump and down other side of peak
      while (iPepStartScan > 0
             && pdFiltered[iPepStartScan] <= pdFiltered[iPepStartScan - 1])
      {
         (iPepStartScan)--;
      }
      while (iPepStartScan > 0
             && pdFiltered[iPepStartScan] >= pdFiltered[iPepStartScan - 1])
      {
         (iPepStartScan)--;
      }

      // Backtrack if the end point is at a flat plateau
      while (pdFiltered[iPepStartScan] == pdFiltered[iPepStartScan + 1])
         (iPepStartScan)++;

      // Walk down to the right
      while (iPepEndScan + 1 < iAnalysisLastScan
             && pdFiltered[iPepEndScan] >= pdFiltered[iPepEndScan + 1])
      {
         (iPepEndScan)++;
      }

      // Backtrack if the end point is at a flat plateau
      while (pdFiltered[iPepEndScan] == pdFiltered[iPepEndScan - 1])
         (iPepEndScan)--;
   }

   if (iPepStartScan < iAnalysisFirstScan)
      iPepStartScan = iAnalysisFirstScan;
   if (iPepStartScan > iAnalysisLastScan)
      iPepStartScan = iAnalysisLastScan;

   if (iPepEndScan < iAnalysisFirstScan)
      iPepEndScan = iAnalysisFirstScan;
   if (iPepEndScan > iAnalysisLastScan)
      iPepEndScan = iAnalysisLastScan;

   if (iPepStartScan >= iPepEndScan)
      iPepStartScan = iPepEndScan - 1;

}                               // FIND_ENDPOINTS


// Return fixed scan # from apex
void XPressPeptideParser::FIND_ENDPOINTS_FIX(
   Array<double> & pdFiltered,
   int &iPepStartScan,
   int &iPepEndScan,
   int iAnalysisFirstScan,
   int iAnalysisLastScan)
{
   int i, bLeftFalling, bRightFalling;

   // Find peak apex and go return fixed amount (5) scans each way

   // check to see if pWhichFilterd is rising or falling
   // as iPepStartScan is decreased.
   i = 1;
   while (pdFiltered[iPepStartScan] == pdFiltered[iPepStartScan - i]
          && iPepStartScan - i > 0)
      i++;
   if (pdFiltered[iPepStartScan] > pdFiltered[iPepStartScan - i])
      bLeftFalling = TRUE;
   else
      bLeftFalling = FALSE;

   i = 1;
   while (pdFiltered[iPepEndScan] == pdFiltered[iPepEndScan + i]
          && iPepStartScan + i < iAnalysisLastScan)
      i++;
   if (pdFiltered[iPepEndScan] > pdFiltered[iPepEndScan + i])
      bRightFalling = TRUE;
   else
      bRightFalling = FALSE;

   if (bLeftFalling == FALSE && bRightFalling == FALSE)
   {
      // error in valley of 2 peaks - just sum up from start to end scan
   }
   else if (bLeftFalling == TRUE && bRightFalling == TRUE)
   {
      // at a peak ... continue down left side to find end
      while (iPepStartScan - 1 > iAnalysisFirstScan
             && pdFiltered[iPepStartScan] >= pdFiltered[iPepStartScan - 1])
         (iPepStartScan)--;


      iPepEndScan = iPepStartScan + pInput_.iFixedScanRange;
      iPepStartScan -= pInput_.iFixedScanRange;
   }
   else if (bLeftFalling == TRUE && bRightFalling == FALSE)
   {
      // Need to walk over right hump and down other side of peak
      while (iPepEndScan + 1 < iAnalysisLastScan
             && pdFiltered[iPepEndScan] <= pdFiltered[iPepEndScan + 1])
         (iPepEndScan)++;

      iPepStartScan = iPepEndScan - pInput_.iFixedScanRange;
      iPepEndScan += pInput_.iFixedScanRange;
   }
   else if (bLeftFalling == FALSE && bRightFalling == TRUE)
   {
      // Need to walk over left hump and down other side of peak
      while (iPepStartScan > 0
             && pdFiltered[iPepStartScan] <= pdFiltered[iPepStartScan - 1])
         (iPepStartScan)--;

      iPepEndScan = iPepStartScan + pInput_.iFixedScanRange;
      iPepStartScan -= pInput_.iFixedScanRange;
   }

   if (iPepStartScan < iAnalysisFirstScan)
      iPepStartScan = iAnalysisFirstScan;
   if (iPepStartScan > iAnalysisLastScan)
      iPepStartScan = iAnalysisLastScan;

   if (iPepEndScan < iAnalysisFirstScan)
      iPepEndScan = iAnalysisFirstScan;
   if (iPepEndScan > iAnalysisLastScan)
      iPepEndScan = iAnalysisLastScan;

   if (iPepStartScan >= iPepEndScan)
      iPepStartScan = iPepEndScan - 1;

}                               // FIND_ENDPOINTS_FIX


// Use my standard filtering routine
void XPressPeptideParser::FILTER_MS(const double *dOrigMS, Array<double> & dFilteredMS, int iStart, int iEnd)
{
   int i;

   // 3rd order butterworth, 0.025 cutoff frequency where 1.0 corresonds to half the sample rate
   static const double a[FILTER_SIZE]={1.0,-2.8430, 2.6980, -0.8546};
   static const double b[FILTER_SIZE]={0.0000561, 0.0001682, 0.0001682, 0.0000561};

   // Add padding on either end to make sure the range we care about is
   // filtered correctly.
   int iOrigStart = iStart;
   int iOrigEnd   = iEnd;
   iStart = iStart - 50;
   if (iStart < 0)
       iStart = 0;
   iEnd = iEnd + 50;
   if (iEnd > pInput_.iScanCount)
       iEnd = pInput_.iScanCount;

   for (i = iStart-FILTER_SIZE; i < iEnd; i++)
   {
      if (i>=0)
      {
	dFilteredMS[i]=0.0;
	if (i>=iOrigStart && i<=iOrigEnd)
	  (*dTmpFilter_)[i] = dOrigMS[i];
      }
   }

   // Pass MS profile through IIR low pass filter:
   // y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
   //      - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
   for (i = iStart; i < iEnd; i++)
   {
      int ii;
      int end = FILTER_SIZE;
      if (end > i + 1)
          end = i + 1;

      dFilteredMS[i] = b[0] * (*dTmpFilter_)[i];

      for (ii = 1; ii < end; ii++)
      {
        dFilteredMS[i] += b[ii] * (*dTmpFilter_)[i - ii];
        dFilteredMS[i] -= a[ii] * dFilteredMS[i - ii];
      }
   }

   // Filtered sequence is reversed and re-filtered resulting
   // in zero-phase distortion and double the filter order.
   for (i = iStart; i < iEnd; i++)
      (*dTmpFilter_)[i] = dFilteredMS[iEnd - 1 - (i - iStart)];

   for (i = iStart; i < iEnd; i++)
   {
      int ii;
      int end = FILTER_SIZE;
      if (end > i + 1)
          end = i + 1;

      dFilteredMS[i] = b[0] * (*dTmpFilter_)[i];
      for (ii = 1; ii < end; ii++)
      {
        dFilteredMS[i] += b[ii] * (*dTmpFilter_)[i - ii];
        dFilteredMS[i] -= a[ii] * dFilteredMS[i - ii];
      }
   }

   // Filtered sequence is reversed again
   for (i = iStart; i < iEnd; i++)
     (*dTmpFilter_)[i] = dFilteredMS[iEnd - 1 - (i - iStart)];
   
   for (i = iStart; i < iEnd; i++)
     dFilteredMS[i] = (*dTmpFilter_)[i];
}                               // FILTER_MS
