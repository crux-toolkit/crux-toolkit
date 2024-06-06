#include "Out2XML.h"

using namespace mzParser;

static double dProtonMass = 1.007276;

void Out2XML::init() {
	pScanIndex = NULL;
}

Out2XML::Out2XML() {
	init();
}

Out2XML::Out2XML(char* path, int topHits, char** argv, int argc) {
  init();
  struct stat fileStat;
#ifdef _MSC_VER
  struct _finddata_t c_file;
  long hFile;
  char szFileMask[SIZE_BUF];
#else
  DIR *dirStream;
  struct dirent *dirEntry;
#endif
  char SEP = '/';
#define OPATHLEN 4086  
  char * pathName = new char [OPATHLEN];
  char * fileName = new char [OPATHLEN];
  char * paramsFile = new char [OPATHLEN];
  pepXmlPath_ = new char [OPATHLEN];
  mzXmlPath_ = new char [OPATHLEN];
  baseName_ = new char [OPATHLEN];
  baseDir_ = new char [OPATHLEN];
  pi_calc_ = NULL;
  maldi_ = False;
  write_all_ = False;
  sequestParams_ = NULL;
  int mass_type_ = 0; // DDS TODO: make enum type if this 0 for default sequest.params value, 1 for average, 2 for monoisotopic
  char* sample_enzyme = new char [128];
  strcpy(sample_enzyme, "trypsin");
  for(int k = 3; k < argc; k++) {
    if(strlen(argv[k]) > 2 && argv[k][0] == '-' && argv[k][1] == 'E') // enzyme
      strcpy(sample_enzyme, argv[k]+2);
    else if(strlen(argv[k]) > 2 && argv[k][0] == '-' && argv[k][1] == 'P') // parameters
      sequestParams_ = new SequestParams(argv[k]+2);
// DEPRECATED OPTION INACTIVE in Sequest2XML else if(strlen(argv[k]) > 2 && argv[k][0] == '-' && argv[k][1] == 'S') // parameters
// strcpy(massspec, argv[k]+2);
// DEPRECATED OPTION INACTIVE in Sequest2XML else if(strlen(argv[k]) > 2 && argv[k][0] == '-' && argv[k][1] == 'V') // schema
// strcpy(schema, argv[k]+2);
    else if(! strcmp(argv[k], "-pI"))
      pi_calc_ = new pICalculator();
    else if(! strcmp(argv[k], "-M"))
      maldi_ = True;
    else if(! strcmp(argv[k], "-m"))
       mass_type_ = 2;
    else if(! strcmp(argv[k], "-a"))
       mass_type_ = 1;
    else if(! strcmp(argv[k], "-all"))
       write_all_ = True;

  }

  numHitsReport_ = topHits; 
  if (numHitsReport_ > 10) 
    numHitsReport_ = 10;
  
  //TODO: Have to support multiple enzymes
  enzyme_ = (new ProteolyticEnzymeFactory())->getProteolyticEnzyme(sample_enzyme); 

  if (!stat(path, &fileStat))
    {
      if (S_ISDIR(fileStat.st_mode))
	{
	  // Since this is a directory name and we
	  // will be appending filenames to it do we
	  // have a separator at the end?
	  if (!isPathSeperator(*(path + strlen(path)-1)))
            {
	      // Append separator
	      sprintf(pathName, "%s%c", path, SEP);
            }
	  else
            {
	      strcpy(pathName, path);
            }
	  
	  sprintf(baseName_, "%s", pathName);
	  sprintf(baseDir_, "%s", pathName);
	  
	  baseName_[strlen(baseName_)-1] = '\0';
	
	  sprintf(pepXmlPath_, "%s%s", baseName_,get_pepxml_dot_ext());
	  rampConstructInputFileName(mzXmlPath_, OPATHLEN, baseName_);

	  //TODO: need an option for this 
	  if (sequestParams_ == NULL) {
	    sprintf(paramsFile, "%ssequest.params", baseDir_); 
	    sequestParams_ = new SequestParams(paramsFile);
	  }

	  if( (pepXmlFile_ = ogzfopen(pepXmlPath_)) == NULL) { // will write as gzip file if name ends in .gz
	    cerr << "cannot open pepXML output file " << pepXmlPath_ << endl;
	    exit(1);
	  }
	  const char *fname;
#ifdef _MSC_VER
      sprintf(szFileMask,"%s*.out",pathName);
      hFile = _findfirst(szFileMask, &c_file);
      if (hFile!=-1) {
        do {
			  fname = c_file.name;
#else
      dirStream = (DIR *) opendir(pathName);
      // Process all files in the directory
      while ((dirEntry = readdir((DIR *) dirStream)))
      {
        // If it contains ".out" process it
	      if (strstr(fname=dirEntry->d_name, ".out"))
        {
#endif
		  //struct SequestOutStruct* tmp = (SequestOutStruct *)malloc(sizeof(SequestOutStruct));
		  SequestOut* tmp = new SequestOut();
		  struct HeaderStruct * hdr = (HeaderStruct *)malloc(sizeof(HeaderStruct));
                  sprintf(fileName, "%s%s", pathName, fname);

		  struct stat s;
		  stat(fileName, &s);
		  if (s.st_size <= 0) {
		    cerr << "ERROR: Out2XML found an empty .out file: " <<  fileName << " this probably means the SEQUEST search didn't finish. Out2XML converter aborting!!!" << endl;
		    exit(1);
		  }

		  this->readOutFile(fileName, tmp, hdr, sequestParams_);
		  outFiles_.insertAtEnd(tmp);
#ifndef _MSC_VER
               }
            }
#else
               } while( _findnext( hFile, &c_file ) == 0 );
            }
#endif
	}
      else 
	{
	  //Do nothing right now
	  //TODO: do something meaningful
	}
      
    }
  else
    {
      fprintf( stderr, "unable to open directory '%s'\n", path );
      exit( 1 );
    }
  delete [] pathName;
  delete [] fileName;
  delete [] paramsFile;
}

Out2XML::~Out2XML() {
  //TODO There's more to clean up
  delete [] pepXmlPath_;
  delete [] mzXmlPath_;
  delete [] baseName_;
  delete [] baseDir_;
  delete enzyme_;
}

void Out2XML::processData() {
  writeXMLHeader();
  writeOutData();
  writeXMLFooter();
}

void Out2XML::writeOutData() {
  for (int i=0; i<outFiles_.size(); i++) {
    fprintf(pepXmlFile_, "  <spectrum_query spectrum=\"%s\"", outFiles_[i]->szBaseFileName);
    
        
    int spec_divisors[] = {-1, -1, -1};
    int index = 2;
    int k = (int)strlen(outFiles_[i]->szBaseFileName)-1;
    if(k < 6) {
      cerr << " error: length of " << outFiles_[i]->szBaseFileName << " less than 7" << endl;
      exit(1);
    }
    
    while((index+1) && k) {
      if(outFiles_[i]->szBaseFileName[k] == '.')
	spec_divisors[index--] = k;
      k--;
    }


    fprintf(pepXmlFile_, " start_scan=\"");
    Boolean scan_start = False;
    for(k = spec_divisors[0]+1; k < spec_divisors[1]; k++) {
      if(! scan_start && outFiles_[i]->szBaseFileName[k] != '0')
	scan_start = True;
      if(scan_start)
	fprintf(pepXmlFile_, "%c", outFiles_[i]->szBaseFileName[k]);
    }
    scan_start = False;
    fprintf(pepXmlFile_, "\" end_scan=\"");
    for(k = spec_divisors[1]+1; k < spec_divisors[2]; k++) {
      if(! scan_start && outFiles_[i]->szBaseFileName[k] != '0')
	scan_start = True;
      if(scan_start)
	fprintf(pepXmlFile_, "%c", outFiles_[i]->szBaseFileName[k]);
    }

    char *chg=(char*)malloc(10 * sizeof(char)); // assume a string of less than 9 characters will represent the charge state
    strncpy(chg, "1", 10); // zero-fill the rest of the array

    if (!maldi_) {
      char* newptr = strrchr(outFiles_[i]->szBaseFileName, '.')+1;
      strncpy(chg, newptr, 10);
    }

    fprintf(pepXmlFile_, "\" precursor_neutral_mass=\"%0.4f\" assumed_charge=\"%s\" index=\"%d\"", 
	    outFiles_[i]->dAMass-dProtonMass, chg, i+1);

    free(chg);
    
    /* add RT from mzXML file */
    if (pFI != NULL)
    {
       int iStartScan = -1;
       char *pStr;
       char szBaseName[256];
 
       strcpy(szBaseName, outFiles_[i]->szBaseFileName);
 
       /* read scan number encoded in .dta name */
       if ((pStr=strrchr(szBaseName, '.'))!=NULL)
       {
          *pStr = '\0';
          if ((pStr=strrchr(szBaseName, '.'))!=NULL)
          {
             *pStr = '\0';
             if ((pStr=strrchr(szBaseName, '.'))!=NULL)
                sscanf(pStr+1, "%d", &iStartScan);
          }
       }
 
       if(iStartScan != -1 && iStartScan <= iAnalysisLastScan)
       {
          struct ScanHeaderStruct scanHeader;
          readHeader(pFI, pScanIndex[iStartScan],&scanHeader);
          fprintf(pepXmlFile_, " retention_time_sec=\"%0.1f\"", scanHeader.retentionTime);
       }
    }
    fprintf(pepXmlFile_, ">\n");

    fprintf(pepXmlFile_, "    <search_result>\n");    
    
    SequestHit* tmphit = outFiles_[i]->getHitByIndex(0);
    int totHits = outFiles_[i]->getNumHits() < numHitsReport_ ? outFiles_[i]->getNumHits() : numHitsReport_;
    for (int schIdx = 0; tmphit != NULL && schIdx < totHits; schIdx++) {

      char matched[4];
      char totions[4];
      sprintf(matched, "%d", outFiles_[i]->getHitByIndex(schIdx)->iIon);
      sprintf(totions, "%d", outFiles_[i]->getHitByIndex(schIdx)->iTot);
      
      Array<Tag*>* modinfo_tags = sequestParams_->getModificationInfoTags(outFiles_[i]->getHitByIndex(schIdx)->szSubPep);
      ModificationInfo* modinfo = NULL;
      if(modinfo_tags != NULL) {
	modinfo = new ModificationInfo(modinfo_tags);
      }

      char * stripdPep = enzyme_->strip(outFiles_[i]->getHitByIndex(schIdx)->szSubPep, 1);
      char * pepWithEnds = new char [strlen(stripdPep)+5];
      sprintf(pepWithEnds, "%c.%s.%c", outFiles_[i]->getHitByIndex(schIdx)->cAA1, stripdPep, outFiles_[i]->getHitByIndex(schIdx)->cAA2);
      
      //TODO: Monoisotopic only for now need to make this an option
      double pepmass = 0;
      if (!mass_type_) {
	pepmass = ModifiedResidueMass::getModifiedPeptideMass(modinfo, stripdPep, sequestParams_->getPrecursorMonoisotopic()); 
      }
      else {
	pepmass = ModifiedResidueMass::getModifiedPeptideMass(modinfo, stripdPep, mass_type_==2); 
      }

      double massdiff = outFiles_[i]->dAMass - outFiles_[i]->getHitByIndex(schIdx)->dMass;
      
      int protCount = 1;
      if (strlen(outFiles_[i]->getHitByIndex(schIdx)->szDup) > 0)
	protCount += atoi(outFiles_[i]->getHitByIndex(schIdx)->szDup);
      if (massdiff < 0.0) {
	fprintf(pepXmlFile_, "      <search_hit hit_rank=\"%d\" peptide=\"%s\" peptide_prev_aa=\"%c\" peptide_next_aa=\"%c\" protein=\"%s\" num_tot_proteins=\"%d\" num_matched_ions=\"%s\" tot_num_ions=\"%s\" calc_neutral_pep_mass=\"%0.4f\" massdiff=\"%f\" num_tol_term=\"%d\" num_missed_cleavages=\"%d\"",  schIdx+1, stripdPep, outFiles_[i]->getHitByIndex(schIdx)->cAA1, outFiles_[i]->getHitByIndex(schIdx)->cAA2, XMLEscape(std::string(outFiles_[i]->getHitByIndex(schIdx)->szProt)).c_str(), protCount, matched, totions, outFiles_[i]->getHitByIndex(schIdx)->dMass - dProtonMass, massdiff, enzyme_->getNumTolTerm(outFiles_[i]->getHitByIndex(schIdx)->cAA1, stripdPep, outFiles_[i]->getHitByIndex(schIdx)->cAA2), enzyme_->getNumMissedCleavages(pepWithEnds));
      }
      else {
		  fprintf(pepXmlFile_, "      <search_hit hit_rank=\"%d\" peptide=\"%s\" peptide_prev_aa=\"%c\" peptide_next_aa=\"%c\" protein=\"%s\" num_tot_proteins=\"%d\" num_matched_ions=\"%s\" tot_num_ions=\"%s\" calc_neutral_pep_mass=\"%0.4f\" massdiff=\"+%f\" num_tol_term=\"%d\" num_missed_cleavages=\"%d\"", schIdx+1 , stripdPep, outFiles_[i]->getHitByIndex(schIdx)->cAA1, outFiles_[i]->getHitByIndex(schIdx)->cAA2, XMLEscape(std::string(outFiles_[i]->getHitByIndex(schIdx)->szProt)).c_str(), protCount, matched, totions, outFiles_[i]->getHitByIndex(schIdx)->dMass - dProtonMass, massdiff, enzyme_->getNumTolTerm(outFiles_[i]->getHitByIndex(schIdx)->cAA1, stripdPep, outFiles_[i]->getHitByIndex(schIdx)->cAA2), enzyme_->getNumMissedCleavages(pepWithEnds));
      }
      delete [] pepWithEnds;
      
      if (pi_calc_ != NULL) {
	fprintf(pepXmlFile_, " calc_pI=\"%0.2f\"", pi_calc_->Peptide_pI(outFiles_[i]->getHitByIndex(schIdx)->szSubPep, modinfo));
      }

      fprintf(pepXmlFile_, " is_rejected=\"%d\">\n", 0);

      if(modinfo_tags != NULL) {  
	for(int k = 0; k < modinfo_tags->length(); k++) 
	  if((*modinfo_tags)[k] != NULL) {
	    fprintf(pepXmlFile_, "        ");
	    (*modinfo_tags)[k]->writeTraditional(pepXmlFile_);
	    delete (*modinfo_tags)[k];
	  }
	delete modinfo_tags;
      }
      if(modinfo != NULL)
	delete modinfo;
      
      
      fprintf(pepXmlFile_, "        <search_score name=\"xcorr\" value=\"%0.3f\"/>\n", outFiles_[i]->getHitByIndex(schIdx)->dXC);
      if (schIdx < (outFiles_[i]->sequestHits_.size() - 1)) { //This value is not known for the last reported hit
	fprintf(pepXmlFile_, "        <search_score name=\"deltacn\" value=\"%0.3f\"/>\n", outFiles_[i]->getHitByIndex(schIdx)->dDeltCn);
	fprintf(pepXmlFile_, "        <search_score name=\"deltacnstar\" value=\"%0.3f\"/>\n", outFiles_[i]->getHitByIndex(schIdx)->dSpecialDeltCn);
      }
      else {
	fprintf(pepXmlFile_, "        <search_score name=\"deltacn\" value=\"%0.3f\"/>\n", 1.0);
	fprintf(pepXmlFile_, "        <search_score name=\"deltacnstar\" value=\"%0.3f\"/>\n", 1.0);
      }
      
      fprintf(pepXmlFile_, "        <search_score name=\"spscore\" value=\"%0.3e\"/>\n", outFiles_[i]->getHitByIndex(schIdx)->dSp);
      fprintf(pepXmlFile_, "        <search_score name=\"sprank\" value=\"%d\"/>\n", outFiles_[i]->getHitByIndex(schIdx)->iRankSp);
      
      fprintf(pepXmlFile_, "      </search_hit>\n");
    }
    fprintf(pepXmlFile_, "    </search_result>\n");
    fprintf(pepXmlFile_, "  </spectrum_query>\n");
  }

}


void Out2XML::writeXMLHeader() {
  fprintf(pepXmlFile_, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");

  fprintf(pepXmlFile_, "<?xml-stylesheet type=\"text/xsl\" href=\"%s%s\"?>\n", PEPXML_NAMESPACE, "_std.xsl");

  // for now, make summary.xml same as xml
#ifndef USING_RELATIVE_WEBSERVER_PATH 
  fprintf(pepXmlFile_, "<msms_pipeline_analysis date=\"%s\" xmlns=\"%s\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"%s %s%s\" summary_xml=\"%s\">\n", Parser::getDateTime(), PEPXML_NAMESPACE, PEPXML_NAMESPACE, getPepXML_std_xsl_web_path(), PEPXML_SCHEMA, pepXmlPath_);
#else
  char* szWebserverRoot = new char[OPATHLEN];
  char* szCommand = new char[OPATHLEN];
  char* szBuf = new char[OPATHLEN];
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
  fprintf(pepXmlFile_, "<msms_pipeline_analysis date=\"%s\" xmlns=\"%s\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"%s %s%s\" summary_xml=\"%s\">\n", 
	  Parser::getDateTime(), PEPXML_NAMESPACE, PEPXML_NAMESPACE, szWebserverRoot, 
	  PEPXML_SCHEMA, pepXmlPath_);
  delete [] szWebserverRoot;
  delete [] szBuf;  
  delete [] szCommand;
#endif
 
  InstrumentStruct* ppstruct = NULL;
  ramp_fileoffset_t indexOffset;
  if ( (pFI=rampOpenFile(mzXmlPath_))==NULL) {
    cout << " warning: cannot open mzXML file \"" << mzXmlPath_  << "\" for reading MS instrument info and scan times." << endl;
    //exit(1);
  }
  else {
    ppstruct = getInstrumentStruct(pFI);
    indexOffset = getIndexOffset(pFI);
    pScanIndex = readIndex( pFI , indexOffset, &iAnalysisLastScan );
  }

  if(ppstruct != NULL) {
    fprintf(pepXmlFile_, "<msms_run_summary base_name=\"%s\" msManufacturer=\"%s\" msModel=\"%s\" msIonization=\"%s\" msMassAnalyzer=\"%s\" msDetector=\"%s\" raw_data_type=\"%s\" raw_data=\"%s\">\n", 
	    baseName_, ppstruct->manufacturer, ppstruct->model, ppstruct->ionisation, \
	    ppstruct->analyzer, ppstruct->detector, "raw", strrchr(mzXmlPath_,'.'));
    free(ppstruct);
  }
  else // leave out mass spec info
    fprintf(pepXmlFile_, "<msms_run_summary base_name=\"%s\" raw_data_type=\"%s\" raw_data=\"%s\">\n", 
	    baseName_, "raw", strrchr(mzXmlPath_,'.'));
  
  
  enzyme_->writeTraditionalPepXMLTags(pepXmlFile_);
  
  
  Array<Tag*>* tags = sequestParams_->getSearchParamTags(baseName_, "SEQUEST");
  for(int k = 0; k < tags->length(); k++) 
    if((*tags)[k] != NULL) {
      (*tags)[k]->writeTraditional(pepXmlFile_);
      delete (*tags)[k];
    }
  delete tags;
  
}

void Out2XML::writeXMLFooter() {
  fprintf(pepXmlFile_, "</msms_run_summary>\n");
  fprintf(pepXmlFile_, "</msms_pipeline_analysis>\n");
  fclose(pepXmlFile_);

  if (pFI != NULL)
     rampCloseFile(pFI);
  if (pScanIndex != NULL) {
     free(pScanIndex);
	 pScanIndex = NULL;
  }
}


void Out2XML::readOutFile(char* fileName, SequestOut* data, struct HeaderStruct * hdr, SequestParams* seqParams) {
  FILE *ppIn;
  
  // Initialize Data Structure
  if ((ppIn = fopen(fileName, "r")) == NULL)
    {
      printf(" Error - can't read .out file %s\n\n", fileName);
      exit(1);
    }

  // Store a copy of the full filename
  strcpy(data->szFileName, fileName);
  char *pStr;

   // Initialize Data Structure

    // Create a base filename from the full filename
   int iLen = (int)strlen(fileName);
   int i;
   for (i = iLen - 1; i > 0; i--)
      if (fileName[i] == '/' || fileName[i] == '\\')
         break;

   if ((i == 0 && (fileName[0] == '/' || fileName[0] == '\\')) || i > 0)
      strcpy(data->szBaseFileName, fileName + i + 1);
   else
      strcpy(data->szBaseFileName, fileName);

   pStr = strstr(data->szBaseFileName, ".out");
   *(pStr) = '\0';
  readOutFile(ppIn, data, hdr, seqParams);
  fclose(ppIn);
}

 
void Out2XML::readOutFile(FILE* ppIn, SequestOut* data, struct HeaderStruct * hdr, SequestParams* seqParams) {

   int i,
       ii,
       iLenPep,
       iRankXC,
       iId,
       iLen,
       iSLen,
       iDiffInSequence = 0,
       iMinLength = 0,
       bTurboSequest = 0,
       bHasSF = 0,  /* has SF column, bioworks 3.2 */
       bHasP = 0,   /* has P column, bioworks 3.2 */
       bNoDeltCnYet = 1,
       bIncludeMod = 0;

   char szBuf[SIZE_BUF],
        szMod[SIZE_BUF],
        szTmp[SIZE_BUF],
        szPep[200], 
        *pStr;
   char *fgot; // make GCC4.3 happy by noting fgets return value

   // Eat lines until time is found or file is exhausted
   while (true) {
     // Grab time/date
     if (!fgets(szBuf, SIZE_BUF, ppIn)) {
       puts("Out2XML: unexpected file format, quitting");
       exit(1);
     }
     if (!strncmp(szBuf, " (M+H)+ mass =", 14)) {
       break;
     }
   }

   // Parse mass type
// fgot=fgets(szBuf, SIZE_BUF, ppIn);
   if (strstr(szBuf, "/MONO"))
      data->iMassType = 1;
   else if (strstr(szBuf, "/AVG"))
      data->iMassType = 0;

   // Grab mass type string
   iLen = (int)strlen(szBuf);
   iSLen = 0;
   for (i = iLen - 1; i > 0; i--)
   {
      if (szBuf[i] == ' ')
         break;
      if (szBuf[i] != '\n' && szBuf[i] != '\r')
         iSLen++;
   }
   if (iSLen > SIZE_PEP - 1)
      iSLen = SIZE_PEP - 1;
   strncpy(hdr->szMassType, szBuf + i + 1, iSLen);
   hdr->szMassType[iSLen] = '\0';

   // Grab aquired mass
   for (i = 1; i < iLen; i++)
   {
      if (szBuf[i] == '=')
      {
         sscanf(szBuf + i + 1, "%lf", &(data->dAMass));
         break;
      }
   }

   // Eat me
   fgot=fgets(szBuf, SIZE_BUF, ppIn);

   // Get sequence database
   fgot=fgets(szBuf, SIZE_BUF, ppIn);
   pStr = strchr(strchr(szBuf, ',') + 1, ',') + 1;
   sscanf(pStr, "%s", data->szDatabase);
   // indexed db searches will have a comma appended to this scanned word
   if (data->szDatabase[strlen(data->szDatabase) - 1] == ',')
      data->szDatabase[strlen(data->szDatabase) - 1] = '\0';

   // Why is this relevant?
   if (strstr(szBuf, "# bases ="))
      data->bNucDb = 1;

   // Eat Ion Series
   fgot=fgets(szBuf, SIZE_BUF, ppIn);

   // Eat the next line.  Looks something like:
   //
   //     display top 10/5, ion % = 0.0, CODE = 001020
   //
   fgot=fgets(szBuf, SIZE_BUF, ppIn);


   //
   // Modification lines
   //
   //     (M# +16.000) (C@ +8.000) C=545.339  {C}
   //
   //
   fgot=fgets(szBuf, SIZE_BUF, ppIn);
   if (strstr(szBuf, "match")) {
     // eat another line
     fgot=fgets(szBuf, SIZE_BUF, ppIn);
   }
   szMod[0] = '\0';

   if (strlen(szBuf) > 2)
   {
      bIncludeMod = 1;

      /* check for modifications */
      if ((pStr = strchr(szBuf, '*')))
         sscanf(pStr + 1, "%lf", &(data->dMass1));
      if ((pStr = strchr(szBuf, '#')))
         sscanf(pStr + 1, "%lf", &(data->dMass2));
      if ((pStr = strchr(szBuf, '@')))
         sscanf(pStr + 1, "%lf", &(data->dMass3));
      if ((pStr = strchr(szBuf, '^')))
         sscanf(pStr + 1, "%lf", &(data->dMass4));
      if ((pStr = strchr(szBuf, '~')))
         sscanf(pStr + 1, "%lf", &(data->dMass5));
      if ((pStr = strchr(szBuf, '$')))
         sscanf(pStr + 1, "%lf", &(data->dMass6));
      if ((pStr = strchr(szBuf, '[')))
         sscanf(pStr + 1, "%lf", &(data->dMassCT));
      if ((pStr = strchr(szBuf, ']')))
         sscanf(pStr + 1, "%lf", &(data->dMassNT));

      pStr = strrchr(szBuf, '=');
      if (pStr && (pStr = strchr(pStr, '(')))
         *pStr='\0';
      pStr = strrchr(szBuf, '=');

      /* look for static modifications */
      if (pStr)
      {
         int  i, iLen;

         if ((pStr = strrchr(szBuf, ')')))
            strcpy(szBuf, pStr + 2);

         iLen = (int)strlen(szBuf);

         for (i = 1; i < iLen; i++)
         {
            if (szBuf[i] == '=' && szBuf[i-1]!='m')  /* make sure not terminal mod */
            {
               double dMass = 0.0;

               sscanf(szBuf + i + 1, "%lf", &(data->dMass));
               sprintf(data->szMod + strlen(data->szMod), "&amp;Mass%c=%0.4f", szBuf[i - 1], data->dMass);
            }
         }
      }

      /* look for N-terminal mods */
      if ((pStr = strstr(szBuf, "+N-term=")))
      {
         double dMass = 0.0;
         sscanf(pStr+8, "%lf", &dMass);
         sprintf(data->szMod + strlen(data->szMod), "&amp;Nterm=%0.4f", dMass);
      }

      /* look for C-terminal mods */
      if ((pStr = strstr(szBuf, "+C-term=")))
      {
         double dMass = 0.0;
         sscanf(pStr+8, "%lf", &dMass);
         sprintf(data->szMod + strlen(data->szMod), "&amp;Cterm=%0.4f", dMass);
      }

      // Eat blank line preceding header
      fgot=fgets(szBuf, SIZE_BUF, ppIn);
   }

   // Load header description line
   fgot=fgets(szBuf, SIZE_BUF, ppIn);
   if (strstr(szBuf, "Id#"))
   {
      bTurboSequest = 1;
      if (strstr(szBuf, " Sf "))
         bHasSF = 1;
      if (strstr(szBuf, " P "))
         bHasP = 1;
   }


   // Eat header underline line
   fgot=fgets(szBuf, SIZE_BUF, ppIn);
   
   for (int hitCount = 0; hitCount < 10; hitCount++) {
      // Now we have the real data
      if (fgets(szBuf, SIZE_BUF, ppIn) == NULL) {
         break;
      }

      iLen = (int)strlen(szBuf);

      if (iLen>5)
      {
         /*
          * if line begins with many spaces, it is a duplicate protein
          * line that should be ignored in parsing.
          */
         if (strncmp(szBuf, "     ", 5))
         {
            // only substitute first two '/' for rank and ion columns
            int iSlashCount=0;

            for (i = 0; i < iLen; i++)
            {
               if (szBuf[i] == '/' && iSlashCount<2)
               {
                  szTmp[i] = ' ';
                  iSlashCount++;
               }
               else
                  szTmp[i] = szBuf[i];
            }
            SequestHit* hitData = new SequestHit();
 
            if (bTurboSequest)
            {
               char szSF_P[500];

               if (bHasSF && bHasP)
               {
                  sscanf(szTmp, "%d. %d %d %d %lf %lf %lf %lf %s %s %d %d %s %s %s",
                        &(hitData->iRank), &iRankXC, &(hitData->iRankSp), &iId,
                        &(hitData->dMass), &(hitData->dDeltCn), &(hitData->dXC), &(hitData->dSp), szSF_P, szSF_P,
                       &(hitData->iIon), &(hitData->iTot), hitData->szProt, hitData->szDup, szPep);
               }
               else if (bHasSF || bHasP)
               {
                  sscanf(szTmp, "%d. %d %d %d %lf %lf %lf %lf %s %d %d %s %s %s",
                       &(hitData->iRank), &iRankXC, &(hitData->iRankSp), &iId,
                       &(hitData->dMass), &(hitData->dDeltCn), &(hitData->dXC), &(hitData->dSp), szSF_P,
                       &(hitData->iIon), &(hitData->iTot), hitData->szProt, hitData->szDup, szPep);
               }
               else
               {
                  sscanf(szTmp, "%d. %d %d %d %lf %lf %lf %lf %d %d %s %s %s",
                       &(hitData->iRank), &iRankXC, &(hitData->iRankSp), &iId,
                       &(hitData->dMass), &(hitData->dDeltCn), &(hitData->dXC), &(hitData->dSp),
                       &(hitData->iIon), &(hitData->iTot), hitData->szProt, hitData->szDup, szPep);
               }
            }
            else
            {
               sscanf(szTmp, "%d. %d %d %lf %lf %lf %lf %d %d %s %s %s",
                     &(hitData->iRank), &iRankXC, &(hitData->iRankSp),
                     &(hitData->dMass), &(hitData->dDeltCn), &(hitData->dXC), &(hitData->dSp),
                     &(hitData->iIon), &(hitData->iTot), hitData->szProt, hitData->szDup, szPep);
            }
      
            if (hitData->szDup[0] != '+')
            {
               strcpy(szPep, hitData->szDup);
               hitData->szDup[0] = '\0';
            }

            /*
             * check variant of n/c-term mod where terminal mod is indicated in separating char
             * nterm: -]VATM*PAPVR.A
             * cterm: K.VATM*PAPVR[A
             */

            if (szPep[1]==']' || szPep[strlen(szPep)-2]=='[')
            {
               int bNtermMod=0;
               int bCtermMod=0;

               if (szPep[1]==']')
                  bNtermMod=1;
               if (szPep[strlen(szPep)-2]=='[')
                  bCtermMod=1;

               ii=0;
               for (i=2; i <(int)strlen(szPep)-2; i++)
               {
                  hitData->szSubPep[ii] = szPep[i];
                  if (i==2 && bNtermMod)
                  {
                     ii++;
                     hitData->szSubPep[ii] = ']';
                  }
                  if (i==strlen(szPep)-3 && bCtermMod)
                  {
                     ii++;
                     hitData->szSubPep[ii] = '[';
                  }
                  ii++;
               }
               hitData->szSubPep[ii] = '\0';
            }
            else
            {
               strcpy(hitData->szSubPep, szPep + 2);
               hitData->szSubPep[strlen(hitData->szSubPep) - 2] = '\0';
            }

            hitData->cAA1 = szPep[0];
            hitData->cAA2 = szPep[strlen(szPep) - 1];

            iLenPep = (int)strlen(hitData->szSubPep);
            ii = 0;
            for (i = 0; i < iLenPep; i++)
            {
               if (!isalpha(hitData->szSubPep[i]))
               {
                  if (hitData->szSubPep[i] == '*')
                     hitData->szDSite[ii - 1] = '1';
                  else if (hitData->szSubPep[i] == '#')
                     hitData->szDSite[ii - 1] = '2';
                  else if (hitData->szSubPep[i] == '@')
                     hitData->szDSite[ii - 1] = '3';
                  else if (hitData->szSubPep[i] == '^')
                     hitData->szDSite[ii - 1] = '4';
                  else if (hitData->szSubPep[i] == '~')
                     hitData->szDSite[ii - 1] = '5';
                  else if (hitData->szSubPep[i] == '$')
                     hitData->szDSite[ii - 1] = '6';
                  else if (hitData->szSubPep[i] == '[')
                     hitData->szDSite[ii - 1] = 'c';
                  else if (hitData->szSubPep[i] == ']')
                     hitData->szDSite[ii - 1] = 'n';
               }
               else
               {
                  hitData->szPlainPep[ii] = hitData->szSubPep[i];
                  hitData->szDSite[ii] = '0';
                  ii++;
               }
            }
            hitData->szPlainPep[ii] = '\0';
            hitData->szDSite[ii] = '\0';
            if (!write_all_ && strchr(hitData->szPlainPep, 'X') != NULL)
            {
               cerr << "WARNING: Peptide " << hitData->szPlainPep << " contains an X, skipping (use -all option to override) ..." << endl;
               delete hitData;
            }
            else
            {
               data->insertNextHit(hitData);
            }
         }
         else
         {
            hitCount--;
         }
      }
      else
      {
         break;
      }
   }
   data->getDeltaCn();
}
