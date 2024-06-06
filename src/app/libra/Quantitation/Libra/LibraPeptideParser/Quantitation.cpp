#include "Common/sysdepend.h"
#include "Quantitation.hpp"

using namespace mzParser;

// Explicit typedef for VC++ 6.0
typedef map<int, int> mapII;

/** constructor needs existing condition object pointer
 * and name of an mzML file.  It sets write to outfile and
 * print to std out as true for defaults
 * @param pointer to condition object
 * @param name of mzXML file
 */
Quantitation::Quantitation(LibraConditionHandler* pCond, char* mzXMLFileName, RAMPFILE* FI) {
  if (dbug) cout << "Quantitation::Quantitation(LibraConditionHandler*,char*)" << endl;

  reset();
  pLibraConditionHandler = pCond;
  mzXMLFile = mzXMLFileName;
  writeToOutFile = true;
  printToStdOut = true;
  pFI = FI;
  // open mzXML file
  //if ( (pFI = rampOpenFile(mzXMLFile)) == NULL) {
  //  printf("could not open input file %s\n", mzXMLFile);
  //  exit(1);
  //}
}


Quantitation::~Quantitation() {
  //rampCloseFile(pFI);
}


void Quantitation::setWriteToOutFile(bool b) {
  writeToOutFile = b;
}


void Quantitation::setPrintToStdOut(bool b) {
  printToStdOut = b;
}


void Quantitation::setRT(double RT) {
  m_RT = RT;
}


vff Quantitation::getMaxima() {
  return m_maxima;
}


int Quantitation::getCentroidingIteration() {
  return pLibraConditionHandler->getNumCentroidingIterations();
}


/********************************************************
 The centroid m/z position is computed using standard
 arithmetical averaging
*********************************************************/
int Quantitation::centroid(int peaksNum) {
  bool isFirst = true;
  int n;
  int flag = 0;
  int peaksCount = 0;
  double delta = 0.;
  double maxDelta = 0.5;
  double tempCentroidMz = 0.;
  double centroidInt = 0.;

  double tolerance = pLibraConditionHandler->getTolerance();
  vff mass_intensity_vff =  m_centroidPeaks;
  m_centroidPeaks.clear();

  int zpOffsetIntensity = 0;
  mapII intMap;
  mapII::iterator iter;

  // get mode of intensities in m/z range of interest (estimate of bias level in spectrum)
  for(n = 0 ; (n < peaksNum) &&
	(mass_intensity_vff[n].first < (pLibraConditionHandler->m_mass.back() + (tolerance + 1)) );
      n++) {
    int count = 0;

    iter = intMap.find(n);

    if (iter != intMap.end())
      count = iter->second;

    count++;

    intMap.insert(make_pair(n,count));

    if (count > zpOffsetIntensity)
      zpOffsetIntensity = count;
  }

  vf::iterator massPos;

  // find reagent profiles: ...zpOffsetIntensity|profile|zpOffsetIntensity...
  for(massPos = pLibraConditionHandler->m_mass.begin();
      massPos < pLibraConditionHandler->m_mass.end(); massPos++) {

    int withinReagentProfile = 0;
    tempCentroidMz = 0.;

    for(n = 0 ; (n < peaksNum) &&
	  (mass_intensity_vff[n].first < (pLibraConditionHandler->m_mass.back() + (tolerance + 1)) );
	n++) {     
      // Until mz of (last isotope + 1) or end of peaklist
      // n+2 since we start from zero and we are going to read n+1

      double mz = mass_intensity_vff[n].first;
      double intensity = mass_intensity_vff[n].second;
      double diffFromReagentLine = mz - *massPos;

      //what's the abs() function in C++?
      if (diffFromReagentLine < 0.)
	diffFromReagentLine = -1.*diffFromReagentLine;

      // if current m/z is within tolerance of a reagent line, and has an intensity
      // above zpOffsetIntensity:
      if ((diffFromReagentLine < tolerance) && (intensity > zpOffsetIntensity)) {
	tempCentroidMz += mz;
	centroidInt += intensity;
	withinReagentProfile = 1;
	peaksCount++;
      }
      else {
	withinReagentProfile = 0;
	tempCentroidMz = 0.;
	centroidInt = 0.;
	peaksCount = 0;
      }

      // if tempCentroidMz is not zero, store values
      if (withinReagentProfile == 1) {
	m_massInt.first = tempCentroidMz / peaksCount;
	m_massInt.second = centroidInt;
	m_centroidPeaks.push_back(m_massInt);
      }
    }
  }

/*
  for( pos = m_centroidPeaks.begin() ; pos < m_centroidPeaks.end() ; pos++ )
  {
  double mz = pos->first;
  double intensity = pos->second;
  cout << "--- mz: " << mz << "  inten: " << intensity << endl;
  }
*/

  return 0;
}


/********************************************************
 This function will find the local maxima in the
 specified m/z region (plus minus TOLERANCE)
 m_nrReagent must be set before calling the function
 Return: maxima are stored in the m_maxima vector as 
         mz-Int pairs.
*********************************************************/
int Quantitation::findMaxima() {
  int    matchedMass = 0;
  double intensity = 0;
  double mz;
  double tempMaxIntensity = 0;
  double tempMaxMz = 0;
  vf::iterator massPos;

  for(massPos = pLibraConditionHandler->m_mass.begin() ; 
      massPos < pLibraConditionHandler->m_mass.end(); massPos++) {
    // Find local maxima around each isotopic mass entered by the user
    for(pos = m_centroidPeaks.begin() ; pos < m_centroidPeaks.end() ; pos++) {
      mz = pos->first;
      intensity = pos->second;

      if( (mz > *massPos - pLibraConditionHandler->m_tolerance) &&
	  (mz < *massPos + pLibraConditionHandler->m_tolerance) ) {
	if(intensity > tempMaxIntensity) {
	  tempMaxIntensity = intensity;
	  tempMaxMz = mz;
	}
      }
    }

    if(tempMaxIntensity > 0) {
      matchedMass++;
    }
 
    m_massInt.first = tempMaxMz > 0 ? tempMaxMz : *massPos;
    m_massInt.second = tempMaxIntensity + pLibraConditionHandler->m_pseudoCount;
    m_maxima.push_back(m_massInt);
    m_maxima_absolute.push_back(m_massInt);
    tempMaxIntensity = 0;
    tempMaxMz = 0;
  }

  return matchedMass;
}


/**
* Intensity weighted means of m/z values within tolerance difference 
* of reagent masses are determined.  The m/z reagent profiles are
* integrated.  The m/z and intensity pairs are stored in 
* m_centroidPeaks.
* @param peaksNum the number of peaks in the spectrum
*/
int Quantitation::centroidWeighted(int peaksNum) {
  bool isFirst = true;
  int n;
  int flag = 0;
  int peaksCount = 0;
  double delta = 0;
  double maxDelta = 0.5;
  double tempCentroidMz = 0;
  double centroidInt = 0;
  double tolerance = pLibraConditionHandler->getTolerance();
  vff mass_intensity_vff =  m_centroidPeaks;
  m_centroidPeaks.clear();

  int zpOffsetIntensity = 0;
  mapII intMap;
  mapII::iterator iter;

  // get mode of intensities in m/z range of interest (estimate of bias level in spectrum)
  for(n = 0 ; (n < peaksNum) && 
	(mass_intensity_vff[n].first < (pLibraConditionHandler->m_mass.back() + (tolerance + 1)) ); 
      n++) {
    int count = 0;

    iter = intMap.find(n);

    if (iter != intMap.end())
      count = iter->second;

    count++;

    intMap.insert(make_pair(n, count));

    if (count > zpOffsetIntensity)
      zpOffsetIntensity = count;
  }

  vf::iterator massPos;

  // find reagent profiles: ...zpOffsetIntensity|profile|zpOffsetIntensity...
  for(massPos = pLibraConditionHandler->m_mass.begin();
      massPos < pLibraConditionHandler->m_mass.end(); massPos++) {

    int withinReagentProfile = 0;
    tempCentroidMz = 0.;

    for(n = 0 ; (n < peaksNum) && 
	  (mass_intensity_vff[n].first < (pLibraConditionHandler->m_mass.back() + (tolerance + 1)) ); 
        n++) {
      // Until mz of (last isotope + 1) or end of peaklist
      // n+2 since we start from zero and we are going to read n+1

      double mz = mass_intensity_vff[n].first;
      double intensity = mass_intensity_vff[n].second;
      double diffFromReagentLine = mz - *massPos;

      //what's the abs() function in C++?
      if (diffFromReagentLine < 0.)
	diffFromReagentLine = -1.*diffFromReagentLine;

      // if current m/z is within tolerance of a reagent line, and has an intensity
      // above zpOffsetIntensity:
      if ((diffFromReagentLine < tolerance) && (intensity > zpOffsetIntensity)) {
	tempCentroidMz += (mz * intensity);
	centroidInt += intensity;
	withinReagentProfile = 1;
      }
      else {
	withinReagentProfile = 0;
	tempCentroidMz = 0.;
	centroidInt = 0.;
      }
      // if tempCentroidMz is not zero, store values
      if (withinReagentProfile == 1) {
	m_massInt.first = tempCentroidMz / centroidInt;
	m_massInt.second = centroidInt;
	m_centroidPeaks.push_back(m_massInt);
      }
    }
  }

  /*
    for( pos = m_centroidPeaks.begin() ; pos < m_centroidPeaks.end() ; pos++ )
    {
    double mz = pos->first;
    double intensity = pos->second;
    cout << "--- mz: " << mz << "  inten: " << intensity << endl;
    }
  */

  return 0;
}


/********************************************************
 This function save the values of the profile and 
 centroid data in two files. It than calls gnuplot and
 generate a superimposed image of the two.
*********************************************************/
int Quantitation::plotCentroid(int nPeaks) {
  int n;
  std::ofstream fout;
  char trash[10];

  fout.open("centroid.txt");
  if(!fout.good()) {
    cerr << "Could not open output file\n";
    exit(2);
  }

  for(pos = m_centroidPeaks.begin() ; pos < m_centroidPeaks.end() ; pos++) {
    fout << pos->first << " " << pos->second << endl;
  }
  
  fout.close();

  fout.open("profile.txt");
  if(!fout.good()) {
    cerr << "Could not open output file\n";
    exit(2);
  }
  
  for(n = 0 ; n < nPeaks ; n++) {
    fout << m_centroidPeaks[n].first << " " << m_centroidPeaks[n].second << endl;
  }
  
  fout.close();
 
  verified_system(GNUPLOT_BINARY" -persist plot");

  cin.read(trash,1);
  return 0;
}



/********************************************************
 This function applies an isotope correction to the
 intensities in the m_maxima vector
*********************************************************/
int Quantitation::isotopeCorrection() {
  int n = 0;
  int m;
  double isotopicCorrection = 0.;

  vff::iterator posAbsolute_iter =  m_maxima_absolute.begin();
  
  for(pos = m_maxima.begin() ; pos < m_maxima.end() ; pos++) {
    // Calculate dot product
    m = 0;
    for(pos2 = m_maxima.begin() ; pos2 < m_maxima.end() ;  pos2++ , m++) {
      isotopicCorrection += (pos2->second * pLibraConditionHandler->m_massIsotopes[m][n]);
    }

    if (dbug) cout << "channel " << n + 1 << " (m/z = " << pos->first
		   << "),  intensity before: " << pos->second ;

    if (pos->second - isotopicCorrection > pLibraConditionHandler->m_pseudoCount)
      pos->second = pos->second - isotopicCorrection;

    if (dbug) cout <<",  intensity after: "<< pos->second << endl;


    if (posAbsolute_iter->second - isotopicCorrection > pLibraConditionHandler->m_pseudoCount)
      posAbsolute_iter->second = posAbsolute_iter->second - isotopicCorrection;

    // null out negative values, as the lines were not detected
    if (pos->second < 0.)
      pos->second = 0.;

    if (posAbsolute_iter->second < 0.)
      posAbsolute_iter->second = 0.;

    posAbsolute_iter++;

    isotopicCorrection = 0.;
    n++;
  }

  return 0;
}


/********************************************************
 This function prints the value in the m_maxima vector
*********************************************************/
int Quantitation::printMaxima() {
  for(pos = m_maxima.begin() ; pos < m_maxima.end() ; pos++) {
    cout << pos->first << " " << pos->second << "\t";
  }
  cout << endl;

  return 0;
}


/**
* Normalize intensities of the reagent mass lines using
* value set in condition file.  
*   -2: for normalize against TIC (not recommended)
*   -1: for normalize against the sum of the reagent lines
*    0: for normalize against most intense peak (not recommended)
*    1: for normalize against 1st listed reagent mass
* ...n: for normalize against nth listed reagent mass
* Default is no normalization.
* Note that the normalization in LibraProtein code , in other words,
* the final protein quantitation, does not use the normalization here.
* This normalization is just for user viewing in the interact.xml file.
*/
int Quantitation::normalizeRatio() {
  double highestIntensity = 0.;
  double sumChannels = 0.;

  switch(pLibraConditionHandler->getNormalPosition()) {
  case 0:
    // Normalize against most intense
    for(pos = m_maxima.begin() ; pos < m_maxima.end() ; pos++) {
      if(pos->second > highestIntensity)
	highestIntensity = pos->second;
    }
    for(pos = m_maxima.begin() ; pos < m_maxima.end() ; pos++) {
      pos->second = pos->second / highestIntensity;
    }
    break;

  case -1:
    // Normalize against sum of reagent channels
    for (pos = m_maxima.begin();  pos < m_maxima.end(); pos++) {
      sumChannels += (pos->second);
    }
    for (pos = m_maxima.begin();  pos < m_maxima.end(); pos++) {
      pos->second = pos->second / sumChannels;
    }
    break;

  case -2:
    // Normalize against TIC
    for(pos = m_maxima.begin() ; pos < m_maxima.end() ; pos++) {
      pos->second = (double) (pos->second / m_tic);
    }
    break;
      
  default:
    // Normalize on one isotope
    int normalizationChannel = pLibraConditionHandler->getNormalPosition()-1;

    // if < 0, no normalization
    if (normalizationChannel >= 0) {
      pos = (m_maxima.begin()) + normalizationChannel;

      double nr = pos->second;

      for (pos = m_maxima.begin();  pos < m_maxima.end(); pos++) {
	if (nr > 0)
	  pos->second = pos->second / nr;
	else
	  pos->second = 0;
      }
    }
    break;
  }
  
  return 0;
}


/** 
 * open an existing file with name based on mzXML file name root 
 * (note changes here should be carried to similar method in
 * LibraConditionHandler)
 * @param mzXML file name
 * @return 0
 */
int Quantitation::openExistingOutFile(char* inFileName) {
  int  nameLen;
  char *test, *nameStart;
  char outFileName[250];

  // Skip the path
  nameStart = inFileName;

  while( (test = strstr(nameStart, "/")) != NULL) {
    nameStart = test + 1;
  }

  // Get rid of mzXML
  if( (nameLen = (int)strlen(nameStart) ) < 250) {
    strncpy(outFileName,nameStart,(nameLen - 5));
  }
  else {
    strncpy(outFileName,nameStart,250);
  }

  strcpy(outFileName + (nameLen - 5),"aq");
  m_fout.open(outFileName, ios::app);

  return 0;
}


/********************************************************
 This function prints the ratios btw the iTRAQ isotopes
*********************************************************/
int Quantitation::printRatio(int scNum) {
  openExistingOutFile(mzXMLFile);

  if(pLibraConditionHandler->m_outputPrefs == 1) {
    m_fout << std::setw(10) << scNum << std::setw(1) << " ";
  }
  else if(pLibraConditionHandler->m_outputPrefs == 2) {
    m_fout << std::setw(10) << m_RT << std::setw(1) << " ";
  }
  else {
    cerr << "Unsupported output format. Check condition file.\n";
    exit(1);
  }
  
  for(pos = m_maxima.begin() ; pos < m_maxima.end() ; pos++) {
    m_fout.width(10);
    m_fout.precision(6);
    m_fout.flags(std::ios::fixed);  
    //m_fout << pos->first << " " << pos->second << " ";
    m_fout << pos->second << " ";
  }
  m_fout << endl;
  m_fout.close();

  return 0;
}


/********************************************************
 This function gets everything ready for the quantification
 of the next scan.
 Remember to call it.
*********************************************************/
int Quantitation::reset() {
  m_maxima.clear();
  m_maxima_absolute.clear();
  m_centroidPeaks.clear();

  return 0;
}


int Quantitation::calculateIntensities(int startScan, int endScan) {
  if (dbug) cout << "Quantitation::calculateIntensities(int,int)"
		 << " where startScan, endScan = " << startScan << " "
		 << endScan << endl;

  // set class attributes:
  startScanNum = startScan;
  endScanNum = endScan;

  int iAnalysisFirstScan;
  int iAnalysisLastScan;
  // ramp_fileoffset_t are 64-bit file offsets
  ramp_fileoffset_t  indexOffset;
  ramp_fileoffset_t  *pScanIndex;

  // Read the offset of the index using pointer to mzXML file
  indexOffset = getIndexOffset(pFI);

  // Read the scan index into a vector, and get LastScan
  pScanIndex = readIndex(pFI, indexOffset, &iAnalysisLastScan);

  iAnalysisFirstScan = 1;

  if (startScan < iAnalysisFirstScan)
    startScanNum = iAnalysisFirstScan;

  if (endScan > iAnalysisLastScan)
    endScanNum = iAnalysisLastScan;

  calculateIntensities(pFI, pScanIndex, startScanNum, endScanNum, iAnalysisLastScan);
  free(pScanIndex);

  return 0;
}


int Quantitation::calculateIntensities() {
  int iAnalysisFirstScan;
  int iAnalysisLastScan;

  // ramp_fileoffset_t are 64-bit file offsets
  ramp_fileoffset_t  indexOffset;
  ramp_fileoffset_t  *pScanIndex;

  // Read the offset of the index using pointer to mzXML file
  indexOffset = getIndexOffset(pFI);

  // Read the scan index into a vector, and get LastScan
  pScanIndex = readIndex(pFI , indexOffset, &iAnalysisLastScan);

  iAnalysisFirstScan = 1;

  if (dbug) cout<< "\nfirst scan is " << iAnalysisFirstScan
		<< "\nlast scan is "  << iAnalysisLastScan 
		<< endl;

  calculateIntensities(pFI , pScanIndex,  iAnalysisFirstScan, iAnalysisLastScan, iAnalysisLastScan);
  free(pScanIndex);

  return 0;
}


int Quantitation::calculateIntensities(RAMPFILE*  pFI, ramp_fileoffset_t* pScanIndex,
      int startScan, int endScan, int iAnalysisLastScan) {

  if (dbug) cout << "Quantitation::calculateIntensities(RAMPFILE*, ramp_fileoffset_t*, int, int) "
		 << "  where startScan, endScan = " << startScan << " " <<  endScan << endl;

  long  scanNum;
  int   numProcessedScan = 0;
  int   numNonLabeledScan = 0;

  // write conditions to outfile if flag is true
  if (writeToOutFile)
    pLibraConditionHandler->writeLibraConditionHandlerToOutFile(mzXMLFile);

  for(scanNum = startScan ; scanNum <= endScan ; scanNum++) {
    int   iCount=0;
    int   n=0;
    RAMPREAL *pPeaks;
    struct ScanHeaderStruct scanHeader;

    readHeader(pFI, pScanIndex[scanNum], &scanHeader);

    setRT( (double)scanHeader.retentionTime + 2.0);  // Skip 2 for the PT

    {
      numProcessedScan++;
      int iReporterIonScan = scanNum;

      // Following code in 'if' block was implemented for Thermo's Synchronous Precursor Selection
      // where the reporter ions are acquired in an MS3 scan.  And the MS3 scan isn't 
      // neccessarily in any sequential order.  So find the right MS3 scan to extract
      // reporter ion peaks from by looking at the filter line.
      // If no MS3 scan is found, it will revert back to MS/MS scan to read reporter peaks.
      if (pLibraConditionHandler->getReporterFromMS3()) {
        char *pStr;
        char szMatchString[512];
        char szPrec[512];

        // MS2 scan:  ITMS + c NSI r d Full ms2 425.2244@cid35.00 [112.0000-860.0000]
        // MS3 scan:  FTMS + p NSI sps d Full ms3 425.2244@cid35.00 474.3074@hcd65.00 [100.0000-500.0000]

        if ( (pStr = strstr(scanHeader.filterLine, "Full ms2 "))==NULL) {
          printf(" Error ms2 filter line: %s\n\n", scanHeader.filterLine);
          exit(1);
        }

        // encode the corresponding MS3 scan filter line text to search for
        sscanf(pStr, "%*s %*s %s", szPrec);
        sprintf(szMatchString, "Full ms3 %s", szPrec);

        for (int iTmp=scanNum ; iTmp<=iAnalysisLastScan; iTmp++) {
          readHeader(pFI, pScanIndex[iTmp], &scanHeader);

          if (scanHeader.msLevel==3 && strstr(scanHeader.filterLine, szMatchString)) {
            iReporterIonScan = iTmp; // ms/ms scan "ctScan" has corresponding ms3 scan at "iTmp"
            break;
          }

          if (iTmp == scanNum + 20) // only look at a small range of scans after current MS/MS scan
            break;
          }
        }


        // MH: This code changes the behavior of SPS/MS3 quant in Libra. Instead of defaulting back to the MS2 scan
        // (which will NOT work for TMTpro if such scans are ion trap), just report nothing - no quant was done for that scan.
        if(iReporterIonScan==0){
          numNonLabeledScan++;

          vf::iterator  massPos;

          for (massPos = pLibraConditionHandler->m_mass.begin();
            massPos < pLibraConditionHandler->m_mass.end(); massPos++) {
            m_massInt.first = *massPos;
            m_massInt.second = pLibraConditionHandler->m_pseudoCount;
            m_maxima.push_back(m_massInt);
            m_maxima_absolute.push_back(m_massInt);
          }

          if (writeToOutFile)
            printRatio(scanNum);

          goto finishingStep;
        }


        pPeaks = readPeaks(pFI, pScanIndex[iReporterIonScan]);
 
        if (dbug) cout << "scan number: "<<scanNum<< "  number of peaks: "
          <<  readPeaksCount(pFI, pScanIndex[scanNum])
          << endl;

        if (readPeaksCount(pFI, pScanIndex[scanNum]) < 1) {
          cout << "ERROR: scan number " << scanNum << " contains no peaks! "
            << "Check that the scan numbers are correct in the source file(s)." << endl
            << "Exiting..." << endl;
          exit(1);
        }

        // store scan's mass and intensity peaks in m_centroidPeaks:
        while (pPeaks[n] != -1) {
          RAMPREAL fMass;
          RAMPREAL fInten;

          fMass = pPeaks[n];
          n++;
          fInten = pPeaks[n];
          n++;

          iCount += 1;

          // copy all peaks of scan into m_centroidPeaks:
          m_massInt.first = fMass;
          m_massInt.second = fInten;
          m_centroidPeaks.push_back(m_massInt);
        }
        free(pPeaks);

        if (m_centroidPeaks[0].first < (pLibraConditionHandler->m_mass.back() + pLibraConditionHandler->getTolerance() + 1) ) {
          if (printToStdOut)
            cout << "Finished copying " << iCount << " peaks\n";
    
          if(pLibraConditionHandler->getCentroidingPref() == 1) {
            for(n = 0 ; n < getCentroidingIteration() ; n++) {
              centroid(scanHeader.peaksCount);
            }
          }
          else if(pLibraConditionHandler->getCentroidingPref() == 2) {
            centroidWeighted(scanHeader.peaksCount);
          }

          // stores m_centroidPeaks in m_maxima and m_maxima_absolute:
          if( !findMaxima() ) {
            if (printToStdOut)
              cout << "Scan " << scanNum << " is not labeled\n";

            numNonLabeledScan++;
          }
          else {
            //plotCentroid(scanHeader.peaksCount);
            if(pLibraConditionHandler->getUseIsotopicCorrection())
              isotopeCorrection();

            // normalize m_maxima if requested:
            if(pLibraConditionHandler->getIsToNormalize()) {
              if (pLibraConditionHandler->getNormalPosition() == -2) {
                m_tic = scanHeader.totIonCurrent;

              if (m_tic > 0) {
                normalizeRatio();
              }
              else {
                cerr << "Can't normalize by TIC as it isn't "
                  << " present in the scan header"<<endl;
              }
            }
            else {
              normalizeRatio();
            }
          } // end normalization
        }

        if (writeToOutFile)
          printRatio(scanNum);

        // There are 2 data structures in Quantitation.cpp
        // --> m_maxima normalized by user choice
        // --> m_maxima_absolute isn't normalized
    
	}
      else {
	if (printToStdOut)
          cout << "Scan " << scanNum << " is not labeled\n";
    
        numNonLabeledScan++;

        vf::iterator  massPos;

        for(massPos = pLibraConditionHandler->m_mass.begin(); massPos < pLibraConditionHandler->m_mass.end(); massPos++) {
          m_massInt.first = *massPos;
          m_massInt.second = pLibraConditionHandler->m_pseudoCount;
          m_maxima.push_back(m_massInt);
          m_maxima_absolute.push_back(m_massInt);
        } 

        if (writeToOutFile)
          printRatio(scanNum);
      }

finishingStep:
      int numNoSignal = 0;
      for(pos = m_maxima_absolute.begin() ; pos < m_maxima_absolute.end() ; pos++) {
	if(pos->second == pLibraConditionHandler->m_pseudoCount)
	  numNoSignal++;
      }
      if (numNoSignal > pLibraConditionHandler->m_maxPseudoChannels) {
	pos2 = m_maxima.begin();
	for(pos = m_maxima_absolute.begin() ; pos < m_maxima_absolute.end() ; pos++) {
	  pos->second = 0;
	  pos2->second = 0;
	  pos2++;
	}
      }
	
      // need to store into maps for pipeline:
      // key = scan number,   value = vector of target masses
      // key = scan number,   value = vector of absolute intensities
      vector<double> t1, t2, t3;
      vff::iterator iter, iter2;

      iter2 = m_maxima.begin();

      for (iter = m_maxima_absolute.begin(); iter < m_maxima_absolute.end(); iter++) {
        t1.push_back(iter->first);
        t2.push_back(iter->second);
        t3.push_back(iter2->second);

        iter2++;
      }

      // store into maps
      mapScanNumAndTargetMasses.insert(make_pair (scanNum, t1));
      mapScanNumAndAbsoluteIntensities.insert(make_pair (scanNum, t2));
      mapScanNumAndNormalizedIntensities.insert(make_pair (scanNum, t3));

      reset();
    }
  }

  if (printToStdOut)
    cout << numNonLabeledScan << " of " << numProcessedScan 
      << " processed scans were not labeled.\n";

  //xxxxxxx
  /*
    std::map<int, std::vector<double> >::iterator ir;
    std::vector<double> tmp;
    ir = mapScanNumAndAbsoluteIntensities.find( startScan );
    if (ir != mapScanNumAndAbsoluteIntensities.end())
    tmp = ir->second;
    for (int i = 0; i < tmp.size(); i++)
    {
    cout << "scan: " << startScan << " " << tmp[i] << endl;
    }
  */

  return 0;
}


void Quantitation::printVFFToStdout(vff m_data) {
  vff::iterator iter;

  for (iter = m_data.begin(); iter < m_data.end(); iter++) {
    cout << "FIRST: " << iter->first << "  SECOND: " << iter->second << endl;
  }
}


vector<double> Quantitation::getTargetMasses(int scanNumber) {
  std::vector<double> tmp;

  std::map<int, std::vector<double> >::iterator iter = 
    mapScanNumAndTargetMasses.find(scanNumber);

  if (iter != mapScanNumAndTargetMasses.end())
    tmp = iter->second;

  return tmp;
}


int Quantitation::getNumberOfChannels(int scanNumber) {
  std::map<int, std::vector<double> >::iterator iter;
  std::vector<double> tmp;

  iter = mapScanNumAndTargetMasses.find(scanNumber);

  if (iter != mapScanNumAndTargetMasses.end())
    tmp = iter->second;

  return (int)tmp.size();
}


vector<double> Quantitation::getAbsoluteIntensities(int scanNumber) {
  std::map<int, std::vector<double> >::iterator iter;
  std::vector<double> tmp;

  iter = mapScanNumAndAbsoluteIntensities.find(scanNumber);

  if (iter != mapScanNumAndAbsoluteIntensities.end())
    tmp = iter->second;

  return tmp;
}


vector<double> Quantitation::getNormalizedIntensities(int scanNumber) {
  std::map<int, std::vector<double> >::iterator iter;
  std::vector<double> tmp;

  iter = mapScanNumAndNormalizedIntensities.find(scanNumber);

  if (iter != mapScanNumAndNormalizedIntensities.end())
    tmp = iter->second;

  return tmp;
}


int Quantitation::getStartScanNumber() {
  return startScanNum;
}

int Quantitation::getEndScanNumber() {
  return endScanNum;
}

