/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Library or "Lesser" General Public      *
 *   License (LGPL) as published by the Free Software Foundation;          *
 *   either version 2 of the License, or (at your option) any later        *
 *   version.                                                              *
 *                                                                         *
 ***************************************************************************/

/****
 
       Date:  01/22/2004
 
    Purpose:  Use RAMP to create .dta, .mgf files from input mzXML file.
 
  Copyright:  Jimmy Eng, Copyright (c) Institute for Systems Biology, 2004.  All rights reserved.
 
   Requires:  ramp_base64.cpp ramp_base64.h ramp.cpp ramp.h from http://sourceforge.net/projects/sashimi
 
****/ 

// 9/27/06 (brendanx) - added to TPP, and added support for MGF files.
// 4/07/09 (natalie) - clairify usage; more function-oriented rewrite; add support for possibleCharges in mzXML 3.1
// 4/09/09 (natalie) - adding peak threshold option


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
//#include <iostream>
//#include <vector>
#include "Common/sysdepend.h"
#include "Common/util.h"
#include "Parsers/mzParser/mzParser.h"
#include "Common/TPPVersion.h" // contains version number, name, revision

#define MINMASS     600.0
#define MAXMASS     5000.0
#define MINPEAKS    5
#define MININTENSITY 0.01
#define MAXOUTPUTCHARGE 7
#define PRECISION_MASS 4
#define PRECISION_INT  0

// tracking charges in an internal array
#define MAXCHARGES  127

#define FIELDSEPARATOR " "

using namespace std;
using namespace mzParser;

static double dProtonMass = 1.00727646688;

 
typedef enum {
  FORMAT_UNDEF = 0,
  DTA,		// SEQUEST DTA format
  XDTA,		// X! Tandem single file DTA
  ODTA,		// OMSSA merged single file DTA
  MGF,		// Mascot Generic Format
  PKL,		// Micromass PKL
  MS2		// SEQUEST MS2 format
} FormatType;
std::string getFormatDescription(const FormatType format);
std::string getFormatExtension(const FormatType format);
FormatType getFormat(const std::string& formatStr);

  
struct OptionsStruct;

struct PeaksStruct
{
   RAMPREAL fMass;
   RAMPREAL fInten;
};


void usage(char *arg);

bool setOption(char *exe, char *arg, struct OptionsStruct *pOptions);

void convert(char *szXMLFile, const OptionsStruct& options);

double calculatePeptideMass(double precursorMZ, int charge) {
  return (precursorMZ * charge) - (charge - 1) * dProtonMass;
}

void writePeaksTSV(FILE *pfOut, const RAMPREAL *pPeaks, const OptionsStruct& options, const double precursorMz, const int precursorCharge = 0);

void outputDTA(FILE *pfOut, const std::string& strBaseName, const std::string& strOutputPath, const OptionsStruct& options, const ScanHeaderStruct& scanHeader, long curScanNum, const bool* chargeArray, const RAMPREAL *pPeaks);
void outputXDTA(FILE *pfOut, const std::string& strBaseName, const OptionsStruct& options, const ScanHeaderStruct& scanHeader, long curScanNum, const bool* chargeArray, const RAMPREAL *pPeaks);
void outputODTA(FILE *pfOut, const std::string& strBaseName, const OptionsStruct& options, const ScanHeaderStruct& scanHeader, long curScanNum, const bool* chargeArray, const RAMPREAL *pPeaks, long& lScanCount);
void outputMGF(FILE *pfOut, const std::string& strBaseName, const OptionsStruct& options, const ScanHeaderStruct& scanHeader, long curScanNum, const bool* chargeArray, const RAMPREAL *pPeaks);
void outputPKL(FILE *pfOut, const OptionsStruct& options, const ScanHeaderStruct& scanHeader, long curScanNum, const bool* chargeArray, const RAMPREAL *pPeaks);
void outputMS2(FILE *pfOut, const OptionsStruct& options, const ScanHeaderStruct& scanHeader, long curScanNum, const bool* chargeArray, const RAMPREAL *pPeaks);

int struct_cmp_mass(const void *a, const void *b);
int struct_cmp_inten(const void *a, const void *b);

bool isNearChargeReducedPrecursors(double mz, double precursorMz, int precursorCharge);
bool isITRAQReporterPeaks(double mz);
bool isTMTReporterPeaks(double mz);

struct OptionsStruct
{
   FormatType  format;
   int    iFirstScan;
   unsigned long    iLastScan;
   int    iMinPeakCount;
   int    iCharge;
   int    iChargeLast;
   int    iMSLevel;
   int    iMSLevelLast;
   int    iPrecisionMass;
   int    iPrecisionInt;
   int    iMaxPrintPeaks;
   double dMinMass;
   double dMaxMass;
   double dMinPeakThreshold;
   double dMinPeakDtaFilter;
   double dMaxOutputCharge;
   char*  outputPath;
   char*  activationMethod;
   bool   useChargeOverrides;
   bool   removeChargeReducedPrecursors;
   bool   removeITRAQReporterPeaks;
   bool   removeTMTReporterPeaks;

   OptionsStruct()
   {
      // defaults
      format = DTA;
      iFirstScan = 0;
      iLastScan = 999999999;
      iCharge = iChargeLast = 0;
      iMSLevel = 2;
      iMSLevelLast = 2;
      iMinPeakCount = MINPEAKS;
      iPrecisionMass = PRECISION_MASS;
      iPrecisionInt = PRECISION_INT;
      iMaxPrintPeaks = 0;
      dMinMass = MINMASS;
      dMaxMass = MAXMASS;
      dMinPeakThreshold = MININTENSITY;
      dMaxOutputCharge = MAXOUTPUTCHARGE;
      dMinPeakDtaFilter = 0;
      outputPath = NULL;
      activationMethod = NULL;
      useChargeOverrides = false;
      removeChargeReducedPrecursors = false;
      removeITRAQReporterPeaks = false;
      removeTMTReporterPeaks = false;
   }

   bool acceptMSLevel(int level) const
   {
      return (iMSLevel <= level && level <= iMSLevelLast);
   }

   bool acceptActivationMethod(char* method) const
   {
      // HENRY -- changed interpretation of -A option
      // When a non-empty -A is given, scans with no activationMethod will NOT be accepted
      // When an empty -A is given (nothing following the -A), ONLY scans with no activationMethod will be accepted
      // [Previous behavior was: all scans with no activationMethod will be accepted regardless of -A specification
      // Previous behavior can be obtained by simply NOT specifying any -A option at the first place!]
      return (activationMethod==NULL || !strcmp(method, activationMethod));
   }

   bool acceptPeaksCount(int peaksCount) const
   {
      return (peaksCount>=iMinPeakCount);
   }
 
   bool acceptPeaksThreshold(const RAMPREAL *pPeaks) const {
      for (int i = 0; pPeaks[i] != -1; i += 2) {
         RAMPREAL fInten = pPeaks[i+1];       
         if (fInten > dMinPeakThreshold) {
            return true;
         }
      }
      return false;
   }
};


int main(int argc, char **argv)
{
   hooks_tpp handler(argc,argv); // set up install paths etc
   printf("\n");

   struct OptionsStruct options;

   int i = 1;
   for ( ; i < argc; i++)
   {
      if (*argv[i] != '-')
          break;

      if (!setOption(argv[0], argv[i], &options)) {
         usage(argv[0]);
       }
   }

   if (i >= argc)
      usage(argv[0]);

   printf(" MzXML2Search - %s\n\n", getFormatDescription(options.format).c_str());
   for ( ; i < argc; i++)
   {
      convert(argv[i], options);
   }

   return 0;
} /*main*/


void convert(char *szXMLFile, const OptionsStruct& options)
{
    FILE *pfOut;

    /*
    * check that ends in .mzXML or other supported extension
    */
    const char *xmlFileExt = rampValidFileType(szXMLFile);
    if (!xmlFileExt)
    {
       const char **mzXMLFileTypes = rampListSupportedFileTypes();
       printf(" Error: Unexpected input file type.  Expected ");
       for (int i=0;mzXMLFileTypes[i];i++) {
          if (i) {
             printf(", ");
          }
          if (!mzXMLFileTypes[i+1]) {
             printf(" or ");
          }
          printf("%s",mzXMLFileTypes[i]);
       }
       printf(" .\n\n");
       exit(EXIT_FAILURE);
    }

#if 0
    //setRampOption(OPTION_ALL_SCANS);
    //setRampOption(OPTION_ORIGIN_SCANS);
    //setRampOption(OPTION_AVERAGE_SCANS);
#endif
    RAMPFILE *pFI = rampOpenFile(szXMLFile);
    if (pFI == NULL)
    {
        printf(" Could not open input file %s\n", szXMLFile);
        exit(EXIT_FAILURE);
    }

    printf(" Reading %s\n", szXMLFile);

    /*
    * Read the offset of the index
    */
    printf(" Getting the index offset\n");
    ramp_fileoffset_t indexOffset = getIndexOffset(pFI);
    if (!indexOffset)
    {
        printf(" Could not find index offset.\n");
        printf(" File may be corrupted.\n");
        exit(EXIT_FAILURE);
    }

    /*
    * Read the scan index into a vector, get LastScan
    */
    printf(" Reading the index\n");
    int iAnalysisFirstScan = 1;
    int iAnalysisLastScan = 0;
    ramp_fileoffset_t* pScanIndices = readIndex(pFI , indexOffset, &iAnalysisLastScan);

    /*
    * If no output path, use a path located beside
    * the .mzXML file, with file's base-name.
    */
    std::string strBasePath = szXMLFile;
    strBasePath.erase(strBasePath.length()-strlen(xmlFileExt));
    const char* szBaseName = findRightmostPathSeperator_const(strBasePath.data());
    if (szBaseName == NULL)
        szBaseName = strBasePath.data();
    else
        szBaseName++;
    std::string strBaseName = szBaseName;

    std::string strOutputPath;
    if (options.outputPath != NULL)
        strOutputPath = options.outputPath;
    else
    {
        strOutputPath = strBasePath;
        if (options.format == DTA)
        {
            std::string strCommand = "mkdir ";
            strCommand.append("\"").append(strBasePath).append("\"");
            if (system(strCommand.data()) == -1)
            {
                printf(" Error - failed to create directory %s\n\n", strBasePath.data());
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            strOutputPath.append(".");
            strOutputPath.append(getFormatExtension(options.format));
        }
    }

    // except for the single file-per-scan DTA output format,
    // open the single output file for writing and write some header info
    if (options.format != DTA)
    {
        if ( (pfOut = fopen(strOutputPath.data(), "w")) == NULL)
        {
            printf(" Error - failed to open %s\n\n", strOutputPath.data());
            exit(EXIT_FAILURE);
        }
        if (options.format == MGF)
        {
            fprintf(pfOut, "COM=Conversion of %s to mascot generic\n", szXMLFile);
            if (options.iCharge == 0)
            {
                // any ion block without a CHARGE attribute will be treated as 2+ and 3+
                fprintf(pfOut, "CHARGE=2+ and 3+\n");
            }
            else if (options.iCharge == options.iChargeLast)
            {
                // all blocks will have the same charge
                fprintf(pfOut, "CHARGE=%d+\n", options.iCharge);
            }
        }
        else if (options.format == MS2)
        {
            fprintf(pfOut, "H\tExtractor\tMzXML2Search\n");
        }
    }
    // end open single file and write header

    printf(" scan: ");

    if (iAnalysisFirstScan < options.iFirstScan)
        iAnalysisFirstScan = options.iFirstScan;
    if (iAnalysisLastScan > (int)options.iLastScan)
        iAnalysisLastScan = options.iLastScan;

    long lScanCount=1;
    for (long ctScan = iAnalysisFirstScan; ctScan <= iAnalysisLastScan; ctScan++)
    {
        /*
        * read scan header
        */
        struct ScanHeaderStruct scanHeader;
        char szTmp[500];
        int x;
        readHeader(pFI, pScanIndices[ctScan], &scanHeader);

        sprintf(szTmp, "%5ld  %.3d%%", ctScan, (int)(100.0*ctScan/iAnalysisLastScan));
        printf("%s", szTmp);
        fflush(stdout);
        for (x=0; x<(int)strlen(szTmp); x++)
           printf("\b");

	// cerr << "Scan #" << scanHeader.acquisitionNum << " is " << scanHeader.activationMethod << endl;
	
	// check if the msLevel, peaksCount, and activationMethod pass
        if (options.acceptMSLevel(scanHeader.msLevel) &&
            options.acceptPeaksCount(scanHeader.peaksCount) &&
            options.acceptActivationMethod(scanHeader.activationMethod))
        {
                /*
                * read scan peaks
                */
                RAMPREAL *pPeaks = readPeaks(pFI, pScanIndices[ctScan]);
		
		// check if we have peaks above the min. peak intensity threshold
		if (options.acceptPeaksThreshold(pPeaks)) {

		  /*
		   * deterimine what charge(s) we'll use when outputing this scan.
		   */

		
		  // the array of possible charges actually to be used for this specific scan's output

		  bool mzXMLChargeArray[MAXCHARGES];
		  bool userChargeArray[MAXCHARGES];
		  bool guessedChargeArray[MAXCHARGES]; // algorithmically guessed
		  bool haveMzXMLCharges = false;
		  bool haveUserCharges = false;
		  bool haveGuessedCharges = false;

		  for (int i=0; i<MAXCHARGES; i++) {
		    mzXMLChargeArray[i] = false;
		    userChargeArray[i] = false;
		    guessedChargeArray[i] = false;
		  }


		  // do we have any charges from the mzXML header,
		  // either precursorMZ or a possibleCharges list?
		
		  if ( (scanHeader.precursorCharge > 0 && 
                        scanHeader.precursorCharge < MAXCHARGES)
		       ||
		       (scanHeader.numPossibleCharges > 0)) {
		    haveMzXMLCharges = true;
		    
		    // copy all charges from mzXML header to our local array
		    if (scanHeader.precursorCharge > 0) {
		      mzXMLChargeArray[scanHeader.precursorCharge] = true;
		    }
		    if (scanHeader.numPossibleCharges > 0) {
		      for (int i=0; i<MAXCHARGES; i++) {
			if (scanHeader.possibleChargesArray[i]) {
			  mzXMLChargeArray[i] = true;
			}
		      }
		    }
		  }
		

		  // do we have a user charge?
		  if (options.iCharge > 0) {
		  
		    // should use something dynamic...
		    if (options.iCharge > MAXCHARGES) {
		      printf("error, cannot handle precursor charges > %d (got %d)\n", MAXCHARGES, options.iCharge);
		      exit(-1);
		    }
		    userChargeArray[options.iCharge] = true;
		    haveUserCharges = true;
		  }

		  // do we have a user charge range?
		  if (options.iChargeLast > options.iCharge) {
		    // should use something dynamic...
		    if (options.iChargeLast > MAXCHARGES) {
		      printf("error, cannot handle precursor charges > %d (got %d)\n", MAXCHARGES, options.iChargeLast);
		      exit(-1);
		    }
		    for (int i=options.iCharge; i <= options.iChargeLast; ++i) {
		      userChargeArray[i] = true;
		    }
		    haveUserCharges = true;
		  }
		

		  bool useUserCharges = false;
		  bool useMzXMLCharges = false;
		  bool useGuessedCharges = false;


		  /********************************************************
		  determine which charge(s) to actually use for this scan
		  */
 
		  // if the user specified a non-zero charge as a charge to force output as, then use user-specified charges
		  if (options.useChargeOverrides && haveUserCharges)
		    {
		      // use the user-specified charges as overrides.
		      useMzXMLCharges = false;
		      useUserCharges = true;
		      useGuessedCharges = false;
		    }

		  // otherwise, if there existing mzXML charge(s) exist, we should use that:
		  else if (haveMzXMLCharges /* && user did NOT choose 'force', we know*/)
		    {
		      useMzXMLCharges = true;
		      useUserCharges = false;
		      useGuessedCharges = false;
		    }

		  // otherwise, the user is *suggesting* charge(s) for otherwise unspecified charges, 
		  // use user-specified charges
		  else if (!options.useChargeOverrides && haveUserCharges /* && we don't have mzXML charges, we know*/ )
		    {
		      // use the user-specified charges for this previously undetermined scan
		      useMzXMLCharges = false;
		      useUserCharges = true;
		      useGuessedCharges = false;
		    }


		  // otherwise try to algorithmically determine the charge
		  // NOTE: +1, +2, +3 only!
		  // TODO: expand to higher charge states?
		  else
		    {
		      useMzXMLCharges = false;
		      useUserCharges = false;
		      useGuessedCharges = true;

		      /*
		       * simple charge state determination
		       */
		      RAMPREAL fSumBelow = 0.0;
		      RAMPREAL fSumTotal = 0.0;
		      for(int i = 0; pPeaks[i] != -1; i += 2 )
			{
			  RAMPREAL fMass = pPeaks[i];
			  RAMPREAL fInten = pPeaks[i+1];

			  fSumTotal += fInten;
			  if (fMass<scanHeader.precursorMZ)
                            fSumBelow += fInten;
			}

		      /*
		       * If greater than 95% signal is below precursor, then
		       * it looks like singly charged peptide.
		       */
		      if (fSumTotal == 0.0 || fSumBelow/fSumTotal>0.95)
			guessedChargeArray[1] = true;
		      else
			{
			  guessedChargeArray[2] = true;
			  guessedChargeArray[3] = true;
			}
		    }



		  // sanity check
		  if (!(useMzXMLCharges || useUserCharges || useGuessedCharges)) {
		    printf("logic error! no charges determined for this scan!\n");
		    printf("exiting with error\n");
		    exit(-1);
		  }

		  int count=0;
		  if (useMzXMLCharges) ++count;
		  if (useUserCharges) ++count;
		  if (useGuessedCharges) ++count;
		  if (count != 1) {
		    printf("logic error! multiple charge methods determined for this scan!\n");
		    printf("exiting with error\n");
		    exit(-1);
		  }

		  //printf("scan %d:", (int)ctScan);
		  bool* chargeArray=NULL;
		  if (useMzXMLCharges) {
		    //printf("using useMzXMLCharges\n");
		    chargeArray=mzXMLChargeArray;
		  }
		  else if (useUserCharges) {
		    //printf("using useUserCharges\n");
		    chargeArray=userChargeArray;
		  }
		  else if (useGuessedCharges) {
		    //printf("using useGuessedCharges\n");
		    chargeArray=guessedChargeArray;
		  }


		  // output the scan(s) for this input scan note: there
		  // may be more than one output scan if the output
		  // format cannot represent multiple charges per
		  // spectra, for example in DTA where we output a +2
		  // and +3 scan for a 2,3 input.


		  switch(options.format) {
		  case MS2:
		    outputMS2(pfOut, options, scanHeader, ctScan, chargeArray, pPeaks);
		    break;
		  case ODTA:
		    outputODTA(pfOut, strBaseName, options, scanHeader, ctScan, chargeArray, pPeaks, lScanCount);
		    break;
		  case XDTA:
		    outputXDTA(pfOut, strBaseName, options, scanHeader, ctScan, chargeArray, pPeaks);
		    break;
		  case PKL:
		    outputPKL(pfOut, options, scanHeader, ctScan, chargeArray, pPeaks);
		    break;
		  case DTA:
		    outputDTA(pfOut, strBaseName, strOutputPath, options, scanHeader, ctScan, chargeArray, pPeaks);
		    break;
		  case MGF:
		    outputMGF(pfOut, strBaseName, options, scanHeader, ctScan, chargeArray, pPeaks);
		    break;
		    
		  default:
		    // shouldn't get here
		    printf("exiting with error: unhandled output format\n");
		    exit(-1);
		    break;
		  }
	
		  //free(pPeaks);
		} // end: if passes peak threshold
        } // end: if passes user filters (activationMethod, etc)
    } // end: looping through scans

    free(pScanIndices);
    rampCloseFile(pFI);

    if (options.format != DTA)
    {
        fclose(pfOut);
    }

    printf("\n");
    printf(" Done.\n\n");

} /*convert*/


void usage(char *arg)
{
    char* exe = findRightmostPathSeperator(arg);
    if (exe == NULL)
        exe = arg;
    else
        exe++;

    printf(" %s (%s)\n", exe, szTPPVersionInfo);
    printf(" Usage:  %s [options] *.mzXML or *.mzML\n", exe);
    printf("\n");
    printf("     options = -dta or -mgf or -pkl or -xdta or -odta or -ms2 output format (default dta)\n");
    printf("               -F<num>      where num is an int specifying the first scan\n");
    printf("               -L<num>      where num is an int specifying the last scan\n");

    printf("               -C<n1>[-<n2>]     \"force charge(s)\": where n1 is an integer\n");
    printf("                                 specifying the precursor charge state (or possible\n");
    printf("                                 charge range from n1 to n2 inclusive) to use; this option\n");
    printf("                                 forces input scans to be output with the user-specified\n");
    printf("                                 charge (or charge range)\n");

    printf("               -c<n1>[-<n2>]     \"suggest charge(s)\": for scans which do not have a\n");
    printf("                                 precursor charge (or charge range) already determined in the\n");
    printf("                                 input file, use the user-specified charge (or charge range)\n");
    printf("                                 for those scans.  Input scans which already have defined\n");
    printf("                                 charge (or charge range) are output with their original,\n");
    printf("                                 unchanged values.\n");

    printf("               -B<num>      where num is a float specifying minimum MH+ mass, default=%0.1f Da\n", MINMASS);
    printf("               -T<num>      where num is a float specifying maximum MH+ mass, default=%0.1f Da\n", MAXMASS);
    printf("               -P<num>      where num is an int specifying minimum peak count, default=%d\n", MINPEAKS);
    printf("               -N<num>      where num is an int specifying max peak count using most intense peaks, default=0 to print all peaks\n");
    printf("               -pm<num>     where num is an int specifying mass precision in peaklist, default=%d\n", PRECISION_MASS);
    printf("               -pi<num>     where num is an int specifying intensity precision in peaklist, default=%d\n", PRECISION_INT);
    printf("               -I<num>      where num is a float specifying minimum threshold for peak intensity, default=%0.2f\n", MININTENSITY);
    printf("               -M<n1>[-<n2>]where n1 is an int specifying MS level to export (default=2)\n");
    printf("                            and n2 specifies an optional range of MS levels to export\n");
    printf("               -A<str>      where str is the activation method, \"CID\", \"ETD\", or \"HCD\"\n");
//  printf("                              if str is blank, then only scans with no activation method are filtered\n");
    printf("                              if this option is not specified, then all scans are included\n");
    printf("               -Z<num>      maximum reported charge state for scans that do have a precursor charge;\n");
    printf("                            useful when scan has a high charge that search engines can't handle.\n");
    printf("                            No charge is reported if charge is larger than max value, default=%d.\n", MAXOUTPUTCHARGE);
    printf("               -X           remove charge-reduced precursors from the spectra (suitable for ETD).\n");
    printf("               -Q           remove iTRAQ reporter peaks in the range 112-122 Th.\n");
    printf("               -G           remove TMT reporter peaks in the range 126-132 Th.\n");
    printf("\n");
    exit(EXIT_FAILURE);
} /*usage*/


bool setOption(char *exe, char *arg, OptionsStruct *pOptions)
{
    int  iTmp = 0;
    int  iTmp2 = 0;
    double dTmp = 0.0;

    FormatType format = getFormat(arg+1);
    if (format != FORMAT_UNDEF) {
       pOptions->format = format;
       printf("output mode selected: %s\n", getFormatDescription(pOptions->format).c_str());
       return true;
    }

    switch (arg[1])
    {
    case 'F':
        if (sscanf(&arg[2], "%d", &iTmp) != 1 || iTmp < 1)
            printf("Bad first scan #: '%s' ... ignored\n", &arg[2]);
        else
            pOptions->iFirstScan = iTmp;
        break;
    case 'L':
        if (sscanf(&arg[2], "%d", &iTmp) != 1 || iTmp < 1)
            printf("Bad last scan #: '%s' ... ignored\n", &arg[2]);
        else
            pOptions->iLastScan = iTmp;
        break;
    case 'B':
        if (sscanf(&arg[2], "%lf", &dTmp) != 1 || dTmp < 0.00)
            printf("Bad minimum mass: '%s' ... ignored\n", &arg[2]);
        else
            pOptions->dMinMass = dTmp;
        break;
    case 'T':
        if (sscanf(&arg[2], "%lf", &dTmp) != 1 || dTmp < 0.000)
            printf("Bad maximum mass: '%s' ... ignored\n", &arg[2]);
        else
            pOptions->dMaxMass = dTmp;
        break;
    case 'P':
        if (sscanf(&arg[2], "%d", &iTmp) != 1 || iTmp < 0)
            printf("Bad minimum peak count: '%s' ... ignored\n", &arg[2]);
        else
            pOptions->iMinPeakCount = iTmp;
        break;
    case 'p':
        if (arg[2]=='m')
        {
           if (sscanf(&arg[3], "%d", &iTmp) != 1 || iTmp < 0 || iTmp>9)
              printf("Bad mass precision: '%s' ... ignored\n", &arg[2]);
           pOptions->iPrecisionMass = iTmp;
        }
        else if (arg[2]=='i')
        {
           if (sscanf(&arg[3], "%d", &iTmp) != 1 || iTmp < 0 || iTmp>9)
              printf("Bad intensity precision: '%s' ... ignored\n", &arg[2]);
           pOptions->iPrecisionInt = iTmp;
        }
        break;
    case 'M':
        if (sscanf(&arg[2], "%d-%d", &iTmp, &iTmp2) == 2 && iTmp >= 2 && iTmp2 >= iTmp)
        {
            pOptions->iMSLevel = iTmp;
            pOptions->iMSLevelLast = iTmp2;
        }
        else if (sscanf(&arg[2], "%d", &iTmp) == 1 && iTmp > 0)
        {
            pOptions->iMSLevel = pOptions->iMSLevelLast = iTmp;
        }
        else
            printf("Bad MS Level: '%s' ... ignored\n", &arg[2]);
        break;
    case 'I':
        if (sscanf(&arg[2], "%d", &iTmp) != 1 || iTmp < 0)
            printf("Bad minimum peak threshhold: '%s' ... ignored\n", &arg[2]);
        else
            pOptions->dMinPeakThreshold = iTmp;
            printf("Minimum peak threshold set to %f\n", pOptions->dMinPeakThreshold);
        break;
    case 'C':
        pOptions->useChargeOverrides = true;

        // Fall through.
    case 'c':
        if (sscanf(&arg[2], "%d-%d", &iTmp, &iTmp2) == 2 && iTmp > 0 && iTmp2 >= iTmp)
        {
            pOptions->iCharge = iTmp;
            pOptions->iChargeLast = iTmp2;
        }
        else if (sscanf(&arg[2], "%d", &iTmp) == 1 && iTmp > 0)
        {
            pOptions->iCharge = pOptions->iChargeLast = iTmp;
        }
        else
            printf("Bad charge: '%s' ... ignored\n", &arg[2]);
        break;
    case 'i':
        if (sscanf(&arg[2], "%d", &iTmp) != 1 || iTmp < 0)
            printf("Bad minimum peak intensity for filtering: '%s' ... ignored\n", &arg[2]);
        else
            pOptions->dMinPeakDtaFilter = iTmp;
        break;
    case 'N':
        if (sscanf(&arg[2], "%d", &iTmp) != 1 || iTmp < 0)
            printf("Bad max num peaks to print per scan: '%s' ... ignored\n", &arg[2]);
        else
            pOptions->iMaxPrintPeaks = iTmp;
        break;
    case 'O':
        pOptions->outputPath = &arg[2];
        break;
    case 'A':
        pOptions->activationMethod = &arg[2];
        break;
    case 'Z':
        if (sscanf(&arg[2], "%d", &iTmp) != 1 || iTmp < 0 || iTmp>=MAXCHARGES)
            printf("Bad maximum output charge: '%s' ... ignored\n", &arg[2]);
        else
            pOptions->dMaxOutputCharge = iTmp;
            printf("Maximum output charge set to %f\n", pOptions->dMaxOutputCharge);
        break;
    case 'X':
        pOptions->removeChargeReducedPrecursors = true;
	break;
    case 'Q':
        pOptions->removeITRAQReporterPeaks = true;
	break;
    case 'G':
        pOptions->removeTMTReporterPeaks = true;
	break;
      
    default:
        printf("ERROR - option not recognized: %s\n\n", arg);
        return false;
        break;
    }
    arg[0] = '\0';

    return true;
} /*setOptions*/


std::string getFormatDescription(const FormatType format) {
  std::string desc;

  switch (format) {
  case DTA:
    desc = "SEQUEST DTA format";
    break;
  case XDTA:
    desc = "X! Tandem single file DTA";
    break;
  case ODTA:
    desc = "OMSSA merged single file DTA";
    break;
  case MGF:
    desc = "Mascot Generic Format";
    break;
  case PKL:
    desc = "Micromass PKL";
    break;
  case MS2:
    desc = "SEQUEST MS2 format";
    break;
  default:
    desc = "unknown format!";
    break;
  }
  return desc;
}


std::string getFormatExtension(const FormatType format) {
  std::string ext;

  switch (format) {
  case DTA:
    ext = "dta";
    break;
  case XDTA:
    ext = "xdta";
    break;
  case ODTA:
    ext = "odta";
    break;
  case MGF:
    ext = "mgf";
    break;
  case PKL:
    ext = "pkl";
    break;
  case MS2:
    ext = "ms2";
    break;
  default:
    ext = "???";
    break;
  }
  return ext;
}


FormatType getFormat(const std::string& formatStr) {
  if (formatStr == "dta") {
    return DTA;
  }
  else if (formatStr == "mgf") {
    return MGF;
  }
  else if (formatStr == "pkl") {
    return PKL;
  }
  else if (formatStr == "ms2") {
    return MS2;
  }
  else if (formatStr == "odta") {
    return ODTA;
  }
  else if (formatStr == "xdta") {
    return XDTA;
  }
  else {
    return FORMAT_UNDEF;
  }
}  


void outputMS2(FILE *pfOut, const OptionsStruct& options,
      const ScanHeaderStruct& scanHeader,
      long curScanNum,
      const bool* chargeArray,
      const RAMPREAL *pPeaks)
{
   fprintf(pfOut, "S\t%ld\t%ld\t%0.2f\n", curScanNum, curScanNum, scanHeader.precursorMZ);
   fprintf(pfOut, "I\tRTime\t%0.4f\n", scanHeader.retentionTime/60.0);
   fprintf(pfOut, "I\tMS1Intensity\t%0.2f\n", scanHeader.precursorIntensity);

   // output charges and prec. mass for this scan
   for (int curCharge = 0; curCharge <= options.dMaxOutputCharge; ++curCharge) {
      if (chargeArray[curCharge]) {
         double dPepMass = calculatePeptideMass(scanHeader.precursorMZ, curCharge); 
         if ( (options.dMinMass <= dPepMass) && (options.dMaxMass >= dPepMass) ) {
            fprintf(pfOut, "Z\t%d\t%0.6f\n", curCharge, dPepMass);
         }
      }
   }

   // print ions
   writePeaksTSV(pfOut, pPeaks, options, scanHeader.precursorMZ);
}


void outputODTA(FILE *pfOut,
      const std::string& strBaseName,
      const OptionsStruct& options,
      const ScanHeaderStruct& scanHeader,
      long curScanNum,
      const bool* chargeArray,
      const RAMPREAL *pPeaks,
      long& lScanCount)
{
   for (int curCharge = 0; curCharge <= options.dMaxOutputCharge; ++curCharge) {
      if (chargeArray[curCharge]) {
         double dPepMass = calculatePeptideMass(scanHeader.precursorMZ, curCharge);
         if ( (options.dMinMass <= dPepMass) && (options.dMaxMass >= dPepMass) ) {
            fprintf(pfOut, "<dta id=\"%ld\" name=\"%s.%.5ld.%.5ld.%d.dta\">\n",
                  lScanCount, strBaseName.data(), curScanNum, curScanNum, curCharge);
            lScanCount++;
            fprintf(pfOut, "%0.6f%s%d\n", dPepMass, FIELDSEPARATOR, curCharge);

            // print ions
            writePeaksTSV(pfOut, pPeaks, options, scanHeader.precursorMZ, curCharge);

            // close the scan entry
            fprintf(pfOut, "</dta>\n\n");
         }
      }
   }
}


void outputXDTA(FILE *pfOut,
      const std::string& strBaseName,
      const OptionsStruct& options,
      const ScanHeaderStruct& scanHeader,
      long curScanNum,
      const bool* chargeArray,
      const RAMPREAL *pPeaks)
{
   for (int curCharge = 0; curCharge <= options.dMaxOutputCharge; ++curCharge) {
      if (chargeArray[curCharge]) {
         double dPepMass = calculatePeptideMass(scanHeader.precursorMZ, curCharge);
         if ( (options.dMinMass <= dPepMass) && (options.dMaxMass >= dPepMass) ) {
            int startScanNum, endScanNum;
            getScanSpanRange(&scanHeader, &startScanNum, &endScanNum);

            fprintf(pfOut, "%0.6f%s%d%s", dPepMass, FIELDSEPARATOR, curCharge, FIELDSEPARATOR);
            fprintf(pfOut, "%s.%.4d.%.4d.%d", strBaseName.data(), startScanNum, endScanNum, curCharge);
            fprintf(pfOut, "\n");

            // print ions
            writePeaksTSV(pfOut, pPeaks, options, scanHeader.precursorMZ, curCharge);

            fprintf(pfOut, "\n");
         }
      }
   }
}


void outputPKL(FILE *pfOut,
      const OptionsStruct& options,
      const ScanHeaderStruct& scanHeader,
      long curScanNum,
      const bool* chargeArray,
      const RAMPREAL *pPeaks)
{
   for (int curCharge = 0; curCharge <= options.dMaxOutputCharge; ++curCharge) {
      if (chargeArray[curCharge]) {
         double dPepMass = calculatePeptideMass(scanHeader.precursorMZ, curCharge);
         if ( (options.dMinMass <= dPepMass) && (options.dMaxMass >= dPepMass) ) {
            int startScanNum, endScanNum;
            getScanSpanRange(&scanHeader, &startScanNum, &endScanNum);
            fprintf(pfOut, "%0.6f%s%0.4f%s%d\n",
                  scanHeader.precursorMZ, FIELDSEPARATOR, scanHeader.precursorIntensity, FIELDSEPARATOR, curCharge);

            // print ions
            writePeaksTSV(pfOut, pPeaks, options, scanHeader.precursorMZ, curCharge);
            fprintf(pfOut, "\n");
         }
      }
   }
}

  
void outputDTA(FILE *pfOut,
      const std::string& strBaseName,
      const std::string& strOutputPath,
      const OptionsStruct& options,
      const ScanHeaderStruct& scanHeader,
      long curScanNum,
      const bool* chargeArray,
      const RAMPREAL *pPeaks)
{
   for (int curCharge = 0; curCharge <= options.dMaxOutputCharge; ++curCharge) {
      if (chargeArray[curCharge]) {
         double dPepMass = calculatePeptideMass(scanHeader.precursorMZ, curCharge);
         if ( (options.dMinMass <= dPepMass) && (options.dMaxMass >= dPepMass) ) {
            int startScanNum, endScanNum;
            char szOut[16384];

            getScanSpanRange(&scanHeader, &startScanNum, &endScanNum);
            sprintf(szOut, "%s/%s.%.5d.%.5d.%d.dta", strOutputPath.data(),
                  strBaseName.data(), startScanNum, endScanNum, curCharge);

            // open single-scan dta file
            if ((pfOut = fopen(szOut, "w")) == NULL) {
               printf(" Cannot open file %s\n", szOut);
               exit(EXIT_FAILURE);
            }

            fprintf(pfOut, "%0.6f%s%d\n", dPepMass, FIELDSEPARATOR, curCharge);

            // print ions
            writePeaksTSV(pfOut, pPeaks, options, scanHeader.precursorMZ, curCharge);

            // close single-scan dta file
            fclose(pfOut);
         }
      }
   }
}


void outputMGF(FILE *pfOut,
      const std::string& strBaseName,
      const OptionsStruct& options,
      const ScanHeaderStruct& scanHeader,
      long curScanNum,
      const bool* chargeArray,
      const RAMPREAL *pPeaks)
{ 

  // since we may have more than 8 possible charges, we use charge
  // state for our TITLE convension, and PEPMASS is calcated from the
  // assumed charge, rather than outputting a list of (less than 8)
  // possible charges for one scan, I output one scan per possible
  // charge state.


  /*
  std::vector<int> chargeVec;
  for (int curCharge = 0; curCharge <= options.dMaxOutputCharge; ++curCharge) {
    if (chargeArray[curCharge]) {
      chargeVec.push_back(curCharge);
    }
  }

  if (chargeVec.size() > 8) {
    printf("exiting with error: MGF output can not specify more than 8 possible charges per scan\n");
    exit(-1);
  }
  // if (!chargeVec.empty()) {
  */
  
   for (int curCharge = 0; curCharge <= options.dMaxOutputCharge; ++curCharge) {
      if (chargeArray[curCharge]) {
         double dPepMass = calculatePeptideMass(scanHeader.precursorMZ, curCharge);
         if ( (options.dMinMass <= dPepMass) && (options.dMaxMass >= dPepMass) ) {
            int startScanNum=0, endScanNum=0;
            getScanSpanRange(&scanHeader, &startScanNum, &endScanNum); 
            
            fprintf(pfOut, "BEGIN IONS\n");
            if (scanHeader.precursorIntensity == 0.0) {
               fprintf(pfOut, "PEPMASS=%0.6f\n", scanHeader.precursorMZ);
            }
            else {
               fprintf(pfOut, "PEPMASS=%0.6f %0.4f\n", scanHeader.precursorMZ, scanHeader.precursorIntensity);
            }

            fprintf(pfOut, "CHARGE=%d+\n", curCharge);
	/*
	int lastChargeIndex = chargeVec.size() - 1;
	for (int i=0; i<lastChargeIndex; i++) {
	  fprintf(pfOut, "%d+", chargeVec[i]);
	  if (i != (lastChargeIndex-1)) {
	    fprintf(pfOut, ",");
	  }
	}
	fprintf(pfOut, "\n");
	*/
            if (startScanNum == endScanNum) {
	      fprintf(pfOut, "SCANS=%d\n", startScanNum);
	    } else {
	      fprintf(pfOut, "SCANS=%d-%d\n", startScanNum, endScanNum);
	    }

            fprintf(pfOut, "RTINSECONDS=%0.2f\n", scanHeader.retentionTime);
            fprintf(pfOut, "TITLE=%s.%.5d.%.5d.%d\n", strBaseName.data(), startScanNum, endScanNum, curCharge);

            // print ions
            writePeaksTSV(pfOut, pPeaks, options, scanHeader.precursorMZ, curCharge);

            // close this entry
            fprintf(pfOut, "END IONS\n");
         }
      }
   }
}


void writePeaksTSV(FILE *pfOut,
      const RAMPREAL *pPeaks,
      const OptionsStruct& options,
      const double precursorMZ,
      const int precursorCharge)
{
   char format[1024];
   int iNumPeaks=0;

   format[0]=0;
   sprintf(format, "%%0.%df%%s%%0.%df\n", options.iPrecisionMass, options.iPrecisionInt);

   if (options.iMaxPrintPeaks > 0)
   {
      // find # ions
      for (int i = 0; pPeaks[i] != -1; i += 2)
         iNumPeaks++;
   }

   if (options.iMaxPrintPeaks==0 || options.iMaxPrintPeaks >= iNumPeaks) // print all ions
   {
      // for iMinPeakCount filter to work, must count peaks that pass intensity threshold
      iNumPeaks=0;
      for (int i = 0; pPeaks[i] != -1; i += 2) {
        RAMPREAL fMass = pPeaks[i];
	RAMPREAL fInten = pPeaks[i+1];

         if (fInten >= options.dMinPeakThreshold &&
	     (!(options.removeChargeReducedPrecursors && isNearChargeReducedPrecursors(fMass, precursorMZ, precursorCharge))) &&
	     (!(options.removeITRAQReporterPeaks && isITRAQReporterPeaks(fMass))) &&
	     (!(options.removeTMTReporterPeaks && isTMTReporterPeaks(fMass))))
            iNumPeaks++;
      }

      if (iNumPeaks >= options.iMinPeakCount) {
         for (int i = 0; pPeaks[i] != -1; i += 2) {
            RAMPREAL fMass = pPeaks[i];
            RAMPREAL fInten = pPeaks[i+1];

            if (fInten >= options.dMinPeakThreshold && 
	        (!(options.removeChargeReducedPrecursors && isNearChargeReducedPrecursors(fMass, precursorMZ, precursorCharge))) &&
                (!(options.removeITRAQReporterPeaks && isITRAQReporterPeaks(fMass))) &&
	        (!(options.removeTMTReporterPeaks && isTMTReporterPeaks(fMass))))
	      
	      fprintf(pfOut, format, fMass, FIELDSEPARATOR, fInten);      
         }
      }
      else {
         // no peaks passed filters ... what to do???  print empty peak 
         fprintf(pfOut, format, 0.0, FIELDSEPARATOR, 0.0);      
      }
   }
   else  // just print out top iMaxPrintPeaks most intense peaks
   {
	   struct PeaksStruct *pTmpPeaks = new struct PeaksStruct [iNumPeaks];

      // assign ions to struct
      iNumPeaks=0;

      for (int i = 0; pPeaks[i] != -1; i += 2) {
         pTmpPeaks[iNumPeaks].fMass = pPeaks[i];
         pTmpPeaks[iNumPeaks].fInten = pPeaks[i+1];

         if (pTmpPeaks[iNumPeaks].fInten >= options.dMinPeakThreshold &&
	     (!(options.removeChargeReducedPrecursors && isNearChargeReducedPrecursors(pTmpPeaks[iNumPeaks].fMass, precursorMZ, precursorCharge))))

            iNumPeaks++;
      }

      // NOTE: pTmpPeaks will only contain peaks above dMinPeakThreshold and not near charge-reduced precursors (if that option is set)

      if (iNumPeaks >= options.iMinPeakCount) {
         int iPrintPeaks = options.iMaxPrintPeaks;

         if (iNumPeaks < options.iMaxPrintPeaks)
            iPrintPeaks = iNumPeaks;

         // sort struct by intensity to get biggest peaks
         qsort(pTmpPeaks, iNumPeaks, sizeof(struct PeaksStruct), struct_cmp_inten);

         // sort N most intense peaks by mass
         qsort(pTmpPeaks, iPrintPeaks, sizeof(struct PeaksStruct), struct_cmp_mass);
      
         for (int i=0; i<iPrintPeaks; i++)
            fprintf(pfOut, format, pTmpPeaks[i].fMass, FIELDSEPARATOR, pTmpPeaks[i].fInten);      
      }
      else {
         // no peaks passed filters ... what to do???  print empty peak 
         fprintf(pfOut, format, 0.0, FIELDSEPARATOR, 0.0);      
      }
	  delete[] pTmpPeaks;
   }
}


int struct_cmp_mass(const void *a, const void *b)
{
   struct PeaksStruct *ia = (struct PeaksStruct *)a;
   struct PeaksStruct *ib = (struct PeaksStruct *)b;

   return (int)(100.0*ia->fMass - 100.0*ib->fMass);
} 


int struct_cmp_inten(const void *a, const void *b)
{
   struct PeaksStruct *ia = (struct PeaksStruct *)a;
   struct PeaksStruct *ib = (struct PeaksStruct *)b;

   return (int)(100.0*ib->fInten - 100.0*ia->fInten);
} 

bool isNearChargeReducedPrecursors(double mz, double precursorMz, int precursorCharge) {
  
  if (precursorMz < 0.0001) return (false); // m_precursorMz somehow not set, can't determine
  
  int lowCharge = precursorCharge;
  int highCharge = precursorCharge;
  
  if (precursorCharge == 0) {
    // charge unknown: remove all possible charge-reduced precursors for precursor charge up to 6.
    lowCharge = 3;
    highCharge = 6;
  }

  for (int possiblePrecursorCharge = lowCharge; possiblePrecursorCharge <= highCharge; possiblePrecursorCharge++) {
    double precursorMass = precursorMz * possiblePrecursorCharge;

    // assuming all possible charge-reduced precursors from 1 up to precursorCharge
    for (int charge = possiblePrecursorCharge; charge >= 1; charge--) {
      // Wipe out -50 to +6. Will remove the water/ammonia loss and the higher isotopes. 
      // Will miss the greater neutral losses -- but we don't want to wipe out too much to catch these less common ones
      if (mz >= (precursorMass - 20.0) / (double)charge && mz <= (precursorMass + 6.0) / (double)charge) {
	return (true);
	
      }
    }
	
  }
  return (false);
  
}

bool isITRAQReporterPeaks(double mz) {
  
  return (mz >= 112.0 && mz <= 122.0);
    
}

bool isTMTReporterPeaks(double mz) {
  
  return (mz >= 126.0 && mz <= 132.0);
  
}
