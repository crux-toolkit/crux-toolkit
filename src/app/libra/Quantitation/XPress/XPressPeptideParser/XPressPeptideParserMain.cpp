/*

Program       : XPressPeptideParser
Author        : J.Eng and Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN Info      : $Id: XPressPeptideParserMain.cpp 8418 2021-03-22 01:05:09Z real_procopio $


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
#include "Common/constants.h"

#include "Common/TPPVersion.h"         // contains version number, name, revision
#include "Parsers/Parser/TagListComparator.h" // for REGRESSION_TEST_CMDLINE_ARG defn

#include <string>
#include <iostream>
#include <iomanip>


void USAGE(char *szCommand,
      InputStruct options);

int main(int argc, char **argv)
{
   hooks_tpp(argc,argv); // handle installdir issues, etc
   InputStruct options;
   int    iNumIsotopicPairs = 0;

   char  *testMode = NULL;      // regression test stuff - bpratt Insilicos LLC, Nov 2005

   // default settings
   options.bUseFixedScanRange = 0;    // 0 = normal, 1=use fix scan range, 2=use fix scan range from peak apex
   options.iFixedScanRange = 5;   // for fix range above, default (set below) is +- 5 scans from determined scan
   options.iMetabolicLabeling = 0;      // 0 = false; 1 = N15; 2 = C13

   options.bUseSameScanRange = TRUE;    // -a, -b => (FALSE)
   options.bIntensityRatio = FALSE;
   options.bLabelFreeMode = FALSE;
   options.bForceMS1Endpoints = FALSE;
   options.bXpressLight1 = 0;   // 0 = unused, 1= light (-L), 2=heavy (-H)
   strcpy(options.szXpressResidues1, "C");      // -n
   strcpy(options.szXpressResidues, "C");       // -n
   options.szXpressResidues2[0] = '\0';
   options.szXpressResidues3[0] = '\0';
   options.dXpressMassDiff1 = 9.0;              // -n
   options.dXpressMassDiff2 = 0.0;
   options.dXpressMassDiff3 = 0.0;
   options.dMassTol = 0.5;                      // -m
   options.iMinNumChromatogramPoints = 5;       // -c
   options.iMinNumIsotopePeaks = 1;             // -p
   options.iChargeState = -1;
   options.szMzXMLDir[0] = '\0';
   options.dMinPprob= 0.5;
   options.bPpmMassTol = FALSE;   // Dalton tolerance by default

   if (argc < 2)
   {
      USAGE(argv[0], options);
   }

   // now look for command line

   for (int k = 2; k < argc; k++)
   {
      if (strlen(argv[k]) > 1 && argv[k][0] == '-')
      {
         if (argv[k][1] == 'b')
         {
            options.bUseSameScanRange = FALSE;
         }
         else if (argv[k][1] == 'i')
         {
            options.bIntensityRatio = TRUE;
         }
         else if (argv[k][1] == 'l')
         {
            options.bLabelFreeMode = TRUE;
         }
         else if (argv[k][1] == 'a')
         {
            options.bPpmMassTol = TRUE;
         }
         else if (argv[k][1] == 'E')
         {
            options.bForceMS1Endpoints = TRUE;
         }
         else if (argv[k][1] == 'f')
         {
            options.bUseFixedScanRange = 1;
            if (strlen(argv[k]) > 2)
               sscanf(argv[k] + 2, "%d", &(options.iFixedScanRange));
            else
               options.iFixedScanRange = 5;
         }
         else if (argv[k][1] == 'F')
         {
            options.bUseFixedScanRange = 2;
            if (strlen(argv[k]) > 2)
               sscanf(argv[k] + 2, "%d", &(options.iFixedScanRange));
            else
               options.iFixedScanRange = 5;
         }
         else if (argv[k][1] == 'L')
         {
            options.bXpressLight1 = 1;
         }
         else if (argv[k][1] == 'H')
         {
            options.bXpressLight1 = 2;
         }
         else if (argv[k][1] == 'M')
         {
            options.iMetabolicLabeling = 1;              // N15 metabolic labeling
            strcpy(options.szXpressResidues, "N15");     // need to set some value for labeled residue
            strcpy(options.szXpressResidues1, "N15");
            options.dXpressMassDiff1 = 0.997034893;      // variable per residue
         }
         else if (argv[k][1] == 'O')
         {
            options.iMetabolicLabeling = 2;              // C13 metabolic labeling
            strcpy(options.szXpressResidues, "C13");     // need to set some value for labeled residue
            strcpy(options.szXpressResidues1, "C13");
            options.dXpressMassDiff1 = 1.0033548378;     // variable per residue
         }
         else if (argv[k][1] == 'm' && strlen(argv[k]) > 2)
         {
            sscanf((char *) (argv[k] + sizeof(char) * 2), "%lf", &options.dMassTol);
         }
         else if (argv[k][1] == 'c' && strlen(argv[k]) > 2)
         {
            sscanf((char *) (argv[k] + sizeof(char) * 2), "%d", &options.iMinNumChromatogramPoints);
         }
         else if (argv[k][1] == 'p' && strlen(argv[k]) > 2)
         {
            sscanf((char *) (argv[k] + sizeof(char) * 2), "%d", &options.iMinNumIsotopePeaks);
         }
         else if (argv[k][1] == 'q' && strlen(argv[k]) > 2)
         {
            sscanf((char *) (argv[k] + sizeof(char) * 2), "%lf", &options.dMinPprob);
         }
         else if (argv[k][1] == 'n' && strlen(argv[k]) > 2)
         {
            char  *pStr;

            if (!(pStr = strchr(argv[k] + 2, ',')))
            {
               fprintf(stderr, "Error - XPRESS residue option ... must specify residues,mass\n");
               fprintf(stderr, "For example:  -nC,9.0\n");
               exit(1);
            }

            *pStr = ' ';

            if (iNumIsotopicPairs == 0)
            {
               options.szXpressResidues[0] = 0;
               sscanf((char *) (argv[k] + sizeof(char) * 2), "%s %lf", options.szXpressResidues1, &options.dXpressMassDiff1);
               strcpy(options.szXpressResidues, options.szXpressResidues1);
            }
            else if (iNumIsotopicPairs == 1)
            {
               sscanf((char *) (argv[k] + sizeof(char) * 2), "%s %lf", options.szXpressResidues2, &options.dXpressMassDiff2);
               strcat(options.szXpressResidues, options.szXpressResidues2);
            }
            else if (iNumIsotopicPairs == 2)
            {
               sscanf((char *) (argv[k] + sizeof(char) * 2), "%s %lf", options.szXpressResidues3, &options.dXpressMassDiff3);
               strcat(options.szXpressResidues, options.szXpressResidues3);
            }
            iNumIsotopicPairs++;
         }
         else if (argv[k][1] == 'd' && strlen(argv[k]) > 2)
         {
            strcpy(options.szMzXMLDir, argv[k] + 2);
         }
         else if (!strncmp (argv[k], REGRESSION_TEST_CMDLINE_ARG, strlen(REGRESSION_TEST_CMDLINE_ARG)))
         {
            testMode = argv[k]; // regression test stuff - bpratt Insilicos LLC, Nov 2005
         }
         else if (argv[k][1] == 'h')
         {                      // display help
            USAGE(argv[0], options);
         }
      }
   }

   //  XPressPeptideParser*  new XPressPeptideParser(argv[1], options);
   XPressPeptideParser * p = new XPressPeptideParser(argv[1], options, testMode);

   delete p;

   return 0;
}


void USAGE(char *szCommand,
      InputStruct options)
{
   cout << fixed << std::setprecision(4) ;
   cout << " " << szCommand << " (" << szTPPVersionInfo << ")" << endl;
   cout << "USAGE:    XPressPeptideParser [xmlfile] [options]" << endl;
   cout << "Options:  -m<num>    change XPRESS mass tolerance (default=" << options.dMassTol << ")" << endl;
   cout << "          -a         tolerance specified by -m is in ppm (default=Daltons)" << endl;
   cout << "          -n<str>,<num>   when specifying multiple isotopic labels, use this option e.g. -nK,6.0 -nR,10.0" << endl;
   cout << "                            also use 'n' for labeled n-terminus and 'c' for c-terminus" << endl;
   cout << "                            (default residue '" << options.szXpressResidues << "' and mass '" << options.dXpressMassDiff1 << "')" << endl;
   cout << "          -b         heavy labeled peptide elutes before light labeled partner" << endl;
   cout << "          -f<num>    fix elution peak as +-<num> scans from start scan (default=5)" << endl;
   cout << "          -F<num>    fix elution peak as +-<num> scans from identified peak apex (default=5)" << endl;
   cout << "          -L         for ratio, set/fix light to 1, vary heavy" << endl;
   cout << "          -H         for ratio, set/fix heavy to 1, vary light" << endl;
   cout << "          -M         for 15N metabolic labeling" << endl;
   cout << "          -O         for 13C metabolic labeling" << endl;
   cout << "          -c<num>    minimum number of chromatogram points needed for quantitation (default=" << options.iMinNumChromatogramPoints << ")" << endl;
   cout << "          -p<num>    number of 13C isotopic peaks to add to precursor chromatogram (default=" << options.iMinNumIsotopePeaks << ")" << endl;
   cout << "          -q<num>    where <num> is minimum probability needed to quantify a peptide (default=" << options.dMinPprob << ")" << endl;
   cout << "          -i         also export intensities and intensity based ratio" << endl;
   cout << "          -l         label free mode: stats on precursor ions only, no ratios" << endl;
   cout << "                       only relevant label-free parameters are -m, -a, -c, and -p" << endl;
   cout << "          -d<path>   path to mzXML file(s), if not in pepXML directory" << endl;
   cout << "          -E         force reported endpoints to correspond to MS1 scans" << endl << endl;
   cout << "Example:  XPressPeptideParser interact"<<get_pepxml_dot_ext()<<" -L" << endl;
   exit(1);
}
