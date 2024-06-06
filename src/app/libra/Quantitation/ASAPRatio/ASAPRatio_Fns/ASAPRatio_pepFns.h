/*
Program       : ASAPRatioPeptideParser
Author        : X.Li and Andrew Keller <akeller@systemsbiology.org>
Date          : 11.27.02
SVN Info      : $Id: ASAPRatio_pepFns.h 8025 2020-02-14 01:05:59Z mhoopmann $


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

#ifndef _ASAPRATIO_PEPFNS_H_
#define _ASAPRATIO_PEPFNS_H_


#include "ASAPRatio_numFns.h"
#include "ASAPRatio_txtFns.h"

//
// Constants
//
#define _MZXML_READER_MXSCANNUM_ 100000
#define _ASAPRATIO_HM_ 1.0079407539 // average mass of H
#define _ASAPRATIO_M_ 1.0033548378 // mass diff. bewteen isotopes
#define _ASAPRATIO_MZBOUND_ 0.5 // bound of an isotopic peak
#define _ASAPRATIO_ISONUM_ 3 // number of isotopes covered
#define _ASAPRATIO_EXPLCRANGE_ 50 // selected range of a LC spectrum
#define NATIVEAANUM 22 // number of native amino acids + N-end + C-end
#define MODAANUM 20 // total number of modified amino acids

#include "Parsers/mzParser/mzParser.h"
#include "Common/spectStrct.h"


//
// Structures
//

// structure of XML index
typedef struct {
  mzParser::RAMPFILE *file;
  mzParser::ramp_fileoffset_t *scanIndx;
  int totScan; // 1 .. totScan
} xmlIndxStrct;

// structure of LC spectra
class lcSpectStrct;


// structure of amino acid residue
typedef struct {
  int sz; // 1 or 2 letter represent
  char rp[3]; // represent: 1 letter for native, 1-2 for modified 
  double ms; // mono-isotopic mass
} residueStrct;


// native amino acid residues
const static residueStrct nativeAA[] = {
  {1, "<", 1.0078250321},
  {1, ">", 17.0027396542},
  {1, "G", 57.0214637236},
  {1, "A", 71.0371137878},
  {1, "V", 99.0684139162},
  {1, "L", 113.0840639804},
  {1, "I", 113.0840639804},
  {1, "S", 87.0320284099},
  {1, "C", 103.0091844778},
  {1, "T", 101.0476784741},
  {1, "M", 131.0404846062},
  {1, "P", 97.0527638520},
  {1, "F", 147.0684139162},
  {1, "Y", 163.0633285383},
  {1, "W", 186.0793129535},
  {1, "H", 137.0589118624},
  {1, "K", 128.0949630177},
  {1, "R", 156.1011110281},
  {1, "D", 115.0269430320},
  {1, "E", 129.0425930962},
  {1, "N", 114.0429274472},
  {1, "Q", 128.0585775114}
};


// structure of pairing amino acid
typedef struct {
  char prtnA[3];
  char prtnB[3];
} pairStrct;


//
// Functions
//

// This function evaluates the light:heavy ratio of a peptide.
void getPepDataStrct(pepDataStrct *data, char *xmlFile,
		     char *pngFileBase, int cgiIndx);
void getPepDataStrct(pepDataStrct *data, char *xmlFile,
		     char *pngFileBase, int cgiIndx, 
		     Boolean quantHighBG, Boolean zeroBG, 
		     double mzBound, bool wavelet);
void getPepDataStrct(pepDataStrct *data, char *xmlFile,
		     char *pngFileBase, int cgiIndx, 
		     Boolean quantHighBG, Boolean zeroBG, 
		     double mzBound, int smoothItrs, 
		     double smoothRTwindow, bool wavelet);


#endif /* _ASAPRATIO_PEPFNS_H_ */

