/*
Program       : XPressPeptideUpdateParser
Author        : J.Eng and Andrew Keller <akeller@systemsbiology.org>, Robert Hubley, and
                open source code
Date          : 11.27.02
Version       : $Id: XPressPeptideUpdateParserMain.cpp 8886 2023-03-14 06:45:24Z real_procopio $


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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "Parsers/mzParser/mzParser.h"

/*
 * Include gd library to create png files
 * see http://www.boutell.com/gd/
 */
#include "gd.h"
#include "gdfonts.h"

#include "XPressPeptideUpdateParser.h"
#include "Parsers/Parser/Tag.h"

#include "Common/TPPVersion.h" // contains version number, name, revision
#include "Common/util.h"

#define szVERSION            "XPRESS"
#define szAUTHOR             "by J.Eng &copy; Institute for Systems Biology, 2000. All rights reserved."

using namespace mzParser;

int  bXpressLight1,
     bLabelFree,  // if set, ignores all Heavy* values, and displays simpler curation interface
     iCidScan;
char szXMLFile[SIZE_FILE],
     szOutputFile[SIZE_FILE],
     szInteractBaseName[SIZE_FILE],
     szInteractDir[SIZE_FILE];
double dProton = 1.00727646688;

struct QuanStruct {
  int  iChargeState;
  int  iLightFirstScan;
  int  iLightLastScan;
  int  iHeavyFirstScan;
  int  iHeavyLastScan;
  int  bPpmMassTol;
  int  bNormalize;
  int  iNumIsotopePeaks;
  double dLightPeptideMass;
  double dHeavyPeptideMass;
  double dLightQuanValue;
  double dHeavyQuanValue;
  double dMassTol;
  // LF
  double dLightFirstScanRT;
  double dLightLastScanRT;
  double dLightIntensity;
  double dLightIntensityRT;
  int iLightIntensityScan;

  char szNewQuan[20];
} pQuan;

void PRINT_FORM(char *szArgv0,
		struct QuanStruct *pQuan,
		char *szXMLFile,
		char *szOutputFile);
void PRINT_LF_FORM(char *szArgv0,
		struct QuanStruct *pQuan,
		char *szXMLFile,
		char *szOutputFile);
void EXTRACT_QUERY_STRING(struct QuanStruct *pQuan,
			  char *szXMLFile,
			  char *szOutputFile,
			  char *szInteractDir,
			  char *szInteractBaseName,
			  int  *bXpressLight1,
			  int  *bLabelFree,
			  int* index,
			  char* xmlfile,
			  int* overwrite,
			  int* quantitating,
			  char* zeroarea,
			  char* quan);
void GET_QUANTITATION(char *szXMLFile,
		      struct QuanStruct *pQuan);
void MAKE_PLOT(int iPlotStartScan,
	       int iPlotEndScan,
	       int iLightStartScan,
	       int iLightEndScan,
	       int iHeavyStartScan,
	       int iHeavyEndScan, 
	       double dMaxLightInten,
	       double dMaxHeavyInten,
	       double *pdLight,
	       double *pdHeavy,
	       double *dLightFiltered,
	       double *dHeavyFiltered,
	       char *szXMLFile,
	       struct QuanStruct *pQuan);
void FILTER_MS(double *dOrigMS,
	       double *dFilteredMS,
	       double *dTmpFilter,
	       double dMaxInten,
	       int    iNumMSScans,
	       int    iPlotStartScan,
	       int    iPlotEndScan); 
void UPDATE_QUAN(void);
void SPECIAL_QUAN();
void flipRatio(char* ratio, char* flipped);
static int GET_SCANNUM_FROM_DTA_NAME(const char* specName);

int update_index;
char update_xmlfile[1024];
int overwrite;
int quantitating;
char zero_area[10];
char update_quan[100];

int main(int argc, char **argv)
{
   hooks_tpp handler(argc,argv); // set up install paths etc

   // Print html header
   printf("Content-type: text/html\n\n");
   printf("<html>\n<head>\n");
   printf("<title>XPRESSPeptide: %s</title>\n", szAUTHOR);
   printf("<link rel='stylesheet' type='text/css' href='%scss/tpp.css'>\n", getHtmlUrl());
   printf("<script type='text/javascript' src='%sjs/tpp.js'></script>\n", getHtmlUrl());
   printf("<script language='JavaScript'>\n\
  var messages = [];\n\
  function showmsgs() {\n\
    var msgs = '';\n\
    for (var msg of messages) {\n\
      msgs += msg + '<br>';\n\
    }\n\
    if (msgs) tpp_showAlerts(msgs);\n\
  }\n");
   printf("</script>\n<head>\n");
   printf("<body onload='self.focus();showmsgs();'>\n<div id='tppWrapper'>\n");
   printf("<div id='tppFiller'><b>XPressPeptideUpdateParser</b> (%s) :: <b>loading...</b>\n",szTPPVersionInfo);

   szXMLFile[0]='\0'; 
   szOutputFile[0]='\0'; 
   szInteractDir[0]='\0'; 
   szInteractBaseName[0]='\0'; 
   bXpressLight1=0;
   bLabelFree=0;
   iCidScan=0;

   update_xmlfile[0] = 0;
   zero_area[0] = 0;
   update_quan[0] = 0;
   overwrite = 0;
   quantitating = 0;

   EXTRACT_QUERY_STRING(&pQuan, szXMLFile, szOutputFile, szInteractDir,
			szInteractBaseName, &bXpressLight1, &bLabelFree,
			&update_index, update_xmlfile, &overwrite, &quantitating, zero_area, update_quan);

   iCidScan=GET_SCANNUM_FROM_DTA_NAME(szOutputFile);

   printf("</div>\n"); // filler
   if (bLabelFree) {
     printf("<div class='tppbanner' banner-bg-text='XPRESS LF'>&nbsp;&nbsp;XPRESS Label Free :: <b class='tppOrange'>%s</b>\n</div>\n",szOutputFile);
     pQuan.iHeavyFirstScan = pQuan.iLightFirstScan; // for x-axis image scale
     pQuan.iHeavyLastScan  = pQuan.iLightLastScan; //
     PRINT_LF_FORM(argv[0], &pQuan, szXMLFile, szOutputFile);
   }
   else {
     printf("<div class='tppbanner' banner-bg-text='XPRESS Peptide'>&nbsp;&nbsp;XPRESS Peptide Ratio :: <b class='tppOrange'>%s</b>\n</div>\n",szOutputFile);
     PRINT_FORM(argv[0], &pQuan, szXMLFile, szOutputFile);
   }
   fflush(stdout);

   GET_QUANTITATION(szXMLFile, &pQuan);

   if (!bLabelFree) {
     printf("<table width='100%%'><tr width='100%%'>\n");

     printf("<td width='33%%' align='right'>");
     SPECIAL_QUAN();
     printf("</td><td width='33%%'>");

     printf("<center>");
     printf("<table><tr><td align='right'>Light</td><td>:</td><td align='left'>Heavy</td></tr>\n");

     if (pQuan.dLightQuanValue!=0.0)
       printf("<tr><td align='right'>1</td><td>:</td><td align='left'>%0.2f</td></tr>\n",
	      pQuan.dHeavyQuanValue / pQuan.dLightQuanValue);
     else
       printf("<tr><td align='right'>1</td><td>:</td><td align='left'>NaN</td></tr>\n");

     if (pQuan.dHeavyQuanValue!=0.0)
       printf("<tr><td align='right'>%0.2f</td><td>:</td><td align='left'>1</td></tr>\n",
	      pQuan.dLightQuanValue /  pQuan.dHeavyQuanValue);
     else
       printf("<tr><td align='right'>NaN</td><td>:</td><td align='left'>1</td></tr>\n");

     if (bXpressLight1==1) {
       if (pQuan.dLightQuanValue == 0.0)
	 sprintf(pQuan.szNewQuan, "1:INF");
       else
	 sprintf(pQuan.szNewQuan, "1:%0.2f", pQuan.dHeavyQuanValue / pQuan.dLightQuanValue);
     }
     else if (bXpressLight1==2) {
       if (pQuan.dHeavyQuanValue == 0.0)
	 sprintf(pQuan.szNewQuan, "INF:1");
       else
	 sprintf(pQuan.szNewQuan, "%0.2f:1", pQuan.dLightQuanValue / pQuan.dHeavyQuanValue);
     }
     else {
       if (pQuan.dLightQuanValue==0.0 && pQuan.dHeavyQuanValue==0.0)
	 sprintf(pQuan.szNewQuan, "?");
       else if (pQuan.dLightQuanValue > pQuan.dHeavyQuanValue)
	 sprintf(pQuan.szNewQuan, "1:%0.2f", pQuan.dHeavyQuanValue / pQuan.dLightQuanValue);
       else
	 sprintf(pQuan.szNewQuan, "%0.2f:1", pQuan.dLightQuanValue / pQuan.dHeavyQuanValue);
     }

     printf("</table></td><td width='33%%' align='left'>");
     UPDATE_QUAN();
     printf("</td></tr></table>\n");
   }
   else {
     printf("<center><table>");
     printf("<tr><td align='right'><b>Peak Area</b> :</td><td align='left'>%0.2e</td></tr>\n",pQuan.dLightQuanValue);
     printf("<tr><td align='right'><b>RT Start-End</b> :</td><td align='left'>%0.0lf -- %0.0lf</td></tr>\n",pQuan.dLightFirstScanRT,pQuan.dLightLastScanRT);
     printf("<tr><td align='right'><b>Peak Intensity</b> :</td><td align='left'>%0.2e</td></tr>\n",pQuan.dLightIntensity);
     printf("<tr><td align='right'><b>Peak RT (scan)</b> :</td><td align='left'>%0.0lf (%d)</td></tr>\n",pQuan.dLightIntensityRT,pQuan.iLightIntensityScan);
     printf("</table><br>");
     UPDATE_QUAN();
     printf("</center>");
   }

   printf("<br style='clear:both'>\n<br/><br/><br/><br/>");
   printf("<footer id='tppPageFooter'>XPressPeptideUpdateParser :: <b>%s</b><br/>\n", szOutputFile);
   printf("<b>%s</b><br>\n", szXMLFile);
   printf("%s<br/></footer>\n", szTPPVersionInfo);
   printf("</div></body>\n</html>");

   return(0);
} /*main*/


void PRINT_LF_FORM(char *szArgv0,
		struct QuanStruct *pQuan,
		char *szXMLFile,
		char *szOutputFile)
{
  printf("<form method='GET' action='%s'>", getenv("SCRIPT_NAME"));
  printf("<input type='hidden' NAME='LabelFree' VALUE='1'>\n"); // bLabelFree==1 got us here
  printf("<input type='hidden' name='InteractDir' VALUE='%s'>\n", szInteractDir);
  printf("<input type='hidden' name='InteractBaseName' VALUE='%s'>\n", szInteractBaseName);
  printf("<input type='hidden' name='bXpressLight1' VALUE='%d'>\n", bXpressLight1);
  printf("<input type='hidden' name='xmlfile' VALUE='%s'>", update_xmlfile);
  printf("<input type='hidden' name='index' VALUE='%d'>", update_index);
  printf("<input type='hidden' name='XMLFile' VALUE='%s'>\n", szXMLFile);
  printf("<input type='hidden' name='OutFile' VALUE='%s'>", szOutputFile);
  printf("<input type='hidden' name='ChargeState' VALUE='%d'>", pQuan->iChargeState);
  printf("<input type='hidden' NAME='LightMass' VALUE='%0.4f'>\n", pQuan->dLightPeptideMass);
  printf("<input type='hidden' name='MassTol' VALUE='%0.4f'>", pQuan->dMassTol);
  printf("<input type='hidden' name='PpmTol' VALUE='%d'>", pQuan->bPpmMassTol);
  printf("<input type='hidden' name='NumIsotopePeaks' VALUE='%d'>", pQuan->iNumIsotopePeaks);

  printf("<center><table class='tppsimpletable' style='min-width:700px;'>\n");
  printf("<tr><th class='subhead' style='position:initial;' colspan='3'>Label-Free Analysis Parameters</th></tr>\n");
  printf("<tr><td class='name'>Start scan</td><td><input type='textarea' name='LightFirstScan' value='%d' size='8'>\n",pQuan->iLightFirstScan);
  printf("<td rowspan='2'><input style='float:right;' type='submit' name='quantitating' value='Quantitate'></td></tr>\n");
  printf("<tr><td class='name'>End scan</td><td><input type='textarea' name='LightLastScan' value='%d' size='8'></td></tr>\n",pQuan->iLightLastScan);
  printf("<tr><td class='name'>Mass / Charge</td><td colspan='2'>%0.4f / +%d</td></tr>\n",pQuan->dLightPeptideMass, pQuan->iChargeState);
  printf("<tr><td class='name'>Tolerance</td><td colspan='2'>%0.4f %s</td></tr>\n",pQuan->dMassTol, pQuan->bPpmMassTol ? "ppm" : "Da");
  printf("<tr><td class='name'>#C13 peaks</td><td colspan='2'>%d</td></tr>\n",pQuan->iNumIsotopePeaks);
  printf("<tr><td class='name'>File</td><td colspan='2'>%s</td></tr>\n",szXMLFile);
  printf("<tr><td class='name'>Spectrum</td><td colspan='2'>%s</td></tr>\n",szOutputFile);
  printf("</table></center><br>\n");
  printf("</form>\n");

}

void PRINT_FORM(char *szArgv0,
		struct QuanStruct *pQuan,
		char *szXMLFile,
		char *szOutputFile)
{

   printf("<FORM METHOD='GET' ACTION='%s'>", getenv("SCRIPT_NAME"));
   printf("<INPUT TYPE='hidden' NAME='InteractDir' VALUE='%s'>\n", szInteractDir);
   printf("<INPUT TYPE='hidden' NAME='InteractBaseName' VALUE='%s'>\n", szInteractBaseName);
   printf("<INPUT TYPE='hidden' NAME='bXpressLight1' VALUE='%d'>\n", bXpressLight1);
   printf("<input TYPE='hidden' NAME='xmlfile' VALUE='%s'>", update_xmlfile);
   printf("<input TYPE='hidden' NAME='index' VALUE='%d'>", update_index);

   printf("<CENTER><TABLE BORDER=0 CELLPADDING=5>\n");

   printf("<TR VALIGN=TOP><TD BGCOLOR='#EEEEEE'>\n");
   printf("<TT>Light:");
   printf(" <INPUT TYPE='textarea' NAME='LightFirstScan' VALUE='%d' SIZE=8>\n", pQuan->iLightFirstScan);
   printf(" <INPUT TYPE='textarea' NAME='LightLastScan'  VALUE='%d' SIZE=8>\n", pQuan->iLightLastScan);
   printf(" mass:<INPUT TYPE='textarea' NAME='LightMass' VALUE='%0.4f' SIZE=8>\n", pQuan->dLightPeptideMass);
   printf(" tol: <INPUT TYPE='textarea' NAME='MassTol' VALUE='%0.4f' SIZE=3>", pQuan->dMassTol);
   printf(" &nbsp;<label>ppm:<INPUT TYPE='checkbox' NAME='PpmTol' VALUE='1' %s></label>\n", (pQuan->bPpmMassTol?"checked":""));
   printf(" &nbsp;#C13:<select name='NumIsotopePeaks'>\n");
   printf("<option value='0'%s>0</option>\n", (pQuan->iNumIsotopePeaks==0?" selected":""));
   printf("<option value='1'%s>1</option>\n", (pQuan->iNumIsotopePeaks==1?" selected":""));
   printf("<option value='2'%s>2</option>\n", (pQuan->iNumIsotopePeaks==2?" selected":""));
   printf("</select>\n");

   printf("<BR>Heavy:");
   printf(" <INPUT TYPE='textarea' NAME='HeavyFirstScan' VALUE='%d' SIZE=8>\n", pQuan->iHeavyFirstScan);
   printf(" <INPUT TYPE='textarea' NAME='HeavyLastScan'  VALUE='%d' SIZE=8>\n", pQuan->iHeavyLastScan);
   printf(" mass:<INPUT TYPE='text' NAME='HeavyMass' VALUE='%0.4f' SIZE=8>\n", pQuan->dHeavyPeptideMass);
   printf(" &nbsp;&nbsp;Z: <INPUT TYPE='textarea' NAME='ChargeState' VALUE='%d' SIZE=3>", pQuan->iChargeState);
   printf(" <label>norm:<INPUT TYPE='checkbox' NAME='Norm' VALUE='1' %s></label>\n", (pQuan->bNormalize?"checked":""));

   printf("<P> file:  <input NAME='XMLFile' type=textarea VALUE='%s' SIZE='65'>\n", szXMLFile);
   printf("<BR> spec:  <input NAME='OutFile' type=textarea VALUE='%s' SIZE='65'>", szOutputFile);


   printf("</tt></td>\n");
   printf("<td bgcolor='#eeeeee' valign='middle'>\n<input type='submit' name='quantitating' value='Quantitate'>");
   printf("</TD>\n");

   printf("</TR>\n");
   printf("</TABLE></CENTER>\n");
   printf("</FORM>\n");

}

void EXTRACT_QUERY_STRING(struct QuanStruct *pQuan,
			  char *szXMLFile,
			  char *szOutputFile,
			  char *szInteractDir,
			  char *szInteractBaseName,
			  int  *bXpressLight1,
			  int  *bLabelFree,
			  int* index,
			  char* xmlfile,
			  int* overwrite,
			  int* quantitating,
			  char* zeroarea,
			  char* quan)
{
   const char *pStr = getenv("REQUEST_METHOD");

   if (pStr==NULL)  {
     printf(" Error - this is a CGI program that cannot be\n");
     printf(" run from the command line.\n\n");
     exit(EXIT_FAILURE);
   }
   else if (!strcmp(pStr, "GET")) {
     int  i;
     char *szQuery,
           szWord[1024];

     // Old files will not have PpmTol attribute,
     // therefore we initialize it to 0 (Dalton by default)
     pQuan->bPpmMassTol = 0;
     pQuan->iNumIsotopePeaks = 0;

     // Get:
     //       ChargeState - charge state of precursor
     szQuery = getenv("QUERY_STRING");

     if (szQuery == NULL) {
       printf("<P>No query information to decode.\n");
       printf("</BODY>\n</HTML>\n");
       exit(EXIT_FAILURE);
     }

     for (i=0; szQuery[0] != '\0'; i++) {
       getword(szWord, szQuery, '=');
       plustospace(szWord);
       unescape_url(szWord);

       if (!strcmp(szWord, "LightFirstScan")) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%d", &(pQuan->iLightFirstScan));
       }
       else if (!strcmp(szWord, "LightLastScan") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%d", &(pQuan->iLightLastScan));
       }
       else if (!strcmp(szWord, "HeavyFirstScan") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%d", &(pQuan->iHeavyFirstScan));
       }
       else if (!strcmp(szWord, "HeavyLastScan") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%d", &(pQuan->iHeavyLastScan));
       }
       else if (!strcmp(szWord, "XMLFile") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%s", szXMLFile);
	 fixPath(szXMLFile,0); // tidy up path separators etc - expect existence
       }
       else if (!strcmp(szWord, "OutFile") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%s", szOutputFile);
	 fixPath(szOutputFile,0); // tidy up path separators etc
       }
       else if (!strcmp(szWord, "InteractDir") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%s", szInteractDir);
	 fixPath(szInteractDir,1); // tidy up path separators etc - expect existence
       }
       else if (!strcmp(szWord, "InteractBaseName") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%s", szInteractBaseName);
	 fixPath(szInteractBaseName,0); // tidy up path separators etc
       }
       else if (!strcmp(szWord, "ChargeState") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%d", &(pQuan->iChargeState));
       }
       else if (!strcmp(szWord, "LightMass") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%lf", &(pQuan->dLightPeptideMass));
       }
       else if (!strcmp(szWord, "HeavyMass") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%lf", &(pQuan->dHeavyPeptideMass));
       }
       else if (!strcmp(szWord, "MassTol") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%lf", &(pQuan->dMassTol));
       }
       else if (!strcmp(szWord, "PpmTol") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%d", &(pQuan->bPpmMassTol));
       }
       else if (!strcmp(szWord, "NumIsotopePeaks") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%d", &(pQuan->iNumIsotopePeaks));
       }
       else if (!strcmp(szWord, "bXpressLight1") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%d", bXpressLight1);
       }
       else if (!strcmp(szWord, "LabelFree") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%d", bLabelFree);
       }
       else if (!strcmp(szWord, "Norm") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%d", &(pQuan->bNormalize));
       }
       else if (!strcmp(szWord, "index") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%d", index);
       }
       else if (!strcmp(szWord, "xmlfile") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 strcpy(xmlfile, szWord);
	 fixPath(xmlfile,1); // tidy up path separators etc - expect existence
	 //sscanf(szWord, "%d", index);
       }
       else if (!strcmp(szWord, "overwrite") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%d", overwrite);
       }
       else if (!strcmp(szWord, "quantitating") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 *quantitating=1;
       }
       else if (!strcmp(szWord, "ZeroArea") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 strcpy(zero_area, szWord);
       }
       else if (!strcmp(szWord, "NewQuan") ) {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
	 strcpy(quan, szWord);
       }
       else {
	 getword(szWord, szQuery, '&'); plustospace(szWord); unescape_url(szWord);
       }
     }
   }
   else {
     printf(" Error not called with GET method\n");
     exit(EXIT_FAILURE);
   }
   // make sure we can do this
   rampValidateOrDeriveInputFilename(szXMLFile,SIZE_FILE,szOutputFile);
} /*EXTRACT_QUERY_STRING*/



// Reads .mzXML files and get quantitation numbers
void GET_QUANTITATION(char *szXMLFile,
        struct QuanStruct *pQuan)
{
   int    i,
          ctScan,
          iAnalysisFirstScan,
          iAnalysisLastScan,
          iLightStartScan,
          iLightEndScan,
          iHeavyStartScan,
          iHeavyEndScan,
          iPlotStartScan,
          iPlotEndScan,
          iStart,
          iEnd;

#define MAX_ISOTOPES 2

   double dLightMass,
          dHeavyMass,
          *dLightMS,
          *dHeavyMS,
          *dLightFilteredMS,
          *dHeavyFilteredMS,
          *dTmpMS,
          *dIsotopeLight[MAX_ISOTOPES+1],
          *dIsotopeHeavy[MAX_ISOTOPES+1],
          dMaxLightInten,
          dMaxHeavyInten,
          dCurrTime;

   struct ScanHeaderStruct scanHeader;

   ramp_fileoffset_t   indexOffset;
   ramp_fileoffset_t  *pScanIndex;
   RAMPFILE   *pFI;

   // confirm presence & open corresponding mzXML file

   if ( (pFI=rampOpenFile(szXMLFile))==NULL) {
     char szXMLFile2[SIZE_FILE];

     // Not ideal but works. XPressCGIProteinDisplayParser.cgi current does not
     // have mzXML/mzML file extension passed to it at all.  So xpress link
     // that it generates always has a .mzXML extension.  Simple fix here is
     // to try given XML file input and, if not found, the corresponding one.
     if (strstr(szXMLFile, ".mzXML")) {
       strcpy(szXMLFile2, szXMLFile);
       szXMLFile2[strlen(szXMLFile2)-5]='\0';
       strcat(szXMLFile2, "mzML");
       if ( (pFI=rampOpenFile(szXMLFile2))==NULL) {
	 printf(" Error - cannot open file \"%s\".\n\n", szXMLFile);
	 exit(0);
       }
     }
     else if (strstr(szXMLFile, ".mzML")) {
       strcpy(szXMLFile2, szXMLFile);
       szXMLFile2[strlen(szXMLFile2)-4]='\0';
       strcat(szXMLFile2, "mzXML");
       if ( (pFI=rampOpenFile(szXMLFile2))==NULL) {
	 printf(" Error - cannot open file \"%s\".\n\n", szXMLFile);
	 exit(0);
       }
     }
     else {
       printf(" Error - cannot open file \"%s\".\n\n", szXMLFile);
       exit(0);
     }
   }

   // Get the analysis start, end scans

   // Read the offset of the index
   indexOffset = getIndexOffset( pFI );
   
   // Read the scan index into a vector, get LastScan
   pScanIndex = readIndex( pFI , indexOffset, &iAnalysisLastScan );
   iAnalysisFirstScan = 1;

   if ( (dLightMS=(double *)malloc(sizeof(double)*(iAnalysisLastScan+5)))==NULL) {
     printf(" Error, cannot malloc dLightMS[%d]\n\n", iAnalysisLastScan+5);
     exit(EXIT_FAILURE);
   }
   if ( (dHeavyMS=(double *)malloc(sizeof(double)*(iAnalysisLastScan+5)))==NULL) {
     printf(" Error, cannot malloc dHeavyMS[%d]\n\n", iAnalysisLastScan+5);
     exit(EXIT_FAILURE);
   }
   if ( (dLightFilteredMS=(double *)malloc(sizeof(double)*(iAnalysisLastScan+5)))==NULL) {
     printf(" Error, cannot malloc dLightFilteredMS[%d]\n\n", iAnalysisLastScan+5);
     exit(EXIT_FAILURE);
   }
   if ( (dHeavyFilteredMS=(double *)malloc(sizeof(double)*(iAnalysisLastScan+5)))==NULL) {
     printf(" Error, cannot malloc dHeavyFilteredMS[%d]\n\n", iAnalysisLastScan+5);
     exit(EXIT_FAILURE);
   }
   if ( (dTmpMS=(double *)malloc(sizeof(double)*(iAnalysisLastScan+5)))==NULL) {
     printf(" Error, cannot malloc dTmpMS[%d]\n\n", iAnalysisLastScan+5);
     exit(EXIT_FAILURE);
   }

   for (int n=0 ; n<=MAX_ISOTOPES; n++) {
     if ((dIsotopeLight[n] = new double [iAnalysisLastScan+5])==NULL) {
       printf("Error malloc dIsotopeLight[%d][%d]\n", n , iAnalysisLastScan+5);
       exit(1);
     }

     if ((dIsotopeHeavy[n] = new double [iAnalysisLastScan+5])==NULL) {
       printf("Error malloc dIsotopeHeavy[%d][%d]\n", n , iAnalysisLastScan+5);
       exit(1);
     }
   }

   // Calculate the precursor mass
   if (pQuan->iChargeState < 1) {
     printf(" Error, charge state = %d\n\n", pQuan->iChargeState);
     exit(EXIT_FAILURE);
   }
   else {
     dLightMass = (dProton * (pQuan->iChargeState - 1) + pQuan->dLightPeptideMass) / pQuan->iChargeState;
     dHeavyMass = (dProton * (pQuan->iChargeState - 1) + pQuan->dHeavyPeptideMass) / pQuan->iChargeState;
   }

   // Clear all values
   memset(dLightMS, 0, sizeof(double)*(iAnalysisLastScan+5));
   memset(dHeavyMS, 0, sizeof(double)*(iAnalysisLastScan+5));
   memset(dLightFilteredMS, 0, sizeof(double)*(iAnalysisLastScan+5));
   memset(dHeavyFilteredMS, 0, sizeof(double)*(iAnalysisLastScan+5));

   for (int n=0 ; n<=pQuan->iNumIsotopePeaks; n++) {
     memset(dIsotopeLight[n], 0, (iAnalysisLastScan+5)*sizeof(double));
     memset(dIsotopeHeavy[n], 0, (iAnalysisLastScan+5)*sizeof(double));
   }

   iStart = pQuan->iLightFirstScan; 
   readHeader( pFI, pScanIndex[iStart], &scanHeader );
   dCurrTime = scanHeader.retentionTime;
   pQuan->dLightFirstScanRT = dCurrTime;
   while (iStart>iAnalysisFirstScan && dCurrTime - scanHeader.retentionTime < 560.0 ) {
     iStart--;
     readHeader( pFI, pScanIndex[iStart], &scanHeader );
   }

   iEnd = pQuan->iLightLastScan;
   readHeader( pFI, pScanIndex[iEnd], &scanHeader );
   dCurrTime = scanHeader.retentionTime;
   pQuan->dLightLastScanRT = dCurrTime;
   while (iEnd<iAnalysisLastScan  && scanHeader.retentionTime-dCurrTime < 560.0 ) {
     iEnd++;
     readHeader( pFI, pScanIndex[iEnd], &scanHeader );
   }

   // Read all MS scan values
   double dLightMassTol;
   double dHeavyMassTol;

   if( pQuan->bPpmMassTol ) {
     dLightMassTol = dLightMass * pQuan->dMassTol / 1000000;
     dHeavyMassTol = dHeavyMass * pQuan->dMassTol / 1000000;       
   }
   else {
     dLightMassTol = pQuan->dMassTol;
     dHeavyMassTol = pQuan->dMassTol;
   }

   double dC13diff = 1.00335483;

   for (ctScan=iStart; ctScan<=iEnd; ctScan++) {
     readHeader( pFI, pScanIndex[ctScan], &scanHeader );

     if (scanHeader.msLevel == 1) {
       RAMPREAL *pPeaks;
       int n = 0;

       // Open a scan
       pPeaks = readPeaks( pFI, pScanIndex[ctScan] );
   
       for (int x=0; x< scanHeader.peaksCount; x++ ) {
	 RAMPREAL fMass;
	 RAMPREAL fInten;

	 fMass=pPeaks[n];
	 n++;
	 fInten=pPeaks[n];
	 n++;

	 for (int ii=0; ii<=pQuan->iNumIsotopePeaks; ii++) {
	   if ( fabs(dLightMass + ii*dC13diff - fMass) <= dLightMassTol) {
	     if (fInten > dIsotopeLight[ii][ctScan])
	       dIsotopeLight[ii][ctScan] = (double) fInten;
	   }

	   if ( fabs(dHeavyMass + ii*dC13diff - fMass) <= dHeavyMassTol) {
	     if (fInten > dIsotopeHeavy[ii][ctScan])
	       dIsotopeHeavy[ii][ctScan] = (double) fInten;
	   }
	 }

       }

       for (int ii=0; ii<=pQuan->iNumIsotopePeaks; ii++) {
	 dLightMS[ctScan] += dIsotopeLight[ii][ctScan];
	 dHeavyMS[ctScan] += dIsotopeHeavy[ii][ctScan];
       }

       //free(pPeaks);
     }
   }

   // Starting from the start and end scans,
   // need to see the real start/end scan of eluting
   // peptide by looking at smoothed/filtered MS profile.

   // Get peptide start & end scans
   iLightStartScan = pQuan->iLightFirstScan;
   iLightEndScan = pQuan->iLightLastScan;
   iHeavyStartScan = pQuan->iHeavyFirstScan;
   iHeavyEndScan = pQuan->iHeavyLastScan;

   // Print out data for plotting
   iPlotStartScan = iLightStartScan;
   iPlotEndScan = iLightEndScan;

   if (iHeavyStartScan < iPlotStartScan)
      iPlotStartScan = iHeavyStartScan;
   if (iHeavyEndScan > iPlotEndScan)
      iPlotEndScan = iHeavyEndScan;

   readHeader( pFI, pScanIndex[iPlotStartScan], &scanHeader );
   dCurrTime = scanHeader.retentionTime;
   while (iPlotStartScan>iAnalysisFirstScan && dCurrTime - scanHeader.retentionTime < 45.0 ) {
     iPlotStartScan--;
     readHeader( pFI, pScanIndex[iPlotStartScan], &scanHeader );
   }

   readHeader( pFI, pScanIndex[iPlotEndScan], &scanHeader );
   dCurrTime = scanHeader.retentionTime;
   while (iPlotEndScan<iAnalysisLastScan  && scanHeader.retentionTime-dCurrTime < 45.0 ) {
     iPlotEndScan++;
     readHeader( pFI, pScanIndex[iPlotEndScan], &scanHeader );
   }

   dMaxLightInten=0.0;
   dMaxHeavyInten=0.0;
   for (i=iPlotStartScan; i<=iPlotEndScan; i++) {
     if (dLightMS[i]>dMaxLightInten)
       dMaxLightInten=dLightMS[i];
     if (dHeavyMS[i]>dMaxHeavyInten)
       dMaxHeavyInten=dHeavyMS[i];
   }

   // Sum up intensity values across chromatrogram
   for (i=iHeavyStartScan; i<=iHeavyEndScan; i++)
     pQuan->dHeavyQuanValue += dHeavyMS[i];
   for (i=iLightStartScan; i<=iLightEndScan; i++) {
     pQuan->dLightQuanValue += dLightMS[i];
     // for label-free:
     if (dLightMS[i] > pQuan->dLightIntensity) {
       pQuan->dLightIntensity = dLightMS[i];
       pQuan->iLightIntensityScan = i;
     }
   }
   readHeader( pFI, pScanIndex[pQuan->iLightIntensityScan], &scanHeader );
   pQuan->dLightIntensityRT = scanHeader.retentionTime;


   memset(dTmpMS, 0, sizeof(double)*(iAnalysisLastScan+5));
   FILTER_MS(dLightMS, dLightFilteredMS, dTmpMS, dMaxLightInten, iAnalysisLastScan+5, iPlotStartScan, iPlotEndScan);

   memset(dTmpMS, 0, sizeof(double)*(iAnalysisLastScan+5));
   FILTER_MS(dHeavyMS, dHeavyFilteredMS, dTmpMS, dMaxHeavyInten, iAnalysisLastScan+5, iPlotStartScan, iPlotEndScan); 

   if (dLightMS[iLightStartScan]!=0.0)
      iLightStartScan--;
   if (dLightMS[iLightEndScan]!=0.0)
      iLightEndScan++;

   if (dHeavyMS[iHeavyStartScan]!=0.0)
      iHeavyStartScan--;
   if (dHeavyMS[iHeavyEndScan]!=0.0)
      iHeavyEndScan++;
   
   MAKE_PLOT(iPlotStartScan, iPlotEndScan, iLightStartScan, iLightEndScan,
      iHeavyStartScan, iHeavyEndScan, dMaxLightInten, dMaxHeavyInten,
      dLightMS, dHeavyMS, dLightFilteredMS, dHeavyFilteredMS, szXMLFile, pQuan);

   free(pScanIndex);
   rampCloseFile(pFI);

   free(dLightMS);
   free(dHeavyMS);
   free(dLightFilteredMS);
   free(dHeavyFilteredMS);
   free(dTmpMS);

} /*GET_QUANTITATION*/


void MAKE_PLOT(int iPlotStartScan,
        int iPlotEndScan,
        int iLightStartScan,
        int iLightEndScan,
        int iHeavyStartScan,
        int iHeavyEndScan, 
        double dMaxLightInten,
        double dMaxHeavyInten,
        double *pdLight,
        double *pdHeavy,
        double *pdLightFiltered,
        double *pdHeavyFiltered,
        char *szXMLFile,
        struct QuanStruct *pQuan)
{
   int  i,
        iImageWidth=700,
        iImageHeight=250,
        iBottomBorder=20,
        iTopBorder=10,
        iGrey,
        iGrey2,
        iWhite,
        iWhite2,
        iBlack,
        iBlack2,
        iRed,
        iRed2,
        iBlue,
        iBlue2,
        rand1,
        rand2;
   char szImageFile[SIZE_FILE],
        szImageFile2[SIZE_FILE],
        szImageDir[SIZE_FILE],
        szImageLink[SIZE_FILE],
        szLabel[SIZE_FILE];
   const char *pStr;
   FILE *fp;
   time_t tStartTime;

   gdImagePtr gdImageLight,
              gdImageHeavy;

   if (bLabelFree)
     iImageHeight*=1.7;  // no Heavy image

   // first color allocated defines background color
   gdImageLight = gdImageCreate(iImageWidth, iImageHeight);   
   iWhite =  gdImageColorAllocate(gdImageLight, 255, 255, 255);
   iGrey  =  gdImageColorAllocate(gdImageLight, 150, 150, 150);
   iBlack =  gdImageColorAllocate(gdImageLight, 0, 0, 0);
   iRed   =  gdImageColorAllocate(gdImageLight, 255, 0, 0);
   iBlue  =  gdImageColorAllocate(gdImageLight, 0, 0, 255);

   gdImageHeavy = gdImageCreate(iImageWidth, iImageHeight);   
   iWhite2 = gdImageColorAllocate(gdImageHeavy, 255, 255, 255);
   iGrey2  = gdImageColorAllocate(gdImageHeavy, 150, 150, 150);
   iBlack2 = gdImageColorAllocate(gdImageHeavy, 0, 0, 0);
   iRed2   = gdImageColorAllocate(gdImageHeavy, 255, 0, 0);
   iBlue2  = gdImageColorAllocate(gdImageHeavy, 0, 0, 255);

   if (bLabelFree)
     sprintf(szLabel, "Prec     : %0.3f m/z, %+d",
	     (pQuan->dLightPeptideMass + dProton*(pQuan->iChargeState -1)) /(pQuan->iChargeState), pQuan->iChargeState);
   else
     sprintf(szLabel, "Light    : %0.3f m/z, %+d",
	     (pQuan->dLightPeptideMass + dProton*(pQuan->iChargeState -1)) /(pQuan->iChargeState), pQuan->iChargeState);
   gdImageString(gdImageLight, gdFontSmall, 3, 3, (unsigned char *)szLabel, iBlue);
   sprintf(szLabel, "area     : %0.2E", pQuan->dLightQuanValue);
   gdImageString(gdImageLight, gdFontSmall, 3, 15, (unsigned char *)szLabel, iBlue);
   sprintf(szLabel, "max inten: %0.2E", dMaxLightInten);
   gdImageString(gdImageLight, gdFontSmall, 3, 27, (unsigned char *)szLabel, iBlue);
 
   sprintf(szLabel, "Heavy    : %0.3f m/z, %+d",
         (pQuan->dHeavyPeptideMass + dProton*(pQuan->iChargeState -1)) /(pQuan->iChargeState), pQuan->iChargeState);
   gdImageString(gdImageHeavy, gdFontSmall, 3, 3, (unsigned char *)szLabel, iBlue2);
   sprintf(szLabel,"area     : %0.2E", pQuan->dHeavyQuanValue);
   gdImageString(gdImageHeavy, gdFontSmall, 3, 15, (unsigned char *)szLabel, iBlue2);
   sprintf(szLabel,"max inten: %0.2E", dMaxHeavyInten);
   gdImageString(gdImageHeavy, gdFontSmall, 3, 27, (unsigned char *)szLabel, iBlue2); 

   if (iCidScan > 0) {
     sprintf(szLabel, "CID");
     gdImageString(gdImageLight, gdFontSmall, 16, 43, (unsigned char *)szLabel, iBlue);
     gdImageFilledEllipse(gdImageLight, 7, 50, 8, 8, iRed);
     gdImageString(gdImageHeavy, gdFontSmall, 16, 43, (unsigned char *)szLabel, iBlue2);
     gdImageFilledEllipse(gdImageHeavy, 7, 50, 8, 8, iRed);
   }


   if (pQuan->bNormalize) {
     if (dMaxLightInten > dMaxHeavyInten)
       dMaxHeavyInten = dMaxLightInten;
     else if (dMaxLightInten < dMaxHeavyInten)
       dMaxLightInten = dMaxHeavyInten;
   }

   // Plot out spectra
   for (i=iPlotStartScan; i<=iPlotEndScan; i++) {
     int iX1Pos, iY1PosLight, iY1PosHeavy;

     iX1Pos= (int)( (i-iPlotStartScan)*iImageWidth/(iPlotEndScan-iPlotStartScan));

     if (dMaxLightInten>0.0)
       iY1PosLight = iImageHeight - iBottomBorder -
	 (int)((pdLight[i]/dMaxLightInten)*(iImageHeight-iBottomBorder-iTopBorder));
     else
       iY1PosLight = iImageHeight - iBottomBorder;

     if (dMaxHeavyInten>0.0)
       iY1PosHeavy = iImageHeight - iBottomBorder -
	 (int)((pdHeavy[i]/dMaxHeavyInten)*(iImageHeight-iBottomBorder-iTopBorder));
     else
       iY1PosHeavy = iImageHeight - iBottomBorder;

     // plot light chromatogram points
     if (i>=iLightStartScan && i<iLightEndScan)
       gdImageLine(gdImageLight, iX1Pos, iImageHeight-iBottomBorder, iX1Pos, iY1PosLight, iRed);
     else
       gdImageLine(gdImageLight, iX1Pos, iImageHeight-iBottomBorder, iX1Pos, iY1PosLight, iGrey);

     // plot heavy chromatogram points
     if (i>=iHeavyStartScan && i<iHeavyEndScan)
       gdImageLine(gdImageHeavy, iX1Pos, iImageHeight-iBottomBorder, iX1Pos, iY1PosHeavy, iRed2);
     else
       gdImageLine(gdImageHeavy, iX1Pos, iImageHeight-iBottomBorder, iX1Pos, iY1PosHeavy, iGrey2);

     // Plot out smoothed trace
     if (dMaxLightInten>0.0) {
       iY1PosLight = iImageHeight - iBottomBorder -
	 (int)((pdLightFiltered[i]/dMaxLightInten)*(iImageHeight-iBottomBorder-iTopBorder));
     }
     else {
       iY1PosLight= iImageHeight - iBottomBorder;
     }

     if (dMaxHeavyInten>0.0) {
       iY1PosHeavy = iImageHeight - iBottomBorder -
	 (int)((pdHeavyFiltered[i]/dMaxHeavyInten)*(iImageHeight-iBottomBorder-iTopBorder));
     }
     else {
       iY1PosHeavy = iImageHeight - iBottomBorder;
     }

     if (i>=iLightStartScan && i<iLightEndScan)
       gdImageSetPixel(gdImageLight, iX1Pos, iY1PosLight, iBlue);
     else
       gdImageSetPixel(gdImageLight, iX1Pos, iY1PosLight, iRed);

     if (i>=iHeavyStartScan && i<iHeavyEndScan)
       gdImageSetPixel(gdImageHeavy, iX1Pos, iY1PosHeavy,  iBlue2);
     else
       gdImageSetPixel(gdImageHeavy, iX1Pos, iY1PosHeavy,  iRed2); 


     sprintf(szLabel, "%d", i);

     // x-label, tick marks
     if (iPlotEndScan-iPlotStartScan<150) {
       if ( !(i%10)) {
	 gdImageString(gdImageLight, gdFontSmall, iX1Pos-10, iImageHeight-13, (unsigned char *)szLabel, iBlack);
	 gdImageString(gdImageHeavy, gdFontSmall, iX1Pos-10, iImageHeight-13, (unsigned char *)szLabel, iBlack2);
       }
       if ( !(i%5)) {
	 // big tick marks
	 gdImageLine(gdImageLight, iX1Pos, iImageHeight-iBottomBorder, iX1Pos, iImageHeight-iBottomBorder+5, iBlack);
	 gdImageLine(gdImageHeavy, iX1Pos, iImageHeight-iBottomBorder, iX1Pos, iImageHeight-iBottomBorder+5, iBlack2);
       }

       // tick marks
       gdImageLine(gdImageLight, iX1Pos, iImageHeight-iBottomBorder, iX1Pos, iImageHeight-iBottomBorder+2, iBlack);
       gdImageLine(gdImageHeavy, iX1Pos, iImageHeight-iBottomBorder, iX1Pos, iImageHeight-iBottomBorder+2, iBlack2);
     }
     else if  (iPlotEndScan-iPlotStartScan<500) {
       if ( !(i%50)) {
	 gdImageString(gdImageLight, gdFontSmall, iX1Pos-10, iImageHeight-13, (unsigned char *)szLabel, iBlack);
	 gdImageString(gdImageHeavy, gdFontSmall, iX1Pos-10, iImageHeight-13, (unsigned char *)szLabel, iBlack2);
       }
       if ( !(i%10)) {
	 // big tick marks
	 gdImageLine(gdImageLight, iX1Pos, iImageHeight-iBottomBorder, iX1Pos, iImageHeight-iBottomBorder+5, iBlack);
	 gdImageLine(gdImageHeavy, iX1Pos, iImageHeight-iBottomBorder, iX1Pos, iImageHeight-iBottomBorder+5, iBlack2);
       }
       if ( !(i%5)) {
	 // tick marks
	 gdImageLine(gdImageLight, iX1Pos, iImageHeight-iBottomBorder, iX1Pos, iImageHeight-iBottomBorder+2, iBlack);
	 gdImageLine(gdImageHeavy, iX1Pos, iImageHeight-iBottomBorder, iX1Pos, iImageHeight-iBottomBorder+2, iBlack2);
       }
     } 
     else {
       if ( !(i%100)) {
	 gdImageString(gdImageLight, gdFontSmall, iX1Pos-10, iImageHeight-13, (unsigned char *)szLabel, iBlack);
	 gdImageString(gdImageHeavy, gdFontSmall, iX1Pos-10, iImageHeight-13, (unsigned char *)szLabel, iBlack2);
       }
       if ( !(i%50)) {
	 // big tick marks
	 gdImageLine(gdImageLight, iX1Pos, iImageHeight-iBottomBorder, iX1Pos, iImageHeight-iBottomBorder+5, iBlack);
	 gdImageLine(gdImageHeavy, iX1Pos, iImageHeight-iBottomBorder, iX1Pos, iImageHeight-iBottomBorder+5, iBlack2);
       }
       if ( !(i%10)) {
	 // tick marks
	 gdImageLine(gdImageLight, iX1Pos, iImageHeight-iBottomBorder, iX1Pos, iImageHeight-iBottomBorder+2, iBlack);
	 gdImageLine(gdImageHeavy, iX1Pos, iImageHeight-iBottomBorder, iX1Pos, iImageHeight-iBottomBorder+2, iBlack2);
       }
     }
   }

   // Draw axis
   gdImageLine(gdImageLight, 0, iImageHeight-iBottomBorder, iImageWidth-1, iImageHeight-iBottomBorder, iBlack);
   gdImageLine(gdImageHeavy, 0, iImageHeight-iBottomBorder, iImageWidth-1, iImageHeight-iBottomBorder, iBlack2);

   // Draw box around image
   gdImageRectangle(gdImageLight, 0, 0, iImageWidth-1, iImageHeight-1, iBlack);
   gdImageRectangle(gdImageHeavy, 0, 0, iImageWidth-1, iImageHeight-1, iBlack2);

   if (iCidScan > 0) {
     gdImageFilledEllipse(gdImageLight, (int)( (iCidScan-iPlotStartScan)*iImageWidth/(iPlotEndScan-iPlotStartScan)), iImageHeight-iBottomBorder, 8, 8, iRed);
     gdImageFilledEllipse(gdImageHeavy, (int)( (iCidScan-iPlotStartScan)*iImageWidth/(iPlotEndScan-iPlotStartScan)), iImageHeight-iBottomBorder, 8, 8, iRed);
   }

   // Create the image file ... image_buffer needs to point to data

   //DDS:
   strcpy(szImageDir, szXMLFile);
   replace_path_with_webserver_tmp(szImageDir,sizeof(szImageDir)); // write these in designated tmp area
   char* tmpStr = findRightmostPathSeperator(szImageDir); 
   if (tmpStr++)
     *tmpStr = '\0';
   strcpy(szImageLink, szImageDir);
   translate_absolute_filesystem_path_to_relative_webserver_root_path(szImageLink);

   // Create the image file ... make unique name for each image file
   srandom(pQuan->iLightFirstScan+pQuan->iHeavyFirstScan+ (int)(pQuan->dLightPeptideMass*5.0));
   tStartTime=time((time_t *)NULL);

   rand1 = random(); rand2 = random();

   strcpy(szImageFile, szImageDir);
   strcpy(szImageFile2, szImageDir);
   sprintf(szImageFile+strlen(szImageFile), "x%ld-%d.png", (long)tStartTime, rand1);
   sprintf(szImageFile2+strlen(szImageFile2), "x%ld-%db.png", (long)tStartTime, rand2);

   gdImageInterlace(gdImageLight, 1);
   if ( (fp=fopen(szImageFile, "wb"))!=NULL) {
     // seem to need the memory based version for MSVC to avoid a crash
     int sz;
     char *img = (char *)gdImagePngPtr(gdImageLight, &sz); 
     if (img) {
       size_t fwrote = fwrite(img, 1, sz, fp);
       free(img);
     }
     else {
       printf(" Error - cannot create file %s<BR>\n", szImageFile);
     }
     fclose(fp);
     gdImageDestroy(gdImageLight);  
   }
   else
     printf(" Error - cannot create file %s<BR>\n", szImageFile);

   strcpy(szImageFile, szImageLink);
   sprintf(szImageFile+strlen(szImageFile), "x%ld-%d.png", (long)tStartTime, rand1);

   if (bLabelFree) {
     printf("<CENTER><img src='%s'></CENTER>",
            makeTmpPNGFileSrcRef(szImageFile).c_str());
   }
   else {
     gdImageInterlace(gdImageHeavy, 1);
     if ( (fp=fopen(szImageFile2, "wb"))!=NULL) {
       // seem to need the memory based version for MSVC to avoid a crash
       int sz;
       char *img = (char *)gdImagePngPtr(gdImageHeavy, &sz); 
       if (img) {
	 size_t fwrote = fwrite(img, 1, sz, fp);
	 free(img);
       }
       else {
	 printf(" Error - cannot create file %s<BR>\n", szImageFile2);
       }
       fclose(fp);
       gdImageDestroy(gdImageHeavy);
     }
     else
       printf(" Error - cannot create file %s<BR>\n", szImageFile2);

     strcpy(szImageFile2, szImageLink);
     sprintf(szImageFile2+strlen(szImageFile2), "x%ld-%db.png", (long)tStartTime, rand2);

     printf("<CENTER><img src='%s'><BR><img src='%s'></CENTER>", 
	    makeTmpPNGFileSrcRef(szImageFile).c_str(), 
	    makeTmpPNGFileSrcRef(szImageFile2).c_str());
   }

} /*MAKE_PLOT*/


#define FILTER_SIZE 4
// Use my standard filtering routine
void FILTER_MS(double *dOrigMS,
        double *dFilteredMS,
        double *dTmpFilter,
        double dMaxInten,
        int    iNumMSScans,
        int    iPlotStartScan,
        int    iPlotEndScan)
{
   int  i,
        iArraySize=iNumMSScans*sizeof(double);
   double dTmpMax;
 
/*
   Defines 5th order butterworth filter w/ cut-off frequency
   of 0.25 where 1.0 corresponse to half the sample rate.

   5th order, 0.010
   double a[FILTER_SIZE]={1.0000, -4.8983, 9.5985, -9.4053, 4.6085, -0.9033},
          b[FILTER_SIZE]={0.000000000909, 0.000000004546, 0.000000009093, 0.000000009093, 0.000000004546, 0.000000000909};

   5th order, 0.075
   double a[FILTER_SIZE]={1.0000, -4.2380, 7.2344, -6.2125, 2.5821, -0.4655},
          b[FILTER_SIZE]={0.0158, 0.0792, 0.1585, 0.1585, 0.0792, 0.0158};

   5th order, 0.10
   double a[FILTER_SIZE]={1.0000, -3.9845, 6.4349, -5.2536, 2.1651, -0.3599},
          b[FILTER_SIZE]={0.0000598, 0.0002990, 0.0005980, 0.0005980, 0.0002990, 0.0000598};

   5th order, 0.15
   double a[FILTER_SIZE]={1.0000, -3.4789, 5.0098, -3.6995, 1.3942, -0.2138},
          b[FILTER_SIZE]={0.0004, 0.0018, 0.0037, 0.0037, 0.0018, 0.0004};

   5th order, 0.20
   double a[FILTER_SIZE]={1.0000, -2.9754, 3.8060, -2.5453, 0.8811, -0.1254},
          b[FILTER_SIZE]={0.0013, 0.0064, 0.0128, 0.0128, 0.0064, 0.0013};

   5th order, 0.25
   double a[FILTER_SIZE]={1.0, -2.4744, 2.8110, -1.7038, 0.5444, -0.0723},
          b[FILTER_SIZE]={0.0033, 0.0164, 0.0328, 0.0328, 0.0164, 0.0033};
*/

   // 3rd order butterworth, 0.025 cutoff frequency
   double a[FILTER_SIZE]={1.0,-2.8430, 2.6980, -0.8546},
          b[FILTER_SIZE]={0.0000561, 0.0001682, 0.0001682, 0.0000561};
 
   memset(dFilteredMS, 0, iArraySize);
   memcpy(dTmpFilter, dOrigMS, iNumMSScans*sizeof(double));
 
   // Pass MS profile through IIR low pass filter:
   // y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
   //      - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
   for (i=0; i<iNumMSScans; i++) {
     int ii;

     dFilteredMS[i]=b[0]*dTmpFilter[i];
     for (ii=1;ii<FILTER_SIZE;ii++) {
       if ((i-ii)>=0) {
	 dFilteredMS[i] += b[ii]*dTmpFilter[i-ii];
	 dFilteredMS[i] -= a[ii]*dFilteredMS[i-ii];
       }
     }
   }
 
   // Filtered sequence is reversed and re-filtered resulting
   // in zero-phase distortion and double the filter order.
   for (i=0; i<iNumMSScans; i++)
     dTmpFilter[i]=dFilteredMS[iNumMSScans-1-i];
 
   memset(dFilteredMS, 0, iArraySize);
   for (i=0; i<iNumMSScans; i++) {
     int ii;
 
     dFilteredMS[i]=b[0]*dTmpFilter[i];
     for (ii=1;ii<FILTER_SIZE;ii++) {
       if ((i-ii)>=0) {
	 dFilteredMS[i] += b[ii]*dTmpFilter[i-ii];
	 dFilteredMS[i] -= a[ii]*dFilteredMS[i-ii];
       }
     }
   }
 
   // Filtered sequence is reversed again
   dTmpMax=0.001;
   for (i=0; i<iNumMSScans; i++) {
     dTmpFilter[i]=dFilteredMS[iNumMSScans-1-i];

     if (i>=iPlotStartScan && i<=iPlotEndScan)
       if (dTmpFilter[i]>dTmpMax)
	 dTmpMax=dTmpFilter[i];
   }

   if (dMaxInten>0.0) {
     for (i=iPlotStartScan; i<=iPlotEndScan; i++)
       dTmpFilter[i] = dTmpFilter[i] * dMaxInten / dTmpMax;
   }

   memcpy(dFilteredMS, dTmpFilter, iArraySize);

} /*FILTER_MS*/


void SPECIAL_QUAN() {
   printf("<TABLE><TR><TD>");
   printf("<FORM METHOD=GET ACTION='%s'>", getenv("SCRIPT_NAME"));
   printf("<INPUT TYPE='hidden' NAME='InteractDir' VALUE='%s'>\n", szInteractDir);
   printf("<INPUT TYPE='hidden' NAME='InteractBaseName' VALUE='%s'>\n", szInteractBaseName);
   printf("<INPUT TYPE='hidden' NAME='bXpressLight1' VALUE='%d'>\n", bXpressLight1);

   printf("<INPUT TYPE='hidden' NAME='LightFirstScan' VALUE='%d'>\n", pQuan.iLightFirstScan);
   printf("<INPUT TYPE='hidden' NAME='LightLastScan'  VALUE='%d'>\n", pQuan.iLightLastScan);
   printf("<INPUT TYPE='hidden' NAME='LightMass' VALUE='%0.3f'>\n", pQuan.dLightPeptideMass);
   printf("<INPUT TYPE='hidden' NAME='MassTol' VALUE='%0.2f'>", pQuan.dMassTol);
   printf("<INPUT TYPE='hidden' NAME='PpmTol' VALUE='%d'>", pQuan.bPpmMassTol);
   printf("<INPUT TYPE='hidden' NAME='NumIsotopePeaks' VALUE='%d'>", pQuan.iNumIsotopePeaks);
   printf("<INPUT TYPE='hidden' NAME='Norm' VALUE='%d'>", pQuan.bNormalize);
   printf("<INPUT TYPE='hidden' NAME='HeavyFirstScan' VALUE='%d'>\n", pQuan.iHeavyFirstScan);
   printf("<INPUT TYPE='hidden' NAME='HeavyLastScan'  VALUE='%d'>\n", 0);
   printf("<INPUT TYPE='hidden' NAME='HeavyMass' VALUE='%0.3f'>\n", pQuan.dHeavyPeptideMass);
   printf("<input TYPE='hidden' NAME='XMLFile' VALUE='%s'>", szXMLFile);
   printf("<input TYPE='hidden' NAME='OutFile' VALUE='%s'>", szOutputFile);
   printf("<input TYPE='hidden' NAME='xmlfile' VALUE='%s'>", update_xmlfile);
   printf("<input TYPE='hidden' NAME='index' VALUE='%d'>", update_index);
   printf("<input TYPE='hidden' NAME='overwrite' VALUE='1'>");
   printf("<INPUT TYPE='hidden' NAME='ChargeState' VALUE='%d'>", pQuan.iChargeState);
   printf("<INPUT TYPE='hidden' NAME='NewQuan' VALUE='%s'>\n", "1:0.00*");
   printf("<INPUT TYPE='hidden' NAME='ZeroArea' VALUE='H'>\n");
   printf("<INPUT TYPE='submit' VALUE='1:0.00'>");

   printf("</FORM>");
   printf("</TD><TD>");
   printf("<FORM METHOD=GET ACTION='%s'>", getenv("SCRIPT_NAME"));

   printf("<INPUT TYPE='hidden' NAME='InteractDir' VALUE='%s'>\n", szInteractDir);
   printf("<INPUT TYPE='hidden' NAME='InteractBaseName' VALUE='%s'>\n", szInteractBaseName);
   printf("<INPUT TYPE='hidden' NAME='bXpressLight1' VALUE='%d'>\n", bXpressLight1);

   printf("<INPUT TYPE='hidden' NAME='LightFirstScan' VALUE='%d'>\n", pQuan.iLightFirstScan);
   printf("<INPUT TYPE='hidden' NAME='LightLastScan'  VALUE='%d'>\n", 0);
   printf("<INPUT TYPE='hidden' NAME='LightMass' VALUE='%0.3f'>\n", pQuan.dLightPeptideMass);
   printf("<INPUT TYPE='hidden' NAME='MassTol' VALUE='%0.2f'>", pQuan.dMassTol);
   printf("<INPUT TYPE='hidden' NAME='PpmTol' VALUE='%d'>", pQuan.bPpmMassTol);
   printf("<INPUT TYPE='hidden' NAME='NumIsotopePeaks' VALUE='%d'>", pQuan.iNumIsotopePeaks);
   printf("<INPUT TYPE='hidden' NAME='Norm' VALUE='%d'>", pQuan.bNormalize);
   printf("<INPUT TYPE='hidden' NAME='HeavyFirstScan' VALUE='%d'>\n", pQuan.iHeavyFirstScan);
   printf("<INPUT TYPE='hidden' NAME='HeavyLastScan'  VALUE='%d'>\n", pQuan.iHeavyLastScan);
   printf("<INPUT TYPE='hidden' NAME='HeavyMass' VALUE='%0.3f'>\n", pQuan.dHeavyPeptideMass);

   printf("<input TYPE='hidden' NAME='XMLFile' VALUE='%s'>", szXMLFile);
   printf("<input TYPE='hidden' NAME='OutFile' VALUE='%s'>", szOutputFile);
   printf("<input TYPE='hidden' NAME='xmlfile' VALUE='%s'>", update_xmlfile);
   printf("<input TYPE='hidden' NAME='index' VALUE='%d'>", update_index);
   printf("<input TYPE='hidden' NAME='overwrite' VALUE='1'>");
   printf("<INPUT TYPE='hidden' NAME='ChargeState' VALUE='%d'>", pQuan.iChargeState);
   printf("<INPUT TYPE='hidden' NAME='NewQuan' VALUE='%s'>\n", "0.00:1*");
   printf("<INPUT TYPE='hidden' NAME='ZeroArea' VALUE='L'>\n");
   printf("<INPUT TYPE='submit' VALUE='0.00:1'>");

   printf("</FORM>");

   printf("</TD><TD>");

   printf("<FORM METHOD=GET ACTION='%s'>", getenv("SCRIPT_NAME"));

   printf("<INPUT TYPE='hidden' NAME='InteractDir' VALUE='%s'>\n", szInteractDir);
   printf("<INPUT TYPE='hidden' NAME='InteractBaseName' VALUE='%s'>\n", szInteractBaseName);
   printf("<INPUT TYPE='hidden' NAME='bXpressLight1' VALUE='%d'>\n", bXpressLight1);

   printf("<INPUT TYPE='hidden' NAME='LightFirstScan' VALUE='%d'>\n", pQuan.iLightFirstScan);
   printf("<INPUT TYPE='hidden' NAME='LightLastScan'  VALUE='%d'>\n", pQuan.iLightFirstScan - 1);
   printf("<INPUT TYPE='hidden' NAME='LightMass' VALUE='%0.3f'>\n", pQuan.dLightPeptideMass);
   printf("<INPUT TYPE='hidden' NAME='MassTol' VALUE='%0.2f'>", pQuan.dMassTol);
   printf("<INPUT TYPE='hidden' NAME='PpmTol' VALUE='%d'>", pQuan.bPpmMassTol);
   printf("<INPUT TYPE='hidden' NAME='NumIsotopePeaks' VALUE='%d'>", pQuan.iNumIsotopePeaks);
   printf("<INPUT TYPE='hidden' NAME='Norm' VALUE='%d'>", pQuan.bNormalize);
   printf("<INPUT TYPE='hidden' NAME='HeavyFirstScan' VALUE='%d'>\n", pQuan.iHeavyFirstScan);
   printf("<INPUT TYPE='hidden' NAME='HeavyLastScan'  VALUE='%d'>\n", pQuan.iHeavyFirstScan - 1); //pQuan.iHeavyLastScan);
   printf("<INPUT TYPE='hidden' NAME='HeavyMass' VALUE='%0.3f'>\n", pQuan.dHeavyPeptideMass);

   printf("<input TYPE='hidden' NAME='XMLFile' VALUE='%s'>", szXMLFile);
   printf("<input TYPE='hidden' NAME='OutFile' VALUE='%s'>", szOutputFile);
   printf("<input TYPE='hidden' NAME='xmlfile' VALUE='%s'>", update_xmlfile);
   printf("<input TYPE='hidden' NAME='index' VALUE='%d'>", update_index);
   printf("<input TYPE='hidden' NAME='overwrite' VALUE='1'>");
   printf("<INPUT TYPE='hidden' NAME='ChargeState' VALUE='%d'>", pQuan.iChargeState);
   printf("<INPUT TYPE='hidden' NAME='NewQuan' VALUE='%s*'>\n", pQuan.szNewQuan);
   printf("<INPUT TYPE='hidden' NAME='NewQuan' VALUE='%s'>\n", "?*");
   printf("<INPUT TYPE='hidden' NAME='ZeroArea' VALUE='LH'>\n");

   printf("<INPUT TYPE='submit' VALUE='?'>");

   printf("</FORM>");
   printf("</TD></TR></TABLE>");
}


void UPDATE_QUAN() {
   printf("<FORM METHOD=GET ACTION='%s'>", getenv("SCRIPT_NAME"));
   printf("<input type='hidden' NAME='LabelFree' VALUE='%d'>\n", bLabelFree);
   printf("<INPUT TYPE='hidden' NAME='InteractDir' VALUE='%s'>\n", szInteractDir);
   printf("<INPUT TYPE='hidden' NAME='InteractBaseName' VALUE='%s'>\n", szInteractBaseName);
   printf("<INPUT TYPE='hidden' NAME='bXpressLight1' VALUE='%d'>\n", bXpressLight1);

   printf("<INPUT TYPE='hidden' NAME='LightFirstScan' VALUE='%d'>\n", pQuan.iLightFirstScan);
   printf("<INPUT TYPE='hidden' NAME='LightLastScan'  VALUE='%d'>\n", pQuan.iLightLastScan);
   printf("<INPUT TYPE='hidden' NAME='LightMass' VALUE='%0.3f'>\n", pQuan.dLightPeptideMass);
   printf("<INPUT TYPE='hidden' NAME='MassTol' VALUE='%0.2f'>", pQuan.dMassTol);
   printf("<INPUT TYPE='hidden' NAME='PpmTol' VALUE='%d'>", pQuan.bPpmMassTol);
   printf("<INPUT TYPE='hidden' NAME='NumIsotopePeaks' VALUE='%d'>", pQuan.iNumIsotopePeaks);
   printf("<INPUT TYPE='hidden' NAME='HeavyFirstScan' VALUE='%d'>\n", pQuan.iHeavyFirstScan);
   printf("<INPUT TYPE='hidden' NAME='HeavyLastScan'  VALUE='%d'>\n", pQuan.iHeavyLastScan);
   printf("<INPUT TYPE='hidden' NAME='HeavyMass' VALUE='%0.3f'>\n", pQuan.dHeavyPeptideMass);

   printf("<input TYPE='hidden' NAME='XMLFile' VALUE='%s'>", szXMLFile);
   printf("<input TYPE='hidden' NAME='OutFile' VALUE='%s'>", szOutputFile);
   printf("<input TYPE='hidden' NAME='xmlfile' VALUE='%s'>", update_xmlfile);
   printf("<input TYPE='hidden' NAME='index' VALUE='%d'>", update_index);
   printf("<input TYPE='hidden' NAME='overwrite' VALUE='1'>");
   printf("<INPUT TYPE='hidden' NAME='ChargeState' VALUE='%d'>", pQuan.iChargeState);
   printf("<input type='submit' ");
   if (!quantitating)
     printf("disabled='1' title='first click on Quantitate to evaluate changes before Saving' ");
   if (bLabelFree)
     printf("value='Save Updated Area'>");
   else
     printf("value='Save Updated Ratio'>");

   printf("</FORM>");

   if(overwrite) { 
     if(strlen(zero_area) > 0) {
       if(strchr(zero_area, 'L') != NULL)
	 pQuan.dLightQuanValue = 0.0;
       if(strchr(zero_area, 'H') != NULL)
	 pQuan.dHeavyQuanValue = 0.0;
     }
     if(strlen(update_quan) > 0)
       strcpy(pQuan.szNewQuan, update_quan);

     char szBuf[4096];
     Tag* replacement;

     if (bLabelFree) {
       replacement = new Tag("xpresslabelfree_result", True, True);

       sprintf(szBuf, "%d", pQuan.iChargeState);
       replacement->setAttributeValue("charge", szBuf);
       sprintf(szBuf, "%d", pQuan.iLightFirstScan);
       replacement->setAttributeValue("first_scan", szBuf);
       sprintf(szBuf, "%d", pQuan.iLightLastScan);
       replacement->setAttributeValue("last_scan", szBuf);
       sprintf(szBuf, "%0.0lf", pQuan.dLightFirstScanRT);
       replacement->setAttributeValue("first_scan_RT_seconds", szBuf);
       sprintf(szBuf, "%0.0lf", pQuan.dLightLastScanRT);
       replacement->setAttributeValue("last_scan_RT_seconds", szBuf);
       //sprintf(szBuf, "%0.4f", (pQuan.dLightPeptideMass + dProton*(pQuan.iChargeState -1)) /(pQuan.iChargeState));
       sprintf(szBuf, "%0.4f", (pQuan.dLightPeptideMass + 1.00727646688*(pQuan.iChargeState -1)) /(pQuan.iChargeState));
       replacement->setAttributeValue("precursor_mz", szBuf);
       sprintf(szBuf, "%0.2e", pQuan.dLightQuanValue);
       replacement->setAttributeValue("peak_area", szBuf);
       sprintf(szBuf, "%0.2e", pQuan.dLightIntensity);
       replacement->setAttributeValue("peak_intensity", szBuf);
       sprintf(szBuf, "%0.0lf", pQuan.dLightIntensityRT);
       replacement->setAttributeValue("peak_intensity_RT_seconds", szBuf);
       sprintf(szBuf, "%d", pQuan.iLightIntensityScan);
       replacement->setAttributeValue("peak_intensity_scan", szBuf);

     }
     else {
       replacement = new Tag("xpressratio_result", True, True);
  
       sprintf(szBuf, "%d", pQuan.iLightFirstScan);
       replacement->setAttributeValue("light_firstscan", szBuf);
       sprintf(szBuf, "%d", pQuan.iLightLastScan);
       replacement->setAttributeValue("light_lastscan", szBuf);

       sprintf(szBuf, "%0.3f", pQuan.dLightPeptideMass);
       replacement->setAttributeValue("light_mass", szBuf);
       sprintf(szBuf, "%d", pQuan.iHeavyFirstScan);
       replacement->setAttributeValue("heavy_firstscan", szBuf);
       sprintf(szBuf, "%d", pQuan.iHeavyLastScan);
       replacement->setAttributeValue("heavy_lastscan", szBuf);
       sprintf(szBuf, "%0.3f", pQuan.dHeavyPeptideMass);
       replacement->setAttributeValue("heavy_mass", szBuf);
       sprintf(szBuf, "%0.3f", pQuan.dMassTol);
       replacement->setAttributeValue("mass_tol", szBuf);
       sprintf(szBuf, "%d", pQuan.bPpmMassTol);
       replacement->setAttributeValue("ppm_tol", szBuf);
       sprintf(szBuf, "%d", pQuan.iNumIsotopePeaks);
       replacement->setAttributeValue("min_num_isotope_peaks", szBuf);

       sprintf(szBuf, "%s", pQuan.szNewQuan);
       replacement->setAttributeValue("ratio", szBuf);
       // now the flipped guy....
       flipRatio(szBuf, szBuf);
       replacement->setAttributeValue("heavy2light_ratio", szBuf);

       sprintf(szBuf, "%0.2e", pQuan.dLightQuanValue);
       replacement->setAttributeValue("light_area", szBuf);
       sprintf(szBuf, "%0.2e", pQuan.dHeavyQuanValue);
       replacement->setAttributeValue("heavy_area", szBuf);

       if(pQuan.dHeavyQuanValue == 0.0) {
	 if (pQuan.dLightQuanValue > 0.0)
	   sprintf(szBuf, "%0.1f", 999.0);
	 else 
	   sprintf(szBuf, "%0.1f", -1.0);
       }
       else
         sprintf(szBuf, "%0.2f", pQuan.dLightQuanValue/pQuan.dHeavyQuanValue);

       replacement->setAttributeValue("decimal_ratio", szBuf);
     }

     XPressPeptideUpdateParser* parser = new XPressPeptideUpdateParser(update_xmlfile, update_index, replacement);

     if (parser->update()) {
       //printf("</TD></TR><TR><TD COLSPAN='3' ALIGN=CENTER> ratio updated, refresh browser to view change"); // ok
       printf("<script language='JavaScript'>messages.push(\"Ratio updated in %s\");</script>\n", update_xmlfile);
       printf("<script language='JavaScript'>messages.push(\"Refresh browser to view change\");</script>\n");
     }
     else {
       //printf("</TD></TR><TR><TD COLSPAN='3' ALIGN=CENTER><font color='red'><b>Error: ratio not updated</b></font>"); // ok
       printf("<script language='JavaScript'>messages.push(\"Error: ratio not updated\");</script>\n");
     }
   }
} /*UPDATE_QUAN*/


void flipRatio(char* ratio, char* flipped) {
  // find the : and 1
  double left, right;
  if (!strcmp(ratio, "?*")) {
    strcpy(flipped, ratio);
  }
  else {
    sscanf(ratio, "%lf:%lf", &left, &right);
    if (left == 1.0)
      ;
    else if (left == 0.0)
      left = 999;
    else if (left >= 999.)
      left = 0.0;
    else
      left = 1.0 / left;

    if(right == 1.0)
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
}

// extract the scan number from the dta file name, which is expected be of the format
// xxxxxxxxxxxxxx.scanstart.scanend.charge
int GET_SCANNUM_FROM_DTA_NAME(const char* specName) {
  int scanNum = -1;
  char szTmp[500];
  char *pStr;

  strcpy(szTmp, specName);

  int fieldCount=3;
  while (fieldCount != 0) {
    pStr = strrchr(szTmp, '.');
    if (pStr == NULL) {
      printf("<script language='JavaScript'>messages.push(\"Error - cannot get scan number from input file <b>%s</b>; unexpected dta name format\");</script>\n", specName);
      return scanNum;
    }
    *pStr = '\0';
    --fieldCount;
  }
  // szTemp is now xxxxxxx.scanstartNULLxxxxxx
  pStr++;
  sscanf(pStr, "%d", &scanNum);

  return scanNum;
}
