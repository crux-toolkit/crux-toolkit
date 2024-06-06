/****************************************************************************
 *
 * Program       : dta-xml.c
 * Author        : Luis Mendoza <l mendoza at systems biology dot org>
 * Date          : 12.06.06
 * Version       : $Id: dta-xml.cpp 8388 2021-02-09 23:10:04Z dshteyn $
 * Purpose       : Extract dta from existing file, tgz, or mzXML
 *                 Returns xml
 * Inputs        : Called via the GET cgi method
 *       * Dta   : full path to location of dta file (does not have to exist)
 *         File  : same as Dta
 *         COMET : if present, assumes files were searched by Comet
 *
 *   Most of the code comes from COMETPLOT by Jimmy Eng (plot-msms.c)
 *
 * Copyright (C) 2006 Luis Mendoza
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * Luis Mendoza
 * Institute for Systems Biology
 * 401 Terry Avenue North
 * Seattle, WA  98109  USA
 * l mendoza at systems biology dot org
 *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "Common/constants.h"
#include "Common/TPPVersion.h" // contains version number, name, revision
#include "Parsers/mzParser/mzParser.h"

using namespace mzParser;

struct EnvironmentStruct {
  int  bCOMET;    /* comet results */   
  int  bDebugging;    /* unused at the moment */

  char szInputFile[SIZE_FILE];
  char szFullPathInputFile[SIZE_FILE];
  char szTarFile[SIZE_FILE];
  char szMZXMLFile[SIZE_FILE];

} pEnvironment;

void EXTRACT_CGI_QUERY(void);
void INITIALIZE(void);
void REPORT_ERROR(const char *szMessage, int iExitCode);
void PRINT_DATA(FILE *ppIn,
		RAMPFILE *fp_,
		int  iFileType,
		int iCharge);

#include "Common/util.h"

int main(int argc, char **argv) {
  FILE *ppIn;

  RAMPFILE *fp_ = NULL;

  int  i,
    iFileType,  /* 0 = .dta file, 1 = read from .tgz file, 2 = read from mzXML file */
    iCharge=0,
    iLen;

  struct stat statbuf;

  char szBaseName[SIZE_FILE],
    szCommand[SIZE_BUF],
    szWebserverRoot[SIZE_BUF],
    szBuf[SIZE_BUF],
    *pStr;
  char* result = NULL;

  hooks_tpp(argc,argv);  // installdir issues, etc

  /*
   * Start the xml document
   */
  printf("Content-type: text/xml\n\n");
  printf("<XML>\n");
  printf("<TPP Version=\"%s\" />\n", szTPPVersionInfo);

  INITIALIZE();

  pStr=getenv("WEBSERVER_ROOT");
  if (pStr==NULL) {
    sprintf(szBuf, "Environment variable WEBSERVER_ROOT does not exist.\n\n");
#ifdef __MINGW__
    strcat(szBuf, " For Windows users, you can set this environment variable\n");
    strcat(szBuf, " through the Advanced tab under System Properties when you\n");
    strcat(szBuf, " right-mouse-click on your My Computer icon.\n\n");
    strcat(szBuf, " Set this environment variable to your webserver's document\n");
    strcat(szBuf, " root directory such as c:\\inetpub\\wwwroot for IIS or\n");
    strcat(szBuf, " c:\\website\\htdocs or WebSite Pro.\n\n");
#endif
    strcat(szBuf, " Exiting.\n");
    REPORT_ERROR(szBuf, 0);
  }
  else {
    strcpy(szWebserverRoot, pStr);
  }

  /*
   * Check if szWebserverRoot is present
   */
  if (access(szWebserverRoot, F_OK)) {
    sprintf(szBuf," Cannot access the webserver's root directory:\n    %s\n", szWebserverRoot);
    strcat(szBuf, " This was set as the environment variable WEBSERVER_ROOT\n\n");
    strcat(szBuf, " For Windows users, you can check this environment variable\n");
    strcat(szBuf, " through the Advanced tab under System Properties when you\n");
    strcat(szBuf, " right-mouse-click on your My Computer icon.\n\n");
    strcat(szBuf, " Exiting.\n");
    REPORT_ERROR(szBuf, 1);
  }

  EXTRACT_CGI_QUERY();

  sscanf(pEnvironment.szInputFile+strlen(pEnvironment.szInputFile)-5, "%d", &iCharge);
  /* small safeguard in case above parsing of charge state from .dta name fails */
  if (iCharge<0 || iCharge>9)
    iCharge=1;


  /*
   * pEnvironment.szInputFile has input .dta file ... need to get this
   * from the tar/gzipped archive now
   *
   * So, first modify szTarFile to be the gzipped tar file name
   */
  strcpy(pEnvironment.szTarFile, pEnvironment.szInputFile);
  iLen=strlen(pEnvironment.szTarFile);
  for (i=iLen-1; i>2; i--) {
    if (pEnvironment.szTarFile[i]=='/' && pEnvironment.szTarFile[i-1]=='.' && pEnvironment.szTarFile[i-2]=='/') {
      pEnvironment.szTarFile[i-2]='\0';

      if (pEnvironment.bCOMET)
	strcat(pEnvironment.szTarFile, ".cmt.tar.gz");
      else
	strcat(pEnvironment.szTarFile, ".tgz");
      break;
    }
    else if (pEnvironment.szTarFile[i]=='/') { /* this case is needed for IE on the Mac ... /./ changes to just / */
      pEnvironment.szTarFile[i]='\0';

      if (pEnvironment.bCOMET)
	strcat(pEnvironment.szTarFile, ".cmt.tar.gz");
      else
	strcat(pEnvironment.szTarFile, ".tgz");
      break;
    }
  }

  iLen=strlen(pEnvironment.szTarFile);
  strcpy(szBaseName, pEnvironment.szTarFile);
  if (pEnvironment.bCOMET)
    szBaseName[strlen(szBaseName)-11]='\0';
  else
    szBaseName[strlen(szBaseName)-4]='\0';

  rampConstructInputFileName(pEnvironment.szMZXMLFile, sizeof(pEnvironment.szMZXMLFile), szBaseName);

  if (!strstr(pEnvironment.szTarFile+iLen-4, ".tgz")
      && !strstr(pEnvironment.szTarFile+iLen-11, ".cmt.tar.gz")) {
    sprintf(szBuf,
	    " Error converting %s to tar file.\nszBaseName=%s\nszTarFile=%s\niLen=%d\n",
	    pEnvironment.szInputFile, szBaseName, pEnvironment.szTarFile, iLen);

    REPORT_ERROR(szBuf, EXIT_FAILURE);
  }

  /*
   * Next, modify pEnvironment.szInputFile to remove full path
   */
  strcpy(pEnvironment.szFullPathInputFile, pEnvironment.szInputFile);
  strcpy(pEnvironment.szInputFile, pEnvironment.szInputFile+i+1);
  sprintf(szCommand, "*%s", pEnvironment.szInputFile);
  strcpy(pEnvironment.szInputFile, szCommand);

  /*
   * JENG: if szInput file starts with pps_ then assume this is a toftof run
   * from internal ISB runtoftof script using the mascot2dta program to create
   * dta files.  There will be an extra label after basename but before
   * scan #s in the .dta file names that needs to be removed from the
   * tar file name.  This next if statement should be irrelevant outside of ISB.
   */
  if (!strncmp(pEnvironment.szInputFile, "*pps_", 5)) {
    /* remove the .tgz */
    pEnvironment.szTarFile[strlen(pEnvironment.szTarFile)-4]='\0';
      
    /* remove the label */
    if ((pStr = strrchr(pEnvironment.szTarFile, '.')))
      *pStr='\0';

    strcat(pEnvironment.szTarFile, ".tgz");
  }

  /*
   * handle searches which weren't done by runsearch, and so don't have tgz'd .out and .dta
   */
  if (stat(pEnvironment.szTarFile,&statbuf) && /* does the tgz file exist? */
      stat(pEnvironment.szFullPathInputFile,&statbuf)) /* does the .dta file exist as named in the command? */
    { /* no, try to remove that middle foo in wwwroot/foo/foo/foo.0001.0001.2.dta */
      char szUntarredDTAfile[SIZE_BUF];
      strncpy(szUntarredDTAfile,pEnvironment.szFullPathInputFile,sizeof(szUntarredDTAfile));
      char *slash = strrchr(szUntarredDTAfile,'/');
      if (slash) {
	char *slash2;
	*slash = 0;
	slash2 = strrchr(szUntarredDTAfile,'/');
	if (slash2)
	  strcpy(slash2+1,slash+1);
      }
      if (!stat(szUntarredDTAfile,&statbuf)) { 
	/* use this as the filename */
	strncpy(pEnvironment.szFullPathInputFile,szUntarredDTAfile,sizeof(pEnvironment.szFullPathInputFile));
      }
    }


  /*
   * (jmt) The goal is to end up with a filepointer that will be
   * reading in dta-format data.  This filepointer may be a direct
   * fopen from an existing dta filename, a pipe-open from a
   * uncompressing a tgz file, or the result of the ramp parser
   * operating on an mzXML file.
   *
   * iFileType == 0: direct from dta file
   * iFileType == 1: from tgz archive
   * iFileType == 2: from mzXML via ramp
   *
   */
  iFileType = 0;

  /* try to open the dta file directly, from the filesystem */
  if ((ppIn=fopen(pEnvironment.szFullPathInputFile, "r"))==NULL) {
    /* if can't open .dta directly, try from .tgz file */
    iFileType = 1;

    /* test if the archive contains the specified dta */
    sprintf(szCommand, "tar --wildcards -tzf %s \"%s\" > /dev/null", pEnvironment.szTarFile, pEnvironment.szInputFile);

    int archiveResult = tpplib_system(szCommand); // like system(), but handles possible tar issues for win32
    if (archiveResult == 0) {
      /* found in archive; open the tar extraction as a pipe */
      sprintf(szCommand, "tar --wildcards  -xzOf %s \"%s\"", 
	      pEnvironment.szTarFile, pEnvironment.szInputFile);
      ppIn=tpplib_popen(szCommand, "r"); // handles possible tar issues for win32
    }

    else { /* not found in archive */
      /* try reading spectrum directly from mzXML file */
      fp_ = NULL;

      if ((fp_ = rampOpenFile(pEnvironment.szMZXMLFile)) == NULL) {
	sprintf(szBuf, " Error - can't read mzXML file %s\n", pEnvironment.szMZXMLFile);
	REPORT_ERROR(szBuf, EXIT_FAILURE);
      }
      else {
	iFileType = 2;   /* success ... read from mzXML file */
      }
    }

  }
  else { /* if here, the dta file actually existed in the filesystem (iFileType == 0) */
    strcpy(szBaseName, pEnvironment.szFullPathInputFile);
    szBaseName[strlen(szBaseName)-4]='\0';   /* remove .dta extension */
  }

  PRINT_DATA(ppIn, fp_, iFileType, iCharge);

  printf("<SOURCE>\n<FILE>");
  if (iFileType == 0) {
    printf("%s",pEnvironment.szFullPathInputFile);
    fclose(ppIn);
  }
  else if (iFileType == 1) {
    printf("%s",pEnvironment.szTarFile);
    pclose(ppIn);
  }
  else if (iFileType == 2) {
    printf("%s",pEnvironment.szMZXMLFile);
    if (fp_ != NULL)
      rampCloseFile(fp_);
  }
  printf("</FILE>\n");
  printf("<TYPE>%d</TYPE>\n", iFileType);
  printf("<SPECTRUM>%s</SPECTRUM>\n", pEnvironment.szInputFile+1); // no asterisk
  printf("</SOURCE>\n");

  printf("<STATUS>OK</STATUS>\n");
  printf("</XML>\n");

  fflush(stdout);

  return(EXIT_SUCCESS);

}


void EXTRACT_CGI_QUERY(void) {
  char *pRequestType,
    *pQS,
    szWord[1024];
  int  i;

  pRequestType=getenv("REQUEST_METHOD");
  if(pRequestType==NULL)
    REPORT_ERROR("This program needs to be called with CGI GET method.\n", EXIT_FAILURE);
  else if (strcmp(pRequestType, "GET"))
    REPORT_ERROR("This program needs to be called with CGI GET method.\n", EXIT_FAILURE);

  // Decode GET method
  pQS = getenv("QUERY_STRING");
  if (pQS == NULL)
    REPORT_ERROR("GET query string empty.\n", EXIT_FAILURE);

  for (i=0; pQS[0]!='\0';i++) {
    getword(szWord, pQS, '=');
    plustospace(szWord);
    unescape_url(szWord);

    //    printf("Word=%s<BR>\n", szWord);

    if (!strcmp(szWord, "File")) {
      getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
      sscanf(szWord, "%s", pEnvironment.szInputFile);
    }
    else if (!strcmp(szWord, "Dta")) {
      getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
      strcpy(pEnvironment.szInputFile, szWord);
    }
    else if (!strcasecmp(szWord, "COMET")) {
      getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
      sscanf(szWord, "%d", &(pEnvironment.bCOMET));
    }
    else if (!strcmp(szWord, "Debug")) {
      getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
      pEnvironment.bDebugging=1;
    }
    else
      getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
  }

  unCygwinify(pEnvironment.szInputFile); // get rid of cygdrive hoohah (no effect in cygwin builds)

}


void INITIALIZE(void) {
  pEnvironment.bCOMET=0;
  pEnvironment.bDebugging=0;
  pEnvironment.szInputFile[0]='\0';

}


void REPORT_ERROR(const char *szMessage, int iExitCode) {
  printf("<STATUS>ERROR</STATUS>\n");
  printf("<MESSAGE>%s</MESSAGE>\n</XML>", szMessage);
  exit(iExitCode);
}


void PRINT_DATA(FILE *ppIn,
		RAMPFILE *fp_,
		int  iFileType,
		int iCharge) {
  int    iScanNum;
  char   szBuf[SIZE_BUF];
  double dPepMass;

  struct ScanHeaderStruct scanHeader;
  int iAnalysisLastScan;
  ramp_fileoffset_t *index_;


  if (iFileType == 2) {
    // mzXML file: read file index and parse scan # from encoded .dta file
    char szTmp[500];
    char *pStr;

    if ( (pStr = strchr(pEnvironment.szInputFile, '.'))==NULL) {
      sprintf(szBuf, "Error - cannot get scan number from input file %s\n", pEnvironment.szInputFile);
      REPORT_ERROR(szBuf, 1);
    }

    strcpy(szTmp, pStr+1);
    pStr = strchr(szTmp, '.');
    *pStr = '\0';

    sscanf(szTmp, "%d", &iScanNum);

    // read mzXML file offset scan index
    index_ = readIndex(fp_, getIndexOffset(fp_), &(iAnalysisLastScan));

    scanHeader.highMZ = 0.0;

    readHeader(fp_, index_[iScanNum], &scanHeader);

    // get precursor mass from mzXML header
    dPepMass = scanHeader.precursorMZ;

    // report mass as expected by dta format
    dPepMass = iCharge*(dPepMass - 1.0) + 1.0;

  }
  else { // is not mzXML
    char *fgot=fgets(szBuf, SIZE_BUF, ppIn);  // read first line of .dta file
    sscanf(szBuf, "%lf %d", &dPepMass, &iCharge);   // override previous charge from .dta name
  }

  printf("<PRECURSOR Mass=\"%f\" Charge=\"%d\" />\n", dPepMass, iCharge);

  // parse through mass/intensity pairs
  printf("<DTA>");
  printf("%f %d\n", dPepMass, iCharge);

  if (iFileType != 2) {
    while (fgets(szBuf, SIZE_BUF, ppIn))
      printf("%s", szBuf);
  }
  else {
    scanHeader.msLevel = readMsLevel(fp_, index_[iScanNum]);

    if (scanHeader.msLevel == 2) {
      RAMPREAL *pPeaks;
      int n = 0;

      // Open a scan
      pPeaks = readPeaks(fp_, index_[iScanNum]);

      while (pPeaks != NULL && pPeaks[n] != -1) {
	RAMPREAL fMass;
	RAMPREAL fInten;

	fMass = pPeaks[n];
	n++;
	fInten = pPeaks[n];
	n++;

	printf("%f %f\n", fMass, fInten);

      }
      //if (pPeaks != NULL)
      //	free(pPeaks);
    }
    else {
      printf("</DTA>\n");
      sprintf(szBuf, " Error - scan %d is not an MS/MS scan in the mzXML file %s\n", iScanNum, pEnvironment.szMZXMLFile);
      REPORT_ERROR(szBuf, 1);
    }

    if (index_ != NULL)  /* done with index */
      delete index_;
  }

  printf("</DTA>\n");

} /*PRINT_DATA*/

