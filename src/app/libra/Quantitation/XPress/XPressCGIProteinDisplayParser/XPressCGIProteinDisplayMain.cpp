/*
Program       : XPressProteinDisplay                                                  
Author        : J.Eng and Andrew Keller <akeller@systemsbiology.org>, Robert Hubley, and 
                open source code                                                       
Date          : 11.27.02 
SVN Info      : $Id: XPressCGIProteinDisplayMain.cpp 8136 2020-05-23 05:21:51Z real_procopio $


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


#include "XPressCGIProteinDisplay.h"

#include "Common/TPPVersion.h" // contains version number, name, revision
#include "Common/util.h"


#define SIZE_BUF      8192
#define SIZE_PEPTIDE   128
#define SIZE_FILE     1024
#define LF 10
#define CR 13


static struct OptionsStruct
{
   char   *szPeptideString;
   char   szInteractFiles[SIZE_FILE];   /* this should be multiple files!! */
   char   szXmlFile[SIZE_FILE];
   char   szProtein[SIZE_FILE];
   double dRatio;
   double dStdDev;
   double dMinProbability;
   int    iNumEntries;
   int    bHeavy2Light;
   char szCgiHome[SIZE_FILE]; 
   char xslt[SIZE_FILE];
   char szMarkAAs[100];
   Boolean bGlyc;

} pOptions;


static void EXTRACT_QUERY_STRING(struct OptionsStruct *pOptions)
{
   int iLen=0,
       i;
   char *pRequestType;
   char *pQS;
   char szWord[SIZE_BUF];

   pRequestType=getenv("REQUEST_METHOD");
   if(pRequestType==NULL)
   {
      printf(" This program needs to be called with CGI GET method.\n");
      exit(EXIT_FAILURE);
   }
   else if (strcmp(pRequestType, "GET"))
   {
      printf(" This program not called with GET method!\n");
      exit(EXIT_FAILURE);
   }

   /*
    * Decode GET method
    */
   pQS = getenv("QUERY_STRING");
   if (pQS == NULL)
   {
      printf("GET query string empty.\n");
      exit(EXIT_FAILURE);
   }

   for (i=0; pQS[0]!='\0';i++)
   {
      getword(szWord, pQS, '=');
      plustospace(szWord);
      unescape_url(szWord);

      if (!strcmp(szWord, "ratio"))
      {
         getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
         sscanf(szWord, "%lf", &(pOptions->dRatio));
//       printf("Ratio:  %lf\n", pOptions->dRatio);
      }
      else if (!strcmp(szWord, "protein"))
      {
         getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
         strcpy(pOptions->szProtein, szWord);
      }
      else if (!strcmp(szWord, "stddev"))
      {
         getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
         sscanf(szWord, "%lf", &(pOptions->dStdDev));
//       printf("StdDev:  %lf\n", pOptions->dStdDev);
      }
      else if (!strcmp(szWord, "num"))
      {
         getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
         sscanf(szWord, "%d", &(pOptions->iNumEntries));
//       printf("NumEntries:  %d\n", pOptions->iNumEntries);
      }
      else if (!strcmp(szWord, "min_pep_prob"))
      {
         getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
         sscanf(szWord, "%lf", &(pOptions->dMinProbability));
//       printf("MinProb:  %lf\n", pOptions->dMinProbability);
      }
      else if (!strcmp(szWord, "xmlfile"))
      {
         getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
         strcpy(pOptions->szXmlFile, szWord);
//       printf("XmlFile:  %s\n", pOptions->szXmlFile);
      }
      else if (!strcmp(szWord, "source_files"))
      {
	getword(szWord, pQS, '&'); // plustospace(szWord); unescape_url(szWord);
         strcpy(pOptions->szInteractFiles, szWord);
	 plustospace(pOptions->szInteractFiles); unescape_url(pOptions->szInteractFiles);
//       printf("InteractFile:  %s\n", pOptions->szInteractFiles);
      }
      else if (!strcmp(szWord, "heavy2light"))
      {
         getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
	 sscanf(szWord, "%d", &(pOptions->bHeavy2Light));
      }
      else if (!strcmp(szWord, "cgihome"))
      {
         getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
	 strcpy(pOptions->szCgiHome, szWord);
      }
      //else if (!strcmp(szWord, "xslt"))
      //{
      //  getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
	 //strcpy(pOptions->xslt, szWord);
      //}
      else if (!strcmp(szWord, "mark_aa"))
      {
         getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
	 strcpy(pOptions->szMarkAAs, szWord);
      }
      else if (!strcmp(szWord, "glyc"))
      {
	pOptions->bGlyc = True;
      }
      else if (!strcmp(szWord, "peptide_string"))
      {
	getword(szWord, pQS, '&'); // plustospace(szWord); unescape_url(szWord);

         pOptions->szPeptideString = (char*)malloc((4+strlen(szWord))*sizeof(char));
         if (pOptions->szPeptideString == NULL)
         {
            printf(" Error - can't malloc szPeptideString\n\n");
            exit(0);
         }

         sprintf(pOptions->szPeptideString, "%s", szWord);
	 //plustospace(pOptions->szPeptideString); 
	 unescape_url(pOptions->szPeptideString);

//       printf("PepString:  '%s'\n", pOptions->szPeptideString);
      }

      else
      {
//       printf("variable:  %s\n", szWord);
         getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
//       printf("   value:  %s\n\n", szWord);
      }
   }

} /*EXTRACT_QUERY_STRING*/



int main(int argc, char** argv) {
  hooks_tpp handler(argc,argv); // set up install paths etc

  pOptions.dStdDev=0.0;
  pOptions.dRatio=0.0;
  pOptions.dMinProbability=0.0;
  pOptions.iNumEntries=0;
  pOptions.bHeavy2Light=0;
  pOptions.szMarkAAs[0] = 0;
  pOptions.bGlyc = False;

  char xslt_proc[1000];
  xslt_proc[0] = 0;

#ifdef XSTL_PROC
  strcat(xslt_proc, XSTL_PROC);
#else
  #ifdef __LINUX__
  strcat(xslt_proc, "/usr/bin/xsltproc");
  #else
  strcat(xslt_proc, getBinPath());
  strcat(xslt_proc, "xsltproc.exe");
  #endif
#endif

  strcpy(pOptions.xslt,xslt_proc);	
  EXTRACT_QUERY_STRING(&pOptions);
  cout << "Content-type: text/html" << endl << endl;
  //cout << pOptions.xslt << endl;
  //cout << "<-- " << szTPPVersionInfo << " -->" << endl << endl;

  XPressCGIProteinDisplay *p = new XPressCGIProteinDisplay(pOptions.szInteractFiles, pOptions.szPeptideString, pOptions.szProtein,
							   pOptions.szCgiHome, pOptions.szXmlFile, (double)pOptions.dMinProbability, 
							   pOptions.bHeavy2Light, pOptions.xslt, pOptions.szMarkAAs, pOptions.bGlyc);

  delete p;
  return 0;
}
