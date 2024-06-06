/*
Program       : main for xpress-prophet-update.cgi                                           
Author        : Andrew Keller <akeller@systemsbiology.org>
                Jimmy Eng <jeng@systemsbiology.org> 
Date          : 11.27.02 
SVN Info      : $Id: xpresscgimain.cpp 8136 2020-05-23 05:21:51Z real_procopio $

Overwrites ProteinProphet XML with user modified protein ratio

Copyright (C) 2003 Andrew Keller, Jimmy Eng

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

#include "XPressCGIParser.h"

#include "Common/TPPVersion.h" // contains version number, name, revision
#include "Common/util.h"

#define SIZE_BUF      8192
#define SIZE_PEPTIDE   128
#define SIZE_FILE     1024
#define LF 10
#define CR 13

struct OptionsStruct
{
   char   szXmlFile[SIZE_FILE];
   char   szProtein[SIZE_FILE];
   double h2l_dRatio;
   double h2l_dStdDev;
   double dRatio;
   double dStdDev;
   int    iNumEntries;
   int    iNobody;
} pOptions;


void EXTRACT_QUERY_STRING(struct OptionsStruct *pOptions)
{
   int iContentLen,
       i;
   char *pRequestType;
   char szVal[SIZE_BUF],
        szName[SIZE_BUF];

   pRequestType=getenv("REQUEST_METHOD");
   if(pRequestType==NULL)
   {
      printf(" This program needs to be called with CGI GET method.\n");
      exit(EXIT_FAILURE);
   }
   else if (strcmp(pRequestType, "POST"))
   {
      printf(" This program not called with POST method!\n");
      exit(EXIT_FAILURE);
   }



   iContentLen = atoi(getenv("CONTENT_LENGTH"));
   if (iContentLen == 0)
   {
      fprintf(stderr,"No query information to decode.</body>\n</html>\n");
      exit(1);
   }

   for (i=0; iContentLen && (!feof(stdin)); i++)
   {
      strcpy(szVal, fmakeword(stdin, '&', &iContentLen));
      plustospace(szVal);
      unescape_url(szVal);
      strcpy(szName, makeword(szVal, '='));

      if (!strcmp(szName, "ratio"))
      {
         sscanf(szVal, "%lf", &(pOptions->dRatio));
      }
      else if (!strcmp(szName, "h2l_ratio"))
      {
         sscanf(szVal, "%lf", &(pOptions->h2l_dRatio));
      }
      else if (!strcmp(szName, "protein") )
      {
         strcpy(pOptions->szProtein, szVal);
      }
      else if (!strcmp(szName, "xmlfile") )
      {
         strcpy(pOptions->szXmlFile, szVal);
      }
      else if (!strcmp(szName, "stddev"))
      {
         sscanf(szVal, "%lf", &(pOptions->dStdDev));
      }
      else if (!strcmp(szName, "h2l_stddev"))
      {
         sscanf(szVal, "%lf", &(pOptions->h2l_dStdDev));
      }
      else if (!strcmp(szName, "num"))
      {
         sscanf(szVal, "%d", &(pOptions->iNumEntries));
      }
      else if (!strcmp(szName, "nobody"))
      {
	 sscanf(szVal, "%d", &(pOptions->iNobody));
      }

   }

} /*EXTRACT_QUERY_STRING*/

int main(int argc, char** argv) {
  hooks_tpp handler(argc,argv); // set up install paths etc

  /*
   * Initialize variables
   */
  pOptions.dStdDev=0.0;
  pOptions.dRatio=0.0;
  pOptions.iNumEntries=0;
  pOptions.szProtein[0]='\0';
  pOptions.szXmlFile[0]='\0';
  pOptions.iNobody=0;

  printf("Content-type: text/html\n\n");
  EXTRACT_QUERY_STRING(&pOptions);

  if (pOptions.iNobody == 0) {
    printf("<HTML>\n");
    printf("<HEAD><TITLE>XPRESS-PROPHET-UPDATE (%s)</TITLE></HEAD>\n",szTPPVersionInfo);
    printf("<BODY BGCOLOR=\"#FFFFFF\" OnLoad=\"self.focus()\">\n");
  }
  printf("<PRE>\n<B>XPRESS-PROPHET-UPDATE (%s)</B>\n\n",szTPPVersionInfo);

  printf("<font color=\"blue\">");
  Parser* cgiparser = new XPressCGIParser(pOptions.szXmlFile, pOptions.szProtein, pOptions.dRatio, pOptions.dStdDev, pOptions.h2l_dRatio, pOptions.h2l_dStdDev, pOptions.iNumEntries);
  printf("</font>");

  delete cgiparser;
  return 0;
}
