/*

Program       : AnalysisSummary                                                  
Author        : Andrew Keller <akeller@systemsbiology.org>                                                     
Date          : 11.27.02 

Primary data object holding all mixture distributions for each precursor ion charge

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


#include "AnalysisSummaryParser.h"

#include "Common/TPPVersion.h" // contains version number, name, revision
#include "Common/util.h"
using namespace std;

#define SIZE_BUF      8192
#define SIZE_PEPTIDE   128
#define SIZE_FILE     1024
#define LF 10
#define CR 13

int getline(char *s, int n, FILE *f) {
    register int i=0;

    while(1) {
        s[i] = (char)fgetc(f);

        if(s[i] == CR)
            s[i] = fgetc(f);

        if((s[i] == 0x4) || (s[i] == LF) || (i == (n-1))) {
            s[i] = '\0';
            return (feof(f) ? 1 : 0);
        }
        ++i;
    }
}

void EXTRACT_QUERY_STRING(char* xmlfile, char* analysis, char* timestamp, int* timestamp_only)
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

      if (!strcmp(szWord, "xmlfile"))
      {
         getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
         strcpy(xmlfile, szWord);
      }
      else if (!strcmp(szWord, "analysis"))
      {
         getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
         strcpy(analysis, szWord);
      }
      else if (!strcmp(szWord, "timestamp"))
      {
         getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
         strcpy(timestamp, szWord);
      }
      else if (!strcmp(szWord, "timestamp_only"))
      {
         getword(szWord, pQS, '&'); plustospace(szWord); unescape_url(szWord);
	 *timestamp_only = atoi(szWord);
         //strcpy(timestamp, szWord);
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

  char xmlfile[2000];
  char analysis[1000];
  char timestamp[1000];
  xmlfile[0] = 0;
  analysis[0] = 0;
  timestamp[0] = 0;

  int timestamp_only = 0; // false

   EXTRACT_QUERY_STRING(xmlfile, analysis, timestamp, &timestamp_only);
   cout << "Content-type: text/html" << endl << endl;

   AnalysisSummaryParser* parser = new AnalysisSummaryParser(xmlfile, analysis, timestamp, timestamp_only);

   cout << "<HTML>" << endl;
   cout << "<!-- " << szTPPVersionInfo << "-->" << endl;
   cout <<"   <HEAD>" << endl;
   cout <<"<TITLE>Analysis Summary (" << szTPPVersionInfo << ")</TITLE>" << endl;
   cout <<"   </HEAD>" << endl;
   cout << "<table width=\"100%\"><tr><td align=\"left\">";

   cout << "Analysis Summary: <b><font color=\"red\">" << analysis << "</font></b>";
   cout << "</td><td align=\"right\">";

   //cout << "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;";
   //cout << "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;";
   cout << "Date: <b>" << timestamp << "</b>";
   cout << "</table>";

   cout << "<p/><HR/><p/>";

   cout << "<center>";
   if(parser != NULL)
     parser->print();

   cout << "</center>";
   cout << "</HTML>";

   delete parser;
   return 0;
}
