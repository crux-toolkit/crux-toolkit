/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Library or "Lesser" General Public      *
 *   License (LGPL) as published by the Free Software Foundation;          *
 *   either version 2 of the License, or (at your option) any later        *
 *   version.                                                              *
 *                                                                         *
 ***************************************************************************/

/*---------------------------------------------------------------------------##
##  File:
##      @@(#) out2summary.c
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##      Jimmy Eng jeng@systemsbiology.org
##  Description:
##      Convert Sequest and TurboSequest *.out files into a
##      single HTML-SUMMARY file ready for use with INTERACT.
##
#******************************************************************************
#*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef _MSC_VER
#include <io.h>
#include <direct.h>
#define S_ISDIR(x) (x&_S_IFDIR)
#else
#include <unistd.h>
#include <dirent.h>
#endif
#include <ctype.h>
#include <string.h>

#define COMETDB_CGI   "comet-fastadb.cgi"
#if defined(__CYGWIN__) || defined(_MSC_VER)
#define CGI_DIR   "isb-bin"
#define PLOT_CGI  "sequest-plot.cgi"
#define OUT_CGI   "sequest-out.cgi"
#else
#define CGI_DIR   "cgi-bin"
#define PLOT_CGI  "plot-msms-js.cgi"
#define OUT_CGI   "sequest-tgz-out.cgi"
#endif


// Lame lame lame.  Should be dynamically
// allocated!
#define SIZE_BUF             8192
#define SIZE_FILE            512
#define SIZE_PEP             128
#define INIT_INTERACT_LINES  1000
#define PERCENTAGE           0.75

//#define __SORCERER__

// Lame I know
// For Unix
//#define SEP                  '/'
// For Windows
#define SEP                  '\\'

// General struct to hold .out contents which
// are needed to create an INTERACT summary line.
struct SequestOutStruct
{
   char cAA1;
   char cAA2;
   char szFileName[SIZE_FILE];
   char szBaseFileName[SIZE_FILE];
   char szProt[SIZE_PEP];
   char szPlainPep[SIZE_PEP];
   char szSubPep[SIZE_PEP];
   char szDSite[SIZE_PEP];
   char szMod[SIZE_FILE];
   char szDup[SIZE_PEP];
   char szDatabase[SIZE_FILE];
   double dAMass;
   double dMass;
   double dXC;
   double dDeltCn;
   double dSp;
   double dMass1;
   double dMass2;
   double dMass3;
   double dMass4;
   double dMass5;
   double dMass6;
   int  iRankSp;
   int  iMassType;
   int  iIon;
   int  iTot;
   int  bSpecialDeltCn;
   int  bNucDb;
} *sequestData;


// General struct to hold .out contents which
// are needed to create an INTERACT header.
struct HeaderStruct
{
   char szDate[SIZE_PEP];
   char szTime[SIZE_PEP];
   char szTimeSuffix[SIZE_PEP];
   char szMassType[SIZE_PEP];
} headerData;

int readOutFile(char *szFileName,
                struct SequestOutStruct *data,
                struct HeaderStruct *hdr);  // return 0 on success, 1 on error (empty file)
void printSummary(struct SequestOutStruct *data,
                 struct HeaderStruct *hdr,
                 int lines,
                 char *szCWD,
                 char *szPeptideLink);
#include "Common/hooks_tpp.h"

/** Program entry point
   */
int main(int argc,
         char **argv)
{
   hooks_tpp(argc,argv); // handle installdir issues, etc
   int  i = 0;
   int  iOutFileCount = 0;
   struct stat fileStat;
#ifdef _MSC_VER
   struct _finddata_t c_file;
   long hFile;
   char szFileMask[SIZE_BUF];
#else
   DIR *dirStream;
   struct dirent *dirEntry;
#endif
   char szPathName[SIZE_BUF];
   char szCWD[SIZE_BUF];
   char szFullQName[SIZE_BUF];
   char szPeptideLink[SIZE_BUF];
   char szBuf[SIZE_BUF]; 

   char *pLink;
   int  iEntriesAllocated = INIT_INTERACT_LINES;
   if (argc < 2)
   {
      printf("Usage: out2summary [filename.out|directory]\n");
      exit(1);
   }

   if ((sequestData = (SequestOutStruct *)malloc(iEntriesAllocated * sizeof(struct SequestOutStruct))) == NULL)
   {
      fprintf(stderr, " Error - cannot allocate initial memory.\n");
      exit(1);
   }

   if (getcwd(szCWD, SIZE_BUF)==NULL)
   {
      fprintf(stderr, " Error - cannot get cwd.\n");
      exit(1);
   }

   /*
    * Get NCBI Blast link
    */
   strcpy(szPeptideLink, "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&amp;LAYOUT=TwoWindows&amp;AUTO_FORMAT=Semiauto&amp;ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;CDD_SEARCH=on&amp;CLIENT=web&amp;COMPOSITION_BASED_STATISTICS=on&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=(none)&amp;EXPECT=1000&amp;FILTER=L&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;I_THRESH=0.005&amp;MATRIX_NAME=BLOSUM62&amp;NCBI_GI=on&amp;PAGE=Proteins&amp;PROGRAM=blastp&amp;SERVICE=plain&amp;SET_DEFAULTS.x=41&amp;SET_DEFAULTS.y=5&amp;SHOW_OVERVIEW=on&amp;END_OF_HTTPGET=Yes&amp;SHOW_LINKOUT=yes&amp;QUERY=");
   pLink = getenv("COMETLINKSFILE");
   if (pLink != NULL)
   {
      FILE *fp;

      if ((fp = fopen(pLink, "r")) != NULL)
      {
         while (fgets(szBuf, SIZE_BUF, fp))
         {
            if (!strncmp(szBuf, "PEPTIDELINK=", 12))
            {
               strcpy(szPeptideLink, szBuf + 12);
               szPeptideLink[strlen(szPeptideLink)-1]='\0';
               break;
            }
         }
         fclose(fp);
      }
   }
   for (i = 1; i < argc; i++)
   {
      // Stat the argument name to see if it's a valid file/dir
      // and to obtain information about it.
      if (!stat(argv[i], &fileStat))
      {
         // If it's a directory lets get all the files 
         // contained within and process them if they
         // are .out files.
         if (S_ISDIR(fileStat.st_mode))
         {
            char *fname;
            // Since this is a directory name and we
            // will be appending filenames to it do we
            // have a separator at the end?
            if (*(argv[i] + strlen(argv[i])) != SEP)
            {
               // Append separator
               sprintf(szPathName, "%s%c", argv[i], SEP);
            }
            else
            {
               strcpy(szPathName, argv[i]);
            }
            if (!strcmp(szPathName, "./") || !strcmp(szPathName,".\\"))
               szPathName[0]='\0';

#ifdef _MSC_VER
            sprintf(szFileMask,"%s*.out",szPathName);
            hFile = _findfirst( szFileMask, &c_file );
            if (hFile!=-1) {
               do {
                  fname = c_file.name;
#else
            dirStream = (DIR *) opendir(argv[i]);

            // Process all files in the directory
            while ((dirEntry = readdir((DIR *) dirStream)))
            {
               // If it contains ".out" process it
               if (strstr(fname=dirEntry->d_name, ".out"))
               {
#endif
                  if (iOutFileCount == iEntriesAllocated)
                  {
                     struct SequestOutStruct *pTmp;

                     pTmp = (SequestOutStruct *)realloc(sequestData, (iOutFileCount + 100) * sizeof(*sequestData));
                     if (pTmp)
                     {
                        sequestData = pTmp;
                        iEntriesAllocated = iOutFileCount + 100;
                     }
                     else
                     {
                        fprintf(stdout, "\n Error - cannot realloc memory for more than %d lines.\n\n", iOutFileCount);
                        fprintf(stdout, " Either create multiple summary files, each with less than %d\n", iOutFileCount);
                        fprintf(stdout, " .out file as input or increase memory.\n\n");
                        exit(1);
                     }
                  }
                  sprintf(szFullQName, "%s%s", szPathName, fname);
                  if (readOutFile(szFullQName, &sequestData[iOutFileCount++], &headerData)) {
                     iOutFileCount--; // that last one was bogus
                  }
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
            // Argument is a filename process if it contains
            // ".out" string.
            if (strstr(argv[i], ".out"))
            {
               if (iOutFileCount == iEntriesAllocated)
               {
                  struct SequestOutStruct *pTmp;

                  pTmp = (SequestOutStruct *)realloc(sequestData, (iOutFileCount + 100) * sizeof(*sequestData));
                  if (pTmp)
                  {
                     sequestData = pTmp;
                     iEntriesAllocated = iOutFileCount + 100;
                  }
                  else
                  {
                     fprintf(stdout, "\n Error - cannot realloc memory for more than %d lines.\n\n", iOutFileCount);
                     fprintf(stdout, " Either create multiple summary files, each with less than %d\n", iOutFileCount);
                     fprintf(stdout, " .out file as input or increase memory.\n\n");
                     exit(1);
                  }
               }
               if (readOutFile(argv[i], &sequestData[iOutFileCount++], &headerData)) {
                 iOutFileCount--; // that last one was bogus
               }
            }
         }
      }
      else
      {
         printf("Error %s is not a file or a directory!\n", argv[i]);
      }
   }
   printSummary(sequestData, &headerData, iOutFileCount, szCWD, szPeptideLink);

   return(0);
}


  /** readOutFile parses a Sequest or TurboSequest *.out
   *  file to obtain the data necessary to create an
   *  INTEARCT file.  The products of the parse are placed
   *  in two data structures: SequestOutStruct holds the
   *  data from the first (highest scoring) data line in
   *  the out file, and HeaderStruct which hold the data
   *  used in the header of a INTERACT file.  
   *
   *  NOTE: This is an embarassingly fragile parser.
   *  This should not be used in a mission critical 
   *  setting unless it has first been beef'd up with
   *  some subtantial amount of error checking code.
   *
   */
int readOutFile(char *szFileName,
                struct SequestOutStruct *data,
                struct HeaderStruct *hdr)
{
   FILE *ppIn;

   int i,
       j,
       k,
       ii,
       iLenPep,
       iRankXC,
       iId,
       iNum1,
       iLen,
       iTmp,
       iRank,
       iSLen,
       iDiffInSequence = 0,
       iMinLength = 0,
       bTurboSequest = 0,
       bHasSF = 0,  /* has SF column, bioworks 3.2 */
       bHasP = 0,   /* has P column, bioworks 3.2 */
       bNoDeltCnYet = 1,
       bIncludeMod = 0;

   double dNextdCn, dTmp;

   char szBuf[SIZE_BUF],
        szMod[SIZE_BUF],
        szTmp[SIZE_BUF],
        szNextDup[SIZE_BUF],
        szNextPep[SIZE_BUF],
        szPlainNextPep[SIZE_BUF],
        szPep[200],
        *pStr;

   // Initialize Data Structure
   data->cAA1 = '\0';
   data->cAA2 = '\0';
   data->szFileName[0] = '\0';
   data->szBaseFileName[0] = '\0';
   data->szProt[0] = '\0';
   data->szPlainPep[0] = '\0';
   data->szSubPep[0] = '\0';
   data->szDSite[0] = '\0';
   data->szMod[0] = '\0';
   data->szDup[0] = '\0';
   data->szDatabase[0] = '\0';
   data->dAMass = 0.0;
   data->dMass = 0.0;
   data->dXC = 0.0;
   data->dDeltCn = 0.0;
   data->dSp = 0.0;
   data->dMass1 = 0.0;
   data->dMass2 = 0.0;
   data->dMass3 = 0.0;
   data->dMass4 = 0.0;
   data->dMass5 = 0.0;
   data->dMass6 = 0.0;
   data->iRankSp = 0;
   data->iMassType = 0;
   data->iIon = 0;
   data->iTot = 0;
   data->bSpecialDeltCn = 0;
   data->bNucDb = 0;

   if ((ppIn = fopen(szFileName, "r")) == NULL)
   {
      printf(" Error - can't read .out file %s: %s\n\n", szFileName, strerror(errno));
      exit(1);
   }

   // Store a copy of the full filename
   strcpy(data->szFileName, szFileName);

   // Create a base filename from the full filename
   iLen = strlen(szFileName);
   for (i = iLen - 1; i > 0; i--)
      if (szFileName[i] == '/' || szFileName[i] == '\\')
         break;

   if ((i == 0 && (szFileName[0] == '/' || szFileName[0] == '\\')) || i > 0)
      strcpy(data->szBaseFileName, szFileName + i + 1);
   else
      strcpy(data->szBaseFileName, szFileName);

   // remove the .out from the end of data->szBaseFileName
   data->szBaseFileName[ strlen(data->szBaseFileName) - 4 ] = '\0';

   // Eat a fixed number of lines....(I know it's bogus)
   char *fgot=fgets(szBuf, SIZE_BUF, ppIn);
   fgot=fgets(szBuf, SIZE_BUF, ppIn);
   fgot=fgets(szBuf, SIZE_BUF, ppIn);
   fgot=fgets(szBuf, SIZE_BUF, ppIn);
#ifndef __SORCERER__
   fgot=fgets(szBuf, SIZE_BUF, ppIn);
#endif

   if (feof(ppIn)) { // sometimes SEQUEST produces empty .out files
      fprintf(stderr,"skipping empty file \"%s\"\n",szFileName);
      fclose(ppIn);
      return 1; // error
   }

   // Grab time/date
   fgot=fgets(szBuf, SIZE_BUF, ppIn);
   sscanf(szBuf, "%s %s %s ", hdr->szDate, hdr->szTime, hdr->szTimeSuffix);

   // Parse mass type
   fgot=fgets(szBuf, SIZE_BUF, ppIn);
   if (strstr(szBuf, "/MONO"))
      data->iMassType = 1;
   else if (strstr(szBuf, "/AVG"))
      data->iMassType = 0;

   // Grab mass type string
   iLen = strlen(szBuf);
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

   // Grab acquired mass
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
   //     display top 10/5, ion % = 0.0, CODE = 001020
   fgot=fgets(szBuf, SIZE_BUF, ppIn);

   // Modification lines
   //     (M# +16.000) (C@ +8.000) C=545.339  {C}
   fgot=fgets(szBuf, SIZE_BUF, ppIn);
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

         iLen = strlen(szBuf);

         for (i = 1; i < iLen; i++)
         {
            if (szBuf[i] == '=' && szBuf[i-1]!='m')
            {
               if (szBuf[i-1]!='m' && szBuf[i-1]!='p' && szBuf[i-1]!='t')  /* make sure not terminal mod */
               {
                  double dMass = 0.0;

                  sscanf(szBuf + i + 1, "%lf", &(data->dMass));
                  sprintf(data->szMod + strlen(data->szMod), "&amp;Mass%c=%0.4f", szBuf[i - 1], data->dMass);
               }
            }
         }
      }

      /* look for static N-terminal mods */
      if ((pStr = strstr(szBuf, "+N-term=")))
      {
         double dMass = 0.0;
         sscanf(pStr+8, "%lf", &dMass);
         sprintf(data->szMod + strlen(data->szMod), "&amp;Nterm=%0.4f", dMass);
      }
      if ((pStr = strstr(szBuf, "+Nterm-pep=")))
      {
         double dMass = 0.0;
         sscanf(pStr+11, "%lf", &dMass);
         sprintf(data->szMod + strlen(data->szMod), "&amp;Nterm=%0.4f", dMass);
      }
      if ((pStr = strstr(szBuf, "+Nterm-prot=")))
      {
         double dMass = 0.0;
         sscanf(pStr+12, "%lf", &dMass);
         sprintf(data->szMod + strlen(data->szMod), "&amp;Nterm=%0.4f", dMass);
      }

      /* look for static C-terminal mods */
      if ((pStr = strstr(szBuf, "+C-term=")))
      {
         double dMass = 0.0;
         sscanf(pStr+8, "%lf", &dMass);
         sprintf(data->szMod + strlen(data->szMod), "&amp;Cterm=%0.4f", dMass);
      }
      if ((pStr = strstr(szBuf, "+Cterm-pep=")))
      {
         double dMass = 0.0;
         sscanf(pStr+11, "%lf", &dMass);
         sprintf(data->szMod + strlen(data->szMod), "&amp;Cterm=%0.4f", dMass);
      }
      if ((pStr = strstr(szBuf, "+Cterm-prot=")))
      {
         double dMass = 0.0;
         sscanf(pStr+12, "%lf", &dMass);
         sprintf(data->szMod + strlen(data->szMod), "&amp;Cterm=%0.4f", dMass);
      }

      /* look for variable C-terminal mods */

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

   // Now we have the real data
   fgot=fgets(szBuf, SIZE_BUF, ppIn);
   iLen = strlen(szBuf);

   if (iLen>5)
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

      szTmp[i]='\0';

      if (bTurboSequest)
      {
         char szSF_P[500];

         if (bHasSF && bHasP)
         {
            sscanf(szTmp, "%d. %d %d %d %lf %lf %lf %lf %s %s %d %d %s %s %s",
                   &iNum1, &iRankXC, &(data->iRankSp), &iId,
                   &(data->dMass), &(data->dDeltCn), &(data->dXC), &(data->dSp), szSF_P, szSF_P,
                   &(data->iIon), &(data->iTot), data->szProt, data->szDup, szPep);
         }
         else if (bHasSF || bHasP)
         {
            sscanf(szTmp, "%d. %d %d %d %lf %lf %lf %lf %s %d %d %s %s %s",
                   &iNum1, &iRankXC, &(data->iRankSp), &iId,
                   &(data->dMass), &(data->dDeltCn), &(data->dXC), &(data->dSp), szSF_P,
                   &(data->iIon), &(data->iTot), data->szProt, data->szDup, szPep);
         }
         else
         {
            sscanf(szTmp, "%d. %d %d %d %lf %lf %lf %lf %d %d %s %s %s",
                   &iNum1, &iRankXC, &(data->iRankSp), &iId,
                   &(data->dMass), &(data->dDeltCn), &(data->dXC), &(data->dSp),
                   &(data->iIon), &(data->iTot), data->szProt, data->szDup, szPep);
         }
      }
      else
      {
         sscanf(szTmp, "%d. %d %d %lf %lf %lf %lf %d %d %s %s %s",
                &iNum1, &iRankXC, &(data->iRankSp),
                &(data->dMass), &(data->dDeltCn), &(data->dXC), &(data->dSp),
                &(data->iIon), &(data->iTot), data->szProt, data->szDup, szPep);
      }
   
      if (data->szDup[0] != '+')
      {
         strcpy(szPep, data->szDup);
         data->szDup[0] = '\0';
      }
   
      strcpy(data->szSubPep, szPep + 2);
      data->szSubPep[strlen(data->szSubPep) - 2] = '\0';

      data->cAA1 = szPep[0];
      data->cAA2 = szPep[strlen(szPep) - 1];
   
      iLenPep = strlen(data->szSubPep);
      ii = 0;
      for (i = 0; i < iLenPep; i++)
      {
         if (!isalpha(data->szSubPep[i]))
         {
            if (data->szSubPep[i] == '*')
               data->szDSite[ii - 1] = '1';
            else if (data->szSubPep[i] == '#')
               data->szDSite[ii - 1] = '2';
            else if (data->szSubPep[i] == '@')
               data->szDSite[ii - 1] = '3';
            else if (data->szSubPep[i] == '^')
               data->szDSite[ii - 1] = '4';
            else if (data->szSubPep[i] == '~')
               data->szDSite[ii - 1] = '5';
            else if (data->szSubPep[i] == '$')
               data->szDSite[ii - 1] = '6';
         }
         else
         {
            data->szPlainPep[ii] = data->szSubPep[i];
            data->szDSite[ii] = '0';
            ii++;
         }
      }
      data->szPlainPep[ii] = '\0';
      data->szDSite[ii] = '\0';
   
      // Continue reading the file to find dCn value
      szBuf[0] = '\0';
      while (fgets(szBuf, SIZE_BUF, ppIn) != NULL)
      {
         /*
          * at end of info
          */
         if (strlen(szBuf) <= 5)
         {
            break;
         }
         /*
          * Make sure it's not a duplicate accession number line.
          */
         else if (bNoDeltCnYet && strncmp(szBuf, "      ", 6))
         {
            char szTmpProt[SIZE_BUF];
            int iSlashCount=0;

            iLen = strlen(szBuf);

            // only substitute first two '/' for rank and ion columns
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
            szTmp[i]='\0';

            szNextDup[0]='\0';
            szNextPep[0]='\0';
            if (bTurboSequest)
            {
               char szSF_P[500];

               if (bHasSF && bHasP)
                  sscanf(szTmp, "%d. %d %d %d %lf %lf %lf %lf %s %s %d %d %s %s %s",
                      &iRank, &iTmp, &iTmp, &iTmp, &dTmp, &dNextdCn,
                      &dTmp, &dTmp, szSF_P, szSF_P, &iTmp, &iTmp, szTmpProt, szNextDup, szNextPep);
               else if (bHasSF || bHasP)
                  sscanf(szTmp, "%d. %d %d %d %lf %lf %lf %lf %s %d %d %s %s %s",
                      &iRank, &iTmp, &iTmp, &iTmp, &dTmp, &dNextdCn,
                      &dTmp, &dTmp, szSF_P, &iTmp, &iTmp, szTmpProt, szNextDup, szNextPep);
               else
                  sscanf(szTmp, "%d. %d %d %d %lf %lf %lf %lf %d %d %s %s %s",
                      &iRank, &iTmp, &iTmp, &iTmp, &dTmp, &dNextdCn,
                      &dTmp, &dTmp, &iTmp, &iTmp, szTmpProt, szNextDup, szNextPep);
            }
            else
            {
               sscanf(szTmp, "%d. %d %d %lf %lf %lf %lf %d %d %s %s %s",
                      &iRank, &iTmp, &iTmp, &dTmp, &dNextdCn,
                      &dTmp, &dTmp, &iTmp, &iTmp, szTmpProt, szNextDup, szNextPep);
            }

            if (strlen(szNextPep) == 0)
            {
               strcpy(szNextPep, szNextDup);
            }
            strcpy(szNextPep, szNextPep + 2);
            szNextPep[strlen(szNextPep) - 2] = '\0';

            /*
             * remove modification characters
             */
            k = 0;
            iLen = strlen(szNextPep);
            for (j = 0; j < iLen; j++)
            {
               if (isalpha(szNextPep[j]))
               {
                  szPlainNextPep[k] = szNextPep[j];
                  k++;
               }
            }
            szPlainNextPep[k] = '\0';

            iDiffInSequence = (int) labs(strlen(szPlainNextPep) - strlen(data->szPlainPep));

            iMinLength = strlen(data->szPlainPep);
            if ((int) strlen(szPlainNextPep) < iMinLength)
               iMinLength = strlen(szPlainNextPep);
 
            for (i = 0; i < iMinLength; i++)
            {
               /* K/Q and I/L don't count as differences */
               if (data->szPlainPep[i] != szPlainNextPep[i])
               {
                  if (!
                      ((data->szPlainPep[i] == 'K' || data->szPlainPep[i] == 'Q')
                       && (szPlainNextPep[i] == 'K' || szPlainNextPep[i] == 'Q'))
                      &&
                      !((data->szPlainPep[i] == 'I'
                         || data->szPlainPep[i] == 'L')
                        && (szPlainNextPep[i] == 'I'
                            || szPlainNextPep[i] == 'L')))
                  {
                     iDiffInSequence++;
                  }
               }
            }

            /*
             * Calculate deltCn only if sequences are less than
             * PERCENTAGE similar;  PERCENTAGE=0.75 for no good reason
             */
            if ((double) ((double) (strlen(data->szPlainPep) - 3.0 - iDiffInSequence) /
                 (double) (strlen(data->szPlainPep) - 3.0)) < PERCENTAGE)
            {
               data->dDeltCn = dNextdCn;
               bNoDeltCnYet = 0;
               if (iRank != 2)
                  data->bSpecialDeltCn = 1;
               else
                  data->bSpecialDeltCn = 0;
            }
         }
      }                            // while

      if (bNoDeltCnYet)
      {
         // Special dCn because there wasn't a second ranked score
         data->dDeltCn = 1.0;
         data->bSpecialDeltCn = 1;
      }
   }
   else
   {
      data->szSubPep[0]='\0';
   }

   fclose(ppIn);

   return (0);

}                               /*readOutFile */


  /** printSummary writes to stdout the contents
   *  of the sequestDataStruct and headerStruct in
   *  INTERACT format.
   */
void printSummary(struct SequestOutStruct *data,
                 struct HeaderStruct *hdr,
                 int lines,
                 char *szCWD,
                 char *szPeptideLink)
{
   int i = 0,
       j = 0,
       iMaxFileNameWidth = 0,
       iMaxRefWidth = 0,
       iMaxSeqWidth = 0;

   // Get max column widths
   for (i = 0; i < lines; i++)
   {
      // FileName
      if ((int) strlen(data[i].szBaseFileName) > iMaxFileNameWidth)
         iMaxFileNameWidth = strlen(data[i].szBaseFileName);

      // Ref
      if ((int) strlen(data[i].szProt) + (int) strlen(data[i].szDup) > iMaxRefWidth)
         iMaxRefWidth = strlen(data[i].szProt)+strlen(data[i].szDup);

      // Sequence
      if ((int) strlen(data[i].szSubPep) > iMaxSeqWidth)
         iMaxSeqWidth = strlen(data[i].szSubPep);
   }
   iMaxRefWidth += 1;

   // Print out easy part of the header
   printf("<HTML>\n");
   printf("<HEAD><TITLE>HTML-SUMMARY</TITLE></HEAD>\n");
   printf("<BODY BGCOLOR=\"#FFFFFF\">\n");
   printf("<PRE><FONT COLOR=\"green\">");
   printf("HTML-SUMMARY</FONT>\n");
   printf("Institute for Systems Biology\n");
   printf("Seattle, WA\n");
   printf("%s %s %s %s, %s\n", hdr->szDate, hdr->szTime, hdr->szTimeSuffix, data[0].szDatabase, hdr->szMassType);
   printf("\n");
   printf("<FONT COLOR=\"green\">");

   // Print complex portion of the header
   printf("   #    File");
   for (i = 0; i < iMaxFileNameWidth - 4; i++)
      printf(" ");

   printf("    MH+                  XCorr ");
   printf("   dCn      Sp    RSp    Ions");

   printf("   Ref");
   for (i = 0; i < iMaxRefWidth; i++)
      printf(" ");

   printf("  Sequence</FONT>\n");
   printf("<FONT COLOR=\"green\">");
   printf("  ----  ");
   for (i = 0; i < iMaxFileNameWidth; i++)
      printf("-");
   printf("  -------------------    ------");
   printf("  -----   ------  ---   ------  ");
   for (i = 0; i < iMaxRefWidth; i++)
      printf("-");
   printf("  ");
   for (i = 0; i < iMaxSeqWidth; i++)
      printf("-");
   printf("</FONT>\n");

   for (i = 0; i < lines; i++)
   {

      // Index
      printf(" %5d", i + 1);

      // File
      if (data[i].szFileName[0] != '/')     /* add extra /./ in path before filename for Pep3D compatibility */
         printf ("  <A TARGET=\"Win1\" HREF=\"/%s/%s?OutFile=%s/./%s\">./%s</A>",    
             CGI_DIR, OUT_CGI, szCWD, data[i].szFileName, data[i].szBaseFileName);
      else
         printf ("  <A TARGET=\"Win1\" HREF=\"/%s/%s?OutFile=%s\">%s</A>",
             CGI_DIR, OUT_CGI, data[i].szFileName, data[i].szBaseFileName);
      for (j = 0; j < iMaxFileNameWidth - (int) strlen(data[i].szBaseFileName); j++)
         printf(" ");

      if (data[i].szSubPep[0]!='\0')
      {
         // MH+
         printf("  %9.4f", data[i].dMass);
         // MHError
         printf(" (%+6.4f)", data[i].dAMass - data[i].dMass);
         // Xcorr
         printf("  %0.4f", data[i].dXC);
         // dCn 
         if (data[i].dDeltCn >= 0.2)
         {
            if (data[i].bSpecialDeltCn)
               printf("  <FONT COLOR=\"#DD00DD\">%0.3f*</FONT>", data[i].dDeltCn);
            else
               printf("  <FONT COLOR=\"#DD00DD\">%0.3f </FONT>", data[i].dDeltCn);
         }
         else
         {
            if (data[i].bSpecialDeltCn)
               printf("  %0.3f*", data[i].dDeltCn);
            else
               printf("  %0.3f ", data[i].dDeltCn);
         }
         // Sp
         printf("  %6.1f", data[i].dSp);
         // Rsp
         if (data[i].iRankSp == 1)
         {
            printf("  <FONT COLOR=\"#DD00DD\">%3d</FONT>", data[i].iRankSp);
         }
         else
         {
            printf("  %3d", data[i].iRankSp);
         }
         // Ions TODO: This needs special parameters based on static/variable 
         // mods
         if (data[i].szBaseFileName[0] != '/')
            printf("  <A TARGET=\"Win1\" HREF=\"/%s/%s?Dta=%s/%s.dta&amp;MassType=%d&amp;NumAxis=1",
                CGI_DIR, PLOT_CGI, szCWD, data[i].szBaseFileName, data[i].iMassType);
         else
            printf("  <A TARGET=\"Win1\" HREF=\"/%s/%s?Dta=%s.dta&amp;MassType=%d&amp;NumAxis=1",
                CGI_DIR, PLOT_CGI, data[i].szBaseFileName, data[i].iMassType);
         if (data[i].dMass1 != 0.0)
            printf("&amp;DMass1=%f", data[i].dMass1);
         if (data[i].dMass2 != 0.0)
            printf("&amp;DMass2=%f", data[i].dMass2);
         if (data[i].dMass3 != 0.0)
            printf("&amp;DMass3=%f", data[i].dMass3);
         if (data[i].dMass4 != 0.0)
            printf("&amp;DMass4=%f", data[i].dMass4);
         if (data[i].dMass5 != 0.0)
            printf("&amp;DMass5=%f", data[i].dMass5);
         if (data[i].dMass6 != 0.0)
            printf("&amp;DMass6=%f", data[i].dMass6);
         if (strlen(data[i].szSubPep) != strlen(data[i].szPlainPep))
            printf("&amp;DSite=%s", data[i].szDSite);

         printf("%s&amp;Pep=%s\">", data[i].szMod, data[i].szPlainPep);
         if (data[i].iIon<100)
            printf(" ");
         if (data[i].iIon<10)
            printf(" ");
         printf("%d/%d", data[i].iIon, data[i].iTot);
         printf("</A>");
         if (data[i].iTot<10)
            printf(" ");
         if (data[i].iTot<100)
            printf(" ");

   
         // Ref
         printf("  <A TARGET=\"Win1\" HREF=\"/%s/%s?Ref=%s&amp;Db=%s&amp;NucDb=%d&amp;Pep=%s&amp;MassType=%d\">%s</A>",
             CGI_DIR, COMETDB_CGI, data[i].szProt, data[i].szDatabase, data[i].bNucDb,
             data[i].szPlainPep, data[i].iMassType, data[i].szProt);

         for (j = 0; j < iMaxRefWidth - (int) strlen(data[i].szProt) - (int) strlen(data[i].szDup); j++)
            printf(" ");
   
         // optional +field
         if (strlen(data[i].szDup) > 1)
         {
            printf("<A TARGET=\"Win1\" HREF=\"/%s/%s?Db=%s&amp;NucDb=%d&amp;Pep=%s&amp;MassType=%d\">%s</A>",
                CGI_DIR, COMETDB_CGI, data[i].szDatabase, data[i].bNucDb, data[i].szPlainPep,
                data[i].iMassType, data[i].szDup);
         }
         printf("  ");
   
         // Seq TODO: Breakup the peptide
         printf("%c.<A TARGET=\"Win1\" HREF=\"%s%s\">%s</A>.%c",
                data[i].cAA1,
                szPeptideLink,
                data[i].szPlainPep,
                data[i].szSubPep,
                data[i].cAA2);

         for (j = 0; j < iMaxSeqWidth - (int) strlen(data[i].szSubPep) + 1; j++)
            printf(" ");
         printf("\n");

      }
      else
      {
         if (data[i].szBaseFileName[0] != '/')
            printf("  <A TARGET=\"Win1\" HREF=\"/%s/%s?Dta=%s/%s.dta&amp;NumAxis=1",
               CGI_DIR, PLOT_CGI, szCWD, data[i].szBaseFileName);
         else
            printf("  <A TARGET=\"Win1\" HREF=\"/%s/%s?Dta=%s.dta&amp;NumAxis=1",
               CGI_DIR, PLOT_CGI, data[i].szBaseFileName);

         printf("\">spectrum</A>\n");
      }
   }

   printf("\n</BODY></HTML>\n");

}
