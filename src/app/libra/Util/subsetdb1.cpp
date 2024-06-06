/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Library or "Lesser" General Public      *
 *   License (LGPL) as published by the Free Software Foundation;          *
 *   either version 2 of the License, or (at your option) any later        *
 *   version.                                                              *
 *                                                                         *
 ***************************************************************************/


/*
 * Create subset sequence databases; reverse database; database statistics
 * Jimmy Eng, Copyright 2004, all rights reserved
 *
 * date initiated: 2004/06/17
 * 
 * $Id: subsetdb1.cpp 7817 2018-09-24 23:30:14Z dshteyn $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <ctime>


#define INIT_SEQ_LEN         1000
#define SIZE_BUF             500000
#define SIZE_FILE            1024
#define MAX_FILTER           50000
#define SIZE_STRING          64

int  iAACount[256];

struct OptionsStruct
{
   int  iNumMatch;
   int  iNumNoMatch;
   int  bOutputDir;          // default FALSE=write to current dir; TRUE=write to input db dir
   int  bVerbose;
   int  bCaseSensitive;      // default FALSE=case insensitive, TRUE=case sensitive
   int  bReverseDb;          // default FALSE, TRUE=reverse sequences
   int  bPrintStats;         // default FALSE, TRUE=print db statistics
   int  bRemoveM;
   int  bMutuallyExclusive;  // default FALSE, TRUE= apply pNoMatch separately
   int  bPrintLength;
   int  iRandomCount;
   int  iSubsetSize;
   int  iSeqMaxLength;
   char szReverseText[SIZE_FILE];
};

void SET_OPTION(char *arg,
      char *szOutput,
      char ** pMatch,
      char ** pNoMatch,
      struct OptionsStruct *pOptions);

void SUBSET_DBASE(char *szInput,
      char *szOutput,
      char **pMatch,
      char **pNoMatch,
      struct OptionsStruct pOptions);

int main(int argc, char **argv)
{
   int  i;
   int  iNumArg;
   int  iStartArgc;
   char szInput[SIZE_FILE];
   char szOutput[SIZE_FILE];
   char *arg;
   struct OptionsStruct pOptions;

   char** pMatch = (char**)malloc(MAX_FILTER*sizeof(char*));
   char** pNoMatch = (char**)malloc(MAX_FILTER*sizeof(char*));

   for (i=0; i<MAX_FILTER; i++) {
     pMatch[i] = (char*)malloc(SIZE_STRING*sizeof(char));
     pNoMatch[i] = (char*)malloc(SIZE_STRING*sizeof(char));
   }
   
   printf("\n");
   printf(" SubsetDB by J.Eng\n");
   printf("\n");

   if (argc<2)
   {
      printf(" USAGE:  %s [options] fasta.database\n", argv[0]);
      printf("\n");
      printf(" Options:   -M<str>    strings to match\n");
      printf("            -N<str>    strings to not match\n");
      printf("                       Note: if -M or -N strings point to file names that exist,\n");
      printf("                             the match/nomatch accessions will be read from those\n");
      printf("                             files; one accession per line.\n");
      printf("            -e         by default, -N option only excludes entries that pass the -M\n");
      printf("                             criteria.  Using this option, both -M and -N options\n");
      printf("                             are mutually exclusive.\n");
      printf("            -O<str>    new subset database name (default appends .new to input name)\n");
      printf("            -P         put subset database in same directory as input database\n");
      printf("            -l<num>    maximum line length for sequences (default 80)\n");
      printf("            -C         use case-sensitive comparisons (default is case in-sensitive)\n");
      printf("            -R         reverse output database (\"rev_\" appended to description lines.)\n");
      printf("            -r<num>    pull out 'num' random sequences from database\n");
      printf("            -s<num>    create multiple database files, each with 'num' sequences in them\n");
      printf("            -D<str>    when reversing output, append <string> instead of \"rev_\" to description lines.\n");
      printf("                          using this parameter invokes -R automatically.\n");
      printf("            -S         print database statistics\n");
      printf("            -L         export protein accession and sequence length to output file\n");
      printf("            -V         verbose reporting\n");
      printf("            -m         cleave/remove N-term methionine from each sequence\n");
      printf("\n");
      printf(" A maximum of %d strings each can be specied for match/no match.  Spaces can\n", MAX_FILTER);
      printf(" be substituted in the the string comparison by using the carat (^) character,\n");
      printf(" for example 'bos^taurus'.\n");
      printf("\n");
      printf(" Output database is written to the current directory unless a full path output\n");
      printf(" is specified or the -P option is used.\n");
      printf("\n");
      printf(" If no strings match/no match or statistics options are specified then output\n");
      printf(" database is validated version of input database with non-alphabetic characters\n");
      printf(" stripped out of the sequence.\n");
      printf("\n");
      printf(" If no command line arguments are given, the entire input database will be read\n");
      printf(" in and written out.\n");
      printf("\n");
      exit(EXIT_FAILURE);
   }

   iNumArg=0;
   arg = argv[iNumArg = 1];
   iStartArgc=1;
   szInput[0]='\0';
   szOutput[0]='\0';
   pOptions.bOutputDir=0;
   pOptions.bVerbose=0;
   pOptions.bCaseSensitive=0;
   pOptions.bReverseDb=0;
   pOptions.bRemoveM=0;
   pOptions.bMutuallyExclusive=0;
   pOptions.bPrintStats=0;
   pOptions.bPrintLength=0;
   pOptions.iNumMatch=0;
   pOptions.iNumNoMatch=0;
   pOptions.iRandomCount=0;
   pOptions.iSubsetSize=0;
   pOptions.iSeqMaxLength=80;
   strcpy(pOptions.szReverseText, "rev_");

   while (iNumArg < argc)
   {
      if (arg[0] == '-')
         SET_OPTION(arg, szOutput, pMatch, pNoMatch, &pOptions);
      else
         break;

      iStartArgc++;
      arg = argv[++iNumArg];
   }

   // need to check if pMatch or pNoMatch entries point to a file; if so
   // read accessions from file.
   if (pOptions.iNumMatch>0)
   {
      for (i=0; i<pOptions.iNumMatch; i++)
      {
         FILE *fp;

         if ( (fp=fopen(pMatch[i], "r"))!=NULL)
         {
            char szBuf[SIZE_BUF];

            pOptions.iNumMatch=0;
            printf(" *** match file %s", pMatch[i]);
            fflush(stdout);
            while (fgets(szBuf, SIZE_BUF, fp))
            {
               char szAccession[SIZE_BUF];

               szAccession[0]=0;

               sscanf(szBuf, "%s", szAccession);
               if (strlen(szAccession)>0 && pOptions.iNumMatch<MAX_FILTER)
               {
                  strncpy(pMatch[pOptions.iNumMatch], szAccession, SIZE_STRING);
                  pOptions.iNumMatch++;
               }
            }
            printf(", %d strings to match\n", pOptions.iNumMatch);
            if (pOptions.iNumMatch == MAX_FILTER)
               printf(" *** warning, maximum number of match strings reached\n");
            break;
         }
      }
   }
   if (pOptions.iNumNoMatch>0)
   {
      for (i=0; i<pOptions.iNumNoMatch; i++)
      {
         FILE *fp;

         if ( (fp=fopen(pNoMatch[i], "r"))!=NULL)
         {
            char szBuf[SIZE_BUF];

            pOptions.iNumNoMatch=0;
            printf(" *** no-match file '%s'", pNoMatch[i]);
            fflush(stdout);
            while (fgets(szBuf, SIZE_BUF, fp))
            {
               char szAccession[SIZE_BUF];

               szAccession[0]=0;

               sscanf(szBuf, "%s", szAccession);
               if (strlen(szAccession)>0 && pOptions.iNumNoMatch<MAX_FILTER)
               {
                  strncpy(pNoMatch[pOptions.iNumNoMatch], szAccession, SIZE_STRING);
                  pOptions.iNumNoMatch++;
               }
            }
            printf(", %d strings to not match\n", pOptions.iNumNoMatch);
            if (pOptions.iNumNoMatch == MAX_FILTER)
               printf(" *** warning, maximum number of no-match strings reached\n");
            break;
         }
      }
   }

   if (iStartArgc == argc)
   {
      printf(" Please enter original database file on the command line.\n\n");
      exit(EXIT_FAILURE);
   }
   else
      strcpy(szInput, argv[iStartArgc]);


   if (szOutput[0]=='\0') // output not specified
   {
      if (pOptions.bOutputDir)     // use output full path
      {
         strcpy(szOutput, szInput);
         strcat(szOutput, ".new");
      }
      else
      {
         char *pStr;

         if ((pStr=strrchr(szInput, '/')))
         {
            strcpy(szOutput, pStr+1);
            strcat(szOutput, ".new");
         }
         else
         {
            strcpy(szOutput, szInput);
            strcat(szOutput, ".new");
         }
      }
   }
   else // output specified
   {
      if (pOptions.bOutputDir && isalnum(szOutput[0]))  // see if input path required
      {
         char *pStr;
         char szTmp[SIZE_FILE];

         strcpy(szTmp, szInput);
         if ((pStr=strrchr(szTmp, '/')))
         {
            *(pStr+1) ='\0';
            strcat(szTmp, szOutput);
            strcpy(szOutput, szTmp);
         }
      }
   }

   for (i=0; i<pOptions.iNumMatch; i++)
   {
      int ii;
      
      for (ii=0; ii<(int)strlen(pMatch[i]); ii++)
      {
         if (!pOptions.bCaseSensitive)
	    pMatch[i][ii]=toupper(pMatch[i][ii]);
 
         if (pMatch[i][ii]=='^')
            pMatch[i][ii]=' ';
      }
   }
   for (i=0; i<pOptions.iNumNoMatch; i++)
   {
      int ii;
      
      for (ii=0; ii<(int)strlen(pNoMatch[i]); ii++)
      {
         if (!pOptions.bCaseSensitive)
            pNoMatch[i][ii]=toupper(pNoMatch[i][ii]);

         if (pNoMatch[i][ii]=='^')
            pNoMatch[i][ii]=' ';
      }
   }


   if (pOptions.bVerbose)
   {
      printf(" match string:\n");
      for (i=0; i<pOptions.iNumMatch; i++)
         printf("  %d   %s\n", i+1, pMatch[i]);

      printf(" no match string:\n");
      for (i=0; i<pOptions.iNumNoMatch; i++)
         printf("  %d   %s\n", i+1, pNoMatch[i]);

      printf(" input %s\n", szInput);
      printf(" output %s\n", szOutput);
      if (pOptions.bReverseDb)
         printf(" decoy text: \"%s\"\n", pOptions.szReverseText);
      printf("\n");
   }

   SUBSET_DBASE(szInput, szOutput, pMatch, pNoMatch, pOptions);

   printf("\n Done.\n\n");
   return(EXIT_SUCCESS);

}


void SET_OPTION(char *arg,
      char *szOutput,
      char **pMatch,
      char **pNoMatch,
      struct OptionsStruct *pOptions)
{
   switch (arg[1])
   {
      case 'M':
         if (pOptions->iNumMatch<MAX_FILTER)
         {
            strncpy(pMatch[pOptions->iNumMatch], arg+2, SIZE_STRING);
            pOptions->iNumMatch += 1;
         }
         break;
      case 'N':
         if (pOptions->iNumNoMatch<MAX_FILTER)
         {
            strncpy(pNoMatch[pOptions->iNumNoMatch], arg+2, SIZE_STRING);
            pOptions->iNumNoMatch += 1;
         }
         break;
      case 'O':
         strcpy(szOutput, arg+2);
         break;
      case 'e':
         pOptions->bMutuallyExclusive = 1;
         break;
      case 'm':
         pOptions->bRemoveM = 1;
         break;
      case 'P':
         pOptions->bOutputDir = 1;
         break;
      case 'V':
         pOptions->bVerbose = 1;
         break;
      case 'C':
         pOptions->bCaseSensitive= 1;
         break;
      case 'R':
         pOptions->bReverseDb = 1;
         break;
      case 'r':
         sscanf(arg+2, "%d", &(pOptions->iRandomCount));
         if (pOptions->iRandomCount < 0)
            pOptions->iRandomCount = 0;
         break;
      case 's':
         sscanf(arg+2, "%d", &(pOptions->iSubsetSize));
         if (pOptions->iSubsetSize < 0)
            pOptions->iSubsetSize = 0;
         break;
      case 'l':
         sscanf(arg+2, "%d", &(pOptions->iSeqMaxLength));
         if (pOptions->iSeqMaxLength < 0)
            pOptions->iSeqMaxLength = 80;
         break;
      case 'D':
         strcpy(pOptions->szReverseText, arg+2);
         pOptions->bReverseDb = 1;
         break;
      case 'S':
         pOptions->bPrintStats = 1;
         break;
      case 'L':
         pOptions->bPrintLength = 1;
         break;
      default:
         break;
   }
   arg[0] = '\0';

}



void SUBSET_DBASE(char *szInput,
		  char *szOutput,
		  char **pMatch,
		  char **pNoMatch,
		  struct OptionsStruct pOptions)
{
   int  i;

   long lLenAllocated,
        lLenSeq,
        lMinLength=0,
        lMaxLength=0;

   long lNumMatch=0,
        lNumTot=0,
        lNumResidues=0;

   FILE *fp,
        *fpOut;


   char *pSeq;

   char * szBuf = (char*) malloc(SIZE_BUF*sizeof(char));

   char * szDef = (char*) malloc(SIZE_BUF*sizeof(char));  // protein def line from database
   char * szCaseDef = (char*) malloc(SIZE_BUF*sizeof(char));  // upper case version of above

   char* szSubsetOut = (char*) malloc(SIZE_FILE*sizeof(char));

   int iSubsetSizeCount = 0;
   int iSubsetSizeBatchCount = 1;

   srand(time(0));

   if ( (fp=fopen(szInput, "r"))==NULL)
   {
      printf(" Error - read:  %s\n\n", szInput);
      exit(EXIT_FAILURE);
   }

   if (pOptions.iSubsetSize > 0)
   {
      sprintf(szSubsetOut, "%s.%03d", szOutput, iSubsetSizeBatchCount);
      if ( (fpOut=fopen(szSubsetOut, "w"))==NULL)
      {
         printf(" Error - write:  %s\n\n", szSubsetOut);
         exit(EXIT_FAILURE);
      }
   }
   else
   {
      if ( (fpOut=fopen(szOutput, "w"))==NULL)
      {
         printf(" Error - write:  %s\n\n", szOutput);
         exit(EXIT_FAILURE);
      }
   }

   if ( (pSeq=(char*)malloc(INIT_SEQ_LEN))==NULL)
   {
      printf(" Error cannot malloc pSeq\n\n");
      exit(EXIT_FAILURE);
   }

   lLenAllocated=INIT_SEQ_LEN;

   if (pOptions.bPrintStats)
      memset(iAACount, 0, sizeof(iAACount));

   printf("\n running ...");
   fflush(stdout);

   while (fgets(szBuf, SIZE_BUF, fp))
   {
      if (szBuf[0]=='>') // definition line
      {
         int cResidue;
         int bPass;

         strcpy(szDef, szBuf);
         szDef[strlen(szDef)-1]='\0';

	 strcpy(szCaseDef, szDef);

         if (!pOptions.bCaseSensitive)
         {
            int ii;
            for (ii=0; ii<(int)strlen(szCaseDef); ii++)
	      if (szCaseDef[ii] >= 'a' && szCaseDef[ii] <= 'z')
		szCaseDef[ii] = toupper(szDef[ii]);
            szCaseDef[strlen(szCaseDef)]='\0';
         }

         lLenSeq=0;

         while ((cResidue=fgetc(fp)))
         {
            if (isalpha(cResidue))
            {
               pSeq[lLenSeq]=cResidue;
               lLenSeq++;

               if (pOptions.bPrintStats)
                  iAACount[toupper(cResidue)] += 1;

               if (lLenSeq == lLenAllocated-1)
               {
                  char *pStr;

                  lLenAllocated += 500;
                  pStr=(char*)realloc(pSeq, lLenAllocated);

                  if (pStr==NULL)
                  {
                     printf(" Error - cannot realloc pSeq[%ld]\n\n", lLenAllocated);
                     exit(EXIT_FAILURE);
                  }
                  pSeq=pStr;
               }
            }
            else if (feof(fp) || cResidue=='>')
            {
               int iReturn;

               iReturn=ungetc(cResidue, fp);
               if (iReturn!=cResidue)
               {
                  printf("Error with ungetc.\n\n");
                  fclose(fp);
                  fclose(fpOut);
                  exit(EXIT_FAILURE);
               }
               break;
            }
         }

         pSeq[lLenSeq]='\0';

         if (lNumTot==0)
         {
            lMaxLength = lLenSeq;
            lMinLength = lLenSeq;
         }
         else
         {
            if (lLenSeq > lMaxLength)
               lMaxLength = lLenSeq;
            if (lLenSeq < lMinLength)
               lMinLength = lLenSeq;
         }

         lNumTot++;
         lNumResidues += lLenSeq;

         if (pOptions.bPrintLength)
         {
            char szWord[256];

            sscanf(szDef+1, "%s", szWord);
            fprintf(fpOut, "%s\t%ld\n", szWord, lLenSeq);
         }
         else
         {

            bPass=0;

            if (pOptions.iNumMatch==0 && pOptions.iNumNoMatch==0)
               bPass=1;
            else
            {
               // check match strings
               for (i=0; i<pOptions.iNumMatch; i++)
                  if (strstr(szCaseDef, pMatch[i]))
                     bPass=1;

               // next check no match strings
               if (pOptions.bMutuallyExclusive)
               {
                  // With the -e option, if protein matches *any* of the -N
                  // entries, it will not be printed exported
                  int bTmp=1;

                  for (i=0; i<pOptions.iNumNoMatch; i++)
                  {
                     if (strstr(szCaseDef, pNoMatch[i]))
                     {
                        bTmp=0;
                        break;
                     }
                  }

                  if (bTmp==1)
                     bPass = 1;
               }
               else
               {
                  if (bPass)
                  {
                     for (i=0; i<pOptions.iNumNoMatch; i++)
                     {
                        if (strstr(szCaseDef, pNoMatch[i]))
                        {
                           bPass=0;
                        }
                     }
                  }
               }
            }

            if (pOptions.iRandomCount>0)
            {
               int xRan;
                     
               xRan=rand()%5+1;
               if (xRan != 1)
                  bPass = 0;
            }

            if (bPass)
            {
               lNumMatch++;

               if (pOptions.bReverseDb)
                  fprintf(fpOut, ">%s%s\n", pOptions.szReverseText, szDef+1);
               else
                  fprintf(fpOut, "%s\n", szDef);

               if (pOptions.bReverseDb)
               {
                  for (i=0; i<lLenSeq; i++)
                  {
                     fprintf(fpOut, "%c", pSeq[lLenSeq-i-1]);
                     if (!((i+1) % pOptions.iSeqMaxLength))
                        fprintf(fpOut, "\n");
                  }
               }
               else
               {
                  int iStart;

                  if (pOptions.bRemoveM && pSeq[0]=='M')
                     iStart=1;
                  else
                     iStart=0;
                  for (i=iStart; i<lLenSeq; i++)
                  {
                     fprintf(fpOut, "%c", pSeq[i]);
                     if (!((i+1) % pOptions.iSeqMaxLength))
                        fprintf(fpOut, "\n");
                  }
               }

               if ((i) % pOptions.iSeqMaxLength)
                  fprintf(fpOut, "\n");

               iSubsetSizeCount++;

               if (pOptions.iSubsetSize > 0 && iSubsetSizeCount == pOptions.iSubsetSize)
               {
                  fclose(fpOut);

                  iSubsetSizeBatchCount++;
                  iSubsetSizeCount = 0;
                  sprintf(szSubsetOut, "%s.%03d", szOutput, iSubsetSizeBatchCount);
                  if ( (fpOut=fopen(szSubsetOut, "w"))==NULL)
                  {
                     printf(" Error - write:  %s\n\n", szSubsetOut);
                     exit(EXIT_FAILURE);
                  }
               }

               if (lNumMatch == pOptions.iRandomCount)
                  break;
            }
         }
      }
   }

   free(pSeq);

   fclose(fp);

   if (pOptions.iRandomCount>0)
   {
      printf("\n\n Extracted %ld random entries.\n", lNumMatch);
   }

   if (pOptions.iNumMatch!=0 || pOptions.iNumNoMatch!=0 || pOptions.bReverseDb)
   {
      fclose(fpOut);

      printf("\n\n");
      if (!pOptions.bPrintLength && !pOptions.bPrintStats)
         printf(" %ld entries matched out of %ld total entries.\n", lNumMatch, lNumTot);
      printf(" File created:  %s\n", szOutput);
   }
   else if (pOptions.bPrintStats)
   {
      printf("\n\n");
      printf(" Statistics for file '%s'\n\n", szInput);
      printf(" Tot num of sequence entries: %ld\n", lNumTot);
      printf("         Tot num of residues: %ld\n", lNumResidues);
      printf("         Min sequence length: %ld\n", lMinLength);
      printf("         Max sequence length: %ld\n", lMaxLength);

      printf("         Residue Composition: A %d\n", iAACount['A']);
      if (iAACount['B'] > 0)
         printf("                              B %d\n", iAACount['B']);
      printf("                              C %d\n", iAACount['C']);
      printf("                              D %d\n", iAACount['D']);
      printf("                              E %d\n", iAACount['E']);
      printf("                              F %d\n", iAACount['F']);
      printf("                              G %d\n", iAACount['G']);
      printf("                              H %d\n", iAACount['H']);
      printf("                              I %d\n", iAACount['I']);
      if (iAACount['J'] > 0)
         printf("                              J %d\n", iAACount['J']);
      printf("                              K %d\n", iAACount['K']);
      printf("                              L %d\n", iAACount['L']);
      printf("                              M %d\n", iAACount['M']);
      printf("                              N %d\n", iAACount['N']);
      if (iAACount['O'] > 0)
         printf("                              O %d\n", iAACount['O']);
      printf("                              P %d\n", iAACount['P']);
      printf("                              Q %d\n", iAACount['Q']);
      printf("                              R %d\n", iAACount['R']);
      printf("                              S %d\n", iAACount['S']);
      printf("                              T %d\n", iAACount['T']);
      if (iAACount['U'] > 0)
         printf("                              U %d\n", iAACount['U']);
      printf("                              V %d\n", iAACount['V']);
      printf("                              W %d\n", iAACount['W']);
      if (iAACount['X'] > 0)
         printf("                              X %d\n", iAACount['X']);
      printf("                              Y %d\n", iAACount['Y']);
      if (iAACount['Z'] > 0)
         printf("                              Z %d\n", iAACount['Z']);
   }
    
}
