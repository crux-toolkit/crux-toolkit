/*************************************************************
 
       Date:  01/21/2009
 
    Purpose:  Add pI column to an input text file where peptide is in 1st column
 
*************************************************************/ 


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// calc_pI should really be separated into a header
//#include "calc_pI.c" 

#define SIZE_FILE        256
#define SIZE_BUF         4096

#define TRUE  1
#define FALSE 0


void SET_OPTION(char *arg);
 
void ADD_PI(char *szArg);

double COMPUTE_PI(char *seq,unsigned long seq_length, int charge_increment);


int main(int argc, char **argv)
{
   int  i,
        iStartArgc,
        iNumArg,
        iNumInputs; 
   char *arg;

   printf("\n");
   printf(" ADD pI by Jimmy Eng\n");
   printf(" Copyright (c) Institute for Systems Biology, 2000.  All rights reserved.\n\n");

   if (argc<2)
   {
      printf(" USAGE:  %s inputfile.txt\n\n", argv[0]);
      printf("\n");
      printf(" Where inputfile.txt contains peptide sequences, one on each row.\n");
      printf(" This program will add a pI column to the file.\n\n");
      exit(EXIT_FAILURE);
   }
   

   iNumInputs=argc;
 
   iNumArg=0;
   arg = argv[iNumArg = 1];
   iStartArgc=1;
   while (iNumArg < argc)
   {
      if (arg[0] == '-')
         SET_OPTION(arg);
      else
         break;

      iStartArgc++;
      arg = argv[++iNumArg];
   }
 
   for (i=iStartArgc; i<iNumInputs; i++)
   {
      ADD_PI(argv[i]);
   }

   printf("\n Done.\n\n");

   return(EXIT_SUCCESS);

} /*main*/


void ADD_PI(char *szArg)
{
   int  bSkip,
        iAnalysisFirstScan,
        iAnalysisLastScan;

   char szBuf[SIZE_BUF],
        szInputFile[SIZE_FILE];

   FILE *fpIn;

   bSkip=FALSE;

   strcpy(szInputFile, szArg);

   printf(" Adding peptide pI to file %s\n", szInputFile);

  /*
   * open each html file
   */
   if ( (fpIn=fopen(szInputFile, "r"))==NULL)
   {
      printf(" Error - cannot open file \"%s\".\n\n", szInputFile);
      bSkip=TRUE;
   } 

   if (!bSkip)
   {
      FILE *fpOut;

      int  ii,
           iPrevLen=0,
           iLen,
           bHasMWColumn;

      char szTmpFile[SIZE_FILE],
           szBuf2[SIZE_BUF],
           szTmp[SIZE_BUF];

      strcpy(szTmpFile, szInputFile);
      sprintf(szTmpFile+strlen(szTmpFile)-4,"tmp");

      if ( (fpOut=fopen(szTmpFile, "w"))==NULL)
      {
         printf(" Error - cannot open output file \"%s\".\n\n", szTmpFile);
      }
      else
      {
        /*
         * Read each summary line, pull out info,
         * generate PI & and print out line
         */
         while (fgets(szBuf, SIZE_BUF, fpIn))
         {
            char szPeptide[SIZE_FILE];
            double dPI;
            char szPI[40];
            unsigned long lLen;

            if (strlen(szBuf)<3)
               break;

            sscanf(szBuf, "%s", szPeptide);

            szBuf[strlen(szBuf)-1]=0;


            lLen=strlen(szPeptide);
            dPI=COMPUTE_PI(szPeptide, lLen, 0);

            fprintf(fpOut, "%s\t%0.2f\n", szBuf, dPI);
         }
   
         fclose(fpIn);

         sprintf(szTmp, "mv %s %s", szTmpFile, szInputFile);
         int dummy = system(szTmp);
      } /*else*/

   } /*if*/

} /*ADD_PI*/


void SET_OPTION(char *arg)
{
   int  sortby;
   char szTempPath[1024];
   double dTmpMass;
 
   switch (arg[1])
   {
      default:
         break;
   }
   arg[0] = '\0';
 
} /*SET_OPTION*/
