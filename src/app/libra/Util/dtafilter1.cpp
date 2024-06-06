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
 * DTAFILTER reads all .dta files in a directory and filters spectra
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "Common/util.h"

#define SIZE_BUF 8192
#define MAXARRAY 50000
#define MINMASS  600.0
#define MAXMASS  4200.0
#define MINPEAKS 5

#define TRUE  1
#define FALSE 0

struct OptionsStruct
{
   int  bAnalyzeAll;
   int  bDeleteBad;
   int  bVerboseOutput;
   int  iFirstScan;
   int  iLastScan;
   int  iMinPeakCount;
   int  iCharge;
   double dMinMass;
   double dMaxMass;
   double dMinPeakThreshHold;
   double dMinPeakDtaFilter;
};

void SET_OPTION(char *arg, struct OptionsStruct *pOptions);
int  ANALYZE(char *szDtaFile, struct OptionsStruct pOptions);

int main(int argc, char **argv)
{
   int  i, iNumArg, iStartArgc;
   char *arg, szBuf[SIZE_BUF], szCommand[SIZE_BUF];
   FILE *fp;
   struct OptionsStruct pOptions;

   if (argc < 2)
   {
      printf("\n");
      printf(" DTA-FILTER\n\n");
      printf(" Usage:  %s [options] *.dta\n", argv[0]);
      printf("\n");
      printf("     options = -F<num>   where num is an int specifying the first scan\n");
      printf("               -L<num>   where num is an int specifying the last scan\n");
      printf("               -C<num>   where num is an int specifying the precursor charge state to analyze\n");
      printf("               -B<num>   where num is a float specifying minimum MH+ mass, default=%0.1f Da\n", MINMASS);
      printf("               -T<num>   where num is a float specifying maximum MH+ scan, default=%0.1f Da\n", MAXMASS);
      printf("               -P<num>   where num is an int specifying minimum peak count, default=%d\n", MINPEAKS);
      printf("               -I<num>   where num is an int specifying minimum peak threshhold for peak count, default=1\n");
      printf("               -i<num>   where num is an int specifying minimum peak intensity filter, default=0\n");
      printf("               -A        analyze all .dta files in current directory.\n");
      printf("               -N        do not delete bad files, just append .bad extension to them\n");
      printf("               -V        verbose output\n");
      printf("\n");
      printf(" Minimum peak count is defined as # of non-zero integer normalized signals\n");
      printf("\n");
      exit(0);
   }

   iNumArg = 1;
   iStartArgc = 1;
   pOptions.iFirstScan = 0;
   pOptions.iLastScan = 99999;
   pOptions.iCharge = 0;
   pOptions.iMinPeakCount = MINPEAKS;
   pOptions.dMinMass = MINMASS;
   pOptions.dMaxMass = MAXMASS;
   pOptions.bAnalyzeAll = 0;
   pOptions.bDeleteBad = 1;
   pOptions.bVerboseOutput = 0;
   pOptions.dMinPeakThreshHold = 1;
   pOptions.dMinPeakDtaFilter = 0;

   arg = argv[iNumArg];
   while (iNumArg < argc)
   {
      if (arg[0] == '-')
         SET_OPTION(arg, &pOptions);
      else
         break;

      iStartArgc++;
      arg = argv[++iNumArg];
   }


   if (pOptions.bAnalyzeAll)
   {
      sprintf(szCommand, "find . -maxdepth 1 -depth -name '*.dta' -print");
   
      if ((fp = tpplib_popen(szCommand, "r")) != NULL)
      {
         int iGood=0;
         int iBad=0;
   
         while (fgets(szBuf, SIZE_BUF, fp))
         {
            int  bPass;
   
            szBuf[strlen(szBuf) - 1] = '\0';
   
            bPass = ANALYZE(szBuf, pOptions);

            if (pOptions.bVerboseOutput)
            {
               if (bPass)
                  printf("  PASS\n");
               else
                  printf("  FAIL\n");
            }
   
            if (bPass==1)
            {
               iGood++;
            }
            else
            {
               iBad++;

               if (pOptions.bDeleteBad)
                  unlink(szBuf);
               else
               {
                  sprintf(szCommand, "%s.bad", szBuf);
                  rename(szBuf, szCommand);
               }
            }
         }
   
         if (pOptions.bVerboseOutput)
            printf("Summary:  %d good, %d bad\n", iGood, iBad);
   
         pclose(fp);
      }
   }
   else
   {
      int iGood=0;
      int iBad=0;

      for (i=iStartArgc; i<argc; i++)
      {
         int bPass;
   
         bPass = ANALYZE(argv[i], pOptions);

         if (pOptions.bVerboseOutput)
         {
            if (bPass)
               printf("  PASS\n");
            else
               printf("  FAIL\n");
         }

   
         if (bPass==1)
         {
            iGood++;
         }
         else
         {
            iBad++;

            if (pOptions.bDeleteBad)
               unlink(argv[i]);
            else
            {
               sprintf(szCommand, "%s.bad", argv[i]);
               rename(argv[i], szCommand);
            }
         }
      }

      if (pOptions.bVerboseOutput)
         printf("Summary:  %d good, %d bad\n", iGood, iBad);
   }

   return (0);
} /*main*/


void SET_OPTION(char *arg,
                struct OptionsStruct *pOptions)
{
   int  iTmp = 0;
   double dTmp = 0.0;

   switch (arg[1])
   {
      case 'F':
         if (sscanf(&arg[2], "%d", &iTmp) != 1 || iTmp < 1)
            printf("Bad first scan #: '%s' ... ignored\n", &arg[2]);
         else
            pOptions->iFirstScan = iTmp;
         break;
      case 'L':
         if (sscanf(&arg[2], "%d", &iTmp) != 1 || iTmp < 1)
            printf("Bad last scan #: '%s' ... ignored\n", &arg[2]);
         else
            pOptions->iLastScan = iTmp;
         break;
      case 'B':
         if (sscanf(&arg[2], "%lf", &dTmp) != 1 || dTmp < 0.00)
            printf("Bad minimum mass: '%s' ... ignored\n", &arg[2]);
         else
            pOptions->dMinMass = dTmp;
         break;
      case 'T':
         if (sscanf(&arg[2], "%lf", &dTmp) != 1 || dTmp < 0.000)
            printf("Bad maximum mass: '%s' ... ignored\n", &arg[2]);
         else
            pOptions->dMaxMass = dTmp;
         break;
      case 'P':
         if (sscanf(&arg[2], "%d", &iTmp) != 1 || iTmp < 0)
            printf("Bad minimum peak count: '%s' ... ignored\n", &arg[2]);
         else
            pOptions->iMinPeakCount = iTmp;
         break;
      case 'I':
         if (sscanf(&arg[2], "%d", &iTmp) != 1 || iTmp < 0)
            printf("Bad minimum peak threshhold: '%s' ... ignored\n", &arg[2]);
         else
            pOptions->dMinPeakThreshHold = iTmp;
         break;
      case 'i':
         if (sscanf(&arg[2], "%d", &iTmp) != 1 || iTmp < 0)
            printf("Bad minimum peak intensity for filtering: '%s' ... ignored\n", &arg[2]);
         else
            pOptions->dMinPeakDtaFilter = iTmp;
         break;
      case 'A':
         pOptions->bAnalyzeAll = 1;
         break;
      case 'N':
         pOptions->bDeleteBad = 0;
         break;
      case 'V':
         pOptions->bVerboseOutput = 1;
         break;
      default:
         break;
   }

   arg[0] = '\0';
}


/*
 * return TRUE if spectrum passes filters, FALSE if spectrum fails
 */
int ANALYZE(char *szDtaFile,
            struct OptionsStruct pOptions)
{
   FILE *fp;
   int  iStartScan = 0;
   int  iEndScan = 0;
   char *pStr;
   char szBuf[SIZE_BUF];
   char szNewFile[SIZE_BUF];
   
   strcpy(szBuf, szDtaFile);

   if (pOptions.bVerboseOutput)
      printf("%s", szBuf);

   pStr=strrchr(szBuf, '.');
   *pStr='\0'; /* crop of .dta */

   pStr=strrchr(szBuf, '.');
   *pStr='\0'; /* crop of charge state */

   pStr=strrchr(szBuf, '.');
   sscanf(pStr+1, "%d", &iStartScan);   /* get start scan */
   *pStr='\0'; /* crop of last scan */

   pStr=strrchr(szBuf, '.');
   sscanf(pStr+1, "%d", &iEndScan);     /* get end scan */


   if (iStartScan<pOptions.iFirstScan)
      return (0);
   if (iEndScan>pOptions.iLastScan)
      return (0);

   if ((fp = fopen(szDtaFile, "r")) != NULL)
   {
      int  i;
      int  iCharge = 0;
      int  iMaxMass=0;
      int  iPeakCount=0;
      int  bFilter=TRUE;
      double dMass = 0.0;
      double dPeakArray[MAXARRAY];
      FILE *fp2;

      if (pOptions.dMinPeakDtaFilter > 0.0)
      {
         sprintf(szNewFile, "%s.new", szDtaFile);
         if ((fp2 = fopen(szNewFile, "w")) == NULL)
            bFilter=FALSE;
      }

      char *fgot=fgets(szBuf, SIZE_BUF, fp);

      if (pOptions.dMinPeakDtaFilter > 0.0 && bFilter)
         fprintf(fp2, "%s", szBuf);

      sscanf(szBuf, "%lf %d\n", &dMass, &iCharge);

      if (pOptions.bVerboseOutput)
         printf("  %6.1f  %d", dMass, iCharge);

      if (pOptions.dMinMass > dMass || dMass > pOptions.dMaxMass)
      {
         fclose(fp);
         if (pOptions.dMinPeakDtaFilter > 0.0)
         {
            fclose(fp2);
            unlink(szNewFile);
         }
         return (0);
      }

      if (pOptions.iCharge != 0 && iCharge != pOptions.iCharge)
      {
         fclose(fp);
         if (pOptions.dMinPeakDtaFilter > 0.0)
         {
            fclose(fp2);
            unlink(szNewFile);
         }
         return (0);
      }

      memset(dPeakArray, 0, sizeof(dPeakArray));

      while (fgets(szBuf, SIZE_BUF, fp))
      {
         double dMass;
         double dInten;

         if (pOptions.dMinPeakDtaFilter > 0.0 && bFilter)
            fprintf(fp2, "%s", szBuf);

         sscanf(szBuf, "%lf %lf", &dMass, &dInten);

         if (dInten > pOptions.dMinPeakDtaFilter && (int)(dMass+0.5) < MAXARRAY-1)
         {
            dPeakArray[(int)(dMass+0.5)] += dInten;
            if ((int)(dMass+0.5) > iMaxMass)
               iMaxMass=(int)(dMass+0.5);
         }
      }

      for (i=0; i<iMaxMass; i++)
      {
         if (dPeakArray[i] >= pOptions.dMinPeakThreshHold)
            iPeakCount++;
      }

      fclose(fp);
      if (pOptions.dMinPeakDtaFilter > 0.0)
         fclose(fp2);

      if (iPeakCount < pOptions.iMinPeakCount)
      {
         if (pOptions.dMinPeakDtaFilter > 0.0)
            unlink(szNewFile);
         return(0);
      }
      else if (pOptions.dMinPeakDtaFilter > 0.0)
      {
         char szCommand[512];
         sprintf(szCommand, "mv %s %s", szNewFile, szDtaFile);
         verified_system(szCommand);
      }

   }
   else if (pOptions.bVerboseOutput)
      printf("... can't open %s\n", szDtaFile);

   return (1);

} /*ANALYZE*/
