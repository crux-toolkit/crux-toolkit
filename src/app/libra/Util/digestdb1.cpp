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
 * Date:    11/05/2002
 * Copyright Jimmy Eng, Institute for Systems Biology
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "Common/AminoAcidMasses.h"

#define MIN_PEPTIDE_MASS     800.0
#define MAX_PEPTIDE_MASS     3000.0
#define ALLOWED_MISSED_CLEAVAGES 1

#define MAX_LEN_PEPTIDE      800
#define MAX_LEN_DEFINITION   100
#define INITIAL_SEQ_LEN      1000
#define SIZE_BUF             500000  // PEFF headers can be HUGE

#define dProtonMass          1.00727649

extern double COMPUTE_PI(char *seq,unsigned long seq_length, int charge_increment);


double dMassAA[128];

struct ParamStruct
{
    int iMissedCleavage;   /* num allowed missed cleavages */
    int iMassType;         /* 0=avg, 1=mono (default) */
   double dMinMass;
   double dMaxMass;
   char szDb[256];
   char szBreak[24];       /* which residues to fragment */
   char szNoBreak[24];
   char szUserMods[SIZE_BUF];
} pInput;

static struct DbStruct
{
   int iLenSeq;
   char szDef[MAX_LEN_DEFINITION];
   char *szSeq;
} pSeq;

struct PepStruct
{
   int iStart;
   int iEnd;
   double dPepMass;
} pPep;


void READ_FASTA(void);
void FIND_TRYPTIC_PEPTIDES();
void SET_OPTION(char *arg);


int main(int argc, char **argv)
{
   int iNumArg;
   int iStartArgc;
   int i;
   char *arg;
 
   pInput.iMissedCleavage = ALLOWED_MISSED_CLEAVAGES; /* Number of missed cleavages */
   pInput.dMinMass = MIN_PEPTIDE_MASS;        /* minimum peptide mass */
   pInput.dMaxMass = MAX_PEPTIDE_MASS;        /* maximum peptide mass */
   pInput.iMassType = 1;
   pInput.szUserMods[0]='\0';
   strcpy(pInput.szBreak, "RK");
   strcpy(pInput.szNoBreak, "P");

   if (argc<2)
   {
      printf("\n");
      printf(" DIGEST PROGRAM, J.Eng\n");
      printf("\n");

      printf(" USAGE:   %s [options] database_file\n\n", argv[0]);

      printf(" Command line options:\n");
      printf("    -l<num>     set minimum peptide mass (<num> is a float; default=%0.2lf)\n", pInput.dMinMass);
      printf("    -h<num>     set maximum peptide mass (<num> is a float; default=%0.2lf)\n", pInput.dMaxMass);
      printf("    -m<num>     set number of missed cleavages (<num> is an int; default=%d)\n\n", pInput.iMissedCleavage);

      printf("    -r<str>     residues to cleave after (default=\"%s\" for trypsin)\n", pInput.szBreak);
      printf("    -n<str>     don't cleave if next AA (default=\"%s\" for trypsin)\n", pInput.szNoBreak);
      printf("                ** only c-term cleavages are support right now so there's no AspN digestion.\n");
      printf("                ** use a dash (-) or leave <str> blank for a null character.\n");
      printf("    -M<str>     static modifications, comma separated of form <mass>@<residue>\n");
      printf("                   for example -M57.021@C,15.995@M\n\n");

      printf("    -t<num>     mass type, 0=average, 1=monoisotopic (default)\n\n");

      printf(" For example:  %s ipi.fasta\n", argv[0]);
      printf("               %s -m3 ipi.fasta       (allow up to 3 missed cleavages)\n", argv[0]);
      printf("               %s -l800.0 -h1000.0 ipi.fasta    (set mass range from 800.0 to 1000.0)\n", argv[0]);
      printf("               %s -rDE -nP ipi.fasta  (glu-C protease cuts after D/E but no if followed by P)\n\n", argv[0]);

      printf(" Results are spit out to stdout in a tab-delimited format.\n");
      printf(" Redirect the output if you want them stored in a file.\n\n");
      printf("     For example:  %s ipi.fasta > digest.txt\n\n", argv[0]);

      printf(" To only get keep certain columns, use the 'awk' command.  For example, to print\n");
      printf(" only the 2nd (protein), 3rd (mass), and 5th (peptide) columns, type\n");
      printf("    %s ipi.fasta | awk '{print $2 \"\\t\" $3 \"\\t\" $4}' > digest.txt\n\n", argv[0]);
      printf(" Asterisks (*) in sequence are treated as proper break points\n\n");

      printf(" Output columns:\n");
      printf(" - peptide length\n");
      printf(" - protein reference\n");
      printf(" - peptide mass\n");
      printf(" - previous amino acid before peptide\n");
      printf(" - peptide sequence\n");
      printf(" - next amino acid after peptide\n");
      printf(" - peptide start location\n");
      printf(" - peptide end location\n");
      printf(" - pI\n\n");
      exit(1);
   }

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

   if (pInput.dMaxMass <= pInput.dMinMass)
   {
      printf(" error - mass range incorrect:  min=%0.2f  max=%0.2f\n\n", pInput.dMinMass, pInput.dMaxMass);
      exit(1);
   }

   if (iStartArgc == argc)
   {
      printf(" Please enter original database file on the command line.\n\n");
      exit(EXIT_FAILURE);
   }
   else
      strcpy(pInput.szDb, argv[iStartArgc]);

   /*
   printf("db=%s, range=%f-%f, missed=%d\n", pInput.szDb, pInput.dMinMass, pInput.dMaxMass, pInput.iMissedCleavage);
   exit(1);
   */

   for (i=0; i<128; i++)
      dMassAA[i]=99999.9;
   INITIALIZE_MASS(dMassAA, pInput.iMassType);

   /*
    * extract static mods if specified
    */
   if (strlen(pInput.szUserMods)>0)
   {
      char *tok;

      tok = strtok(pInput.szUserMods, ",");
      while (tok != NULL)
      {
         char *pStr;
         double dMass;
         char szResidue[100];

         pStr=strchr(tok, '@');
         *pStr = ' ';

         sscanf(tok, "%lf %s", &dMass, szResidue);

         /*
          * <mass>@<residue>; assume only a single character is entered for each residue
          */
         if (szResidue[0]>='A' && szResidue[0]<='Z')
         {
//          printf("adding static mod mass %f to residue %c\n", dMass, szResidue[0]);
            dMassAA[toupper(szResidue[0])] += dMass;
         }

         tok = strtok(NULL, ",");
      }
   }


   /* minimum & maximum peptide masses to print out */
   READ_FASTA();

   return(0);

} /*main*/


void SET_OPTION(char *arg)
{
   int iTmp;
   double dTmp;
   char szTmp[SIZE_BUF];

   switch (arg[1])
   {
      case 'l':
         if (sscanf(&arg[2], "%lf", &dTmp) != 1)
            printf("Bad #: '%s' ... ignored\n", &arg[2]);
         else
         {
            if (dTmp >= 0.0)
               pInput.dMinMass = dTmp;
            else
               printf(" error: negative minimum mass specified (%0.2f).\n", dTmp);
         }
         break;
      case 'h':
         if (sscanf(&arg[2], "%lf", &dTmp) != 1)
            printf("Bad #: '%s' ... ignored\n", &arg[2]);
         else
         {
            if (dTmp >= 0.0)
               pInput.dMaxMass = dTmp;
            else
               printf(" error: negative maximum mass specified (%0.2f).\n", dTmp);
         }
         break;
      case 'm':
         if (sscanf(&arg[2], "%d", &iTmp) != 1)
            printf("Bad #: '%s' ... ignored\n", &arg[2]);
         else
            pInput.iMissedCleavage= iTmp;
         break;
      case 'r':
         sscanf(&arg[2], "%s", szTmp);
         strncpy(pInput.szBreak, szTmp, 24);
         break;
      case 't':
         if (sscanf(&arg[2], "%d", &iTmp) != 1)
            printf("Bad #: '%s' ... ignored\n", &arg[2]);
         else
         {
            if (iTmp<0 || iTmp>1)
               iTmp=1;
            pInput.iMassType= iTmp;
         }
         break;
      case 'n':
         sscanf(&arg[2], "%s", szTmp);
         strncpy(pInput.szNoBreak, szTmp, 24);
         break;
      case 'M':
         sscanf(&arg[2], "%s", szTmp);
         strcpy(pInput.szUserMods, szTmp);
         break;
      default:
         break;
   }
   arg[0] = '\0';

} /*SET_OPTION*/


void READ_FASTA(void)
{
   int  iLenAllocSeq;
   char szBuf[SIZE_BUF];
   FILE *fp;
   
   iLenAllocSeq = INITIAL_SEQ_LEN;

   pSeq.szSeq = (char*) malloc(iLenAllocSeq);
   if (pSeq.szSeq == NULL)
   {
      printf(" Error malloc pSeq.szSeq (%d)\n\n", iLenAllocSeq);
      exit(0);
   }

   if ( (fp=fopen(pInput.szDb, "r"))==NULL)
   {
      printf(" Error - cannot database file %s\n\n", pInput.szDb);
      exit(EXIT_FAILURE);
   }

   while (fgets(szBuf, SIZE_BUF, fp))
   {
      if (szBuf[0]=='>')
      {
         int cResidue;

         strncpy(pSeq.szDef, szBuf+1, MAX_LEN_DEFINITION-1);
         pSeq.szDef[MAX_LEN_DEFINITION-1]='\0';
         pSeq.iLenSeq=0;

         while ((cResidue=fgetc(fp)))
         {
            if (isalpha(cResidue) || cResidue=='*')
            {
               pSeq.szSeq[pSeq.iLenSeq]=cResidue;
               pSeq.iLenSeq++;
               if (pSeq.iLenSeq == iLenAllocSeq-1)
               {
                  char *pTmp;

                  iLenAllocSeq += 500;
                  pTmp = (char*) realloc(pSeq.szSeq, iLenAllocSeq);
                  if (pTmp == NULL)
                  {
                     printf(" Error realloc pSeq.szSeq (%d)\n\n", iLenAllocSeq);
                     exit(1);
                  }

                  pSeq.szSeq = pTmp;
               }
            }
            else if (feof(fp) || cResidue=='>')
            {
               int iReturn;

               iReturn=ungetc(cResidue, fp);

               if (iReturn!=cResidue)
               {
                  printf("Error with ungetc.\n\n");
                  exit(EXIT_FAILURE);
               }
               break;
            }
         }

         pSeq.szSeq[pSeq.iLenSeq]='\0';
         FIND_TRYPTIC_PEPTIDES();
      }
   }

   fclose(fp);
    
} /*READ_FASTA*/


void FIND_TRYPTIC_PEPTIDES()
{
   int  iMissed,
        iStartNextPeptide;

   pPep.iStart=0;
   pPep.iEnd=0;
   iMissed=0;
   iStartNextPeptide=0;

   do
   {
      if (( (strchr(pInput.szBreak, pSeq.szSeq[pPep.iEnd])
              && !strchr(pInput.szNoBreak, pSeq.szSeq[pPep.iEnd+1]))
            || pSeq.szSeq[pPep.iEnd]=='*')
         || pPep.iEnd==pSeq.iLenSeq-1)
      {
         int i;
 
         if (pSeq.szSeq[pPep.iEnd]=='*')
            pPep.iEnd--;

         pPep.dPepMass = dMassAA['o']+ dMassAA['h'] + dMassAA['h'] + dProtonMass;

         for (i=pPep.iStart; i<=pPep.iEnd; i++)
         {
            // should array subscript really be int, or is this a typo?
             pPep.dPepMass += dMassAA[(int)pSeq.szSeq[i]]; 
         }

         /*
          * These peptides pass mass filter
          */
         if (    pPep.dPepMass>=pInput.dMinMass
              && pPep.dPepMass<=pInput.dMaxMass
              && pPep.iEnd-pPep.iStart+2<MAX_LEN_PEPTIDE)
         {
            int iTmp;
            char szPeptide[MAX_LEN_PEPTIDE],
                 szDefWord[MAX_LEN_DEFINITION];
            double dPI;

            iTmp=0;
            for (i=pPep.iStart; i<=pPep.iEnd; i++)
            {
               szPeptide[iTmp]=pSeq.szSeq[i];
               iTmp++;
            }
            szPeptide[iTmp]='\0';

            sscanf(pSeq.szDef, "%s", szDefWord);
               
            /* print header & MH+ mass */
            printf("%d\t%s\t%11.6f\t", (int)strlen(szPeptide), (char*)szDefWord, pPep.dPepMass);
   
            /* print previous amino acid */
            if (pPep.iStart-1 >= 0)
               printf("%c\t", pSeq.szSeq[pPep.iStart-1]);
            else
               printf("-\t");
   
            /* print peptide */
            printf("%s", szPeptide);
   
            /* print following amino acid */
            if (pPep.iEnd+1 < pSeq.iLenSeq)
               printf("\t%c", pSeq.szSeq[pPep.iEnd+1]);
            else
               printf("\t-");

            dPI = COMPUTE_PI(szPeptide, strlen(szPeptide), 0);
   
            printf("\t%d\t%d\t%0.2f", pPep.iStart, pPep.iEnd, dPI);
            printf("\n");
         }

         if (pSeq.szSeq[pPep.iEnd+1] == '*')
         {
            iMissed=0;

            if (pPep.iStart==iStartNextPeptide)
            {
               pPep.iStart=pPep.iEnd+2;
               iStartNextPeptide=pPep.iEnd+2;
               pPep.iEnd=pPep.iStart;
            }
            else
            {
               pPep.iStart=iStartNextPeptide;
               pPep.iEnd=pPep.iStart;
            }
         }
         else
         {
            iMissed++;
            if (iMissed==1)  /* first break point is start of next peptide */
               iStartNextPeptide=pPep.iEnd+1;

            if (iMissed <= pInput.iMissedCleavage
                  && pPep.dPepMass<pInput.dMaxMass
                  && pPep.iEnd < pSeq.iLenSeq-1)
            {
               pPep.iEnd++;
            }
            else
            {
               iMissed=0;
               pPep.iStart=iStartNextPeptide;
               pPep.iEnd=pPep.iStart;
            }
         }

      }
      else
      {
         pPep.iEnd++;
      }

   } while (pPep.iStart<pSeq.iLenSeq);

} /*FIND_TRYPTIC_PEPTIDES*/
