/*
    Program: translateDNA2AA-FASTA
    Description: translates nucletoide to AA database in 6 reading frames

    Copyright (C) Jimmy Eng, ISB Seattle

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
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


#define SIZE_DESCR     4096




struct  pDB        /* pDB stores temporary name & sequence info. */
{               /* for each protein in the database to be analyzed */
   int iLen;
   int iMaxLenAllocated;
  char szProt[SIZE_DESCR];
  char *szSeq;
} pDB;

int iAllocatedProtSeqLength,
    iProteinSeqLength;

char *szTranslatedSeq;


void LOAD_DBASE(char *szBuf,
        FILE *fpIn);
void TRANSLATE(int *frame,
        int iForward);
char AA_ENCODE(int i,
        int iForward);


int main (int argc, char **argv)
{
   int i;
   char szBuf[SIZE_DESCR],
        szFileOut[256];
   FILE *fpIn,
        *fpOut;
   
   printf("\n");
   printf(" TRANSLATE\n\n");

   if (argc<2)
   {
      printf(" This program takes in a FASTA formatted nucleotide database as\n");
      printf(" input and creates an output database of translated amino acid\n");
      printf(" sequences.  Sequences are translated in all six reading frames.\n");
      printf(" Stop codons are encoded as '*'.\n\n");

      printf(" Please enter nucleotide db on command line.\n\n");
      printf("    for example:  %s mydb.seq\n\n", argv[0]);
      printf(" Output file has '.new' appended i.e. mydb.seq.new\n\n");
      exit(0);
   }

   if ( (fpIn=fopen(argv[1], "r"))==NULL)
   {
      printf(" Error - cannot open input file %s to read.\n\n", argv[1]);
      exit(0);
   }

   strcpy(szFileOut, argv[1]);
   strcat(szFileOut, ".new");

   if ( (fpOut=fopen(szFileOut, "w"))==NULL)
   {
      printf(" Error - cannot open output file %s to write.\n\n", szFileOut);
      exit(0);
   }

   printf("   Input File:  %s\n", argv[1]);

#define INIT_SEQ_LEN   60000
#define INIT_PROT_LEN  20000

   /*
    * Allocate memory to store database sequence
    */
   pDB.szSeq=(char *)malloc(INIT_SEQ_LEN);
   if (pDB.szSeq== NULL)
   {
      fprintf(stderr, " Error malloc(pDB.szSeq)\n");
      exit(EXIT_FAILURE);
   }
   pDB.iMaxLenAllocated=INIT_SEQ_LEN;

  /*
   * Allocate memory to store translated nucleotide sequence 
   */
   szTranslatedSeq=(char *)malloc(INIT_PROT_LEN);
   if (szTranslatedSeq==NULL)
   {
      fprintf(stderr, " Error malloc(szTranslatedSeq).\n");
      exit(EXIT_FAILURE);
   }
   iAllocatedProtSeqLength=INIT_PROT_LEN;

  /*
   * Load database entry header
   */
   char* dummy = fgets(szBuf, SIZE_DESCR, fpIn);

   /*
    * Loop through entire database
    */
   while(!feof(fpIn))
   {
      char *szTemp;
      int  j;

      /*
       * Expect a '>' for sequence header line
       */
      if (szBuf[0]!='>')
      {
         fprintf(stderr, "\n\n Error in database file.\n");
         exit(EXIT_FAILURE);  /* print out error message and exit */
      }
      memcpy(pDB.szProt, szBuf+1, SIZE_DESCR);
      pDB.szProt[strlen(pDB.szProt)-1]='\0';


      LOAD_DBASE(szBuf, fpIn);

      for (j=0;j<3;j++)
      {
         int k,
             iLen;
         TRANSLATE(&j, 1);
         fprintf(fpOut, ">F+%d_%s\n", j+1, pDB.szProt);

         iLen=strlen(szTranslatedSeq);
         for (k=0; k<iLen; k++)
         {
            fprintf(fpOut, "%c", szTranslatedSeq[k]);
            if ( !((k+1)%50))
               fprintf(fpOut, "\n");
         }
         if ((k)%50)
            fprintf(fpOut, "\n");
      }

     /*
      * Generate complimentary strand
      */
      szTemp=(char *)malloc(pDB.iLen);
      if (szTemp == NULL)
      {
         fprintf(stderr, " Error malloc(szTemp) (length sequence=%d)\n", pDB.iLen);
         exit(EXIT_FAILURE);
      }
 
      memcpy(szTemp, pDB.szSeq, pDB.iLen);
 
      for (j=0;j<pDB.iLen;j++)
      {
         if (szTemp[j]=='G')
            pDB.szSeq[j]='C';
         else if (szTemp[j]=='C')
            pDB.szSeq[j]='G';
         else if (szTemp[j]=='T')
            pDB.szSeq[j]='A';
         else if (szTemp[j]=='A')
            pDB.szSeq[j]='T';
         else
            pDB.szSeq[j]=szTemp[j];
      }
      free(szTemp);

      for (j=0;j<3;j++)
      {
         int k,
             iLen;

         TRANSLATE(&j, -1);
         fprintf(fpOut, ">F-%d_%s\n", j+1, pDB.szProt);

         iLen=strlen(szTranslatedSeq);
         for (k=0; k<iLen; k++)
         {
            fprintf(fpOut, "%c", szTranslatedSeq[k]);
            if ( !((k+1)%50))
               fprintf(fpOut, "\n");
         }
         if ((k)%50)
            fprintf(fpOut, "\n");
      } 

   } /*while*/

   free(pDB.szSeq);
   free(szTranslatedSeq);
   fclose(fpIn);
   fclose(fpOut);

   printf(" File created:  %s\n\n", szFileOut);

   exit(EXIT_SUCCESS);
} /*main*/


/*
 * load database entry
 */
void LOAD_DBASE(char *szBuf,
        FILE *fpIn)
{
   int iAACount=0,
       iTmpCh;

   while ( (iTmpCh=getc(fpIn))!='>' && iTmpCh!=EOF )
   {
      if (isalnum(iTmpCh))
      {
         pDB.szSeq[iAACount]=toupper(iTmpCh);

         iAACount++;

         if (iAACount>=pDB.iMaxLenAllocated)
         {
            char *pTmp;

            pTmp = (char *)realloc(pDB.szSeq, iAACount+100);

            if (pTmp == NULL)
            {
               fprintf(stderr, " Error in realloc pDB.szSeq (%d)\n\n", iAACount);
               exit(EXIT_FAILURE);
            }

            pDB.szSeq = pTmp;
            pDB.iMaxLenAllocated=iAACount+99;
         } 
      } 
   }

   pDB.iLen = iAACount;
   pDB.szSeq[iAACount] = '\0';

   if (iTmpCh!=EOF)
   {
      ungetc(iTmpCh, fpIn);
      char* dummy = fgets(szBuf, SIZE_DESCR, fpIn);
   }

} /*LOAD_DBASE*/


/*
 * For nucleotide search, translate from DNA to amino acid
 */
void TRANSLATE(int *frame,
        int iForward)
{
   int i,
       j=0;
                   
   if (iForward== 1)
   { 
      /*
       * forward reading frame
       */
      i=(*frame);
      while ((i+2)<pDB.iLen)
      {
         if (j>=iAllocatedProtSeqLength)
         {
            char *pTmp;

            pTmp=(char *)realloc(szTranslatedSeq, j+100);
            if (pTmp == NULL)
            {
               fprintf(stderr, " Error realloc(szTranslatedSeq) ... size=%d\n", j);
               exit(EXIT_FAILURE);
            }

            szTranslatedSeq=pTmp;
            iAllocatedProtSeqLength=j+99;
         }

         *(szTranslatedSeq+j) = AA_ENCODE(i, 1);
         i += 3;
         j++;
      }
      iProteinSeqLength = j;
      szTranslatedSeq[j]='\0';
   }
   else
   {
      /*
       * reverse reading frame
       */
      i=pDB.iLen-(*frame)-1;
      while (i>=2)   /* 2,1,0 makes the last AA */
      {
         if (j>=iAllocatedProtSeqLength)
         {
            char *pTmp;

            pTmp=(char *)realloc(szTranslatedSeq, j+100);
            if (pTmp == NULL)
            {
               fprintf(stderr, " Error realloc(szTranslatedSeq) ... size=%d\n", j);
               exit(EXIT_FAILURE);
            }

            szTranslatedSeq=pTmp;
            iAllocatedProtSeqLength=j+99;
         }

         *(szTranslatedSeq+j) = AA_ENCODE(i, -1);
         i -= 3;
         j++;
      }
      iProteinSeqLength = j;
      szTranslatedSeq[j]='\0';
   }

} /*TRANSLATE*/



char AA_ENCODE(int i,
        int iForward)
{
   int iBase2=i+iForward,     /* i=position of 1st base, iBase2=2nd, iBase3=3rd) */
       iBase3=i+iForward*2;

   if (pDB.szSeq[i]=='G')    /* First base G */
   {
      if (pDB.szSeq[iBase2]=='T')
         return ('V');
      else if (pDB.szSeq[iBase2]=='C')
         return ('A');
      else if (pDB.szSeq[iBase2]=='G')
         return ('G');
      else if (pDB.szSeq[iBase2]=='A')
      {
         if (pDB.szSeq[iBase3]=='T' || pDB.szSeq[iBase3]=='C')
            return ('D');
         else if (pDB.szSeq[iBase3]=='A' || pDB.szSeq[iBase3]=='G')
            return ('E');
      }
   }

   else if (pDB.szSeq[i]=='C')    /* First base C */
   {
      if (pDB.szSeq[iBase2]=='T')
         return ('L');
      else if (pDB.szSeq[iBase2]=='C')
         return ('P');
      else if (pDB.szSeq[iBase2]=='G')
         return ('R');
      else if (pDB.szSeq[iBase2]=='A')
      {
         if (pDB.szSeq[iBase3]=='T' || pDB.szSeq[iBase3]=='C')
            return ('H');
         else if (pDB.szSeq[iBase3]=='A' || pDB.szSeq[iBase3]=='G')
            return ('Q');
      }
   } /*else*/

   else if (pDB.szSeq[i]=='T')   /* First base T */
   {
      if (pDB.szSeq[iBase2]=='C')
         return ('S');
      else if (pDB.szSeq[iBase2]=='T')
      {
         if (pDB.szSeq[iBase3]=='T' || pDB.szSeq[iBase3]=='C')
            return ('F');
         else if (pDB.szSeq[iBase3]=='A' || pDB.szSeq[iBase3]=='G')
            return ('L');
      }
      else if (pDB.szSeq[iBase2]=='A')
      {
         if (pDB.szSeq[iBase3]=='T' || pDB.szSeq[iBase3]=='C')
            return ('Y');
         else if (pDB.szSeq[iBase3]=='A' || pDB.szSeq[iBase3]=='G')
            return ('*');
      }
      else if (pDB.szSeq[iBase2]=='G')
      {
         if (pDB.szSeq[iBase3]=='T' || pDB.szSeq[iBase3]=='C')
            return ('C');
         else if (pDB.szSeq[iBase3]=='A')
            return ('*');
         else if (pDB.szSeq[iBase3]=='G')
            return ('W');
      }
   } /*elseif*/

   else if (pDB.szSeq[i]=='A')    /* First base A */
   {
      if (pDB.szSeq[iBase2]=='C')
         return ('T');
      else if (pDB.szSeq[iBase2]=='T')
      {
         if (pDB.szSeq[iBase3]=='T' || pDB.szSeq[iBase3]=='C' || pDB.szSeq[iBase3]=='A')
            return ('I');
         else if (pDB.szSeq[iBase3]=='G')
            return ('M');
      }
      else if (pDB.szSeq[iBase2]=='A')
      {
         if (pDB.szSeq[iBase3]=='T' || pDB.szSeq[iBase3]=='C')
            return ('N');
         else if (pDB.szSeq[iBase3]=='A' || pDB.szSeq[iBase3]=='G')
            return ('K');
      }
      else if (pDB.szSeq[iBase2]=='G')
      {
         if (pDB.szSeq[iBase3]=='T' || pDB.szSeq[iBase3]=='C')
            return ('S');
         else if (pDB.szSeq[iBase3]=='A' || pDB.szSeq[iBase3]=='G')
            return ('R');
      }
   } /*elseif*/
   return ('*');   /* Return * if matches nothing */
} /*AA_ENCODE*/
