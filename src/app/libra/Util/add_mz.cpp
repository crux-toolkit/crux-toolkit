#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "Common/AminoAcidMasses.h"

void INITIALIZE_MASS(double *dAA);
void ANALYZE(double *dAA,
      char *szInputPep,
      int iCharge,
      FILE *fout);

int main(int argc, char **argv)
{
   int i;
   double dAA[128];

   if (argc<2)
   {
      printf("\n");
      printf(" This program takes in a text file containing peptide sequence\n");
      printf(" and charge state as first 2 columns.  Output appends additional\n");
      printf(" columns of precursor m/z & number of NxS/T motifs.  All Asp in NxS/T\n");
      printf(" motifs are considered glycosylated +0.984016.  All Cys are +57.021464.\n");
      printf(" Met are +15.994915 if indicated modified in the sequence.  Every other\n");
      printf(" modification is ignored!  Peptides must contain preceeding & trailing\n");
      printf(" amino acids.\n");
      printf("\n");
      printf(" Jimmy Eng, 20080520\n");
      printf("\n");
      exit(1);
   }

   for (i=0; i<128; i++)
      dAA[i]=99999.9;

   INITIALIZE_MASS(dAA, 1);

   dAA['C'] += 57.021464;

   printf("\n");
   printf("\n add_mz running ...");
   printf("\n");

   for (i=1; i<argc; i++)
   {
      FILE *fin;

      if ( (fin=fopen(argv[i], "r"))==NULL)
      {
         printf(" Error - cannot read %s, skipping\n", argv[i]);
      }
      else
      {
         FILE *fout;
         char szOut[256];
         char szBuf[4096];

         printf("   input: %s\n", argv[i]);
         sprintf(szOut, "%s.new", argv[i]);
         if ( (fout=fopen(szOut, "w"))==NULL)
         {
            printf(" Error - cannot write %s output\n", szOut);
            exit(1);
         }

         /*
          * get & print header line
          */
	 char* errCheckMe = NULL;
         errCheckMe = fgets(szBuf, 4096, fin);
         while (iscntrl(szBuf[strlen(szBuf)-1]))
               szBuf[strlen(szBuf)-1]='\0';
         fprintf(fout, "%s\tm/z\tNxST_count\n", szBuf);

         while (fgets(szBuf, 4096, fin))
         {
            char szInputPep[256];
            int iCharge;
            int iLen;

            iLen=strlen(szBuf);

            sscanf(szBuf, "%s %d", szInputPep, &iCharge);

            if (szInputPep[1]!='.' || szInputPep[strlen(szInputPep)-2]!='.')
            {
               printf("\n");
               printf(" Error - peptide %s is not in expected format.\n", szInputPep);
               printf(" Peptides must be of format X.XXXXXXX.X\n");
               printf("\n");
               exit(1);
            }

            /*
             * Strip out line termination
             */
            while (iscntrl(szBuf[strlen(szBuf)-1]))
               szBuf[strlen(szBuf)-1]='\0';

            fprintf(fout, "%s", szBuf);

            ANALYZE(dAA, szInputPep, iCharge, fout);
         }
         fclose(fout);
         fclose(fin);

         printf(" created: %s\n\n", szOut);
      }
   }
}


/*
 * reads input pep, notes M & N mods, prints m/z
 */
void ANALYZE(double *dAA,
      char *szInputPep,
      int iCharge,
      FILE *fout)
{
   int i;
   int j;
   int iLen;
   int iNumOxMet=0;
   int iNumNtoD=0;
   char szModPep[256];
   double dPeptideMass;
   double dMZ;

   iLen=strlen(szInputPep);

   /*
    * count # modified M
    */
   for (i=2; i<iLen-2; i++)
   {
      if (szInputPep[i]=='M' && i<iLen-1
            && szInputPep[i+1]=='['
            && szInputPep[i+2]=='1'
            && szInputPep[i+3]=='4'
            && szInputPep[i+4]=='7')
         iNumOxMet++;
   }

   /*
    * strip peptide and count # of modified N via NxS/T motif
    */
   j=0;
   for (i=2; i<iLen; i++)
   {
      if (szInputPep[i]=='.')
      {
      }
      else if (szInputPep[i]!='[')
      {
         szModPep[j]=szInputPep[i];
         j++;
      }
      else
      {
         while (szInputPep[i]!=']')
            i++;
      }
   }
   szModPep[j]='\0';

   /*
    * Now count number of modified N
    */
   iLen=strlen(szModPep);
   for (i=0; i<iLen-2; i++)
   {
      if (szModPep[i]=='N' && (szModPep[i+2]=='S' || szModPep[i+2]=='T'))
         iNumNtoD++;
   }

   /*
    * Sum up mass and calculate m/z; do not consider last/trailing AA
    * that's present to check motif but not part of identified peptide
    */
   dPeptideMass = 0;
   for (i=0; i<iLen-1; i++)
   {
      dPeptideMass += dAA[szModPep[i]];
   }

   dPeptideMass += dAA['o'] + dAA['h'] + dAA['h'] ;  /* neutral mass */
   dPeptideMass += iNumNtoD*0.984016 + iNumOxMet*15.994915;  /* modified neutral mass */

   dMZ = (dPeptideMass + iCharge * 1.00727646686) / iCharge;

   fprintf(fout, "\t%0.6f\t%d\n", dMZ, iNumNtoD);

/*
   printf("Input: %s, Strip: %s, #M: %d, #N: %d, mz %f\n\n",
         szInputPep, szModPep, iNumOxMet, iNumNtoD, dMZ);
*/
}
