/*
 * Date:    2004/12/08
 * Author:  Jimmy Eng, FHCRC
 * Purpose: validates that a sequence database doesn't have duplicate headers as
 *          expected by our db search program.  Specifically check that the first
 *          word, up to N number of characters or first white space, is unique.
 *          Also reports if any header lines are longer than MAX_DESCR_LEN.
 *
 * Requires SHA1.cpp and SHA1.h; a copy of these can be found on Sashimi CVS
 * under the ReAdW directory.
 *
 * Compile:  g++ -O3 -o checkdb checkdb1.c SHA1.cpp
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define INIT_HDR_LEN         1000
#define INIT_HDRSTRUCT_LEN   25000
#define SIZE_BUF             4096
#define SIZE_FILE            1024
#define MAX_HEADER_LEN       50    /* maximum length of first word to compare */
#define MAX_DESCR_LEN        8000

#define TITLE               "CHECKDB"

#include "SHA1_TPP.h"

SHA1        m_Csha1;

struct HdrStruct
{
   char szHdrWord[MAX_HEADER_LEN];
   long liLenHdr;
};


void SET_OPTION(char *arg,
        int  *bVerbose,
        int  *bCreateCheckFile);

/*
 * grabs all database headers and makes sure there are no duplicates
 */
void VALIDATE_DB(char *szInput,
        int bVerbose,
        int bCreateCheckFile);

/*
 * if checksum file exists, compares input file checksum to the checksum file
 */
void VALIDATE_CHECKSUM(char *szInput,
        int bVerbose);

int CHECK_FILE_EXISTS(char *szFile); /* return 0 or 1 pending readability of .check file */

int iSortHdr(const void *p0,
      const void *p1);


int main(int argc, char **argv)
{
   int  i;
   int  iNumArg;
   int  iStartArgc;
   int  bVerbose;          /* spits out verbose info */
   int  bCreateCheckFile;      /* run in automated mode; will create .check files */
   char szInput[SIZE_FILE];
   char *arg;

   printf("\n");
   printf(" %s\n", TITLE);
   printf("\n");

   if (argc<2)
   {
      printf(" Please enter one or more fasta files as input on the command line.");
      printf("\n");
      printf(" USAGE:  %s [options] fasta.database\n", argv[0]);
      printf("\n");
      printf(" Options:   -V         verbose reporting\n");
      printf("            -N         do not write .check (SHA1 checksum) file\n");
      printf("\n");
      printf(" This program checks for duplicate accession numbers (actually just the\n");
      printf(" first word in the definition lines of a fasta file).  If there are no\n");
      printf(" duplicates, a checksum is generated and written to a .check file.\n");
      printf(" If a .check file exists, it compares that checksum to the input database\n");
      printf(" to confirm the fasta file has not changed.  Long descriptions are also\n");
      printf(" flagged for any headers greater than %d characters.\n", MAX_DESCR_LEN);
      printf("\n");
      exit(EXIT_FAILURE);
   }

   iNumArg=0;
   arg = argv[iNumArg = 1];
   iStartArgc=1;
   bVerbose=0;
   bCreateCheckFile=1;
   szInput[0]='\0';

   /*
    * Get command line arguments
    */
   while (iNumArg < argc)
   {
      if (arg[0] == '-')
         SET_OPTION(arg, &bVerbose, &bCreateCheckFile);
      else
         break;

      iStartArgc++;
      arg = argv[++iNumArg];
   }

   if (iStartArgc == argc)
   {
      printf(" Please enter database file on the command line.\n");
      printf(" For example:\n");
      printf("    %s human.IPI.fasta\n\n, ", argv[0]);
      exit(EXIT_FAILURE);
   }
   else
      strcpy(szInput, argv[iStartArgc]);


   /*
    * Loop through input files
    */
   for (i=iStartArgc; i<argc; i++)
   {
      printf(" File: %s  ", argv[i]);
      fflush(stdout);

      if (CHECK_FILE_EXISTS(argv[i]))
      {
         /*
          * if check file exists, validate checksum
          */
         VALIDATE_CHECKSUM(argv[i], bVerbose);
      }
      else
      {
         /*
          * validate database
          */
         VALIDATE_DB(argv[i], bVerbose, bCreateCheckFile);
      }

      if (bVerbose)
         printf("\n");
   }

   printf("\n Done.\n\n");
   return(EXIT_SUCCESS);

} /*main*/


void SET_OPTION(char *arg,
      int  *bVerbose,
      int  *bCreateCheckFile)
{

   switch (arg[1])
   {
      case 'V':
         *bVerbose = 1;
         break;
      case 'N':
         *bCreateCheckFile = 0;
         break;
      default:
         break;
   }
   arg[0] = '\0';

} /*SET_OPTION*/


/*
 * grabs all database headers and makes sure there are no duplicates
 */
void VALIDATE_DB(char *szInput,
      int bVerbose,
      int bCreateCheckFile)
{
   int  iNumDuplicates;

   long liLenHdrAllocated,
        liLenHdr,
        liNumHdrStruct,
        liNumHdrStructAllocated,
        liMaxLength=0;
   long li;

   FILE *fp;

   char szBuf[SIZE_BUF],
        szSHA1[1024],
        *pHeader,
        szHdrWord[MAX_HEADER_LEN];  /* first word */

   struct HdrStruct *pHdrStruct;

   if ( (fp=fopen(szInput, "r"))==NULL)
   {
      printf(" Read_Error\n");
   }
   else
   {
      int cAA;

      cAA=getc(fp);

      if (cAA == '>')  /* expect first character in the fasta file to be '>' */
      {
         int bSuccess;

         ungetc(cAA, fp);

         if ( (pHeader=(char *)malloc(INIT_HDR_LEN*sizeof(char)))==NULL)
         {
            printf(" Error cannot malloc pHeader\n\n");
            exit(EXIT_FAILURE);
         }
      
         liLenHdrAllocated=INIT_HDR_LEN;
      
         if ( (pHdrStruct = (struct HdrStruct *)malloc(INIT_HDRSTRUCT_LEN*sizeof(struct HdrStruct)))==NULL)
         {
            printf(" Error cannot malloc pHdrStruct\n\n");
            exit(EXIT_FAILURE);
         }
      
         liNumHdrStructAllocated=INIT_HDRSTRUCT_LEN;
         liNumHdrStruct=0;
      
         fflush(stdout);
      
         /*
          * read through database and grab headers
          */
         while (1)
         {
            cAA=getc(fp);

            if (cAA == '>')
            {
               liLenHdr=0;

               while ((cAA = fgetc(fp)))  /* want to grab whole line in case '>' appears within def line */
               {
                  if (!feof(fp) && cAA!='\n')
                  {
                     pHeader[liLenHdr]=cAA;
                     liLenHdr++;
      
                     if (liLenHdr >= liLenHdrAllocated)
                     {
                        char *pStr;
      
                        liLenHdrAllocated += 500;
                        pStr=(char *)realloc(pHeader, liLenHdrAllocated);
      
                        if (pStr==NULL)
                        {
                           printf(" Error - cannot realloc pHeader[%ld]\n\n", liLenHdrAllocated);
                           exit(EXIT_FAILURE);
                        }
                        pHeader=pStr;
                     }
                  }
                  else
                     break;
               }
      
               /*
                * The description line should begin with an '>' an there should be no
                * space between the '>' and first letter of the identifier.
                */
               bSuccess=1;
               if (pHeader[0]==' ' || !isalnum(pHeader[0]))
               {
                  printf(" Failure: space between '>' and identifier in the description line.\n");
                  printf(" >%s\n", pHeader);
		  char* errCheckMe = NULL;
                  errCheckMe = fgets(szBuf, SIZE_BUF, fp);
                  printf(" %s\n", szBuf);
                  bSuccess = 0;
                  break;
               }

               /* get the first word of each header */
               pHeader[MAX_HEADER_LEN-1]='\0';
               pHeader[liLenHdr]='\0';
               sscanf(pHeader, "%s", szHdrWord);
      
               if (liNumHdrStruct >= liNumHdrStructAllocated)
               {
                  struct HdrStruct *pTmpHdr;

                  liNumHdrStructAllocated += 500;
                  pTmpHdr=(struct HdrStruct *)realloc(pHdrStruct, liNumHdrStructAllocated*sizeof(struct HdrStruct));
      
                  if (pTmpHdr==NULL)
                  {
                     printf(" Error - cannot realloc pHdrStruct[%ld]\n\n", liNumHdrStructAllocated);
                     exit(EXIT_FAILURE);
                  }
                  pHdrStruct=pTmpHdr;
               }
      
               /* store the header word */
               strcpy( (pHdrStruct+liNumHdrStruct)->szHdrWord, szHdrWord);
               (pHdrStruct+liNumHdrStruct)->liLenHdr = liLenHdr;

               liNumHdrStruct += 1;
      
               if (liLenHdr > liMaxLength)
                  liMaxLength = liLenHdr;

            }
            if (feof(fp))
               break;
         }

         free(pHeader);
         fclose(fp);
   
         if (bSuccess)
         {
            /*
             * get SHA1 checksum for input file
             */
            szSHA1[0] = 0;
            m_Csha1.Reset();
            if( !(m_Csha1.HashFile( szInput)) ) {
               printf(" Cannot open file %s for sha-1 calculation\n", szInput);
               exit(-1);     // Cannot open file
            }
            m_Csha1.Final();
            m_Csha1.ReportHash( szSHA1, SHA1::REPORT_HEX );

            /*
             * sort headers
             */
            qsort(pHdrStruct, liNumHdrStruct, sizeof(struct HdrStruct), iSortHdr);

            /*
             * find # of duplicate headers
             */
            iNumDuplicates=0;
            for (li=1; li<liNumHdrStruct; li++)
            {
               if (!strcasecmp( (pHdrStruct+li)->szHdrWord, (pHdrStruct+li-1)->szHdrWord ))
                  iNumDuplicates++;
            }

            bSuccess=0;
            if (iNumDuplicates==0 && liMaxLength < MAX_DESCR_LEN)
            {
               if (bCreateCheckFile)
               {
                  if (strcmp(szSHA1, "sha1-error"))
                  {
                     FILE *fpCheck;
                     char szCheckFile[SIZE_FILE];
      
                     sprintf(szCheckFile, "%s.check", szInput);
      
                     if ( (fpCheck=fopen(szCheckFile, "w"))==NULL)
                     {
                        printf(" Error - cannot write check file %s\n", szCheckFile);
                        exit(EXIT_FAILURE);
                     }
      
                     fprintf(fpCheck, "%s\tsha1 %s\n", szInput, szSHA1);
                     fclose(fpCheck);

                     bSuccess=1;
                     printf(" Success:  file validated\n");
                  }
                  else
                  {
                     printf(" Failure: cannot calculate checksum (%s)\n", szInput);
                  }
               }
               else
               {
                  bSuccess=1;
                  printf(" Success:  file validated, no .check file created\n");
               }
            }
            
            if (bSuccess & bVerbose)
            {
               printf("           MaxHdrLen=%ld, NumHdrs=%ld, Sha1=%s\n",
                     liMaxLength, liNumHdrStruct, szSHA1);
            }
            else if (!bSuccess)
            {
               printf(" Failure: MaxHdrLen=%ld, NumHdrs=%ld, NumDupl=%d\n",
                     liMaxLength, liNumHdrStruct, iNumDuplicates);

            }
      
            /*
             * print the duplicate headers for verbose output
             */
            if (bVerbose)
            {
               if (iNumDuplicates>0)
               {
                  int iCt=0;

                  printf("    Duplicates:\n");
                  for (li=1; li<liNumHdrStruct; li++)
                     if (!strcasecmp( (pHdrStruct+li)->szHdrWord, (pHdrStruct+li-1)->szHdrWord ))
                        printf("    %d.  %s %s\n", ++iCt, (pHdrStruct+li)->szHdrWord, (pHdrStruct+li-1)->szHdrWord);
               }

               if (liMaxLength >= MAX_DESCR_LEN)
               {
                  printf("    Long Headers:\n");
                  for (li=0; li<liNumHdrStruct; li++)
                     if ((pHdrStruct+li)->liLenHdr >= MAX_DESCR_LEN)
                        printf("    %s   HdrLen=%ld\n", (pHdrStruct+li)->szHdrWord, (pHdrStruct+li)->liLenHdr);
               }
            }
         }
      }
      else
      {
         printf(" Failure: not fasta file (%s, first char=%c; expecting '>')\n", szInput, cAA);
         ungetc(cAA, fp);
	 char* errCheckMe=NULL;
         errCheckMe = fgets(szBuf, SIZE_BUF, fp);
         printf(" %s\n", szBuf);
         errCheckMe = fgets(szBuf, SIZE_BUF, fp);
         printf(" %s\n", szBuf);
      }
   }

} /*VALIDATE_DB*/


/*
 * This sorts the sequence headers stored in pHdrStruct in ascending order
 */
int iSortHdr(const void *p0,
      const void *p1)
{
   return (strcasecmp(((struct HdrStruct *) p0)->szHdrWord, ((struct HdrStruct *) p1)->szHdrWord));
}


int CHECK_FILE_EXISTS(char *szFile)
{
   FILE *fp;
   char szCheckFile[SIZE_FILE];

   sprintf(szCheckFile, "%s.check", szFile);

   if ( (fp=fopen(szCheckFile, "r"))==NULL)
      return 0;
   else
   {
      fclose(fp);
      return 1;
   }
}


/*
 * if checksum file exists, compares input file checksum to the checksum file
 */
void VALIDATE_CHECKSUM(char *szInput,
        int bVerbose)
{
   FILE *fp;
   char szBuf[SIZE_BUF],
        szSHA1check[SIZE_BUF], /* checksum read from .check file */
        szSHA1calc[SIZE_BUF],  /* calculated checksum of input file */
        szTmp[SIZE_BUF];
   int  bError=0;

   /*
    * calculate SHA1 check sum for input file
    */
   szSHA1calc[0] = 0;
   m_Csha1.Reset();
   if( !(m_Csha1.HashFile( szInput)) ) {
      printf(" Cannot open file %s for sha-1 calculation\n", szInput);
      bError = 1;
   }
   m_Csha1.Final();
   m_Csha1.ReportHash( szSHA1calc, SHA1::REPORT_HEX );

   if (!bError)
   {
      /*
       * get SHA1 check sum from .check file
       */
      strcpy(szSHA1check, "_error");
      sprintf(szTmp, "%s.check", szInput);
      if ( (fp=fopen(szTmp, "r"))!=NULL)
      {
	 char* errCheckMe = NULL;
         errCheckMe = fgets(szBuf, SIZE_BUF, fp);
         sscanf(szBuf, "%s %s %s\n", szTmp, szTmp, szSHA1check);
         fclose(fp);

         /*
          * compare checksums
          */
         if (!strcmp(szSHA1check, szSHA1calc))
         {
            printf(" Success: checksum validated\n");
         }
         else
         {
            printf(" Failure: checksum_mismatch (%s)\n", szInput);
            printf("        stored sha1: %s\n", szSHA1check);
            printf("    calculated sha1: %s\n", szSHA1calc);
         }
      }
      else
      {
         printf(" Failure: cannot read checksum file %s.check\n", szInput);
      }
   }
}
