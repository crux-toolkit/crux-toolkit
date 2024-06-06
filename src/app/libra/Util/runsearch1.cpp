
#include "Common/util.h"
#include "Parsers/mzParser/mzParser.h"
#include <ctype.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "Common/TPPVersion.h"

using namespace mzParser;

#define szVersion     "RUNSEARCH"
#define szAcknowledge "Written by J.Eng (c) 2000"

#define SIZE_FILE  512
#define SIZE_BUF  1024

#define MAX_SEARCHES 3

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 1
#endif
#ifndef EXIT_FAILURE
#define EXIT_FAILURE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif


void SET_OPTION(char *arg,
        int *iFirstScan,
        int *iLastScan,
        int *iMinNumPeaks,
        int *bRerunSearches,
        int *iWhichCharge,
        int *bStartAtSummary,
        char *szParamsFile);
void USAGE(int failure,
        char *command);


int main(int argc, char **argv)
{
   char *arg,
        szBaseName[SIZE_FILE],
        szInputFile[SIZE_FILE],
        szDatabaseFile[SIZE_FILE],
        szEmailAddress[SIZE_FILE],
        szCommand[SIZE_BUF],
        szPWD[SIZE_FILE],
        szParamsFile[SIZE_FILE],
        szSequest[SIZE_FILE];
   int  i,
        iFirstScan,
        iLastScan,
        iMinNumPeaks,
        iNumArg,
        iStartArgc,
        iWhichCharge,
        bStartAtSummary,
        bRerunSearches;
   FILE *fp;
   char *command;
   hooks_tpp handler(argc,argv); // set up install paths etc


   printf("\n %s\n %s\n %s\n\n", szVersion, szAcknowledge, szTPPVersionInfo);
   command = argv[0];
   if (argc<2)
   {
      USAGE(0, command);
   }

   iNumArg=0;
   iFirstScan=0;
   iLastScan=0;
   iMinNumPeaks=5;
   iStartArgc=1;
   iWhichCharge=0;
   bRerunSearches=FALSE;
   bStartAtSummary=FALSE;
   strcpy(szParamsFile, "sequest.params");
   iNumArg=1;
   arg = argv[iNumArg];


   /* processing arguments */

   while (iNumArg < argc)
   {
      if (arg[0] == '-')
         SET_OPTION(arg, &iFirstScan, &iLastScan, &iMinNumPeaks,
            &bRerunSearches, &iWhichCharge, &bStartAtSummary, szParamsFile);
      else
         break;

      iStartArgc++;
      arg = argv[++iNumArg];
   }

   /*
    * no input files specified, die
    */
   if (argc-iStartArgc == 0)
      exit(0);

   if (iFirstScan>0 && iLastScan>0)
   {
      if (iFirstScan > iLastScan)
      {
         printf(" Error - first scan # (%d) is less than last scan # (%d)\n\n",
            iFirstScan, iLastScan);
         exit(EXIT_FAILURE);
      }
   }

   szEmailAddress[0]='\0';
   szDatabaseFile[0]='\0';

   if ( (fp=fopen(szParamsFile, "r"))==NULL)
   {
      printf(" Cannot read %s.\n\n", szParamsFile);
      exit(EXIT_FAILURE);
   }
   while (fgets(szCommand, SIZE_BUF, fp))
   {
      if (!strncmp(szCommand, "first_database_name = ", 22) ||
          !strncmp(szCommand, "database_name = ", 16))
      {
         char szTmp[SIZE_BUF];

         sscanf(szCommand, "%s %s %s", szTmp, szTmp, szDatabaseFile);
      }
      if (!strncmp(szCommand, "email_address = ", 16))
      {
         sscanf(szCommand+15, "%s", szEmailAddress);
      }
   }
   fclose(fp);

   if (strlen(szDatabaseFile)==0)
   {
      printf(" No database specified in the %s file.\n\n", szParamsFile);
      printf(" Requires entry 'first_database_name = '.\n\n");
      exit(EXIT_FAILURE);
   }

   if ( (fp=fopen(szDatabaseFile, "r"))==NULL)
   {
      printf(" Cannot read the database file '%s'.  Check that this filename is correct.\n\n", szDatabaseFile);
      exit(EXIT_FAILURE);
   }
   fclose(fp);

   /*
    * Check existence of sequest binary.
    */
   strcpy(szSequest, "/cygdrive/c/xcalibur/system/programs/bioworksbrowser/sequest27.exe");
   unCygwinify(szSequest); // no effect in cygwin builds
   if ( (fp=fopen(szSequest, "r"))==NULL)
   {
      strcpy(szSequest, "/cygdrive/c/xcalibur/system/programs/bioworksbrowser/sequest.exe");
	  unCygwinify(szSequest); // no effect in cygwin builds
      if ( (fp=fopen(szSequest, "r"))==NULL)
      {
         printf(" Cannot find Sequest executable in c:\\xcalibur\\system\\programs\\bioworksbrowser\\.");
         exit(EXIT_FAILURE);
      }
      }
     fclose(fp);


   /*
    * Get PWD
    */
   char *ret=getcwd(szPWD, sizeof(szPWD));
   if (strchr(szPWD, ' ')
         && strchr(szPWD, '+')
         && strchr(szPWD, '!')
         && strchr(szPWD, '#')
         && strchr(szPWD, '|')
         && strchr(szPWD, '>')
         && strchr(szPWD, '<')
         && strchr(szPWD, '%')
         && strchr(szPWD, '&')
         && strchr(szPWD, '*')
         && strchr(szPWD, '[')
         && strchr(szPWD, ']'))
   {
      printf(" Current directory '%s' has an illegal character.\n", szPWD);
      exit(EXIT_FAILURE);
   }

   sprintf(szCommand, "chmod g+w .");
   verified_system(szCommand);

   /* all remaining args are input files; loop and process them */
   for (i=iStartArgc; i<argc; i++)
   {
      FILE *fpTmp;

      strcpy(szInputFile, argv[i]);

      if (strchr(szInputFile, ' ')
            && strchr(szInputFile, '+')
            && strchr(szInputFile, '!')
            && strchr(szInputFile, '#')
            && strchr(szInputFile, '|')
            && strchr(szInputFile, '>')
            && strchr(szInputFile, '<')
            && strchr(szInputFile, '%')
            && strchr(szInputFile, '&')
            && strchr(szInputFile, '*')
            && strchr(szInputFile, '[')
            && strchr(szInputFile, ']'))
      {
         printf(" Input '%s' has an illegal character in the filename, skipping ...\n", szInputFile);
      }
      else /* valid filename */
      {
         char szMZXML[512];
         struct stat pFileStat;

         szMZXML[0]='\0';

         if (!stat(szInputFile, &pFileStat))
	 /* yes, input file seems to exist */
         {
            if (!strcasecmp(szInputFile+strlen(szInputFile)-4, ".raw"))   /* finnigan */
            {
	       /* convert raw input file to mzXML if necessary */

               char szCommand[1024];

               strcpy (szMZXML, szInputFile);
               strcpy(szMZXML+strlen(szMZXML)-3, "mzXML");
      
               if (!(!stat(szMZXML, &pFileStat) && S_ISREG(pFileStat.st_mode)))  /* if mzXML does not already exist, */
               {
                  /* yes, convert to mzXML, as no existing mzXML found */
                  sprintf(szCommand, "readw.exe %s", szInputFile);
                  printf(" Converting Finnigan raw file to mzXML:  %s\n", szCommand);  /* micromass */
                  printf(" be patient ...\n");
                  verified_system(szCommand);

                  /* check validity of converted file and then archive original */
               }
            }
			else if (rampValidFileType(szInputFile))  /* mzXML or other ramp-supported mass spec format */
	    /* input file is already mzXML */
            {
               strcpy (szMZXML, szInputFile);
            }


	    /* at this point, we have a good mzXML file (either input, or converted from input).
	       so now we can convert mzXML->dta files, and run the search!
	    */
            if (!stat(szMZXML, &pFileStat) && S_ISREG(pFileStat.st_mode))
            {
               strcpy(szBaseName, szMZXML);
			   //if (!rampTrimBaseName(szBaseName)) { // strip .mzXML, .mzXML.gz, etc
                  *strrchr(szBaseName,'.') = 0; // trim .ext
			   //} MH: just strip any extension.
      
               /*---------------------------------*/
               /*  Create subdir and run searches */
               /*---------------------------------*/
               if (!bStartAtSummary)
               {
                  char szTarFile[SIZE_FILE];
      
                  sprintf(szTarFile, "%s.tgz", szBaseName);
      
                  /*---------------------*/ 
                  /* Create search files */
                  /*---------------------*/ 
                  sprintf(szCommand, "MzXML2Search -dta ");

                  if (iFirstScan!=0)
                     sprintf(szCommand+strlen(szCommand), " -F%d", iFirstScan); 
                  if (iLastScan!=0)
                     sprintf(szCommand+strlen(szCommand), " -L%d", iLastScan); 
                  if (iWhichCharge!=0)
                     sprintf(szCommand+strlen(szCommand), " -C%d", iWhichCharge); 
                  sprintf(szCommand+strlen(szCommand), " -P%d %s", iMinNumPeaks, szMZXML); 
                  printf("%s\n", szCommand);
                  verified_system(szCommand);

                  sprintf(szCommand, "chmod g+w %s; cp %s %s", szBaseName, szParamsFile, szBaseName);
                  printf("%s\n", szCommand);
                  verified_system(szCommand); // like system(), but handles multiple commands for win32
      
                  if ((fpTmp=fopen(szTarFile, "r"))!=NULL)  /* extract out any existing search data */
                  {
                     fclose(fpTmp);
      
                     sprintf(szCommand, "cd %s; tar --wildcards -xzf ../%s.tgz '*.out'", szBaseName, szBaseName);
                     printf("%s\n", szCommand);
                     verified_system(szCommand); // like system(), but handles multiple commands for win32
                  }
         
                  /*--------------*/ 
                  /* Run searches */
                  /*--------------*/ 
                  sprintf(szCommand, "cd %s; find . -maxdepth 1 -depth -name \"%s*.dta\" -print | xargs --max-args=50 %s -P%s",
                     szBaseName, szBaseName, szSequest, szParamsFile);
//                if (!bRerunSearches)
//                   sprintf(szCommand+strlen(szCommand), " -S");
                  printf("%s\n", szCommand);
                  verified_system(szCommand); // like system(), but handles multiple commands for win32
               }
               else  /* we get here if "start at summary set"; deprecated.
			extract out already completed search results for re-doing summary step */
               {
                  char szTarFile[SIZE_FILE];
      
                  sprintf(szTarFile, "%s.tgz", szBaseName);
      
                  if ((fpTmp=fopen(szTarFile, "r"))!=NULL)
                  {
                     fclose(fpTmp);
      
                     sprintf(szCommand, "mkdir %s; chmod g+w %s; cp %s %s/; cd %s; tar xzf ../%s.tgz",
                           szBaseName, szBaseName, szParamsFile, szBaseName, szBaseName, szBaseName);
                     printf("%s\n", szCommand);
                     verified_system(szCommand); // like system(), but handles multiple commands for win32
                  }
               }


	       /* at this point search is complete (either just ran, or was already run),
		  and a directory should have been populated with .out and .dta files */

               /*----------------*/ 
               /* Create pepXML file from .out file search results */
               /*----------------*/
               /* sprintf(szCommand, "cd %s; out2summary . > ../%s.html ", szBaseName, szBaseName);*/
               sprintf(szCommand, "Out2XML %s 1", szBaseName);

               printf("%s\n", szCommand);
               verified_system(szCommand);
      
               if ((fpTmp=fopen(szBaseName, "r"))!=NULL)
               {
                  fclose(fpTmp);
                  sprintf(szCommand, "cd %s ; find . -maxdepth 1 -depth -name \"*\" | tar czf ../%s.tgz --files-from -; cd .. ; rm -rf %s",
                        szBaseName, szBaseName, szBaseName);
                  verified_system(szCommand); // like system(), but handles multiple commands for win32
               }

      
            }
         } /*stat*/
         else
         {
            printf(" Error - %s does not exist as a file or directory.\n", szInputFile);
         }
      }
   }
   
   printf("\n");
   return(EXIT_SUCCESS);

} /*main*/


void SET_OPTION(char *arg,
        int *iFirstScan,
        int *iLastScan,
        int *iMinNumPeaks,
        int *bRerunSearches,
        int *iWhichCharge,
        int *bStartAtSummary,
        char *szParamsFile)
{
   int iTmp=0;

   switch (arg[1])
   {
      case 'F':
         if (sscanf(&arg[2], "%d", &iTmp) != 1 || iTmp<1)
            printf("Bad first scan #: '%s' ... ignored\n", &arg[2]);
         else
            *iFirstScan = iTmp;
         break;
      case 'L':
         if (sscanf(&arg[2], "%d", &iTmp) != 1 || iTmp<1)
             printf("Bad last scan #: '%s' ... ignored\n", &arg[2]);
         else
            *iLastScan = iTmp;
         break;
      case 'P':
         if (sscanf(&arg[2], "%d", &iTmp) != 1 || iTmp<1)
            printf("Bad min num of peaks: '%s' ... ignored\n", &arg[2]);
         else
            *iMinNumPeaks = iTmp;
         break;
      case 'p':
         if (arg[2] == '\0')
         {
            fprintf(stderr, "Error - no alternate params file specified\n");
            exit(0);
         }
         else
            strcpy(szParamsFile, arg+2);
         break;
      case 'R':
         *bRerunSearches = TRUE;
         break;
      case 'S':
         *bStartAtSummary= TRUE;
         break;
      case '1':
         *iWhichCharge = 1;
         break;
      case '2':
         *iWhichCharge = 2;
         break;
      case '3':
         *iWhichCharge = 3;
         break;
      default:
         break;
   }
   arg[0] = '\0';

} /*SET_OPTIONS*/


void USAGE(int failure,
      char *command)
{
   printf(" RUNSEARCH usage:  %s [options] [.mzXML files]\n\n", command);
   printf(" options = -F<num>   where num is an integer specifying the first scan\n");
   printf("           -L<num>   where num is an integer specifying the last scan\n");
// printf("           -S        re-do Summaries only (no searches)\n");
   printf("           -P<num>   minimum # of peaks required to create .dta file (default=15)\n");
   printf("           -p<str>   specify alternate SEQUEST parameters file.\n");
   printf("\n");
   printf("           -1        analyze +1 peptides only\n");
   printf("           -2        analyze +2 peptides only\n");
   printf("           -3        analyze +3 peptides only\n");
   printf("\n");
   printf("\n");

   exit(EXIT_FAILURE);

} /*USAGE*/
