/*

Program       : win32_system.c, helper functions for TPP under win32 
Author        : Brian Pratt, Insilicos LLC                                                       
Date          : 5.1.06 


Copyright (C) 2006 Insilicos LLC

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

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "sysdepend.h"
#include "util.h"
#include <errno.h>
#include <sys/stat.h>

// here is where we actually implement what we redefined elsewhere
#undef system
#undef popen
using namespace std;
#ifdef __CYGWIN__
//
// need to make sure we hit the utilities in cygwin an not any possible
// native windows implementations
//
static char *process_cmd_for_cygwin(const char *in_cmd) {
   char *cmd = strdup(in_cmd);
   char *result = strdup("");
   char *next = cmd;
   while (next && *next) {
      // see if we can specify this as being in /usr/bin
      struct stat info;
      char *space, *full, *newcmd;
	  const char *binpath;	  
      char *copy = strdup(next);
      // locate next command start (as ';')
      next = strchr(next,';');
      if (next) {
         *next = 0; // temp truncate to end of command
      }
      // locate space between cmd and arg
      space = strchr(copy,' ');
      // now try prepending path
      binpath = "/usr/bin/";
      full = (char *)malloc(strlen(copy)+strlen(binpath)+1);
      if (space) {
         *space = 0; // truncate after cmd to get a filename
      }
      strcpy(full,binpath);
      strcat(full,copy);
      // check to see if addition of /usr/bin makes a real filename
      newcmd = stat(full,&info)?copy:full;
      // and add that to the result we're building
      char *tmp = result;
      result = (char *)malloc(strlen(result)+strlen(newcmd)+(space?strlen(space+1):0)+3);
      strcpy(result,tmp);
      strcat(result,newcmd);
      if (space) {
         *space = ' '; // restore
         strcat(result,space); // append args
      }
      if (next) {
         strcat(result,";");
         next++;
      }
      free(copy);
      free(full);
      free(tmp);     
   }
   free(cmd);
   return result;
}

int win32_system(const char *cmd) {  
   char *processed = process_cmd_for_cygwin(cmd);
   int result = system(processed);
   free(processed);
   return result;
}
FILE *win32_popen(const char *cmd, const char *mode) {
   char *processed = process_cmd_for_cygwin(cmd);
   FILE *result = popen(processed,mode);
   free(processed);
   return result;
}

#elif defined(WINDOWS_NATIVE) // MSVC or MinGW

// local helper functions
static char *process_command_for_windows(const char *cmd);
static char *setComSpec() {
	// now swap out cmd.exe for our copy of cmd.exe (installer should copy local cmd.exe
	// to tpp bin
	char *oldComSpec = NULL;
	char *cp;
	struct stat info;
	char *newComSpec = (char *)malloc(strlen("COMSPEC=")+strlen(getBinPath())+strlen("cmd.exe")+3);
	if (cp=getenv("COMSPEC")) {
		oldComSpec=strdup(cp);
	} else if (cp=getenv("ComSpec")) {
		oldComSpec=strdup(cp);
	} 
	sprintf(newComSpec,"%scmd.exe",getBinPath());
	if (oldComSpec && (-1==stat(newComSpec,&info))) {
		copy_file(oldComSpec,newComSpec);
	}
	sprintf(newComSpec,"COMSPEC=%scmd.exe",getBinPath());
	fixpath(newComSpec); // get those slashes pointing left
	putenv(newComSpec);
	sprintf(newComSpec,"ComSpec=%scmd.exe",getBinPath());
	fixpath(newComSpec); // get those slashes pointing left
	putenv(newComSpec);
	free(newComSpec);
	return oldComSpec;
}
static void restoreComSpec(char *oldComSpec) {
	if (oldComSpec) { // restore previous comspec
		char *newComSpec = (char *)malloc(strlen(oldComSpec)+9);
		sprintf(newComSpec,"COMSPEC=%s",oldComSpec);
		putenv(newComSpec);
		sprintf(newComSpec,"ComSpec=%s",oldComSpec);
		putenv(newComSpec);
		free(oldComSpec);
		free(newComSpec);
	}	
}

static void stripquotes(char *fname) {
	if (('\"' == *fname)||('\''==*fname)) {
		int len = (int)strlen(fname);
		if (fname[len-1]==*fname) {
			fname[len-1]=0;
			memmove(fname,fname+1,len);
		}
	}
}

// split up multipart commands, avoid use of cygwin-ese
int win32_system(const char *cmd) {  
   int result = 0;
   const char *next=cmd;
   char curdir[1024];
   char *ret=getcwd(curdir,sizeof(curdir)); // preserve cwd
   while (next && *next) {
      char *copy=strdup(next);
      char *processed;
      char *cp;
      int local=0;
      const char *semi = strchr(next,';');
      if (semi) {
         copy[semi-next]=0;
      }
      processed = process_command_for_windows(copy);
      // quick check on gnuplot issues
      cp = strstr(processed,"pgnuplot");

      if (cp && !strchr(cp,'<')) { // pgnuplot only reads from stdin
         char *oldprocessed = processed;
         processed=(char *)malloc(strlen(oldprocessed)+2);
         strcpy(processed,oldprocessed);
         cp = strchr(processed,' ');
         if (cp) {
           memmove(cp+2,cp,strlen(cp)+1);
           *(cp+1)='<';
           *(cp+2)=' ';
           fixpath(cp+3);
         }
         free(oldprocessed);
      }


      if (!strncmp(copy,"cd ",3)) {
	     stripquotes(copy+3); // get path seps pointing left
         verified_chdir(copy+3);
      } else if ((!strncmp(copy,"mkdir ",6)) && strncmp(copy+6,"-p ",3)) {
	     stripquotes(copy+6); // get path seps pointing left
         mkdir(copy+6);
      } else {
		 char *oldComSpec=NULL;
		 if (getIsInteractiveMode()) { // TPP (web oriented) vs LabKey (headless) usage style
			// now swap out cmd.exe for our copy of cmd.exe (installer should copy local cmd.exe
			// to tpp bin)
			oldComSpec=setComSpec();
		 }
         fflush(stdout);
         fflush(stderr);
         result = system(processed);

         if (result) {
            printf("\ncommand \"%s\" failed: %s\n",processed,strerror(result));
         }
         restoreComSpec(oldComSpec); // restore previous comspec
	   }
      if (local<0) {
         int err=errno;
         result = local;
         if (copy!=processed) {
            printf("cmd %s (%s) failed, %s \n", copy, processed, strerror(err));
         } else {
            printf("cmd %s failed, %s \n",copy, strerror(err));
         }
      } else {
         result |= local;
      }
      if (copy!=processed) {
         free(processed);
      }
      free(copy);
      next = semi;
      if (next) { // had a semicolon, advance to next command
         do {
            next++;
         } while (*next && isspace(*next));
      }
   }
   verified_chdir(curdir); // restore cwd
   return result;
}

FILE *win32_popen(const char *cmd, const char *mode) {
   FILE * result;
   char *copy = strdup(cmd);
   char *processed = process_command_for_windows(copy);
   char *oldComSpec=NULL;

   if (getIsInteractiveMode()) { // TPP (web oriented) vs LabKey (headless) usage style
	   // now swap out cmd.exe for our copy of cmd.exe (installer should copy local cmd.exe
	   // to tpp bin)
	   oldComSpec=setComSpec();
   }
   AllocConsole(); // IIS doesn't furnish a console by default, so stdout pipe fails
   result = _popen(processed,mode);
   if (!result) {
      printf("pipe \"%s\" failed: %s\n",processed,strerror(errno));
   }

   restoreComSpec(oldComSpec); // restore previous comspec

   if (copy!=processed) {
      free(processed);
   }
   free(copy);
   return result;
}

static char *process_command_for_windows(const char *cmd) {
   // handle perl etc
   char *result = (char *)cmd;
   char *buf;
   char *space;
   int needs_cgi_path = 0;
   int l;
   struct stat statbuf;

   const char *cp=strstr(cmd,".pl");
   const char *explicit_perl=strstr(cmd,"perl "); // is perl explicitly mentioned already?
   if (cp && ((!explicit_perl) || (cp < explicit_perl)))
   {
	   cp += 3;
	   if (*cmd == '"' && *cp == '"')
		   cp++;

      if (*cp == ' ' || *cp == 0)
	  {
		  result = (char *)malloc(strlen(cmd)+25); // extra in case you want the profile option
		  strcpy(result,"perl ");
		  // strcpy(result,"perl -d:DProf ");
		  strcat(result,cmd);
	  }
   } 
   if (!strncmp(cmd,"tar ",4)) { // mingw (bsd) tar doesn't need or want "--wildcards"
      const char *arg = strstr(cmd,"--wildcards");
      if (arg) { // mingw tar (like many) doesn't know that switch
		 result = strdup(cmd);
         int skip = (int)strlen("--wildcards")+1;
		 char *rarg = result + (arg-cmd);
         memmove(rarg,rarg+skip,strlen(rarg+skip)+1);
      }
   } 
   buf = (char *)malloc(strlen(result)+strlen(getBinPath())+8);
   // do we need to add the path? or is this a shell cmd like "dir"?
   strcpy(buf,getBinPath());
   fixpath(buf); // get those slashes pointing left
   strcat(buf,result);
   space = strchr(result,' ');
   if (space) { // cut off the arglist
      *(buf+strlen(getBinPath())+(space-result))=0;
   }
   if (!stat(buf,&statbuf)) {
      needs_cgi_path = 1; // adding cgi path gives an actual exe path
   } else {
      strcpy(buf+strlen(buf),".exe"); // missing .ext?
      if (!stat(buf,&statbuf)) {
         needs_cgi_path = 1; // adding cgi path gives an actual exe path
      }
   }
   if (needs_cgi_path) { 
      // enclose the program name in quotes in case of path with spaces
      if (strchr(buf,' ')) {
         memmove(buf+1,buf,strlen(buf)+1);
         *buf='\"';
         strcat(buf,"\"");
      }
      if (space) { // add the args back in
         strcat(buf,space);
      }
   } else {
      strcpy(buf,result);
   }
   if (cmd!=result) {
      free(result); // we did a malloc to insert perl in blah.pl
   }
   // kill trailing spaces
   for (l=(int)strlen(buf);' '==buf[l-1];) {
      buf[--l]=0;
   }
   return buf;
}
#endif

