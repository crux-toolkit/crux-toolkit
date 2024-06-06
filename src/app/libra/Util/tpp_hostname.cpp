//
// tpp_hostname.cpp
//
// DESCRIPTION:
//
// a little exe to give perl access to same hostname selection code as C++
// also can be used for checking existence of mzXML/mzData file for a file basename,
// and for pulling the pepxml and protxml filename extensions - used
// pretty much in place of a proper perl module since those have proven to be
// nightmarish for end users to compole
//
// Copyright (c) 2006,2009 Insilicos, LLC
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// 
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// 
// $Author: bpratt $
//
//
// NOTES:
//
//
// TODO:
//

#include "Common/sysdepend.h"
#include "Common/util.h"
#include "Parsers/mzParser/mzParser.h"
#include "Common/TPPVersion.h" // contains version number, name, revision

using namespace mzParser;

int main(int argc, char *argv[]) {
   hooks_tpp hooks(argc,argv); // installdir issues etc  
   char *p=strdup(getPepXML_std_xsl_web_path()); // like http://foobar/schema/pepxml.xsd
   if (strstr(p, "http://")) { 
     char *h=p+strlen("http://");
     if (h) {
       memmove(p,h,strlen(h)+1); // like foobar/schema/pepxml.xsd
     }
   }
   char *s = strchr(p,'/');
   if (s) {
      *s = 0; // like foobar
   }
   if (argc > 1) { // optional args
	   const char *cmd_uncompress_to_tmpfile = "uncompress_to_tmpfile!";
	   if (!strcmp(argv[1],"GET_PEPXML_EXT!")) {
		   printf("%s",get_pepxml_dot_ext()); // return canonical pepxml .ext
	   } else if (!strcmp(argv[1],"GET_PROTXML_EXT!")) {
		   printf("%s",get_protxml_dot_ext()); // return canonical protxml .ext
	   } else if (!strcmp(argv[1],"hasValidPepXMLFilenameExt!")) {
		   const char *ext = (argc > 2)?hasValidPepXMLFilenameExt(argv[2]):NULL;
		   if (ext) {
			  printf("%s",ext); // return .ext if found
		   }
	   } else if (!strcmp(argv[1],"hasValidProtXMLFilenameExt!")) {
		   const char *ext = (argc > 2)?hasValidProtXMLFilenameExt(argv[2]):NULL;
		   if (ext) {
			  printf("%s",ext); // return .ext if found
		   }
	   } else if (!strcmp(argv[1],"versionInfo!")) {
 		  printf("%s",szTPPVersionInfo); 
	   } else if (!strncmp(argv[1],cmd_uncompress_to_tmpfile,strlen(cmd_uncompress_to_tmpfile))) {
		   // check for optional maxlength arg after the "!"
		   int maxchar = atoi(argv[1]+strlen(cmd_uncompress_to_tmpfile));
		   std::string fname = (argc > 2)?uncompress_to_tmpfile(argv[2],maxchar):"";
		   printf("%s",fname.c_str()); // returns original fname if no uncompress needed
	   } else if (!strcmp(argv[1],"getGnuplotBinary!")) {
		   const char *fname = getGnuplotBinary();
		   if (fname) {
			  printf("%s",fname); // return .ext if found
		   }
	   } else if (rampValidFileType(argv[1])) {
		   // assume it's a check to see if this is a supported mzxml/mzdata/mzml/? type
		   printf("%s",argv[1]);
	   } else if (!strcmp(argv[1],"HomePath!")) {
 		  printf("%s",getHomePath());
	   } else if (!strcmp(argv[1],"TPPPort!")) {
	          printf("%s",getTPPPort()); 
	   } else if (!strcmp(argv[1],"BinPath!")) {
 		  printf("%s",getBinPath()); 
	   } else if (!strcmp(argv[1],"CgiPath!")) {
 		  printf("%s",getCgiPath()); 
	   } else if (!strcmp(argv[1],"ConfPath!")) {
 		  printf("%s",getConfPath()); 
	   } else if (!strcmp(argv[1],"HtmlPath!")) {
 		  printf("%s",getHtmlPath()); 
	   } else if (!strcmp(argv[1],"DataPath!")) {
 		  printf("%s",getDataPath()); 
	   } else if (!strcmp(argv[1],"LogPath!")) {
 		  printf("%s",getLogPath()); 
	   } else if (!strcmp(argv[1],"BaseUrl!")) {
 		  printf("%s",getBaseUrl()); 
	   } else if (!strcmp(argv[1],"DataUrl!")) {
 		  printf("%s",getDataUrl()); 
	   } else if (!strcmp(argv[1],"CgiUrl!")) {
 		  printf("%s",getCgiUrl()); 
	   } else if (!strcmp(argv[1],"HtmlUrl!")) {
 		  printf("%s",getHtmlUrl()); 
	   } else {
		   // assume it's a basename to be expanded to an mzxml/mzdata file if possible
		   int len = 1000+(int)strlen(argv[1]);
		   char *buf = (char *)malloc(len);
		   char *result = rampConstructInputFileName(buf,len,argv[1]);
		   printf("%s",result?result:"no_raw_data_with_this_basename");
		   free(buf);
	   }
   } else {
      printf("%s",p);
   }
   free(p);
   return 0;
}
