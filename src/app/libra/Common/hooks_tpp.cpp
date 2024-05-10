//
// hooks_tpp.cxx
//
// handle install dir issues etc
//
// Copyright (c) 2006 Insilicos, LLC
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
// TODO: smarten this up with a config file or somesuch - only uses default for now
//



#include "util.h"
#include "constants.h"
#include "TPPVersion.h"
#include <iostream>

//
// installdir issues - initialize defaults
//
bool getIsInteractiveMode() { // TPP (web oriented) vs LabKey (headless) usage style
	// env var may be set by hooks_tpp in response to XML_ONLY in argv, or at system level
#ifdef XML_ONLY
	return false;
#else
	return getenv("XML_ONLY")?false:true;	
#endif
}

#ifdef _DEBUG
static void my_setenv(const char *var,const char *val) {
	char *buf = (char *)malloc(strlen(var)+strlen(val)+2);
	strcpy(buf,var);
	strcat(buf,"=");
	strcat(buf,val);
	putenv(buf);
	// don't free(buf); - on some systems that will alter the environment
}
#endif
   
      
//
// instatiate at program start to handle install dir issues etc
//
hooks_tpp::hooks_tpp(int &argc, char *argv[]) {
	for (int i=argc;i--;) {
		if (!strcmp(argv[i],"XML_ONLY")) {
			// this is the sign to skip web stuff like bitmaps in proteinprophet
			// note we use strdup as on some systems the proto specifies non-const arg
			putenv(strdup("XML_ONLY=1")); // for use in getIsInteractiveMode()
			// done with this flag
			for (int n=i+1;n<argc;n++) {
				argv[n-1] = argv[n];
			}
			argc--;
		} else if (!strcmp(argv[i],"-installtest")) {
			std::cout << "testing " << argv[0] << " (" << szTPPVersionInfo << ") installation..." ;
			std::cout << "  OK." << std::endl;
			exit(0);
		}  else if (!strcmp(argv[i],"-versioncheck")) {
			std::cout << szTPPVersionInfo ;
			exit(0);
		}
	}
#ifdef _DEBUG
	// set up for CGI debug - transfer commandline to env
	my_setenv("REQUEST_METHOD","GET");
	my_setenv("SCRIPT_NAME",argv[0]);
	if (argc>1) {
		std::string query(argv[1]);
		for (int arg=2;arg<argc;arg++) {
			query+=argv[arg];
		}
		my_setenv("QUERY_STRING",query.c_str());
	}
#endif
}


//
// let's use the hostname instead of localhost, so we can be a server
//
#ifndef MAXHOSTNAMELEN
#define MAXHOSTNAMELEN 512
#endif
static char pepxml_std_xsl[MAXHOSTNAMELEN+10];
const char *getPepXML_std_xsl() {
   if (!pepxml_std_xsl[0]) {  // first access?
      const char *http_host=getenv("HTTP_HOST"); // provide by webserver
      const char *webserver_url; // TPP style
      if (http_host) { // take your cue from the webserver
         snprintf(pepxml_std_xsl,MAXHOSTNAMELEN+10,"http://%s/",http_host);
      } else if (NULL!=(webserver_url=getenv("WEBSERVER_URL"))) { // take your cue from config env var
		 pepxml_std_xsl[0] = 0; // we'll construct it
         if (strncmp(webserver_url,"http://",7)) {
            strncpy(pepxml_std_xsl,"http://",MAXHOSTNAMELEN+10); // url didn't start with http
         }
         size_t len = strlen(pepxml_std_xsl);
         strncpy(pepxml_std_xsl+len,webserver_url,MAXHOSTNAMELEN+10-len);
         len = strlen(pepxml_std_xsl);
         if (pepxml_std_xsl[len-1] != '/') {
            strncpy(pepxml_std_xsl+len,"/",MAXHOSTNAMELEN+10-len);
         }
      } else {
//         char hostname[MAXHOSTNAMELEN+1];
#ifdef WINDOWS_NATIVE // stinky winsock
        // WSADATA wsaData;
         //bool sockOK = !WSAStartup(MAKEWORD( 2, 0 ),&wsaData);
#endif
	 //DDS: this broke validation against the schema
	 //         if (!gethostname(hostname,MAXHOSTNAMELEN)) {
	 //            sprintf(pepxml_std_xsl,"http://%s/",hostname);
	 //         } else {
#ifdef WINDOWS_NATIVE // stinky winsock
	//  int err = WSAGetLastError();
#endif
	    strcpy(pepxml_std_xsl,DEFAULT_PEPXML_STD_XSL); // default
	//      }
#ifdef WINDOWS_NATIVE // stinky winsock
	    //if (sockOK) {
	   //WSACleanup();
	    //}
#endif
      }
   }
   return pepxml_std_xsl;
}

static char pepxml_std_xsl_web_path[MAXHOSTNAMELEN+20];
const char *getPepXML_std_xsl_web_path() {
   if (!pepxml_std_xsl_web_path[0]) {  // first access?
      strcpy(pepxml_std_xsl_web_path,getPepXML_std_xsl());
#ifdef WINDOWS_NATIVE  // go to a different dir to get non-cygwin version of schemas
      strcat(pepxml_std_xsl_web_path,"schema/");
#endif
   }
   return pepxml_std_xsl_web_path;
}
