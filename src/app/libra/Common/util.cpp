/*
 *
 * Revision 1.13  2006/05/25 22:25:43  brendanx
 * Win32 build for VC++7.1 includes:
 * - All relative paths (to new module win_lib)
 * - Compile time warnings fixed
 * - Build instructions
 * - All char[line_width_] allocations moved off stack to avoid stack overflow
 * - Better use of PATH for XML_ONLY build
 *
 * Also in this check-in (wish they could have been separate):
 * - PeptideProphet support for X! Tandem Native scoring (v0.01)
 *
 * Revision 1.11  2006/03/17 23:27:57  pcbrefugee
 * added read_dta_or_out_from_tgz_file() to assist in Mascot2XML regression tests
 *
 * Revision 1.10  2006/01/06 20:05:23  pcbrefugee
 * There were a few places in the TPP code with the comment "//Convert to same case first for lack of strstri" wherein filenames being checked for inclusion of the webserver_root path got trashed from the *nix case sensitive filename point of view.  The solution was just to implement strstri() (in util.c).
 *
 * Revision 1.9  2005/12/16 20:03:45  pcbrefugee
 * added a bit more sophistication to regression test handling: more crossplatform tolerant, more informative error messages
 *
 * Revision 1.8  2005/12/08 22:05:10  pcbrefugee
 * need to include <errno.h> for linux build
 *
 * Revision 1.7  2005/12/08 21:58:32  pcbrefugee
 * add a utility routine to copy one file by name to another
 *
 * Revision 1.6  2005/11/21 23:48:30  pcbrefugee
 * portability tweak (stricmp vs strcasecmp)
 *
 * Revision 1.5  2005/11/21 23:38:02  pcbrefugee
 * portability tweak (stricmp vs strcasecmp)
 *
 * Revision 1.4  2005/11/21 23:11:36  pcbrefugee
 * 1) Added integrated regression test support for major C++ implemented TPP components called by xinteract.  See TESTING doc for details.
 * 2) Moving toward use of .pepXML instead of .xml as standard pepXML filename extension (much as we do with .mzXML), but still actually using .xml at this time - the actual choice has been virtualized, though.  Still need to modify perl CGIs to handle the switch (or rather, to handle any extension thrown at them, as continued support for .xml extension for pepXML files is required).
 * 3) Some cleanups of unused variables etc that were causing noisy compiles, and portability tweaks.
 *
 * Revision 1.3  2005/10/20 00:27:59  pcbrefugee
 * rename util.c's getline() to get_line(), which won't conflict with stdio.h
 *
 * Revision 1.2  2005/10/19 23:49:12  pcbrefugee
 * TPP components now display release version info and build numbers.  This text from the Makefile explains:
 *
 * #
 * # A note about TPP Version and Build Numbers
 * #
 * # The file src/common/TPPVersion.h controls the version info displayed by TPP components.
 * # This info is furnished to the various exe and cgi files by linking with an object file
 * # that also contains a build number (timestamp, actually).  The source code for this object
 * # file is generated automatically, see target TPPVersionInfo.cxx: below.  For perl scripts
 * # the info is generated as an include file TPPVersionInfo.pl.
 * #
 * # The build number (timestamp) is refreshed when src/common/TPPVersion.h is changed,
 * # or during a "make clean all".  The intent is that the build number be the same across
 * # all TPP components delivered to the end user.  It could of course be refreshed with
 * # every "make all" but that would mean every app would relink every time, which would bog
 * # down the development cycle.  If you wish to force a build number update without doing
 * # a "make clean", use "make new_buildnum".
 * #
 * #
 *
 * Revision 1.1  2005/06/20 19:21:20  dshteyn
 * First commit after reorg
 *
 * Revision 1.2  2005/02/09 20:13:36  adkeller
 * update before TPP branch
 *
 * Revision 1.1  2003/03/11 00:46:43  rhubley
 *   Several changes:
 * 	- Moved the location of util.c since it is now shared by
 * 	  both bin and cgi applications.
 *         - Made sure programs would compile appropriately under
 *           windows and unix.
 *         - Fixed problem with missing type definitions in Jimmy's
 *           recent changes.
 *         - Final checkin before release 6.0
 *
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "util.h"
#include "constants.h"
#include <errno.h>
#include <sys/stat.h>
#include "TPPVersion.h"
#include <iostream>
#include <vector>
#include <algorithm>

#ifndef _MSC_VER
#include <stdarg.h>
#endif

// for quoting PERL_BIN
//
// double macro expansion required:
// see http://gcc.gnu.org/onlinedocs/cpp/Stringification.html#Stringification
#define PREPROCESSOR_QUOTE_2(s) #s
#define PREPROCESSOR_QUOTE(s) PREPROCESSOR_QUOTE_2(s)

#ifdef  PERL_BIN
#define QUOTED_PERL_BIN PREPROCESSOR_QUOTE(PERL_BIN)
#endif

#define LF 10
#define CR 13
#include <string>
using namespace std;

//
// TPP Installation/Configuration information ----------------------------------
//
// (All filesystem directory paths have '/' appended to them)

// TPP "home", the installation directory
static char _tpp_home_path[1024]={0};
const char *getHomePath() {
  if (!_tpp_home_path[0]) { // first access?
    // try TPP_HOME env first then default
    if (!check_env_var("TPP_HOME",_tpp_home_path,sizeof(_tpp_home_path),true)) {
      // no env var set
      strncpy(_tpp_home_path,TPP_HOME,sizeof(_tpp_home_path));
      // pass this on to any child processes
      std::string setenv = "TPP_HOME=";
      setenv += _tpp_home_path;
      putenv(strdup(setenv.c_str())); // intentional leak, freeing may affect env
    }
  }
  return _tpp_home_path;
}

// TPP "port", the connection
static char _tpp_port[6]={0};
const char *getTPPPort() {
  if (!_tpp_port[0]) { // first access?
    // try TPP_HOME env first then default
    if (!check_env_var("TPP_PORT",_tpp_port,sizeof(_tpp_port),true)) {
      // no env var set
      strncpy(_tpp_port,TPP_PORT,sizeof(_tpp_port));
      // pass this on to any child processes
      std::string setenv = "TPP_PORT=";
      setenv += _tpp_port;
      putenv(strdup(setenv.c_str())); // intentional leak, freeing may affect env
    }
  }
  return _tpp_port;
}

// Seq2MS "url", the connection
static char _seq2ms_url[1024]={0};
const char *getSeq2MSUrl() {
  if (!_seq2ms_url[0]) { // first access?
    // try TPP_HOME env first then default
    if (!check_env_var("SEQ2MS_SOURCE_URL",_seq2ms_url,sizeof(_seq2ms_url),true)) {
      // no env var set
      strncpy(_seq2ms_url,"",sizeof(_seq2ms_url));
      // pass this on to any child processes
      std::string setenv = "SEQ2MS_SOURCE_URL=";
      setenv += _seq2ms_url;
      putenv(strdup(setenv.c_str())); // intentional leak, freeing may affect env
    }
  }
  return _seq2ms_url;
}

// Seq2MS "model", the model path
static char _seq2ms_model[1024]={0};
const char *getSeq2MSModel() {
  if (!_seq2ms_model[0]) { // first access?
    // try TPP_HOME env first then default
    if (!check_env_var("SEQ2MS_DEFAULT_MODEL",_seq2ms_model,sizeof(_seq2ms_model),true)) {
      // no env var set
      strncpy(_seq2ms_model,"",sizeof(_seq2ms_model));
      // pass this on to any child processes
      std::string setenv = "SEQ2MS_DEFAULT_MODEL=";
      setenv += _seq2ms_model;
      putenv(strdup(setenv.c_str())); // intentional leak, freeing may affect env
    }
  }
  return _seq2ms_model;
}

// Directory path containg TPP executable programs
static char _tpp_bin_path[1024]={0};
const char *getBinPath() {
  if (!_tpp_bin_path[0]) { // first access?
    std::string bin(getHomePath());
    bin += "bin/";
    strncpy(_tpp_bin_path,bin.c_str(),sizeof(_tpp_bin_path));
  }
  return _tpp_bin_path;
}

// Directory path containing TPP common gateway interface programs
static char _tpp_cgi_path[1024]={0};
const char *getCgiPath() {
  if (!_tpp_cgi_path[0]) { // first access?
    std::string cgi(getHomePath());
    cgi += "cgi-bin/";
    strncpy(_tpp_cgi_path,cgi.c_str(),sizeof(_tpp_cgi_path));
  }
  return _tpp_cgi_path;
}

// Directory path containing TPP configuration files
static char _tpp_conf_path[1024]={0};
const char *getConfPath() {
  if (!_tpp_conf_path[0]) { // first access?
    std::string conf(getHomePath());
    conf += "conf/";
    strncpy(_tpp_conf_path,conf.c_str(),sizeof(_tpp_conf_path));
  }
  return _tpp_conf_path;
}

// Directory path containing TPP static public www files
static char _tpp_html_path[1024]={0};
const char *getHtmlPath() {
  if (!_tpp_html_path[0]) { // first access?
    std::string html(getHomePath());
    html += "html/";
    strncpy(_tpp_html_path,html.c_str(),sizeof(_tpp_html_path));
  }
  return _tpp_html_path;
}

// TPP "data", the directory containing users data
static char _tpp_data_path[1024]={0};
const char *getDataPath() {
  if (!_tpp_data_path[0]) { // first access?
    // try TPP_DATADIR env first then default
    if (!check_env_var("TPP_DATADIR",_tpp_data_path,sizeof(_tpp_data_path),true)) {
      // no env var set
      strncpy(_tpp_data_path,TPP_DATADIR,sizeof(_tpp_data_path));
      // pass this on to any child processes
      std::string setenv = "TPP_DATADIR=";
      setenv += _tpp_data_path;
      putenv(strdup(setenv.c_str())); // intentional leak, freeing may affect env
    }
  }
  return _tpp_data_path;
}

// Directory path containing TPP log files
static char _tpp_log_path[1024]={0};
const char *getLogPath() {
  if (!_tpp_log_path[0]) { // first access?
    std::string log(getHomePath());
    log += "log/";
    strncpy(_tpp_log_path,log.c_str(),sizeof(_tpp_log_path));
  }
  return _tpp_log_path;
}

// Base prefix for all TPP URLs
static char _tpp_base_url[1024]={0};
const char *getBaseUrl() {
  if (!_tpp_base_url[0]) { // first access?
    // try TPP_BASEURL env first then default
    if (!check_env_var("TPP_BASEURL",_tpp_base_url,sizeof(_tpp_base_url),true)) {
      // no env var set
      strncpy(_tpp_base_url,"/" TPP_BASEURL,sizeof(_tpp_base_url));
      // pass this on to any child processes
      std::string setenv = "TPP_BASEURL=";
      setenv += _tpp_base_url;
      putenv(strdup(setenv.c_str())); // intentional leak, freeing may affect env
    }
  }
  return _tpp_base_url;
}

// Base prefix for TPP Data URLs
static char _tpp_data_url[1024]={0};
const char *getDataUrl() {
  if (!_tpp_data_url[0]) { // first access?
    // try TPP_DATAURL env first then default
    if (!check_env_var("TPP_DATAURL",_tpp_data_url,sizeof(_tpp_data_url),true)) {
      // no env var set
      strncpy(_tpp_data_url,"/" TPP_DATAURL,sizeof(_tpp_data_url));
      // pass this on to any child processes
      std::string setenv = "TPP_DATAURL=";
      setenv += _tpp_data_url;
      putenv(strdup(setenv.c_str())); // intentional leak, freeing may affect env
    }
  }
  return _tpp_data_url;
}

// Base prefix for TPP CGI URLs
static char _tpp_cgi_url[1024]={0};
const char *getCgiUrl() {
  std::string url(getBaseUrl());
  url += "cgi-bin/";
  strncpy(_tpp_cgi_url,url.c_str(),sizeof(_tpp_cgi_url));
  return _tpp_cgi_url;
}

// Base prefix for TPP HTML URLs
static char _tpp_html_url[1024]={0};
const char *getHtmlUrl() {
  std::string url(getBaseUrl());
  url += "html/";
  strncpy(_tpp_html_url,url.c_str(),sizeof(_tpp_html_url));
  return _tpp_html_url;
}


//
// TPP Path/URL Conversions ----------------------------------------------------
//
// (Uses above TPP directory functions and are meant to replace the deprecated
// translate_* functions below)
//

//
// Converts a path to a file under the TPP data directory to a TPP url by 
// replacing the TPP base data path with the TPP base data url.  The Url is 
// expected to be large enough to hold results.  Returns 1 if successful.
//
int tppDataPath2Url( const char *path, char *url ) {

  char* rel = strstri( path, getDataPath() );
  if ( NULL != rel ) {	// found data path
    rel += strlen(getDataPath());
    strcpy( url, getDataUrl() );
    strcat( url, rel );
    return 1;
  } else {	
    strcpy( url, path );
    return 0;
  }
}

//
// Converts a relative TPP url to a path to a file under the TPP data directory
// by replacing the base TPP url part with the base TPP data directory part.
// Path is expected to be large enough to hold rsults. Returns 1 if successful.
//
int tppDataUrl2Path( const char *url, char *path ) {

  char* rel = strstri( url, getDataUrl() );
  if ( NULL != rel ) {	// found data url
    rel += strlen(getDataUrl());
    strcpy( path, getDataPath() );
    strcat( path, rel );
    return 1;
  } else {			// just append it
    if ( '/' == *rel ) rel++;
    sprintf( path, "%s/%s", getDataPath(), rel );
    return 0;
  }
}


//
// Miscellaneous --------------------------------------------------------------
//

void getword(char *word, char *line, char stop) {
  int x = 0,y;

  for(x=0;((line[x]) && (line[x] != stop));x++)
    word[x] = line[x];

  word[x] = '\0';
  if(line[x]) ++x;
  y=0;

  while( (line[y++] = line[x++]) );
}

char *makeword(char *line, char stop) {
  int x = 0,y;
  char *word = (char *) malloc(sizeof(char) * (strlen(line) + 1));

  for(x=0;((line[x]) && (line[x] != stop));x++)
    word[x] = line[x];

  word[x] = '\0';
  if(line[x]) ++x;
  y=0;

  while( (line[y++] = line[x++]) );
  return word;
}

char *fmakeword(FILE *f, char stop, int *cl) {
  int wsize;
  char *word;
  int ll;

  wsize = 102400;
  ll=0;
  word = (char *) malloc(sizeof(char) * (wsize + 1));

  while(1) {
    word[ll] = (char)fgetc(f);
    if(ll==wsize) {
      word[ll+1] = '\0';
      wsize+=102400;
      word = (char *)realloc(word,sizeof(char)*(wsize+1));
    }
    --(*cl);
    if((word[ll] == stop) || (feof(f)) || (!(*cl))) {
      if(word[ll] != stop) ll++;
      word[ll] = '\0';
      return word;
    }
    ++ll;
  }
}

char x2c(char *what) {
  register char digit;

  digit = (what[0] >= 'A' ? ((what[0] & 0xdf) - 'A')+10 : (what[0] - '0'));
  digit *= 16;
  digit += (what[1] >= 'A' ? ((what[1] & 0xdf) - 'A')+10 : (what[1] - '0'));
  return(digit);
}

void unescape_url(char *url) {
  register int x,y;

  for(x=0,y=0;url[y];++x,++y) {
    if((url[x] = url[y]) == '%') {
      url[x] = x2c(&url[y+1]);
      y+=2;
    }
  }
  url[x] = '\0';
}

std::string urlencode(const std::string &s)
{
  //RFC 3986 section 2.3 Unreserved Characters (January 2005)
  const std::string unreserved = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_.~";

  std::string escaped="";
  for(size_t i=0; i<s.length(); i++)
    {
      if (unreserved.find_first_of(s[i]) != std::string::npos)
        {
	  escaped.push_back(s[i]);
        }
      else
        {
	  escaped.append("%");
	  char buf[3];
	  sprintf(buf, "%.2X", s[i]);
	  escaped.append(buf);
        }
    }
  return escaped;
}

void plustospace(char *str) {
  register int x;

  for(x=0;str[x];x++) if(str[x] == '+') str[x] = ' ';
}

int rind(char *s, char c) {
  register int x;
  for(x=(int)strlen(s) - 1;x != -1; x--)
    if(s[x] == c) return x;
  return -1;
}

int get_line(char *s, int n, FILE *f) { // "getline" conflicts with stdio.h
  register int i=0;

  while(1) {
    s[i] = (char)fgetc(f);

    if(s[i] == CR)
      s[i] = fgetc(f);

    if((s[i] == 0x4) || (s[i] == LF) || (i == (n-1))) {
      s[i] = '\0';
      return (feof(f) ? 1 : 0);
    }
    ++i;
  }
}

//
// case insensitive strstr
//
char *strstri( const char *str, const char *strCharSet ) {
  char *lstr = strdup(str);
  char *lstrCharSet = strdup(strCharSet);
  char *result;
  char *cp;
  // strlwr
  for (cp=lstr;*cp;cp++) {
    *cp = tolower(*cp);
  }
  for (cp=lstrCharSet;*cp;cp++) {
    *cp = tolower(*cp);
  }
  // search
  result = strstr(lstr,lstrCharSet);
  if (result) {
    // result pointer must be to original data
    result = (char *)str+(result-lstr);
  }
  // no leaks!
  free(lstr);
  free(lstrCharSet);
  return result;
}

//
// case insensitive strstr, finds rightmost match
//
char *strstrir( const char *str, const char *strCharSet ) {
  char *search = strstri( str, strCharSet );
  char *result = search;
  while (search && *search) {
    search = strstri(result+1,strCharSet);
    if (search) {
      result = search;
    }
  }
  return result;
}

// copy one file to another, return 0 on success
int copy_file(const char *fromName, const char *toName) {
  int result = 0;
  FILE *from = fopen(fromName,"rb");
  if (from) {
    FILE *to = fopen(toName,"wb");
    if (to) {
      send_fd(from,to);
      fclose(to);
    }  else {
      result = errno;
    }
    fclose(from);
  } else {
    result = errno;
  }
  return result;
}

void send_fd(FILE *f, FILE *fd)
{
  char c[0xfff];
  int n;
  while ( (n=(int)fread(c,1,sizeof(c),f)) ) {
    if (n!=fwrite(c,1,n,fd)) {
      puts("write error!");
      break;
    }
  }
}

string XMLEscape(const string& s)
{
  string ret;
  for (size_t i = 0; i < s.length(); i++)
    {
      char ch = s.at(i);
      switch (ch)
        {
        case '<':
	  ret.append("&lt;");
	  break;
        case '>':
	  ret.append("&gt;");
	  break;
        case '&':
	  ret.append("&amp;");
	  break;
        case '"':
	  ret.append("&quot;");
	  break;
        case '\'':
	  ret.append("&apos;");
	  break;
        case '\n':
	  ret.append(" ");
	  break;
        case '\r':
	  ret.append(" ");
	  if (i+1 < s.length() && s.at(i+1) == '\n')
	    i++;
	  break;
        default:
	  ret.append(1, ch);
	  break;
        }
    }
  return ret;
}

/** Basic routine to escape spaces in a string with
 *  the "\" character.  This can be used to create
 *  a POSIX path/filename from a string.  The escaped
 *  space path/filename can then be used in a "system()"
 *  call.
 */
char * escape_spaces (char * inStr)
{
  int i;
  int space_count = 0;
  int inStrLen = 0;
  char * pinStr = inStr;
  char * outStr;
  char * poutStr;

  while ( ( pinStr = strchr(pinStr,' ') ) )
    {
      pinStr++;
      space_count++;
    }

  if ( space_count > 0 ) 
    {
      poutStr = outStr = 
        (char *) malloc(sizeof(char) * (strlen(inStr) + 1 + space_count));

      if ( outStr == NULL )
	{
	  printf("Error mallocing memory in escape_spaces()!\n");
	  return(0);  //DCT \0 == 0 stop clang warning
	}

      for ( i = 0; inStr[i]; i++ ) 
	{
	  if ( inStr[i] == ' ' )
	    {
	      *poutStr++ = '\\';
	    }
	  *poutStr++ = inStr[i];
	}
      *poutStr = '\0';

      return(outStr);

    }

  return(inStr);

}


// locate any leading cygwin path stuff and make it windowsy instead
void force_unCygwinify(std::string &filename) {
  char *copy = strdup(filename.c_str());
  force_unCygwinify(copy);
  filename = copy;
  free(copy);
}
void force_unCygwinify(char *filename) {
  char *cp;
  for (cp=filename;cp&&*cp;) {
    char *cp2=strstr(cp,"/usr/bin/");
    cp2=cp2?cp2:strstr(cp,"\\usr\\bin\\"); // windows commandline hijinx
    if (cp2) { 
      memmove(cp2,cp2+9,strlen(cp2+9)+1);
    } else {
      break;
    }
    cp++;
  }
  // drop any single quotes
  for (cp=filename;cp&&*cp;) {
    if ('\''==*cp) {
      memmove(cp,cp+1,strlen(cp));
    } else {
      cp++;
    }
  }
  for (cp=filename;cp&&*cp;) {
    char *cp2=strstr(cp,"/cygdrive/");
    cp2=cp2?cp2:strstr(cp,"\\cygdrive\\"); // windows commandline hijinx
    if (cp2) { // convert /cygdrive/c to c:
      cp=cp2;
      cp[0] = cp[10];
      cp[1] = ':';
      memmove(cp+2,cp+11,strlen(cp+11)+1);
    }
    if (cp[0] && (':'==cp[1])) {
      cp[0] = tolower(cp[0]); // drive letters in lowercase
    }
    cp++;
  }
  for (cp=filename;*cp;cp++) {
    if ('\\'==*cp) { // we want consistent path separators
      *cp = '/';
    }
  }
}

#ifdef WINDOWS_NATIVE // change "/" to "\" 
void fixpath(char *in) {
  // mingw drive spec?
  if (('/'==*in) && isalpha(*(in+1)) && ('/'==*(in+2))) {
    *in = *(in+1);
    *(in+1) = ':';
  }
  while (*in) {
    if ('/'==*in) {
      if ('/'==*(in+1)) { // leave URLs alone - http:// etc
	while (*in && (' '!=*in)) {
	  in++;
	}
      } else {
	*in++='\\';
      }
    } else {
      in++;
    }
  }
}
void fixpath(std::string &in) {
  char *tmp=strdup(in.c_str());
  fixpath(tmp);
  in = tmp;
  free(tmp);
}
#else
#include <glob.h>
#include <stdexcept>
#include <dirent.h>
#endif

// helpful functions for determining filetypes

// check to see if fname ends in .gz
// return ptr to .gz if found, else NULL
const char* isDotGZ(const std::string &fname) {
  return isDotGZ(fname.c_str());
}
const char* isDotGZ(const char *fname) {
  const char *gzExt = ".gz";
  const char *gz = fname?strstrir(fname,gzExt):NULL;
  if (gz && *(gz+strlen(gzExt))) { // found .gz at end of fname?
    gz = NULL; // wasn't the last thing in fname
  }
  return gz;
}

// look for the indicated .ext, or .ext.gz, in the filename
// return pointer to .ext if found, else NULL
const char *hasExtOrExtDotGZ(const char *fname, const char *ext) {
  const char *result = fname?strstrir(fname,ext):NULL; // find rightmost occurance of .ext
  if (result && *(result+strlen(ext))) { // found .ext, but not at end, maybe .ext.gz?
    if (!isDotGZ(fname)) {
      result = NULL;
    }
  }
  return result;
}

#include "Util/RACI/RACI.h"
#include "zlib.h"

// run gzip on this file unless it's already gzipped
int do_gzip(const char *fname) {
  std::string str(fname);
  return do_gzip(str);
}
int do_gzip(std::string &fname) {
  bool err = false;
  bool bNeedsCompression;
  { // scope the constructor/destructor so file handle isn't locked after we use it here
    RACI test(fname.c_str());
    bNeedsCompression = (test.is_open() && !test.getCompressionType());
  }
  if (bNeedsCompression) { // it exists and it's not gzipped
    const char *dotGZ = isDotGZ(fname.c_str());
    if (dotGZ) { // it already has .gz ext, remove for now
      std::string uncompressed_fname(fname);
      uncompressed_fname.erase(fname.length()-strlen(dotGZ));
      struct stat statbuf;
      if (!stat(uncompressed_fname.c_str(),&statbuf)) {
	std::cout << "warning: removing old copy of " << uncompressed_fname << " to create " << fname << std::endl; 
	remove(uncompressed_fname.c_str());
      }
      safe_rename(fname.c_str(),uncompressed_fname.c_str());
      fname = uncompressed_fname;
    } else { // no .gz ext yet, but does one exist?
      std::string compressed_fname(fname);
      compressed_fname += ".gz";
      struct stat statbuf;
      if (!stat(compressed_fname.c_str(),&statbuf)) {
	std::cout << "warning: removing old copy of " << compressed_fname << " so we can compress " << fname << std::endl; 
	remove(compressed_fname.c_str());
      }
    }
    // now do gzip (using zlib to avoid external gzip app dependency)
    std::string compressed_fname(fname);
    compressed_fname += ".gz";
    gzFile out = gzopen(compressed_fname.c_str(),"wb");
    if (out) {
      gzsetparams(out,Z_BEST_COMPRESSION,Z_DEFAULT_STRATEGY);
      FILE *in = fopen(fname.c_str(),"rb");
      if (in) {
#define ZBUFLEN 32768
	char *buf = new char[ZBUFLEN];
	if (buf) {
	  int n;
	  while ((n=(int)fread(buf,1,ZBUFLEN,in))>0) {
	    gzwrite(out,buf,(unsigned) n);
	  }
	  delete[] buf;
	} else {
	  err = true;
	}
	fclose(in);
      } else {
	err = true;
      }
      gzclose(out);
    } else {
      err = true;
    }
    if (err) {
      std::cout << "error: unable to gzip file " << fname << std::endl;
    } else {
      unlink(fname);
    }
  }
  return err?-1:0;
}

// helper func
bool comparelen( const char * elem1, const char * elem2 ){
  return strlen(elem1) > strlen(elem2); // longest first
}
static const char *hasExt(const char *fname, std::vector<const char *> &exts) {
  std::sort(exts.begin(),exts.end(),comparelen); // sort longest to shortest
  for (size_t i=0;i<exts.size();i++) {
    const char *result;
    if ( (result = hasExtOrExtDotGZ(fname,exts[i])) ) {
      return result;
    }
  }
  return NULL;
}

// accept either .xml or .pepXML as valid pepXML filename extensions
// also accept these extensions with .gz appended 
// return NULL if invalid, else return pointer to .ext
const char *hasValidPepXMLFilenameExt(const char *fname) {
  if (hasValidProtXMLFilenameExt(fname)) { // watch for .prot.xml
    return NULL;
  }
  std::vector<const char *>exts; // a list of possible extenstions
  exts.push_back(get_pepxml_dot_ext()); // user config
  exts.push_back(DEFAULT_PEPXML_FILENAME_DOTEXT); // usually .pep.xml
  exts.push_back(".xml"); // pre Jan 2008
  exts.push_back(".pepXML");// never really took off...
  return hasExt(fname,exts);
}

// accept either -prot.xml or .protXML as valid protXML filename extensions,
// also accept these extensions with .gz appended 
// return NULL if invalid, else return pointer to .ext
const char *hasValidProtXMLFilenameExt(const char *fname) {
  std::vector<const char *>exts; // a list of possible extenstions
  exts.push_back(get_protxml_dot_ext()); // user config
  exts.push_back("-prot.xml"); // pre Jan 2008
  exts.push_back(".protXML");// never really took off...
  exts.push_back(DEFAULT_PROTXML_FILENAME_DOTEXT);
  return hasExt(fname,exts);
}

//
//
// is c a path seperator for linux or windows?
//
int isPathSeperator(char c) {
  return (c=='\\')||(c=='/');
}

// convert rel path to abs path - caller must free() result
char *makeFullPath(const char *fname) {
  char *result;
  if (isAbsolutePath(fname)) {
    result = strdup(fname);
  } else {
    char szFile[SIZE_FILE];
    int len;
    safepath_getcwd(szFile, SIZE_FILE);
    len = (int)strlen(szFile);
    strncat(szFile, "/", SIZE_FILE-len++);
    strncat(szFile, fname, SIZE_FILE-len);
    result = strdup(szFile);
  }
  unCygwinify(result); // no effect in cygwin builds
  return result;
}

void makeFullPath(std::string &fname) {
  char *str = makeFullPath(fname.c_str());
  fname = str;
  free(str);
}


const char *findRightmostPathSeperator_const(const char *path) { // return pointer to rightmost / or \ .
  const char *result = path+strlen(path);
  while (result-->path) {
    if (isPathSeperator(*result)) {
      return result;
    }
  }
  return NULL; // no match
}

int findRightmostPathSeperator(const std::string &str) {
  int slashPos = (int)str.find_last_of('/');
  if (slashPos == string::npos) {
    slashPos = (int)str.find_last_of('\\');
  }
  return slashPos;
}

char *findRightmostPathSeperator(char *path) { // return pointer to rightmost / or \ .
  return (char *)findRightmostPathSeperator_const(path);
}

//
// compare the paths, ignoring filename and / vs \
//
int pathcmp(const char *path1, const char *path2) {
  const char *end1 = findRightmostPathSeperator_const(path1);
  const char *end2 = findRightmostPathSeperator_const(path2);
  if ((!end1) && (!end2)) {
    return 0; // no path means same path
  }
  if ((!end1) || (!end2)) {
    return 1; // no path on one but not the other 
  }
  while ((path1<end1) && (path2<end2)) {
    if (*path1!=*path2) {
      if (!(isPathSeperator(*path1) &&  // / or \ are the same to us
	    isPathSeperator(*path2))) {
	return(*path1-*path2);
      }
    }
    path1++;
    path2++;
  }
  return 0;
}

//
// local helper func for cached env read, with path separators cleaned up
//
static char * check_env_var(const char *env,char *buf, int buflen, bool bIsPath, const char *defaultval) {
  if (!buf[0]) { // first access?
    const char *e=getenv(env);
    if (e) {
      strncpy(buf,e,buflen);
    } else {
      buf[0] = '\0'; // mark as failed
    }
    if ((strlen(buf) <= 0) && defaultval) { // no env var value, use default if provided
      strncpy(buf,defaultval,buflen);
    }
    // tidy up path seps if needed
    if (bIsPath && (buf[0]>0)) {
      char *cp;
      for (cp=buf;*cp;cp++) {
	if ('\\'==*cp) {
	  *cp = '/';
	}
      }
      // add trailing path sep if needed
      if ((cp > buf) && ('/'!=*(cp-1))) {
	if ((cp+1-buf) < buflen) {
	  *cp++ = '/';
	  *cp = 0;
	}
      }
    }
  }
  return (buf[0]>0)?buf:NULL; // buf[0]<0 means env var does not exist
}


// Pre TPP5 config/filesystem path stuff
static char true_wsr[1024]={0};
const char *getTrueWebserverRoot() { // when using tmpdir, webserver root may be messed with
  if (!true_wsr[0]) { // first access?
    if (!check_env_var("WEBSERVER_ROOT_TPP",true_wsr,sizeof(true_wsr),true)) { 
      // no env var set
      const char *wsr = getWebserverRoot();
      if (wsr) {
	strncpy(true_wsr,wsr,sizeof(true_wsr));
      }
      // pass this on to any child processes
      std::string twsr = "WEBSERVER_ROOT_TPP=";
      twsr += true_wsr;
      putenv(strdup(twsr.c_str())); // intentional leak, freeing may affect env
    }
  }
  return true_wsr;
}

//
// DEPRECATED: use getDataPath();
//
// get the WEBSERVER_ROOT env var, with path separators cleaned up
//
static char wsroot[1024]={0};
const char *getWebserverRoot() {
  // return check_env_var("WEBSERVER_ROOT",wsroot,sizeof(wsroot),true);
  return getDataPath();
}

//
// set the WEBSERVER_ROOT variable, preserve previous value
// this is done in the case where xinteract wants to do its work in a 
// local tmpdir instead of working across a wire to a webserver root which 
// is actually a network share, or AWS S3 bucket, etc
//
static std::vector<std::string> wsrootstack;
void pushWebserverRoot(const char *path) {
  getTrueWebserverRoot(); // side effect avoidance: cache this value before we mess with webserver root
  wsrootstack.push_back(getWebserverRoot());
  wsroot[0]=0;
  std::string e("WEBSERVER_ROOT=");
  e+=path;
  putenv(strdup(e.c_str())); // intentional leak, otherwise env can change on free
}

//
// restore previous WEBSERVER_ROOT
//
void popWebserverRoot() {
  wsroot[0]=0;
  std::string e("WEBSERVER_ROOT=");
  e+=wsrootstack[wsrootstack.size()-1];
  wsrootstack.pop_back();
  putenv(strdup(e.c_str())); // intentional leak, otherwise env can change on free
}


//
// get the WEBSERVER_TMP env var, with path separators cleaned up
// make sure it exists, too
//
static char wst[1024]={0};
bool bCheckedWebserverTmpPath=false;
const char *getWebserverTmpPath() {
  const char *result = check_env_var("WEBSERVER_TMP",wst,sizeof(wst),true);
  if (result && !bCheckedWebserverTmpPath) { // verify that tmpdir exists, create if needed
    bCheckedWebserverTmpPath = true;
    cygwinify(wst,sizeof(wst)); // no effect in non-cygwin builds
    struct stat statbuf;
    std::string tmpdir(result);
    tmpdir = tmpdir.substr(0,tmpdir.length()-1); // kill the trailing slash
    if (stat(tmpdir.c_str(),&statbuf)) { // maybe it doesn't exist yet?
      std::string cmd("mkdir -p ");
      cmd += tmpdir;
      int ret = tpplib_system(cmd.c_str());
      if (stat(tmpdir.c_str(),&statbuf)) {
	cout << "warning: the webserver temporary directory " << tmpdir << 
	  " specified by WEBSERVER_TMP does not exist and could not be created.  " <<
	  "Tempfiles will created in data directories instead.";
	wst[0]='\0'; // mark it as failed
	result = NULL;
      }
    }		 
  }
  return result;
}

//
// set the WEBSERVER_TMP env var for this process and any children it spawns
//
static std::vector<std::string> wstmpstack;
void pushWebserverTmpPath(const char *path) {
  wstmpstack.push_back(getWebserverTmpPath());
  bCheckedWebserverTmpPath = false;
  wst[0]=0;
  std::string e("WEBSERVER_TMP=");
  e+=path;
  putenv(strdup(e.c_str())); // intentional leak, otherwise env can change on free
  getWebserverTmpPath(); // reset
}

//
// restore previous WEBSERVER_TMP
//
void popWebserverTmp() {
  wst[0]=0;
  bCheckedWebserverTmpPath = false;
  std::string e("WEBSERVER_TMP=");
  e+=wstmpstack[wstmpstack.size()-1];
  wstmpstack.pop_back();
  putenv(strdup(e.c_str())); // intentional leak, otherwise env can change on free
  getWebserverTmpPath(); // reset
}

//
// return a copy of the input filename with a webserver root path
//
std::string translate_relative_webserver_root_path_to_absolute_filesystem_path(const char *path) {
  char *result;
  char *pathcopy = strdup(path);
  char *szWebserverRoot=strdup(getWebserverRoot());
  const char *pStr;
  unCygwinify(pathcopy); // no effect in cygwin build
  unCygwinify(szWebserverRoot); // no effect in cygwin build
  if (!findRightmostPathSeperator(pathcopy)) { // no path info
    char *bbuf = (char *)malloc(1025+strlen(pathcopy));
    safepath_getcwd(bbuf,1023);
    strcat(bbuf,"/");
    strcat(bbuf,pathcopy);
    free(pathcopy);
    pathcopy = bbuf;
  }
  pStr = strstri(pathcopy, szWebserverRoot);
  if (pStr==pathcopy) {
    free(szWebserverRoot);
    free(pathcopy);
    return strdup(path); // already in terms of filesystem
  }
  result = (char *)malloc(strlen(szWebserverRoot)+strlen(pathcopy)+1);
  sprintf(result, "%s%s", szWebserverRoot, pathcopy);
  free(szWebserverRoot);
  free(pathcopy);
  std::string strResult(result);
  free(result);
  return strResult;
}

//
// in-place replacement of path with webserver's tmp dir, if any
// so /blarg/foo becomes /inetpub/wwwroot/tmp/foo
void replace_path_with_webserver_tmp(std::string &path) { // write this in tmpdir if we have one
  const char *szWebserverTmp=getWebserverTmpPath();
  if (szWebserverTmp) { // set up to call into char* implementation
    int pathmaxlen = (int)(strlen(szWebserverTmp)+path.length())+2;
    char *pathbuf = new char[pathmaxlen];
    strcpy(pathbuf,path.c_str());
    replace_path_with_webserver_tmp(pathbuf,pathmaxlen);
    path = pathbuf;
    delete[] pathbuf;
  }
}

//
// in-place replacement of path with webserver's tmp dir, if any
// so /blarg/foo becomes /inetpub/wwwroot/tmp/foo
// do nothing if incoming path is a subdir of webserver's tmp path
void replace_path_with_webserver_tmp(char *path,int pathmaxlen) { // write this in tmpdir if we have one
  const char *szWebserverTmp=getWebserverTmpPath();
  if (szWebserverTmp) {
    if (strncasecmp(path,szWebserverTmp,strlen(szWebserverTmp))) {
      replace_path(path,pathmaxlen,szWebserverTmp);
    }
  }
}


//
// in-place replacement of path with newpath
// so /blarg/foo becomes newpath/foo
//
void replace_path(char *path,int pathmaxlen,const char *newpath) {
  char *fname = findRightmostPathSeperator(path);
  fname = strdup(fname?(fname+1):path);
  strncpy(path, newpath, pathmaxlen);
  strncat(path,fname,pathmaxlen-strlen(path));
  free(fname);
}

void replace_path(std::string &path,const char *newpath) {
  int len;
  char *tmp = (char *)malloc(len = (int)(path.length()+strlen(newpath)+2));
  replace_path(tmp,len,newpath);
  path = tmp;
  free(tmp);
}


//
// remove the webserver root portion of the input path, if any
// so /inetpub/wwwroot/foo/blah becomes /foo/blah
//
void translate_absolute_filesystem_path_to_relative_webserver_root_path(std::string &path) {
  char *p=strdup(path.c_str());
  translate_absolute_filesystem_path_to_relative_webserver_root_path(p);
  path = p;
  free(p);
}
void translate_absolute_filesystem_path_to_relative_webserver_root_path(char *path) {
  const char *szWebserverRoot = getWebserverRoot();
  if (!szWebserverRoot) return;
  char *tmpStr = strstri(path, szWebserverRoot);
  if (tmpStr != NULL) {
    tmpStr += strlen(szWebserverRoot);
    char *relpath = strdup(tmpStr);
    if ('/'!=*relpath) {
      *path++ = '/';
    }
    strcpy(path, relpath);
    free(relpath);
  }
}


//
// return a copy of the input filename with the filesystem webserver root prepended
// caller must free() the result
//
static char *prepend_helper(const char *path,const char *prependPath) {
  char *result;
  char *pathcopy = strdup(path?path:"");
  char *prepend=strdup(prependPath?prependPath:"");
  char *pStr;
  unCygwinify(pathcopy); // no effect in cygwin build
  unCygwinify(prepend); // no effect in cygwin build
  pStr = strstri(pathcopy, prepend);
  if (pStr != NULL) {
    free(prepend);
    free(pathcopy);
    return strdup(path?path:""); // presumably set up already
  }
  result = (char *)malloc(strlen(prepend)+strlen(pathcopy)+1);
  sprintf(result, "%s%s", prepend, pathcopy); 
  free(prepend);
  free(pathcopy);
  for (pStr=result;*pStr;) { // convert // to /
    if (('/'==*pStr)&&(*pStr==*(pStr+1))) {
      memmove(pStr,pStr+1,strlen(pStr));
    } else {
      pStr++;
    }
  }
  return result;
}

//
// return a copy of the input filename with the filesystem webserver root prepended
// caller must free() the result
//
std::string prepend_webserver_root(const char *path) {
  return prepend_helper(path,getWebserverRoot());
}

// do we need to fix up the path at all?
// returns:
//   path of file that can be opened,
//   or just a copy of input path
std::string resolve_root(const char *path) {
  struct stat statbuf;
  if (stat(path,&statbuf)) { // didn't find it
    std::string str = prepend_helper(path,getTrueWebserverRoot()); // is it in the wwwroot area?
    if (!stat(str.c_str(),&statbuf)) {
      // that worked
      return str;
    }
#ifdef WINDOWS_NATIVE
    str = prepend_helper(path,"c:/cygwin/");
    if (!stat(str.c_str(),&statbuf)) {
      // that worked
      return str;
    }
#endif
    // try cwd
    char *buf = strdup(str.c_str());
    char *slash = findRightmostPathSeperator(buf);
    if (slash) {
      memmove(buf,slash+1,strlen(slash)); 
      if (!stat(buf,&statbuf)) {
	// that worked
	char *bbuf = (char *)malloc(1025+strlen(buf));
	safepath_getcwd(bbuf,1023);
	strcat(bbuf,"/");
	strcat(bbuf,buf);
	std::string result(bbuf);
	free(buf);
	free(bbuf);
	return result;
      }
    }
    free(buf);
  }
  return std::string(path);
}

// fix up path's root dir if needed
void resolve_root_dir(char *path,int buflen) {
  unCygwinify(path); // no effect in cygwin builds
  char *copy = strdup(path);
  char *slash = strchr(copy+('/'==*copy),'/');
  if (slash) {
    struct stat s;
    *slash = 0;
    if (stat(copy,&s)) {
      std::string str = resolve_root(copy);
      const char *ccopy = str.c_str();
      if (!stat(ccopy,&s)) {
	if ((int)(strlen(ccopy)+strlen(path)+1) < buflen) {
	  strcpy(path,ccopy);
	  strcat(path,"/");
	  strcat(path,slash+1);
	}
      }
    }
  }
  free(copy);
}


//
// get the WEBSERVER_ROOT env var, with path separators cleaned up
//
static char wsurl[1024]={0};
const char *getWebserverUrl() {
  if (!wsurl[0]) {
    const char *e=getenv("WEBSERVER_URL");
    if (e) {
      int len;
      strncpy(wsurl,e,sizeof(wsurl));
      len = (int)strlen(wsurl);
      if (wsurl[len-1] != '/') {
	wsurl[len] = '/';
	wsurl[len+1] = '\0';
      }
    }
    else {
      e=getenv("SERVER_NAME");	
      if (e) {
	sprintf(wsurl,"%s%s/", "http://", e);
      }
    }
  }
  return wsurl[0]?wsurl:NULL;
}

//
// getcwd with path separators cleaned up
//
char *safepath_getcwd(char *buf, int buflen) {
  char *cp;
  char *result;
  result = getcwd(buf,buflen);
  for (cp=buf;*cp;cp++) {
    if ('\\'==*cp) {
      *cp = '/';
    }
  }
  return result;
}

std::string safepath_getcwd() {
  char buf[1024];
  std::string result = safepath_getcwd(buf,sizeof(buf));
  return result;
}

//
// return nonzero if fname appears to be absolute
// does not check for actual existence of file
//
int isAbsolutePath(const char *fname) {
  return
#ifdef WINDOWS_NATIVE
    (*fname && (':'==fname[1]) && (('\\'==fname[2])||('/'==fname[2]))) ||
#endif
    isPathSeperator(*fname);
}
int isAbsolutePath(const std::string &fname) {
  return isAbsolutePath(fname.c_str());
}
//
// snprintf may not nullterm on overflow, so use this
//
void safe_snprintf(char *buf, int buflen, const char *format, /* args */ ...)
{
  va_list va;
  va_start(va, format);
  if (buflen) {
#ifdef _MSC_VER
    _vsnprintf(buf, buflen, format, va);
#else
    vsnprintf(buf, buflen, format, va);
#endif
    buf[buflen-1] = 0;
  } else {
    *buf = 0;
  }
}

//
// does path end with a path seperator?
//
bool endsWithPathSeperator(const char *path) {
  return path && *path &&isPathSeperator(path[strlen(path)-1]);
}


// get consistent path seperators, optionally check for existence
// attempt to clean up confused pathnames
int fixPath(char *path, int expectExist) {
  struct stat s;
  int statval;
  int success=0;
  if (!(path && *path)) {
    return !expectExist;
  }
  unCygwinify(path);  // no effect in non-cygwin builds
  char *backup = strdup(path);
  if (isPathSeperator(path[strlen(path)-1])) { // MSVC stat gets confused by trailing /
    char *up = strdup(path);
    up[strlen(up)-1] = 0;
    statval = stat(up,&s);
    free(up);
  } else {
    statval = stat(path,&s);
  }

  if (statval && expectExist) {
    // maybe a confused basename, as in some Phenyx data (foo.mzXML.mzXML)
    char *up = strdup(path); 
    char *ext = strrchr(up,'.');
    if (ext) {
      char * ext2;
      *ext++ = 0; // now reads foo.mzXML
      ext2 = strrchr(up,'.');
      if (ext2 && !strcmp(ext,++ext2)) {
	if ((!stat(up,&s)) && (S_IFREG & s.st_mode))  {
	  statval = 0;  // it exists as a regular file foo.mzXML
	  strcpy(path, up);
	}
      }
    }
    free(up);
  }
  success = !statval;
  if (statval && expectExist) {
    char *up = strdup(path); 
    // maybe a confused basename, as from mascot2xml
    // from ISB/bob/F02341/greeble.mzXML try ISB/bob/greeble.mzXML
    char *name = findRightmostPathSeperator(up);
    if (name) {
      char *slash;
      *name++ = 0;
      slash = findRightmostPathSeperator(up);
      if (slash) {
	strcpy(slash+1,name);
	if (!stat(up,&s)) {
	  strcpy(path,up);
	  success = 1;
	}
      }
    }
    free(up);
    if (!success) { 
      // or could be a path
      // from ISB/bob/F02341/ try ISB/bob/
      if (isPathSeperator(path[strlen(path)-1])) {
	char *up = strdup(path);
	char *slash;
	up[strlen(up)-1] = 0;
	slash = findRightmostPathSeperator(up);
	if (slash) {
	  *++slash = 0;
	  if ((!stat(up,&s)) && (S_IFDIR & s.st_mode)) {
	    // now that looks like a path
	    strcpy(path,up);
	    success = 1;
	  }
	}
	free(up);
      }
    }
    if (!success) { 
      // or could be a path
      // from ISB/bob/F02341 try ISB/bob
      char *up = strdup(path);
      char *slash;
      up[strlen(up)-1] = 0;
      slash = findRightmostPathSeperator(up);
      if (slash) {
	*++slash = 0;
	if ((!stat(up,&s)) && (S_IFDIR & s.st_mode))  {
	  // now that looks like a path
	  strcpy(path,up);
	  success = 1;
	}
      }
      free(up);
    }
  }
  if (!expectExist) {
    if (!success) { 
      // or could be a file not yet written - 
      // does its path exist?
      // from ISB/bob/F02341/out try ISB/bob/out
      // but watch for ISB/bob/out 
      char *up = strdup(path);
      char *name = findRightmostPathSeperator(up);
      if (name) {
	*name++ = 0;
	if (!((!stat(up,&s)) && (S_IFDIR & s.st_mode))) {        
	  // example ISB/bob/F02341 is not an existing directory
	  char *slash = findRightmostPathSeperator(up);
	  if (slash) {
	    *slash = 0;
	    if ((!stat(up,&s)) && (S_IFDIR & s.st_mode)) {
	      // example ISB/bob is a directory
	      strcpy(slash,"/");
	      strcpy(slash+1,name);
	      strcpy(path,up);
	      success = 1;
	    }
	  }
	}
      }
      free(up);
    }
  }
  if (expectExist && !success) { // restore (uncigwinified) input
    strcpy(path,backup);
  }
  free(backup);
  return expectExist?(success!=0):0;
}

//
// format a png file reference to either use the
// show_tmp_pngfile.pl CGI or not (using the CGI causes
// the png file to be unlinked after use)
//
// example: 
// pngfile="foo.png";
//
// sprintf(fp, "<IMG SRC=\"%s\"/>", makeTmpPNGFileSrcRef(pngfile);
// output would be <IMG SRC="/tpp-bin/show_tmp_pngfile.pl?foo.png"/>
//
// sprintf(fp, "<IMG SRC=\"%s\"/>", makeTmpPNGFileSrcRef(pngfile,TRUE);
// output would be <IMG SRC="foo.png"/>
//
std::string makeTmpPNGFileSrcRef(const char *filename, bool persist) {
  std::string result;
  if (!persist) { // 
    result = getCgiUrl();
    if (!endsWithPathSeperator(result.c_str())) {
      result += '/';
    }
    result+="show_tmp_pngfile.pl?file=";
  }
  result += filename;
  return result;
}

// get the canonical pepxml filename extension, including period
static char pepxmldotext[1024]={0};
const char *get_pepxml_dot_ext() {
  if (!pepxmldotext[0]) { // first pass
    check_env_var("PEPXML_EXT",pepxmldotext,sizeof(pepxmldotext),false,DEFAULT_PEPXML_FILENAME_DOTEXT);
    if (!strchr(pepxmldotext,'.')) { // user didn't put a dot in there?
      memmove(pepxmldotext+1,pepxmldotext,strlen(pepxmldotext)+1);
      pepxmldotext[0] = '.';
    }
  }
  return pepxmldotext;
}

// get the canonical protxml filename extension, including period
static char protxmldotext[1024]={0};
const char *get_protxml_dot_ext() {
  if (!protxmldotext[0]) { // first pass
    check_env_var("PROTXML_EXT",protxmldotext,sizeof(protxmldotext),false,DEFAULT_PROTXML_FILENAME_DOTEXT);
    if (!strchr(protxmldotext,'.')) { // user didn't put a dot in there?
      memmove(protxmldotext+1,protxmldotext,strlen(protxmldotext)+1);
      protxmldotext[0] = '.';
    }
  }
  return protxmldotext;
}

//
// construct a tempfile name, possibly in the tmp dir if so configured
//
std::string make_tmpfile_name(const char *basis) {
  std::string outfile(basis);
  outfile += ".tmp.XXXXXX";
  replace_path_with_webserver_tmp(outfile); // do this in designated tmpdir, if any
  safe_fclose(FILE_mkstemp(outfile)); // replace XXXXXX with unique string
  return outfile;
}

//
// construct a tempfile name, in the current directory
//
std::string make_temp_name(const char *basis) {
  std::string outfile(basis);
  outfile += ".tmp.XXXXXX";
  //replace_path_with_webserver_tmp(outfile); // do this in designated tmpdir, if any
  safe_fclose(FILE_mkstemp(outfile)); // replace XXXXXX with unique string
  return outfile;
}

//
// like mkstemp, but returns open FILE * instead of file handle
// replaces XXXXXX with some unique string in place
//
FILE *FILE_mkstemp(char *fname_XXXXXX) {
#ifdef WINDOWS_NATIVE		
  mktemp(fname_XXXXXX);
#else
  int fd = mkstemp(fname_XXXXXX);
  close(fd);
#endif
  return fopen(fname_XXXXXX,"wb");
}

FILE *FILE_mkstemp(std::string &fname_XXXXXX) {
  char *tmp = strdup(fname_XXXXXX.c_str());
  FILE *result = FILE_mkstemp(tmp);
  fname_XXXXXX = tmp;
  free(tmp);
  return result;
}

// standard rename() doesn't work across filesystem boundaries, which we may need to do
int safe_rename(const char *from, const char *to) {
  int result;

  // First try rename because some systems don't have mv!!!
  result= rename( from , to );
  if ( result == 0 )
    return result;

  
  const char *quot=getCmdlineQuoteChar(); // get a system-appropriate quote char
  std::string fromto = quot;
  fromto += from;
  fromto += quot;
  fromto += " ";
  fromto += quot;
  fromto += to;
  fromto += quot;
  // note: using system("mv -f from to") instead of library call rename(from, to)
  // as the library call does not work across different filesystems
  std::string cmd = "mv -f ";
  cmd += fromto;
  result = tpplib_system(cmd.c_str());
  return result;
}

// unlink with verbose error check
int verified_unlink(const char *fname) {
  int result = unlink(fname);
  if (result) {
    // don't bark about nonexistent files
    int save_errno = errno;
    struct stat statbuf;
    if (!stat(fname,&statbuf)) { 
      printf("could not delete %s, error %d \"%s\"\n",fname,save_errno,strerror(save_errno));
    }
  }
  return result;
}

// rmdir with verbose error check
int verified_rmdir(const char *fname) {
  int result = rmdir(fname);
  if (result) {
    printf("could not remove dir %s, error %d \"%s\"\n",fname,errno,strerror(errno));
  }
  return result;
}

// chdir with verbose error check
int verified_chdir(const char *fname) {
  int result = chdir(fname);
  if (result) {
    printf("could not chdir to %s, error %d \"%s\"\n",fname,errno,strerror(errno));
  }
  return result;
}

// system() with verbose error check
int verified_system(const char *cmd) {
  int result = tpplib_system(cmd);
  if (result) {
    printf("command \"%s\" failed with exit code %d\n",cmd,result);
  }
  return result;
}

// determine if named file is gzipped, if so unzip to a tempfile
// returns: tmpfile name, or original file name if not gzipped
// optional maxchar arg limits effect of huge files when you just want to see header
std::string uncompress_to_tmpfile(const char *fname,int maxchar) {
  std::string result;
  bool err = false;
  if (!fname || !*fname) {
    result = ""; // empty
  } else if (isDotGZ(fname)) {
    // check to see if it's really compressed or just so named
    { // scope the constructor/destructor so file handle isn't locked after we use it here
      RACI test(fname);
      if (test.is_open() && !test.getCompressionType()) {
	result = fname; // doesn't need unzip
	return result; // done here
      }
    }
    // decompress to a tmpfile
    result = make_tmpfile_name(fname);
    gzFile z = gzopen(fname, "rb");
    if (z) {
      char *buf = new char[ZBUFLEN];
      FILE *tmp = fopen(result.c_str(),"wb");
      if (!tmp || !buf) {
	err = true;
	result = fname; // unable to process
      } else {
	int r;
	int wtotal=0;
	while ((r=gzread(z,buf,ZBUFLEN))>0) {
	  int w = (int)fwrite(buf,1,r,tmp);
	  if (r!=w) {
	    err = true;
	    result = fname; // unable to process
	    break;
	  }
	  if (maxchar>0) { // only output the first maxchar bytes
	    wtotal+=w;
	    if (wtotal > maxchar) {
	      break;
	    }
	  }
	}
	fclose(tmp);
      }
      gzclose(z);
      delete[] buf;
    } else {
      err = true;
      result = fname; // unable to process
    }
  } else {
    result = fname; // doesn't need unzip
  }
  if (err) {
    std::cout << "error: could not unzip file " << fname << std::endl;
  }
  return result;
}

// get the name of the gnuplot binary
const char *getGnuplotBinary() {
  return GNUPLOT_BINARY;
}

// returns either "perl" or the value of PERL_BIN as defined in Makefile.config.incl
const char *getPerlBinary(void) {
#ifdef PERL_BIN
  return QUOTED_PERL_BIN;
#else
  return "perl";
#endif
}

void remove_files_olderthan(const std::string &pathmask, time_t olderThanElapsed) {
  std::vector<std::string> matchingPaths;
  find_files_olderthan(pathmask,olderThanElapsed,matchingPaths);
  for (int i=(int)matchingPaths.size();i--;) {
    remove(matchingPaths[i].c_str());
  }
}
void find_files_olderthan(const std::string &pathmask, time_t olderThanElapsed,
			  std::vector<std::string>& matchingPaths)
{
  time_t now = time(NULL);
  std::string path;
  const char *slash = findRightmostPathSeperator_const(pathmask.c_str());
  if (slash) {
    path = pathmask.substr(0,1+slash-pathmask.c_str());
  }

#ifdef WIN32
  WIN32_FIND_DATA fdata;
  HANDLE srcFile = FindFirstFile(pathmask.c_str(), &fdata);
  if (srcFile == INVALID_HANDLE_VALUE) {
    return; // no matches
  }
  do
    {
      if (strcmp(fdata.cFileName, ".") != 0 &&
	  strcmp(fdata.cFileName, "..") != 0 )
	{
	  struct stat curEntryData;
	  stat(fdata.cFileName, &curEntryData);
	  if ((now-curEntryData.st_atime) > olderThanElapsed) {
	    std::string full(path);
	    full+=fdata.cFileName;
	    matchingPaths.push_back( full );
	  }
	}
    }
  while (FindNextFile(srcFile, &fdata));

  FindClose(srcFile);

#else

  glob_t globbuf;
  int rv = glob(pathmask.c_str(), 0, NULL, &globbuf);
  if(rv > 0 && rv != GLOB_NOMATCH)
    throw std::runtime_error("FindFilesByMask(): glob() error");

  DIR* curDir = opendir(".");
  struct stat curEntryData;

  for (size_t i=0; i < globbuf.gl_pathc; ++i)
    {
      stat(globbuf.gl_pathv[i], &curEntryData);
      if ((S_ISDIR(curEntryData.st_mode) ||
	   S_ISREG(curEntryData.st_mode) ||
	   S_ISLNK(curEntryData.st_mode)) && 
	  now-curEntryData.st_atime > olderThanElapsed ) {
	matchingPaths.push_back(globbuf.gl_pathv[i]);
      }
    }
  closedir(curDir);

  globfree(&globbuf);

#endif
}
