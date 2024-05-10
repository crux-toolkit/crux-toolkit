/*
 *
 * declarations for TPP's util.c
 *
 */

#ifndef UTIL_H_INCLUDED
#define UTIL_H_INCLUDED

#include <string>
#include <vector>
#include "sysdepend.h"

//
// TPP Configuration, Path, and URL management
//
// The following functions return absolute filesystem paths and relative URLs
// for use in TPP.  All return a value with linux style forward '/' slashes
// and with a final '/' appended to the value. Windows in general works with
// this style path format so long as its not used in a DOS command shell.
//
// Values returned are based on defaults set at compile time but can be
// overridden by environment variables (TPP_HOME) at runtime.  These functions
// are meant to replace the previous windows/cygwin/IIE/Apache infrastructure
// that was quite unwieldy.
//

// Filesystem paths (all have '/' appended to them)
const char *getHomePath();  // full path to TPP "home", the install directory
const char *getBinPath();   // full path to directory for program/scripts
const char *getCgiPath();   // full path to directory for cgi program/scripts
const char *getConfPath();  // full path to directory for conf files
const char *getHtmlPath();  // full path to directory for static images/js/css
const char *getDataPath();  // full path to directory containing data
const char *getLogPath();   // full path to directory containing log files 

// HTTP URLs (all have '/' appended to them)
const char *getBaseUrl();   // TPP base for all TPP URLs, e.g /tpp/
const char *getCgiUrl();    // TPP URL for cgi,  e.g. /tpp/cgi-bin/
const char *getHtmlUrl();   // TPP URL for html, e.g. /tpp/html/
const char *getDataUrl();   // TPP URL for data, e.g. /tpp/data/

const char *getTPPPort();   // TPP Port for web

const char *getSeq2MSUrl();   // Seq2MS Url
const char *getSeq2MSModel(); // Seq2MS Model


// Path/URL Conversion Functions
int tppDataPath2Url( const char* path, char* url );
int tppDataUrl2Path( const char* url, char* path );


// 
//
// Legacy filesystem path handling functions (DEPRECATED)
//

// Installdir issues
bool getIsInteractiveMode(); // TPP (web oriented) vs LabKey (headless)
const char *getPepXML_std_xsl();
const char *getPepXML_std_xsl_web_path();

//
/* compare the paths, ignoring filename and \ vs / */
//
int pathcmp(const char *path1, const char *path2);

//
// Get the WEBSERVER_ROOT env var, with path separators cleaned up. Webserver
// root is now actually the getDataPath() directory since the location of TPP's
// data directory isn't necessarily the web server's document root.
//
const char *getWebserverRoot();

//
// Set the WEBSERVER_ROOT variable, preserve previous
//
void pushWebserverRoot(const char *path);

//
// Restore previous WEBSERVER_ROOT
//
void popWebserverRoot();

//
// It's possible to inherit an environment with a webserver root that's actually
// a tmp dir (useful when actual webserver root is on a network share, or 
// AWS S3 etc) this lets you know the real location
//
const char *getTrueWebserverRoot(); 

//
// Get the WEBSERVER_ROOT env var, concatenated with webserver's tmp dir, with path separators cleaned up
//
const char *getWebserverTmpPath();

//
// Set the webserver tmp env var for this process and any children it spawns
//
void pushWebserverTmpPath(const char *path);

//
// Restore previous webserver tmp path
//
void popWebserverTmpPath();

//
// In-place replacement of path with webserver's tmp dir, if any
// so /blarg/foo becomes /inetpub/wwwroot/tmp/foo
//
void replace_path_with_webserver_tmp(char *path,int pathmaxlen);	// write this in tmpdir if we have one
void replace_path_with_webserver_tmp(std::string &path);			// write this in tmpdir if we have one

//
// Construct a tempfile name, possibly in the tmp dir if so configured
//
std::string make_tmpfile_name(const char *basis);

//
// Construct a tempfile name, in the current directory
//
std::string make_temp_name(const char *basis);


//
// Miscellaneous 
//
void getword(char *word, char *line, char stop);

char *makeword(char *line, char stop);

char *fmakeword(FILE *f, char stop, int *cl);

char x2c(char *what);

void unescape_url(char *url);

std::string urlencode(const std::string &s);

void plustospace(char *str);

int rind(char *s, char c);

// copy one file to another, return 0 on success
int copy_file(const char *fromName, const char *toName);

void send_fd(FILE *f, FILE *fd);


int get_line(char *s, int n, FILE *f); // was "getline", which conflicts with stdio.h

/** Basic routine to escape spaces in a string with
 *  the "\" character.  This can be used to create
 *  a POSIX path/filename from a string.  The escaped
 *  space path/filename can then be used in a "system()"
 *  call.
 */
char * escape_spaces (char * inStr);

void force_unCygwinify(char *filename); // strip out cygwin stuff from a path name, even in cygwin builds
void force_unCygwinify(std::string &filename); // strip out cygwin stuff from a path name, even in cygwin builds
#ifdef WINDOWS_NATIVE
#define unCygwinify(a) force_unCygwinify(a); // strip out cygwin stuff from a path name
#else
#define unCygwinify(a) /* nothing */
#endif
#define cygwinify(a,b) /* nothing */
inline bool isCygwinishPath(const char *path) {
  return !strncasecmp(path,"/cygdrive/",10);
}
inline bool isCygwinishPath(const std::string &path) {
  return isCygwinishPath(path.c_str());
}

// TODO: are these used anywhere?  If so are they needed?  In almost all cases
// Window's accepts the Linux form of paths 
#ifdef WINDOWS_NATIVE // change "/" to "\" 
void fixpath(char *in);
void fixpath(std::string &in);
#else
#define fixpath(in) (in)
#endif

//
// Return nonzero if fname appears to be absolute does not check for actual 
// existence of file
//
int isAbsolutePath(const char *fname);
int isAbsolutePath(const std::string &fname);

// Accept either .xml, .pepXML .pep.xml, or user config value 
// as valid pepXML filename extensions,
// or these extensions plus ".gz"
// return NULL if invalid, else return pointer to .ext
const char *hasValidPepXMLFilenameExt(const char *fname);

// accept either -prot.xml, .protXML, .prot.xml, or user config value
// as valid protXML filename extensions
// or these extensions plus ".gz"
// return NULL if invalid, else return pointer to .ext
const char *hasValidProtXMLFilenameExt(const char *fname);

// helpful functions for determining filetypes
const char *isDotGZ(const std::string &fname); // return ptr to ext if fname ends with .gz (case insensitive)
const char *isDotGZ(const char *fname); // return ptr to ext if fname ends with .gz (case insensitive)
// used by hasValidProtXMLFilenameExt etc, case insensitively finds ext or ext.gz at end
// of fname and returns pointer to its occurrance in fname, or NULL
const char *hasExtOrExtDotGZ(const char *fname, const char *ext); 

// determine if named file is gzipped, unzip to a tempfile if so
// returns: tmpfile name, or original file name if not gzipped
// optional maxchar arg limits effect of huge files when you just want to see header
std::string uncompress_to_tmpfile(const char *fname,int maxchar=-1);

// run gzip on file if it hasn't been done already
// may update filename to reflect addition of .gz to name
// returns system code, or 0 if no gzip needed
int do_gzip(std::string &filename);
int do_gzip(const char *filename);

//
// like mkstemp, but returns open FILE * instead of file handle
// replaces XXXXXX with some unique string in place
//
FILE *FILE_mkstemp(char *fname_XXXXXX);
FILE *FILE_mkstemp(std::string &fname_XXXXXX);

//
// in-place replacement of path with newpath
// so /blarg/foo becomes newpath/foo
// 
void replace_path(char *path,int pathmaxlen,const char *newpath);
void replace_path(std::string &path,const char *newpath);

//
// return a copy of the input filename with an absolute path in the webserver_root tree
// caller must free() the result
//
std::string translate_relative_webserver_root_path_to_absolute_filesystem_path(const char *path);

//
// remove the webserver root portion of the input path, if any
// so /inetpub/wwwroot/foo/blah becomes /foo/blah
//
void translate_absolute_filesystem_path_to_relative_webserver_root_path(std::string &path);
void translate_absolute_filesystem_path_to_relative_webserver_root_path(char *path);


//
// return a copy of the input filename with the filesystem webserver root prepended
//
std::string prepend_webserver_root(const char *path);

// fix up path's root dir if needed
void resolve_root_dir(char *path,int buflen);

// do we need to fix up the path at all?
// returns:
//   path of file that can be opened,
//   or just a copy of input path
std::string resolve_root(const char *path);

//
// get the WEBSERVER_URL env var, with path separators cleaned up
//
const char *getWebserverUrl();

// get the canonical ppepxml filename extension, including period
const char *get_pepxml_dot_ext();

// get the canonical protxml filename extension, including period
const char *get_protxml_dot_ext();

//
// case insensitive strstr
//
char *strstri( const char *str, const char *strCharSet );
char *strstrir( const char *str, const char *strCharSet ); // finds rightmost match

//
// snprintf may not nullterm on overflow, so use this
//
void safe_snprintf(char *buf, int buflen, const char *format, /* args */ ...);

//
// getcwd with path separators cleaned up
//
char *safepath_getcwd(char *buf, int buflen);
std::string safepath_getcwd(); // std::string version

// standard rename() doesn't work across filesystem boundaries, which we may need to do
int safe_rename(const char *from, const char *to);

//
// is c a path separator for linux or windows?
//
int isPathSeperator(char c); // return nonzero if c is a path separator / or \ .

//
// does path end with a path separator?
//
bool endsWithPathSeperator(const char *path);

int isAbsolutePath(const char *fname); // return nonzero if fname appears to a full path

char *makeFullPath(const char *fname); // return a full path version of fname - caller must free()
void makeFullPath(std::string &fname); // make a full path version of fname, in place

char *findRightmostPathSeperator(char *path); // return pointer to rightmost / or \ , or NULL.
const char *findRightmostPathSeperator_const(const char *path); // return pointer to rightmost / or \ , or NULL.
int findRightmostPathSeperator(const std::string &str); // return position to rightmost / or \, or npos.

// tidy up path separators, and if expected to exist do a bit
// of searching around in case of confused pathnames etc (may alter the
// input string - won't make it any longer, of course)
int fixPath(char *path, int bExpectExist); // return non0 if expect to exist and does exist

// fclose a potentially NULL fileptr
inline void safe_fclose(FILE *f) {
  if (f) {
    fclose(f);
  }
}

// unlink for std::string
inline void unlink(const std::string &fname) {
  unlink(fname.c_str());
}

// unlink with verbose error check
int verified_unlink(const char *fname);
inline int verified_unlink(const std::string &fname) {
  return verified_unlink(fname.c_str());
}

// rmdir with verbose error check
int verified_rmdir(const char *fname);
inline int verified_rmdir(const std::string &fname) {
  return verified_rmdir(fname.c_str());
}

// chdir with verbose error check
int verified_chdir(const char *fname);
inline int verified_chdir(const std::string &fname) {
  return verified_chdir(fname.c_str());
}

// system() with verbose error check
int verified_system(const char *cmd);

std::string XMLEscape(const std::string& s);

// get the name of the gnulot binary
const char *getGnuplotBinary();

// returns either "perl" or the value of PERL_BIN as defined in Makefile.config.incl
const char *getPerlBinary(void);

//
// format a png file reference to either use the
// show_tmp_pngfile.pl CGI or not (using the CGI causes
// the png file to be unlinked after use)
//
// example: 
// pngfile="foo.png";
//
// sprintf(fp, "<IMG SRC=\"%s\"/>", makeTmpPNGFileSrcRef(pngfile);
// output would be <IMG SRC="show_tmp_pngfile.pl?foo.png"/>
//
// sprintf(fp, "<IMG SRC=\"%s\"/>", makeTmpPNGFileSrcRef(pngfile,TRUE);
// output would be <IMG SRC="foo.png"/>
//
std::string makeTmpPNGFileSrcRef(const char *filename, bool persist=false);

// kill files matching the mask that have been around longer than the given elapsed time
void remove_files_olderthan(const std::string &pathmask, time_t olderThanElapsed);
// find files matching the mask that have been around longer than the given elapsed time
void find_files_olderthan(const std::string &pathmask, time_t olderThanElapsed,
			  std::vector<std::string>& matchingPaths);

//
// local helper func for cached env read, with path separators cleaned up
//
static char *check_env_var(const char *env,char *buf, int buflen, bool bIsPath, const char *defaultval=NULL);

#endif // UTIL_H_INCLUDED
