/*

Program       : tpp_tarball.h            
Purpose       : stuff for creating a tgz tarball, as in Mascot2XML
                and for reading from tgz or zip files

Copyright (C) 2007 Insilicos LLC

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

#if !defined(TPP_TARBALL_H_INCLUDED_)
#define TPP_TARBALL_H_INCLUDED_

#include <string>
class tarball_info; // forward ref for private data

class tpp_tarball {
public:
	tpp_tarball(const char *tarfilename, const char *directoryname);
	~tpp_tarball();
	// add the contents of the string to the tarball as the named file
	void create_memberfile(const char *tarfilename, const char *outfilename, const std::string &output);
	// possibly add the contents of the string, unless create_memberfile gets called with same filename later
	void create_tentative_memberfile(const char *tarfilename, const char *outfilename, const std::string &output);
	// write any tentative members, and close the tarfile
	void commit(const char *tarfilename,const char *directoryname);
private:
	tarball_info *tarball_;
};

//
// return a pipe to read a file from a tar file - caller must close_dta_or_out_from_tgz_file(result)
//
FILE *read_dta_or_out_from_tgz_file(char *szTarFile,  // not const, contents may be altered
                         char *szInputFile, // not const, contents may be altered
                         char *szFullPathInputFile, // not const, contents may be altered
                         int buflen); // assumed buffer size for all three inputs
void close_dta_or_out_from_tgz_file(FILE *fp); // call this to close

//
// general tarball/zipfile reader, as in out2summary
//
class tarinfo;
tarinfo *tarball_read_open(const char *tarfilename, const char *filematch);
FILE * tarball_read_next(tarinfo *tarball, char **filename) ;
void tarball_close_stream(FILE *f);

#endif // !defined(TPP_TARBALL_H_INCLUDED_)
