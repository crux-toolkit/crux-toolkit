/*
Program       : tpp_tarball.cpp            
Purpose       : stuff for creating a tgz tarball, as in Mascot2XML                   

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

#include <sys/stat.h>
#include <iostream>
#include <vector>
#include <errno.h>
#include <fstream>
#include "sysdepend.h"
#include "util.h"
#include <map>

// libarchive stuff
#include <archive.h>
#include <archive_entry.h>

#include "tpp_tarball.h"


class tarball_info { //  private data
public:
	tarball_info() : tgz_(NULL) {
	};
	~tarball_info() {
		clear_commitlist();
	};
	struct archive	*tgz_;
	struct stat tgz_stat_;
	std::vector<char *> committed_filenames;
	std::vector<char *> tentative_filenames;
	std::vector<char *> tentative_file_contents;
	struct archive_entry *allocate_archive_entry_hdr(const char *filename,unsigned filesize);
#define DELETEARRAY(a) {for(int i=(int)(a).size();i--;){free((a)[i]);}(a).clear();}
	void clear_commitlist() {
		DELETEARRAY(committed_filenames); 
		DELETEARRAY(tentative_filenames); 
		DELETEARRAY(tentative_file_contents); 
	}
};

// helper for tgz file creation
struct archive_entry *tarball_info::allocate_archive_entry_hdr(const char *filename,unsigned filesize) {
	struct archive_entry	*hdr = archive_entry_new();
	archive_entry_set_pathname(hdr,filename);
	archive_entry_set_atime(hdr,tgz_stat_.st_atime, 0);
	archive_entry_set_gid(hdr,tgz_stat_.st_gid);
	archive_entry_set_size(hdr,filesize);
	archive_entry_set_mode(hdr,tgz_stat_.st_mode);
	return hdr;
}


tpp_tarball::tpp_tarball(const char *tarfilename, const char *directoryname) {
	FILE *tgztest = fopen(tarfilename,"wb");
	if (tgztest) { // can we overwrite?
		fputs("test",tgztest);
		fclose(tgztest);
	}
	tarball_ = new tarball_info;
	if ((!tgztest) || stat(tarfilename,&(tarball_->tgz_stat_))) {
		printf("cannot create tgz file \"%s\" (error \"%s\"), quitting",
			tarfilename,strerror(errno));
		exit(1);
	} else {
		tarball_->tgz_ = archive_write_new();
		archive_write_set_compression_gzip(tarball_->tgz_);
		archive_write_set_format_ustar(tarball_->tgz_);
		if (archive_write_open_file(tarball_->tgz_,tarfilename)) {
			printf("cannot create tgz file \"%s\" (error \"%s\"), quitting",
				tarfilename,archive_error_string(tarball_->tgz_));
			exit(1);
		}
	}
}

void tpp_tarball::create_memberfile(const char *tarfilename, const char *outfilename, const std::string &output) {
	const char *fname = findRightmostPathSeperator_const(outfilename);
	struct archive_entry	*hdr = tarball_->allocate_archive_entry_hdr(fname?(fname+1):outfilename,(unsigned int)output.length());
	int err;
	tarball_->committed_filenames.push_back(strdup(outfilename));
	if (!(err = archive_write_header(tarball_->tgz_, hdr))) {
		archive_write_data(tarball_->tgz_, output.c_str(), output.length());
	}
	if (err) {
		std::cerr << "error writing file " << outfilename << ": " << archive_error_string(tarball_->tgz_) << std::endl;
		exit(1);
	}
	archive_entry_free(hdr);
}

void tpp_tarball::create_tentative_memberfile(const char *tarfilename, const char *outfilename, const std::string &output) {
	tarball_->tentative_filenames.push_back(strdup(outfilename));
	tarball_->tentative_file_contents.push_back(strdup(output.c_str()));
}

void tpp_tarball::commit(const char *tarfilename, const char *directoryname) {
	for (int i=0;i<(int)tarball_->tentative_filenames.size();i++) {
		int j;
		for (j=(int)tarball_->committed_filenames.size();j--;) {
			if (!strcmp(tarball_->committed_filenames[j],tarball_->tentative_filenames[i])) {
				break;
			}
		}
		if (-1==j) {
			// no other file by that name, go ahead an commit
			std::string str(tarball_->tentative_file_contents[i]);
			create_memberfile(tarfilename, tarball_->tentative_filenames[i], str);
		}
	}
	if (archive_write_close(tarball_->tgz_)) {
		std::cerr << "error writing tgz file " << tarfilename << std::endl;
	}
	archive_write_finish(tarball_->tgz_);
	tarball_->clear_commitlist();
}


// roll our own tmpfile mechanism, to deal with IIS permissions issues

static std::map<FILE *,std::string> tmpfiles;

FILE *tpp_tmpfile() {
   char name[1024];
   safepath_getcwd(name,sizeof(name));
   strcat(name,"/tpp.tmp.XXXXXX");
   replace_path_with_webserver_tmp(name,sizeof(name)); // do this in designate tmp dir if any
#ifdef WINDOWS_NATIVE // MSVC or MinGW
   mktemp(name);
   FILE *result = fopen(name,"wb+");
#else
   int fd = mkstemp(name);
   FILE *result = fdopen(fd,"wb+");
#endif   
   tmpfiles[result] = std::string(name);
   return result;
}

void tpp_fclose(FILE *fp) {
   fclose(fp);
   if (tmpfiles[fp].length()) {
      unlink(tmpfiles[fp].c_str());
   }
   tmpfiles.erase(fp);
}

//
// handle zip files in addition to tgz
//
FILE * read_dta_or_out_from_tgz_file(char *szTarFile,  // not const, contents may be altered
                         char *szInputFile, // not const, contents may be altered
                         char *szFullPathInputFile, // not const, contents may be altered
                         int buflen) {
   struct stat statbuf;
   FILE *result=NULL;
   // tidy up filenames, construct a plausible tgz file name
   char *cp,*copy;
   int len;
   char *slash,*szTmp;

   // deal with WEBSERVER_ROOT relative paths
   resolve_root_dir(szTarFile,buflen);
   resolve_root_dir(szInputFile,buflen);
   resolve_root_dir(szFullPathInputFile,buflen);

   /*
   * check that szInputFile references a real directory
   */
   copy = strdup(szInputFile);
   if ((slash=findRightmostPathSeperator(copy))) {
      *slash = 0;
      if (stat(copy,&statbuf)) { // directory doesn't exist
         /* try to remove that middle foo in wwwroot/foo/foo/foo.0001.0001.2.dta */
         char *dir;
         free(copy);
         copy = strdup(szInputFile);
         slash = findRightmostPathSeperator(copy);
         if (slash) {
            char *slash2;
            *slash = 0;
            slash2 = findRightmostPathSeperator(copy);
            if (slash2) {
               strcpy(slash2+1,slash+1);
            }
         }
         dir = strdup(copy);
         slash = findRightmostPathSeperator(dir);
         if (slash) {
            *slash = 0;
         }
         if (!stat(dir,&statbuf)) {
            /* use this as the filename */
            strncpy(szInputFile,copy,buflen);
         }
         free(dir);
      }
   }
   free(copy);

   /*
   * modify szInputFile to remove full path
   */
   strncpy(szFullPathInputFile, szInputFile, buflen);
   cp = findRightmostPathSeperator(szInputFile);
   len = (int)strlen(cp=(cp?cp+1:szInputFile));
   szTmp = (char *)malloc(len+2);
   sprintf(szTmp, "*%s", cp);
   strncpy(szInputFile, szTmp, buflen);
   free(szTmp);

   //
   // check the tarfile name - on zoom click it may come in malformed
   // seeking "c:/Inetpub/wwwroot/foo/foo0020.0986.0986.2.dta" in
   // "c:/Inetpub/wwwroot/foo/foo.tgz" instead of in
   // "c:/Inetpub/wwwroot/foo/foo0020.tgz"
   if (stat(szTarFile,&statbuf)) {
      char *copy = strdup(szFullPathInputFile);
      // locate start of .scan.scan.charge.ext
      char *dot=NULL;
      int i;
      for (i=4; i-- && (dot = strrchr(copy,'.'));) {
         *dot = 0;
      }
      if (dot) {
         strcat(copy,".tgz");
      }
      if (!stat(copy,&statbuf)) {
         strcpy(szTarFile,copy);
      }
      free(copy);
   }
   // maybe a case of 
   //  "c:/Inetpub/wwwroot/ISB/Test/5396/Thu-Mar-29-07-24-03-2007.000.000.1.dta"
   // where c:/Inetpub/wwwroot/ISB/Test/5396.tgz exists and contains Thu-Mar-29-07-24-03-2007.000.000.1.dta
   if (stat(szTarFile,&statbuf)) {
      char *copy = strdup(szTarFile);
      char *slash = findRightmostPathSeperator(copy);
      if (slash && (strlen(slash)>4)) {
         strcpy(slash,".tgz");
         if (stat(copy,&statbuf)) {
            strcpy(slash,".zip");
         }
         if (!stat(copy,&statbuf)) {
            strcpy(szTarFile,copy);
         }
      }
      free(copy);
   }   

   /*
   * JENG: if szInput file starts with pps_ then assume this is a toftof run
   * from internal ISB runtoftof script using the mascot2dta program to create
   * dta files.  There will be an extra label after basename but before
   * scan #s in the .dta file names that needs to be removed from the
   * tar file name.  This next if statement should be irrelevent outside of ISB.
   */
   if (!strncmp(szInputFile, "*pps_", 5)) {
      char *pStr;
      /* remove the .tgz */
      szTarFile[strlen(szTarFile)-4]='\0';
      
      /* remove the label */
      if ( (pStr = strrchr(szTarFile, '.')) ) {
         *pStr='\0';
      }
      strcat(szTarFile, ".tgz");
   }
   if (stat(szFullPathInputFile,&statbuf)) { // unadorned filename doesn't exist, look for an archive
      if (stat(szTarFile,&statbuf)) { // maybe a zipfile?
         char *zip = strdup(szTarFile);
         char *dot = strstri(zip,".tgz");
         if (dot) {
            *++dot = 'z';
            *++dot = 'i';
            *++dot = 'p';
            if (!stat(zip,&statbuf)) { // maybe a zipfile?
               strcpy(szTarFile,zip);
            }
         }
         free(zip);
      }

      struct archive	*a = archive_read_new(); // allocate an archive
      archive_read_support_compression_all(a);
      archive_read_support_format_all(a);
      // that leading * doesn't really help us here
      const char *target = ('*'==*szInputFile)?szInputFile+1:szInputFile;
      if (!archive_read_open_file(a,szTarFile,0xfff)) {
         struct archive_entry	*entry;
         while (!archive_read_next_header(a, &entry)) {
            if (strstri(archive_entry_pathname(entry),target)) {
               const void *buff;
               size_t size;
               int64_t offset; // weird casting here to deal with mingw dll + VC6 exe off_t differences
               result = tpp_tmpfile(); // need to delete on fclose
               if (!result) {
                  return NULL;
               }
               while (!archive_read_data_block(a, &buff, &size, (off_t *)&offset)) {
				   if (size!=fwrite(buff,1,size,result)) {
		               fprintf(stderr,"problem writing temp file");
				       break;
				   }
               }
               rewind(result);
               break;
            }
            if (archive_read_data_skip(a)) { // move on to next header
               fprintf(stderr,"problem reading %s",szTarFile);
               break;
            }
         }
         archive_read_finish(a);
      }
   } else {  // read it directly
      result = fopen(szFullPathInputFile,"r");
   }
   return result;
}

void close_dta_or_out_from_tgz_file(FILE *fp) {
   tpp_fclose(fp); // if opened with tpp_tmpfile(), this will delete it
}


// stuff for reading from a tgz tarball, as in out2summary
class tarinfo {
public:
   tarinfo(struct archive	*a,const char *filename,const char *filematch):
      m_a(a),m_filename(filename),m_filematch(filematch) 
   {
   };
   ~tarinfo() {
      archive_read_finish(m_a);
   };
   struct archive	*m_a;
   std::string m_filename;
   std::string m_filematch;
};

tarinfo *tarball_read_open(const char *tarfilename, const char *filematch) {
   tarinfo *result = NULL;
   struct archive	*a = archive_read_new(); // allocate an archive
   archive_read_support_compression_all(a);
   archive_read_support_format_all(a);
   if (!archive_read_open_file(a,tarfilename,0xfff)) {
      result = new tarinfo(a,tarfilename,filematch);
   }
   return result;
}


FILE * tarball_read_next(tarinfo *tarball, char **filename) {
   FILE *result = NULL;
   tarinfo *info = (tarinfo *)tarball;
   struct archive_entry	*entry;
   while (!archive_read_next_header(info->m_a, &entry)) {
      if (strstri(archive_entry_pathname(entry),info->m_filematch.c_str())) {
         const void *buff;
         size_t size;
         int64_t offset; // weird casting here to deal with mingw dll + VC6 exe off_t differences
         *filename = strdup(archive_entry_pathname(entry));
         result = tpp_tmpfile(); // needs to be deleted on fclose
         while (!archive_read_data_block(info->m_a, &buff, &size, (off_t *)&offset)) {
			if (size!=fwrite(buff,1,size,result)) {
				fprintf(stderr,"problem writing temp file");
				break;
			}
         }
         rewind(result);
         return(result);
      }
      if (archive_read_data_skip(info->m_a)) { // move on to next header
         fprintf(stderr,"problem reading %s",info->m_filename.c_str());
         break;
      }
   }
   delete info;
   return NULL;
}

void tarball_close_stream(FILE *f) {
   tpp_fclose(f);
}
