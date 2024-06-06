/*

Program       : TagListComparator                                                    
Author        : Brian Pratt, Insilicos LLC                                                       
Date          : 11-10-05 

Class for comparing lists of tags, used in regression tests.  Allows
programs to write uniquely named files with the XML tags they're responsible
for creating, and read them back in subsequent runs for regression tests.
Support for naive comparison is built in, and virtual functions for smarter
handling of differences are provided.



Copyright (C) 2005 Insilicos LLC, All Rights Reserved

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

#ifndef TAGLISTCOMPARATOR_INCL
#define TAGLISTCOMPARATOR_INCL

#include "Tag.h"

// regression test stuff - bpratt Insilicos LLC, Nov 2005
#define REGRESSION_TEST_CMDLINE_ARG "-t"  // as in -tmytest (run) or -t!mytest (learn) or -t#mytest (run, continue on error)
#define REGRESSION_TEST_CMDLINE_ARG_MOD_LEARN '!'
#define REGRESSION_TEST_CMDLINE_ARG_MOD_FORCE '#'

// no test, learn test, run test, run test without exit on failure
enum eTagListFilePurpose {NO_TEST,LEARN_TEST,RUN_TEST,FORCE_TEST}; 

class TagListComparator {
public:
   TagListComparator(const char *progname,eTagListFilePurpose why,const Array<Tag*> &lhs,const Array<Tag*> &rhs); // compare taglists
   TagListComparator(const char *progname,eTagListFilePurpose why,const Array<Tag*> &lhs,const char *rhsFilename); // create file (learn), or read from file then compare
   TagListComparator(const char *progname,eTagListFilePurpose why,const char *lhsFilename,const char *rhsFilename); // copy file (learn), or read from files then compare
   virtual ~TagListComparator();
private:
   const Array<Tag*> &m_lhs;
   const Array<Tag*> &m_rhs;
   Boolean m_bManageLHS;
   Boolean m_bManageRHS;
   const char *compare(const char *progname,eTagListFilePurpose why); // return NULL if "same", otherwise name of attribute that differed
public:
   static int check_difference(const char *attribute, // return nonzero if these values are too different for this attribute
       const char *value1, const char *value2);
protected:
   // allow things like date to differ -  child classes can extend this
   virtual const char *unacceptable_difference(Tag &lhs,Tag &rhs); // return NULL if "same", useful text otherwise
   // this is the one you're most likely to extend, though
   virtual int unacceptable_difference(const char *attribute, // return nonzero if these values are too different for this attribute
       const char *value1, const char *value2);
};

//
// helper functions
//
// construct a filename for regression tests (caller must delete[])
// for example for cmdline "peptideProphetParser ICAT raft0020.xml"
//   constructTagListFilename("raft0020.xml","ICAT","peptideProphet")
// returns
//   "raft0020.xml.ICAT.peptideProphet.tagList"
char *constructTagListFilename(const char *jobFileName, // input data file
                               const char *args, // args to the program
                               const char *programTagName,  // program name
                               eTagListFilePurpose why); // use NONE for silence

// parse commandline test arg for regression test learn or run, then remove -t portion
void checkRegressionTestArgs(char *testModeArg,eTagListFilePurpose &TestType);

// write the indicated tags to the indicated filename
Boolean writeTagListFile(const char *TagListFilename, const Array<Tag*> &tags);

// read the contents of the named file into an array (caller must delete the array)
Array<Tag*> *readTagListFile(const char *TagListFilename);

#endif // TAGLISTCOMPARATOR_INCL
