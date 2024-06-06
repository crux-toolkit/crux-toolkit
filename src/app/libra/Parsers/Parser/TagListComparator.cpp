/*

Program       : TagListComparator
Author        : Brian Pratt, Insilicos LLC
Date          : 11-10-05
SVN Info      : $Id: TagListComparator.cpp 8025 2020-02-14 01:05:59Z mhoopmann $

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

#include "TagListComparator.h"
#include "stdlib.h"
#include <iostream>
#include <fstream>
#include <errno.h>
#include "Common/util.h"
#include "Parsers/mzParser/mzParser.h" 
#include <math.h>
#include <vector>
#include <ctype.h>
#include "Util/RACI/RACI.h"

using namespace mzParser;

TagListComparator::TagListComparator(const char *progname,eTagListFilePurpose why,const Array<Tag*> &lhs,const Array<Tag*> &rhs) :
   m_lhs(lhs),m_rhs(rhs),m_bManageLHS(false),m_bManageRHS(false) {
   	if (LEARN_TEST==why) {
		cerr << progname << " internal error in TagListComparator" << endl;
	} else {
		compare(progname,why);
	}
}

static const Array<Tag*> dummy;

// depending on mode, either create the test file or test against it
TagListComparator::TagListComparator(const char *progname,eTagListFilePurpose why,const Array<Tag*> &lhs,const char *rhsFilename) :
   m_lhs(lhs),m_rhs((LEARN_TEST==why)?dummy:*readTagListFile(rhsFilename)),
   m_bManageLHS(false), m_bManageRHS(LEARN_TEST!=why) {
	if (LEARN_TEST==why) {  // write tags to file
		writeTagListFile(rhsFilename,m_lhs);
	} else {
		compare(progname, why);
	}
}
TagListComparator::TagListComparator(const char *progname,eTagListFilePurpose why,const char *lhsFilename,const char *rhsFilename) :
   m_lhs((LEARN_TEST==why)?dummy:*readTagListFile(lhsFilename)),
   m_rhs((LEARN_TEST==why)?dummy:*readTagListFile(rhsFilename)),
   m_bManageLHS(LEARN_TEST!=why),
   m_bManageRHS(LEARN_TEST!=why) {
   if (LEARN_TEST==why) { // copy file to file
      FILE *in=fopen(lhsFilename,"rb");
	  FILE *out=fopen(rhsFilename,"wb");
	  send_fd(in,out);
	  fclose(in);
	  fclose(out);
   } else {
	  compare(progname,why);
   }
}

TagListComparator::~TagListComparator() {
   if (m_bManageLHS) { // we allocated it
      for(int k = 0; k < m_lhs.size(); k++) {
         delete m_lhs[k];
      }
      delete &m_lhs;
   }
   if (m_bManageRHS) { // we allocated it
      for(int k = 0; k < m_rhs.size(); k++) {
         delete m_rhs[k];
      }
      delete &m_rhs;
   }
}

// return true if value appears to be a timestamp
static bool isdatetime(const char *value) {
    int yy,mo,dd,hh,mn,ss;
    return ((6==sscanf(value,"%d:%d:%d:%d:%d:%d", &yy,&mo,&dd,&hh,&mn,&ss)) ||
            (6==sscanf(value,"%d-%d-%dT%d:%d:%d", &yy,&mo,&dd,&hh,&mn,&ss)));
}

const char * TagListComparator::compare(const char *progname,eTagListFilePurpose why) { // return null if "same", nonnull otherwise
   const Array<Tag*> &lhs = m_lhs;
   const Array<Tag*> &rhs = m_rhs;
   const char *result = NULL;
   //
   // default stoopid implementation, you probably want to subclass this
   //
   if (lhs.size() != rhs.size()) {
      result = "different tag count"; // different size
      cerr << result << ": " << lhs.size() << " vs " << rhs.size() << endl ;
   }
   int count = lhs.size();
   if (count > rhs.size()) {
      count = rhs.size();
   }
   for (int i=0;i<count;i++) {
      if (lhs[i] && rhs[i]) {
         if (!lhs[i]->isIdentical(*rhs[i])) {
            result = unacceptable_difference(*lhs[i],*rhs[i]);
            if (result) {
               (*lhs[i]).write(cerr);
               cerr << "\tdoes not agree with expected value" << endl;
               (*rhs[i]).write(cerr);
               cerr << "\tat \"" << result << "\" in tag #" << i+1 << " at line " << (*rhs[i]).getLineNumber() << endl;
               break;
            }
         }
      } else if (lhs[i] != rhs[i]) {
         // one is null, the other isn't
         result = "???";
         cerr << result << endl << "\t\t";;
         break;
      }
   }
   if (result) {
	  const char *p=findRightmostPathSeperator_const(progname);
	  std::string errfilename(p?p:progname);
	  errfilename += ".bad.tagList";
      cerr << "\tregression test FAILURE in " << progname << " (see file " << errfilename << ")" ;
	  writeTagListFile(errfilename.c_str(),lhs);
	  if (FORCE_TEST != why) {
		  cerr << " - QUIT" << endl;
		  exit(1);
	  }
	  cerr << endl << "\t\t";
   } else {
      cout << "\t" << progname << " regression test success" << endl;
   }

   return result; // NULL==same, nonNULL==different
}

const char * TagListComparator::unacceptable_difference(Tag &lhs,Tag &rhs) { // return 0 if "same", nonzero otherwise
   if (!lhs.getName() && !rhs.getName()) {
      return NULL; // both empty, so same
   }
   if (!lhs.getName() || !rhs.getName()) {
      return lhs.getName()?lhs.getName():rhs.getName(); // one is empty
   }
   if (strcmp(lhs.getName(),rhs.getName())) {
      if ((!strncmp(lhs.getName(),"!--#",4)) &&
          (!strncmp(rhs.getName(),"!--#",4))) {
         return NULL; // just a directive, probably different webserver at creation time
      }
      return lhs.getName(); // different name
   }
   if ((lhs.isStart() != rhs.isStart()) || (lhs.isEnd() != rhs.isEnd())) {
      return lhs.getName();
   }
   if (lhs.getNameSpace()==NULL) {
      if (rhs.getNameSpace()) {
         return rhs.getNameSpace();
      }
   } else if (rhs.getNameSpace()==NULL) {
      return lhs.getNameSpace();
   } else if (strcmp(lhs.getNameSpace(),rhs.getNameSpace())) {
      return lhs.getNameSpace();
   }
   // compare attribute/value pairs
   int laIndex=0;
   int raIndex=0;
   while ((laIndex<lhs.getAttributeCount()) &&
          (raIndex<rhs.getAttributeCount())) {
      const char *attributeL = lhs.getAttribute(laIndex);
      const char *attributeR = rhs.getAttribute(raIndex);
      if (strcmp(attributeL,attributeR)) {
         const char *optional[] = {"win-cyg_","organism","peptide_prev_aa","peptide_next_aa",NULL}; // win-cyg_reference_database, win-cyg_source_files
         const char *optional_if_zero[] = {"num_input_4_spectra","num_input_5_spectra",NULL};
         int o,oz;
         for (o=0;optional[o];o++) {
            if (!strncasecmp(attributeL,optional[o],strlen(optional[o]))) {
            laIndex++;
               break;
            } else  if (!strncasecmp(attributeR,optional[o],strlen(optional[o]))) {
            raIndex++;
               break;
            }
         }
         if (!optional[o]) { // ran off end
            for (oz=0;optional_if_zero[oz];oz++) {
               if (!strncasecmp(attributeL,optional_if_zero[oz],strlen(optional_if_zero[oz]))) {
                  if (!atoi(lhs.getAttributeValue(laIndex))) {
                     laIndex++;
                     break;
                  }
               } else  if (!strncasecmp(attributeR,optional_if_zero[oz],strlen(optional_if_zero[oz]))) {
                  if (!atoi(rhs.getAttributeValue(raIndex))) {
                     raIndex++;
                     break;
                  }
               }
            }
         }
         if (!optional[o]&&!optional_if_zero[oz]) { // hit end of list
            return attributeL;
         }
      } else { // attribute names match, check values
         const char *valueLHS = lhs.getAttributeValue(laIndex);
         const char *valueRHS = rhs.getAttributeValue(raIndex);
         if (!valueRHS) {
            return attributeL; // no attribute by that name
         }
         if (strcmp(valueLHS,valueRHS)) {
             // values differ - is that OK?
             if (unacceptable_difference(attributeL,valueLHS,valueRHS)) {
                 // watch for bogus xml where everything is a "parameter", as in tandem output
                 if (!strcmp(attributeL,"value")) {
                     const char *name = lhs.getAttributeValue("name"); // get a more meaningful name
                     if (name) {
                         if (unacceptable_difference(name,valueLHS,valueRHS)) {
                             return name;
                         }
                         if (!strncmp(name,"timing, ",8)) {
                             cout << "note: " << name << " was " << valueLHS << " is " << valueRHS << endl;
                         }
                     }
                 } else {
                     return attributeL; // unacceptable difference
                 }
            }
         }
         laIndex++;
         raIndex++;
      }
   }
   return NULL; // same
}


static bool isNAN(const char *txt) {
   return (!strncasecmp(txt,"NaN",3)) ||  // MSVC NAN
          (!strncasecmp(txt,"-1.#J",5));  // GCC NAN
}

// return 0 if not a numeric string, else return length
// may modify input string to remove leading paren or comma
static int isnumeric(char *v) {
   if (v && *v) {
      Boolean bIsNum=true;
      Boolean done = false;
      int ndecimals=0;
      int nexps=0;
      char *vv=v;
	  while (isspace(*vv)) {
		  vv++; // ignore leading spaces
	  }
      for (;bIsNum && *vv && !done;vv++) {
         switch (*vv) {
         case '(':
            if (vv==v) {
               memmove(v,v+1,strlen(v+1));
               vv--; // so we don't advance
            } else {
               bIsNum = false;
               break;
            }
            break;
         case ')': // in a list?
            *vv=0;
            done = true;
            break;
         case ',': // in a list?
            if (vv==v) { // start of a list
               memmove(v,v+1,strlen(v+1));
               vv--; // so we don't advance
            } else {
               *vv=0;
               done = true;
            }
            break;
         case 'N':
         case 'n':
            if (isNAN(vv)) { // ,"NaN" MSVC NaN
               vv+=2;
            } else {
               bIsNum = false;
            }
            break;
         case '-':
            if (isNAN(vv)) { // ,"-1.#J",5)) gcc NaN
               vv+=4;
               break;
            }
         case '+':
            bIsNum&=((v==vv)||nexps);
            break;
         case 'e':
         case 'E':
            bIsNum &=!nexps++;
            break;
         case '0':
         case '1':
         case '2':
         case '3':
         case '4':
         case '5':
         case '6':
         case '7':
         case '8':
         case '9':
            break;
         case '.':
            bIsNum &= !ndecimals++;
            break;
         default:
            bIsNum = false;
            break;
         }
      }
      return bIsNum?(vv-v):0;
   } else {
      return 0;
   }
   
}

class cNumericList {
public:
   cNumericList(const char *text) {
      char *copy = strCopy(text);
      for (char *cp=copy;*cp;) {
         int n = isnumeric(cp);
         if (n) {
            m_numbers.push_back(std::string(cp));
            cp+=n;
         } else {
            m_numbers.clear();
            break;
         }
      }
      delete[] copy;
   }

   std::vector<std::string> m_numbers;

   int dissimilar(cNumericList &list2,double toler=0.0) {
      int result = 1; // until proven similar
      if (m_numbers.size()==list2.m_numbers.size()) {
         for (int i=(int)m_numbers.size();i--;) {
            double d1=atof(m_numbers[i].c_str());
            double d2=atof(list2.m_numbers[i].c_str());
            if ((d1==d2)||(isNAN(m_numbers[i].c_str())&&isNAN(list2.m_numbers[i].c_str()))) {
               result=0;
            } else { // not identical - at limit of accuracy?
               const char *dec = strchr(m_numbers[i].c_str(),'.');
               if (dec++) { // got decimal point, count digits to end or to exponent
                  int prec = 0;
                  while (isdigit(*dec)) {
                     prec++;
                     dec++;
                  }
                  if (prec>6) {
                     prec=6; // ridiculous precision
                  }
                  if (tolower(*dec)=='e') {
                     prec -= atoi(dec+1); // handle the exponent
                  }
                  double dprec = pow(10.0,-prec);
                  result=(fabs(d1-d2) > (1.5*dprec));
               } // end if we got a decimal point
                  if (result) {
                     if (d1 && d2) { // check simple % difference
                        double tol = toler?toler:((Max(fabs(d1),fabs(d2)) < .1)?.03:.01);
                        result = ((1.0-fabs(Min(d1,d2)/Max(d1,d2))) > tol);
                        if (result) {
                     break;
                  }
                     }
                  }
            }
         }
      }
      return result;
   }
};



int TagListComparator::unacceptable_difference(const char *attribute, const char *value1, const char *value2) {
   return check_difference(attribute, value1, value2);
}

int TagListComparator::check_difference(const char *attribute, const char *value1, const char *value2) {
   int result;
   if (!strcmp(attribute,"raw_data")) {
	   result = !(rampValidFileType(value1)&&rampValidFileType(value2)); // ignore mzML vs mzXML etc
   } else 
   if (strcmp(attribute,"time") &&
      strncmp(attribute,"timing, ",8) &&
      strcmp(attribute,"date") &&
      strcmp(attribute,"version") &&
      strcmp(attribute,"xsi:schemaLocation") &&
      strncmp(attribute,"win-cyg_",8)) // win-cyg_reference_database, win-cyg_source_files
   {  // not just a date/time/version/schemalocation thing, look closer
      result = 1; // assume failure
      int charge;
      // formatting issue? "1.023" vs "1.022"
      cNumericList v1(value1);
      cNumericList v2(value2);
      if (v1.m_numbers.size() && !v1.dissimilar(v2)) {
         result = 0;
      } else if (strstr(value1,"cygdrive")||strstr(value2,"cygdrive") ||
          (findRightmostPathSeperator_const(value1)&&findRightmostPathSeperator_const(value2))) {
         // could just be a linux/cygwin/win32 thing
         char *v1=strCopy(value1);
         char *v2=strCopy(value2);
         force_unCygwinify(v1); 
         strlwr(v1);
         force_unCygwinify(v2); 
         strlwr(v2);
         result = (0!=strcasecmp(v1,v2));
         if (result) {
            // could also be a http://<hostname> vs $WEBSERVER_ROOT thing
            if (getWebserverRoot()) {
               char *wsr = strCopy(getWebserverRoot());
               strlwr(wsr);
               char *c1=NULL,*c2;
               if (NULL!=(c1=strstr(v1,wsr))) {
                  c2 = v2;
               } else if (NULL!=(c1=strstr(v2,wsr))) {
                  c2 = v1;
               }
               if (c1 && (NULL!=(c2=strstr(c2,"http://")))) {
                  c1+=strlen(wsr);
                  char *slash = strchr(c2+strlen("http://")+1,'/');
                  if (slash) {
                     c2 = slash+1;
	       }
                  result = (0!=strcmp(c1,c2));
               }
               if (result && // still not OK 
		   ((strstr(v1,"http://") && strstr(v2,"http://")) || // possible hostname difference
		    (strstr(v1,wsr)||(strstr(v2,wsr))))) { // possible webserver root difference
                  // work backwards to find start of split, if any commonality call it good
		  bool saw_slash=false;
                  c1 = v1+strlen(v1);
                  c2 = v2+strlen(v2);
                  while ((c1>v1) && (c2>v2) && !strcasecmp(c1,c2)) {
		     if ('/'==*c1) {
                       saw_slash = true;
                     }
                     c1--;
                     c2--;
                  }
                  result = !(saw_slash); // result of 0 means acceptable match
               }
               delete[] wsr;
            }
         }
         if (result) { // could be a root path difference
            result = (0!=strcasecmp(resolve_root(v1).c_str(),resolve_root(v2).c_str()));
         }
         delete[] v1;
         delete[] v2;
      } else if ((!strcmp("num_incorr",attribute))||(!strcmp("num_corr",attribute))||(1==sscanf(attribute,"obs_%d_distr",&charge))) {
         // the calculation for this int value is based on doubles,
         // susceptible to rounding error - also random seed effects
         result = (abs(atoi(value1)-atoi(value2)) > 3);
         if (result) {
            result = v1.dissimilar(v2,.15);
         }
      } else if ((!strcmp("error",attribute))||(!strcmp("sensitivity",attribute))||(1==sscanf(attribute,"obs_%d_distr",&charge))) {
         // susceptible to rounding error and effects of random seeds
         result = v1.dissimilar(v2,.15);
      } else if ((!strcasecmp(value1,value2)) && (!strcasecmp(value1,"unknown"))) {
         result = 0; // acceptable
      } else if (isdatetime(value1) && isdatetime(value2)) {
         result = 0; // acceptable
      } else {
         // could be something like
         // name="MASCOT discrim score -22.35 error: -1.#J">
         char *v1=strCopy(value1);
         char *v2=strCopy(value2);
         char *vv1=v1;
         char *vv2=v2;
         char *vv1end=v1+strlen(v1);
         char *vv2end=v2+strlen(v2);
         result = 0;
         while (!result) {
            while (*vv1 && (*vv1==*vv2)) {
               vv1++;
               vv2++;
            }
            if (*vv1 && *vv2) {
               char *space1 = vv1+strcspn(vv1," \t");
               char *space2 = vv2+strcspn(vv2," \t");
               *space1=0;
               *space2=0;
               cNumericList l1(vv1);
               cNumericList l2(vv2);
               if (l1.m_numbers.size() && !l1.dissimilar(l2)) {
                   vv1=Min(vv1end,space1+1);
                   vv2=Min(vv2end,space2+1);
               } else {
                  result = 1; // not the same
               }
            } else {
               result = (*vv1 || *vv2);
               break; // reached end of one or both
            }
         }
         delete[] v1;
         delete[] v2;
      }
   } else {
      result = 0; // acceptable
   }
   return result;
}

//
// helper functions
//

// construct a filename for regression tests (caller must delete[])
// for example for cmdline "peptideProphetParser ICAT raft0020.xml"
//   constructTagListFilename("raft0020.xml","ICAT","peptideProphet")
// returns
//   "raft0020.xml.ICAT.peptideProphet.tagList"
char *constructTagListFilename(const char *jobFileName,const char *args, 
                               const char *programTagName, 
                               eTagListFilePurpose why) { 
   const char *testFileExt = ".tagList";
   char *name = new char[strlen(jobFileName)+strlen(programTagName)+strlen(args?args:"")+20];
   if (strstr(args,jobFileName) && strstr(args,programTagName) && strstr(args,testFileExt)) {
      strcpy(name,args); // args is actually the filename
   } else {
      strcpy(name,jobFileName);
      strcat(name,".");
      strcat(name,programTagName);
      if (args&&args[0]) {
         strcat(name,".");
         strcat(name,args);
      }
      strcat(name,testFileExt);
   }
   // remove any cygdrive stuff for best crossplatform checking
   char *cygdrive;
   while ( (cygdrive=strstr(name,"/cygdrive/")) ) {
      char *bump=cygdrive+strlen("/cygdrive/"); // copy from here
      cygdrive[0] = *bump; // drive letter
      cygdrive[1] = ':'; // make it look windowsy
      memmove(cygdrive+2,bump+1,strlen(bump));
   }
   // no weird characters in names, please
   for (char *cp=name;*cp;cp++) {
      if (strchr("*?/\\\"' \t:+;<>",*cp)) { 
         *cp = '_';
      }
   }

   // constructed names can be long...
   FILE *test;
   const char *mode=NULL;
   switch (why) {
   default:
      break;
   case LEARN_TEST:
      mode = "wb";
      break;
   case RUN_TEST:
   case FORCE_TEST:
      mode = "rb";
      break;
   }
   char *orig_name = strdup(name);
   while (NULL== (test = fopen(name,mode))) {
      int e = errno;
      const char *ellipses="_._._";
      switch(e) {
      case ENAMETOOLONG:
      case ENOENT:
         if (strlen(name)>3*strlen(ellipses)) { // shorten it up by eating out the middle
            size_t len = strlen(name);
            strcpy(name+(len/2)-2,ellipses);
            len = strlen(name);
            memmove(name+len,name+len+2,strlen(name+len+2)+1);
            break; // try again
         } 
      default:
         cerr << "problem with basis file \"" << orig_name << "\": "<<strerror(e) << endl;
         exit(e);
         break;
      }
      if (!e) {
         break;
      }
   }
   free(orig_name);
   fclose(test);

   switch (why) {
   default:
      break;
   case LEARN_TEST:
      cout << "writing regression test basis to file " << name << endl;
      break;
   case RUN_TEST:
   case FORCE_TEST:
      cout << "running regression test on basis file " << name << endl;
      break;
   }
   return name;
}

Boolean writeTagListFile(const char *TagListFilename,const Array<Tag*> &tags) {
   unlink(TagListFilename); // sometimes old cygwin files won't reopen for write
   ofstream fout(TagListFilename,ios::binary);
   if (!fout) {
      cerr << "cannot write regression test basis to file " << TagListFilename << endl;
      exit(1);
   }
   for(int k = 0; k < tags.size(); k++) {
      if(tags[k] != NULL) {
         tags[k]->write(fout);
      }
   }
   fout.close();
   return true;
}

Array<Tag*> *readTagListFile(const char *TagListFilename) {
   Array<Tag*> *tags;
   tags = new Array<Tag*>;
   RACI fin(TagListFilename); // can read gzipped xml
   
   if(! fin) {
      cerr << "error opening tagList file " << TagListFilename << endl;
      exit(1);
   }
   int linecount = 0;
   char *nextline = new char[LINE_WIDTH];
   while(fin.getline(nextline, LINE_WIDTH)) {
      linecount++;
      for (char *data = nextline; NULL!=(data = strchr(data, '<'));data++) {
         Tag *tag = new Tag(data);
         tag->setLineNumber(linecount);
         tag->trim(); // trim off excess vector capacity for smallest mem footprint
         tags->insertAtEnd(tag);
      }
   }
   
   delete [] nextline;
   return tags;
}

// parse commandline test arg for regression test learn or run, then remove -t portion
void checkRegressionTestArgs(char *testMode,eTagListFilePurpose &testType) {
   testType = NO_TEST;
   if (testMode) {
      // is there a test arg in there?
      char *cp=strstr(testMode,REGRESSION_TEST_CMDLINE_ARG);
	  // make sure this is really a switch "-t..." or "... -t..."
      if (cp && ((testMode==cp) || (' '==*(cp-1)))) { 
         size_t testarglen=strlen(REGRESSION_TEST_CMDLINE_ARG);
         if (REGRESSION_TEST_CMDLINE_ARG_MOD_LEARN==cp[testarglen]) {
            testType = LEARN_TEST;
            testarglen++;
	 } else if (REGRESSION_TEST_CMDLINE_ARG_MOD_FORCE==cp[testarglen]) {
            testType = FORCE_TEST;
            testarglen++;
         } else {
            testType = RUN_TEST;
         }
         // erase the testmode arg
         memmove(testMode,cp+testarglen,(strlen(cp)-testarglen)+1);
      }
   }
}
