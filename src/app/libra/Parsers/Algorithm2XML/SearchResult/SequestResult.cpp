#include "SequestResult.h"

/*

Program       : SequestResult for discr_calc of PeptideProphet 
Author        : Andrew Keller <akeller@systemsbiology.org>                                                       
Date          : 11.27.02 


Copyright (C) 2003 Andrew Keller

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

Andrew Keller
Institute for Systems Biology
401 Terry Avenue North 
Seattle, WA  98109  USA
akeller@systemsbiology.org

*/


/*
format -1 uncertain
format 0 no mw column
format 1 has mw column
*/

SequestResult::SequestResult() { 

  xcorr_ = 0;
  delta_ = 0;
  deltastar_ = 0;
  rank_ = 0;
  sp_score_ = 0;

}

SequestResult::SequestResult(Array<Tag*>* tags) : SearchResult(tags) {
  char score_tag[] = "search_score";
  /*
  char result_tag[] = "search_result";
  char hit_tag[] = "search_hit";
  char pepproph_tag[] = "peptideprophet_result";
  double neutral_prec_mass;
  */
  if(processed_) {
    const int num_nec_fields = 5;
    Boolean found[num_nec_fields];
    int k;
    for(k = 0; k < num_nec_fields; k++)
      found[k] = False;

    Tag* tag;
    expect_ = -1;
    for(k = 0; k < tags->length(); k++) {
      tag = (*tags)[k];

      if(! strcmp(tag->getName(), score_tag) && tag->isStart()) {
	if(! found[0] && (! strcmp(tag->getAttributeValue("name"), "xcorr") || ! strcmp(tag->getAttributeValue("name"), "XCorr"))) {
	  xcorr_ = atof(tag->getAttributeValue("value"));
	  found[0] = True;
	}
	else if(! found[1] && ! strcmp(tag->getAttributeValue("name"), "deltacn")) {
	  delta_ = atof(tag->getAttributeValue("value"));
	  found[1] = True;
	}
	else if(! found[2] && ! strcmp(tag->getAttributeValue("name"), "deltacnstar")) {
	  deltastar_ = atof(tag->getAttributeValue("value"));
	  found[2] = True;
	}
	else if(! found[3] && ! strcmp(tag->getAttributeValue("name"), "sprank")) {
	  rank_ = atoi(tag->getAttributeValue("value"));
	  found[3] = True;
	}
	else if(! found[4] && (! strcmp(tag->getAttributeValue("name"), "spscore") || ! strcmp(tag->getAttributeValue("name"), "SpScore")) ) {
	  sp_score_ = atof(tag->getAttributeValue("value"));
	  found[4] = True;
	}

	if(! strcmp(tag->getAttributeValue("name"), "expect")) {
	  expect_ = atof(tag->getAttributeValue("value"));
	}

      }

/*
    else if(! strcmp(tag->getName(), pepproph_tag) && tag->isStart()) {
      probability_ = atof(tag->getAttributeValue("probability"));
    }
*/

    } // next tag
    for(k = 0; k < num_nec_fields; k++)
      if(! found[k])
	processed_ = False;

    if(xcorr_ == 0.0)
      processed_ = False;
    else 
      processed_ = True;



  } // if main score processed

}

SequestResult::SequestResult(char* szBuf, Boolean preexisting_probs, int format) {
  init(); // initial settings

  int formats[] = {1,0};  // 1 for bHasMWColumn True, 0 for False
  char validator[] = "http://www.ncbi.nlm.nih.gov"; //EXPECT="; //"DATABASE=";
  //char db_tag[] = "&amp;Db=";

  process(szBuf, preexisting_probs, format, formats, sizeof(formats)/sizeof(int), validator);
}

char* SequestResult::extractDatabase(char* html) {
  char start[] = "&amp;Db=";
  char stop[] = "&amp";
  return extractDatabaseWithTags(html, start, stop);
}


void SequestResult::process(char* szBuf, Boolean preexisting_probs, int format, int* formats, int numformats, char* valid) {
  if(strstr(szBuf, valid) == NULL) {
    return;
  }
  if(
     (format == -1 && ! parseUnknown(szBuf, preexisting_probs, formats, numformats)) ||
     (format >= 0 && ! parse(szBuf, preexisting_probs, format))
     ){
    cerr << "error parsing " << szBuf << " with format " << format << endl;
    exit(1);
  }
}

Boolean SequestResult::parseUnknown(const char* szBuf, Boolean preexisting_probs, int* formats, int numformats) {
  // formats dictates order

  Boolean output = False;
  if(formats == NULL)
    return output;
  for(int k = 0; k < numformats; k++) 
    if(! output)
      output = parse(szBuf, preexisting_probs, formats[k]);
  return output;
}

// adapted from Interact by J.Eng
Boolean SequestResult::parse(const char* szBuf, Boolean preexisting_probs, int format) {

   int iCharge,
       iMW,
       iNumDupl;

   int rank;
   // char szNewBuf[SIZE_BUF],    /* string with html tags removed */
   
   char* szNewBuf = NULL;
   char szDup[SIZE_BUF],
        szCn[SIZE_BUF],
        szTmp[SIZE_BUF],
        szFileName[128],
        szProtein[128],
        szPeptide[128];

   double dCn,
          dXC,
          dMH,
          dMHerror;


   double prob = 2.0; // illegal value



   /*
    * Remove HTML tags in szBuf ... put new string in szNewBuf
    */
   szNewBuf = stripHTML(szBuf);
   if(szNewBuf == NULL)
     return False;
   /*
   iLen=strlen(szBuf);
   iBufCt=0;
   iNewBufCt=0;
   while (iBufCt < iLen)
   {
      if (szBuf[iBufCt]=='<')
      {
         while (szBuf[iBufCt++] != '>');
      }
      else
      {
         szNewBuf[iNewBufCt]=szBuf[iBufCt];
         if (iNewBufCt-1 >= 0 && szNewBuf[iNewBufCt]=='/' && szNewBuf[iNewBufCt-1]!='.')
            szNewBuf[iNewBufCt]=' ';

         iBufCt++;
         iNewBufCt++;
      }
   }
   szNewBuf[iNewBufCt]='\0';
   */




   iCharge=0;
   iMW=0;
   iNumDupl=0;

   dCn=0.0;
   dXC=0.0;
   dMH=0.0;
   dMHerror=0.0;

   szProtein[0]='\0';
   szPeptide[0]='\0';
   szDup[0]='\0';
   
   int result = -1;
   char* substring = NULL;
   int min_num;
   int num_matched_ions = -1;
   int tot_num_ions = -1;
   char spscore[100];
   
   format = 0;
   if(preexisting_probs) {

     if(format == 1) {
       min_num = 15;
       result = sscanf(szNewBuf, "%lf %s %s %lf (%lf) %lf %s %s %d %d %d %d %s %s %s",
		&prob, szTmp, szFileName, &dMH, &dMHerror, &dXC, szCn, spscore, &rank, &num_matched_ions, &tot_num_ions, &iMW, szProtein, szDup, szPeptide);
       format_ = format;

       if (szDup[0]!='+')
	 min_num--; // peptide is actually in szDup

       // check for space before protein name (if prot name starts with number, will fool sscanf...)
       substring = strstr(szNewBuf, szProtein);

       if(substring == NULL || strlen(szNewBuf) == strlen(substring) || ((substring-1)[0] != ' ' && (substring-1)[0] != '\t'))
	 result = False; // cannot be true

       if(result < min_num - VAL_UNCERTAINTY)
	 return False;

     }
     else if(format == 0) { // format 0
       min_num = 14;
       result = sscanf(szNewBuf, "%lf %s %s %lf (%lf) %lf %s %s %d %d %d %s %s %s",
		&prob, szTmp, szFileName, &dMH, &dMHerror, &dXC, szCn, spscore, &rank, &num_matched_ions, &tot_num_ions, szProtein, szDup, szPeptide);
       format_ = format;

       if(result < min_num - VAL_UNCERTAINTY)
	 return False;
     }
     else {
       cerr << "cannot interpret format " << format << endl;
       exit(1);
     }
   }
   // no extra probability column
   else {


     if(format == 1) {
       min_num = 14;
       result = sscanf(szNewBuf, "%s %s %lf (%lf) %lf %s %s %d %d %d %d %s %s %s",
		szTmp, szFileName, &dMH, &dMHerror, &dXC, szCn, spscore, &rank, &num_matched_ions, &tot_num_ions, &iMW, szProtein, szDup, szPeptide);
       format_ = format;
       if (szDup[0]!='+')
	 min_num--; // peptide is actually in szDup
       // check for space before protein name (if prot name starts with number, will fool sscanf...)
       substring = strstr(szNewBuf, szProtein);

       if(substring == NULL || strlen(szNewBuf) == strlen(substring) || ((substring-1)[0] != ' ' && (substring-1)[0] != '\t'))
	 return False; // cannot be true

       if(result < min_num - VAL_UNCERTAINTY)
	 return False;
      
     }
     else if(format == 0) { // format 0

       min_num = 13;
         result = sscanf(szNewBuf, "%s %s %lf (%lf) %lf %s %s %d %d %d %s %s %s",
			 szTmp, szFileName, &dMH, &dMHerror, &dXC, szCn, spscore, &rank, &num_matched_ions, &tot_num_ions, szProtein, szDup, szPeptide);
	 format_ = 0;

	 if(result < min_num - VAL_UNCERTAINTY)
	   return False;

     }
     else {
       cerr << "cannot interpret format " << format << endl;
       exit(1);
     }

   } // if not preexisting probs

   /*
    * Some columns have the # of duplicate protein entries while others don't ... need to adjust for this
    */


   if (szDup[0]!='+')
   {
      strcpy(szPeptide, szDup);
      iNumDupl=0;
   }
   else
   {
      sscanf(szDup, "+%d", &iNumDupl);
   }
   
   /*
    * Just take raw peptide sequence, not preceeding and trailing amino acids
    * This sequence still would have the modifications in it
    */

   if(strlen(szCn) > 0 && szCn[strlen(szCn)-1] == '*') {
     deltastar_ = 1.0;
   }
   else {
     deltastar_ = 0.0;
     //cout << "no delta star for " << szCn << endl;
   }
   sscanf(szCn, "%lf", &dCn);   /* Cn sometime has '*' character; this gets rid of it in dCn */
   sscanf(szFileName+strlen(szFileName)-1, "%d", &iCharge);


   // here set values
   
   charge_ = iCharge;

   // snip off preceding useless './'
   if(strlen(szFileName) > 2 && szFileName[0] == '.' && szFileName[1] == '/')
     spectrum_ = strCopy(szFileName+2);
   else
     spectrum_ = strCopy(szFileName);
   massdiff_ = dMHerror;
   xcorr_ = dXC;
   delta_ = dCn;
   rank_ = rank;
   if(preexisting_probs)
     probability_ = prob;

   protein_ = strCopy(szProtein);
   if(protein_ == NULL) {
     cerr << "could not copy protein " << szProtein << endl;
     exit(1);
   }

   peptide_ = strCopy(szPeptide);
   if(peptide_ == NULL) {
     cerr << "could not copy peptide " << szPeptide << endl;
     exit(1);
   }

   maldi_ = isMaldi(spectrum_);

   if(maldi_)
     charge_ = 1; // make the correction
   degen_ = iNumDupl; // > 0);

   neutral_mass_ = dMH - 1.008;
   num_matched_ions_ = num_matched_ions;
   tot_num_ions_ = tot_num_ions;
   sp_score_ = atof(spscore);

   processed_ = True;


   if(charge_ < 1) {

   
     printf("Error: illegal charge for result %s\npeptide=%s, protein=%s, file=%s, charge=%d, XCorr=%f, dCn=%f, MH+=%f, MHerror=%f, iNumDupl=%d\n",
	    szNewBuf,
	    szPeptide,
	    szProtein,
	    szFileName,
	    iCharge,
	    dXC,
	    dCn,
	    dMH,
	    dMHerror,
	    iNumDupl); 

   }
   /*
   if(! maldi_) 
     printf("not ");
   printf("maldi\n");
   */
   if(szNewBuf != NULL)
     delete szNewBuf;
 

   return True;

}

const char* SequestResult::getName() { return "Sequest"; }
