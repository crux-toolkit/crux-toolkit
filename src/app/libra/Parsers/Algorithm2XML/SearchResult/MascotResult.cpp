#include "MascotResult.h"

/*

Program       : MascotResult for discr_calc of PeptideProphet 
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

MascotResult::MascotResult(char* szBuf, Boolean preexisting_probs){
  init();
  char validator[] = "http://www.ncbi.nlm.nih.gov"; //EXPECT="; //"DATABASE=";
  format_ = 0; // set
  ave_identity_ = 0.0; 
  process(szBuf, preexisting_probs, validator);

}

MascotResult::MascotResult(Array<Tag*>* tags) : SearchResult(tags) {
  char score_tag[] = "search_score";
  /*
  char result_tag[] = "search_result";
  char hit_tag[] = "search_hit";
  char pepproph_tag[] = "peptideprophet_result";
  double neutral_prec_mass;
  */
  //for(int k = 0; k < tags->length(); k++) 
  // (*tags)[k]->write(cout);


  if(processed_) {
    const int num_nec_fields = 5;
    Boolean found[num_nec_fields];
    int k;
    for(k = 0; k < num_nec_fields; k++)
      found[k] = False;

    Tag* tag;
    for(k = 0; k < tags->length(); k++) {
      tag = (*tags)[k];
      //tag->write(cout);
      if(! strcmp(tag->getName(), score_tag) && tag->isStart()) {
	if(found[0] == False && ! strcmp(tag->getAttributeValue("name"), "ionscore")) {
	  ionscore_ = atof(tag->getAttributeValue("value"));
	  found[0] = True;
	}
	else if(found[1] == False && ! strcmp(tag->getAttributeValue("name"), "identityscore")) {
	  identity_ = atof(tag->getAttributeValue("value"));
	  found[1] = True;
	}
	else if(found[2] == False && ! strcmp(tag->getAttributeValue("name"), "star")) {
	  star_ = ! strcmp(tag->getAttributeValue("value"), "1");
	  found[2] = True;
	}
	else if(found[3] == False && ! strcmp(tag->getAttributeValue("name"), "expect")) {
	  expect_ = atof(tag->getAttributeValue("value"));
	  found[3] = True;
	}
	else if(! strcmp(tag->getAttributeValue("name"), "homologyscore")) {
	  homology_ = atof(tag->getAttributeValue("value"));
	  found[4] = True;
	}

      } // if score
    } // next tag
    for(k = 0; k < num_nec_fields; k++)
      if(! found[k]) {
	cerr << "WARNING: for spectrum " << spectrum_ << ": couldn't find Mascot field " << (k+1) << endl;
	//processed_ = False;
      }

    //  if(processed_)
    //    cout << "processed" << endl;
    //  else
    //   cout << "NOT Processed" << endl;

  } // if score proc
/*
    else if(! strcmp(tag->getName(), pepproph_tag) && tag->isStart()) {
      probability_ = atof(tag->getAttributeValue("probability"));
    }
*/


}

char* MascotResult::extractDatabase(char* html) {
  char start[] = "&amp;Db=";
  char stop[] = "&amp";
  return extractDatabaseWithTags(html, start, stop);
}


void MascotResult::process(char* szBuf, Boolean preexisting_probs, char* valid) {
  if(strstr(szBuf, valid) == NULL) {
    return;
  }
   int iCharge,
       iMW,
       iNumDupl;

   char* szNewBuf = NULL;
   char szDup[SIZE_BUF],
        szTmp[SIZE_BUF],
        szFileName[2000],
        szProtein[2000],
        szPeptide[128];

   double dCn,
          dXC,
          dMH,
          dMHerror;


   double prob = -1.0;



   /*
    * Remove HTML tags in szBuf ... put new string in szNewBuf
    */
   szNewBuf = stripHTML(szBuf);
   if(szNewBuf == NULL)
     return;



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
   
   char ionscore[100];
   double ident;
   double homol;

   if(preexisting_probs) {

     min_num = 13;
       result = sscanf(szNewBuf, "%lf %s %s %lf (%lf) %s %lf %lf %s %s %s %s %s",
		&prob, szTmp, szFileName, &dMH, &dMHerror, ionscore, &ident, &homol, szTmp, szTmp, szProtein, szDup, szPeptide);


       if (szDup[0]!='+')
	 min_num--; // peptide is actually in szDup

       // check for space before protein name (if prot name starts with number, will fool sscanf...)
       substring = strstr(szNewBuf, szProtein);

       if(substring == NULL || strlen(szNewBuf) == strlen(substring) || ((substring-1)[0] != ' ' && (substring-1)[0] != '\t'))
	 result = 0; // cannot be true

       if(result < min_num - VAL_UNCERTAINTY) {
	 cerr << "returning with result " << result << " for " << szNewBuf << endl;
	 return; // done
       }

   }
   // no extra probability column
   else {

     min_num = 12;
     result = sscanf(szNewBuf, "%s %s %lf (%lf) %s %lf %lf %s %s %s %s %s",
			 szTmp, szFileName, &dMH, &dMHerror, ionscore, &ident, &homol, szTmp, szTmp, szProtein, szDup, szPeptide);

     substring = strstr(szNewBuf, szProtein);

     if(substring == NULL || strlen(szNewBuf) == strlen(substring) || ((substring-1)[0] != ' ' && (substring-1)[0] != '\t'))
       result = 0; // cannot be true
     if(result < min_num - VAL_UNCERTAINTY) {
       cerr << "returning with result " << result << " for " << szNewBuf << endl;
       return;
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

   //cerr << ionscore << endl;
   if(strlen(ionscore) > 0 && ionscore[strlen(ionscore)-1] == '*') {
     star_ = True;
     //cerr << " found a delta star for " << szFileName << endl;
     ionscore[strlen(ionscore)-1] = 0; // truncate
   }
   else 
     star_ = False;

   ionscore_ = atof(ionscore);

   sscanf(szFileName+strlen(szFileName)-1, "%d", &iCharge);


   // here set values
   
   charge_ = iCharge;
   spectrum_ = new char[strlen(szFileName)+1];
   strcpy(spectrum_, szFileName);


   massdiff_ = dMHerror;
   identity_ = ident;
   if(identity_ < 0.0)
     identity_ = 0.0;
   homology_ = homol;
   if(homology_ < 0.0)
     homology_ = 0.0;

   char * tprot = new char[strlen(szProtein)+1];
   strcpy(tprot, szProtein);
   if(tprot == NULL) {
     cerr << "could not copy protein " << szProtein << endl;
     exit(1);
   }
   proteins_->insertAtEnd(tprot);
   protein_ =  (*proteins_)[0];
   peptide_ = new char[strlen(szPeptide)+1];
   strcpy(peptide_, szPeptide);

   if(peptide_ == NULL) {
     cerr << "could not copy peptide " << szPeptide << endl;
     exit(1);
   }

   if(preexisting_probs)
     probability_ = prob;
   maldi_ = isMaldi(spectrum_);

   if(maldi_)
     charge_ = 1; // make the correction
   degen_ = iNumDupl; // > 0);

   processed_ = True;

   if(szNewBuf != NULL)
     delete szNewBuf;

   /*
   printf("prob=%0.2f, peptide=%s, protein=%s, file=%s, charge=%d, IonSc=%lf, Id=%lf, HomSc=%lf, MHerror=%lf, iNumDupl=%d\n",
	  probability_,
        peptide_,
         protein_,
         spectrum_,
         charge_,
         ionscore_,
         identity_,
         homology_,
         dMHerror,
         iNumDupl); 
   if(! maldi_) 
     printf("not ");
   printf("maldi\n");

   */





}

const char* MascotResult::getName() {
  return "Mascot";
}
