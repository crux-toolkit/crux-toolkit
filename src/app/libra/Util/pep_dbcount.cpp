#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>

#define SIZE_BUF 4096

using namespace std;



struct PeptideContainer { // public class
  string peptideName_;
  double mass_;
  vector<string> proteinList_;

  PeptideContainer() : peptideName_(""), mass_(-1) {
  }


  void clear(void) {
    peptideName_ = "";
    mass_ = -1;
    proteinList_.clear();
  }


  void print(void) {    
    // build a space-separated list of protein names for this peptide's proteins
    string proteins = "";
    for (vector<string>::size_type curProtIndx = 0; curProtIndx < proteinList_.size(); ++curProtIndx) {
      if (proteins != "") {
         proteins += " ";
      }
      proteins += proteinList_[curProtIndx];
    }
    
    cout << peptideName_
      << "\t"
      << proteinList_.size()
      << "\t"
      << mass_
      << "\t"
      << proteins
      << endl;
  }

  // only store the protein if it hasn't been seen before
  void addUniqueProtein(const string& proteinString) {
    // horribly lazy lazy O(N) search here.  actually, same as orig. c
    // version; a motivated person could replace with map
     bool found = false;
    for (vector<string>::size_type curProtIndx = 0; curProtIndx < proteinList_.size(); ++curProtIndx) {
      if (proteinString == proteinList_[curProtIndx]) {
        found = true;
        break;
      }
    }
    if (!found) {
      // add the protein
      proteinList_.push_back(proteinString);
    }
  }


};


int main(int argc, char **argv)
{
   FILE *fp;
   char szBuf[SIZE_BUF];
   string digestDBFilename = "";

   if (argc != 2)
   {
     cout << endl
       << " USAGE:  " << argv[0] << " <digest_db output>" << endl
       << endl
       << " This program takes output from digestdb and prints out tab delimited file" << endl
       << " including" << endl
       << " - peptide sequence" << endl
       << " - # proteins the peptide appears in" << endl
       << " - peptide mass" << endl
       << " - comma separated protein list" << endl
       << endl
       << " NOTE: this program assumes that the digest_db output is already sorted by peptide name (field 5)" << endl
       << endl;
     exit(1);
   }


   digestDBFilename = argv[1];

   // try to open digestDB input file
   if ( (fp=fopen(digestDBFilename.c_str(), "r"))==NULL)
   {
      printf("\n");
      printf(" Error - cannot open input file %s\n\n", digestDBFilename.c_str());
      exit(1);
   }


   PeptideContainer curPeptide;

   // NOTE:
   // assuming the peptides are already sorted in the digestDB output (the input file we're looping through)

   char* dummy = fgets(szBuf, SIZE_BUF, fp); // prime the pump
   bool firstTime = true;
   while(!feof(fp))
   {
      char szProtein[500];
      char szPeptide[500];
      char cPrevAA;
      double dMass;
      int iLen;
      sscanf(szBuf, "%d %s %lf %c %s", &iLen, szProtein, &dMass, &cPrevAA, szPeptide); 
      
      //cout << "read line " << szBuf << endl;

      string curPeptideName = szPeptide;
      string curProteinName = szProtein;

      // are we on to the next peptide?
      // NOTE! assuming sorted digest_DB input!
      if (curPeptideName != curPeptide.peptideName_) {
        if (firstTime) {
          firstTime = false;
        }
        else {
          // we're going on to the next peptide
          // print the current peptide
          curPeptide.print();
          // and clear
          curPeptide.clear();
        }

        // we're on to a new peptide. so, save the pep data and the
        // first protein for this peptide
        curPeptide.peptideName_ = curPeptideName;
        curPeptide.mass_ = dMass;
        curPeptide.addUniqueProtein(curProteinName);
      }
      else {
        // we've encountered another protein for the current peptide.
        // we'll add the protein to the current the protein list

        // but first, verify that things haven't changed
         if (dMass != curPeptide.mass_) {
           cerr << "error in digest_db file: multiple masses reported for "
              << curPeptideName 
              << ": " << endl
              << curPeptide.mass_ << " and " << dMass << endl;
           exit(1);
         }
        
         curPeptide.addUniqueProtein(curProteinName);
      }

      dummy = fgets(szBuf, SIZE_BUF, fp);
   }

   fclose(fp);

   // print the last peptide
   curPeptide.print();
}

