#include "CStPeterMatrix.h"
#include "Quantitation/StPeter/ProtXMLParser.h"

#include <string>
#include <iostream>

#define VERSION "2.0.0a2"
#define BDATE "May 1 2020"

using namespace std;

void processCmd(int argc, char* argv[], CStPeterMatrix& matrix, vector<string>& files, string& global);
void usage();

int main(int argc, char* argv[]){
  cout << "StPeter2Matrix " << VERSION << " " << BDATE << endl;
  cout << "Copyright 2020, Michael Hoopmann, Institute for Systems Biology" << endl;

  if(argc<3) {
    usage();
    return 1;
  }

  string global;
  vector<string> files;
  CStPeterMatrix stpMatrix;
  CProteinMasterList masterList;
  int iPercent;
  int iTmp;

  //process command line
  processCmd(argc, argv,stpMatrix,files,global);

  //Set progress meter
  iPercent = 0;
  printf("Loading ProtXML files: %2d%%", iPercent);
  fflush(stdout);

  //read all protXML files and create a masterlist of all observed proteins  
  vector<ProtXMLParser> protXML(files.size());
  for(size_t i=0;i<files.size();i++){
    protXML[i].parse(files[i].c_str());
    masterList.addProtXML(protXML[i],files[i]);
    stpMatrix.addName(files[i]);

    //Update progress meter
    iTmp = (int)((double)i / files.size() * 100);
    if (iTmp>iPercent){
      iPercent = iTmp;
      printf("\b\b\b%2d%%", iPercent);
      fflush(stdout);
    }

  }

  //Finalize progress meter
  printf("\b\b\b100%%");
  cout << endl;
  
  /* For testing purposes
  for(size_t j=0;j<files.size();j++){
    size_t theIndex=-1;
    
    for(size_t i=0;i<masterList.size();i++){
      if(masterList[i].protein.compare("sp|P02748|CO9_HUMAN")==0) {
        //cout << "CO9_HUMAN: " << masterList[i].index[j] << endl;
        //theIndex = masterList[i].index[j];
      }
      if (masterList[i].protein.compare("sp|Q96BY6|DOC10_HUMAN") == 0) {
        if (masterList[i].index[j]>-1 && protXML[j][masterList[i].index[j]].stPeter.dSIn>0){
          cout << files[j] << endl;
          cout << "DOC10_HUMAN: " << protXML[j][masterList[i].index[j]].stPeter.dSIn << endl;
        }
      }
      if (masterList[i].protein.compare("sp|Q5JSL3|DOC11_HUMAN") == 0) {
        //cout << "DOC11_HUMAN: " << masterList[i].index[j] << endl;
        
      }
      if (masterList[i].protein.compare("sp|Q9BZ29|DOCK9_HUMAN") == 0) {
        //if (theIndex != -1){
        //  if (masterList[i].index[j] < theIndex + 2 && masterList[i].index[j] > theIndex - 2) {
        //    cout << files[j] << endl;
        //    cout << "DOC9_HUMAN: " << masterList[i].index[j] << endl;
        //  }
        //}
      }
    }
  }
  exit(2);
  */
  
  //extract quantified proteins into matrix.
  stpMatrix.extractStPeter(protXML,masterList);
  cout << stpMatrix.size() << " total proteins in matrix." << endl;

  //cross-reference with global list of proteins
  if(global.size()>0) {
    cout << "Cross-referencing with " << global << " to control global protein error rate." << endl;
    ProtXMLParser p;
    p.parse(global.c_str());
    stpMatrix.crossFilter(p);
    cout << stpMatrix.size() << " total proteins in matrix." << endl;
  }

  //apply filters (if any)
  stpMatrix.decoyFilter();
  cout << stpMatrix.size() << " total proteins in matrix." << endl;

  //export results
  if(!stpMatrix.exportMatrix()){
    cout << "ERROR: Failed to export matrix." << endl;
    return -2;
  }

  return 0;
}

//Process the command line, divving up the words into options, input, and output.
//Nothing smart here, so don't get cute.
void processCmd(int argc, char* argv[], CStPeterMatrix& matrix, vector<string>& files, string& global){
  int i;
  string str;
  
  global.clear();
  files.clear();
  for(i=1;i<argc-1;i++){

    //process options
    if(strcmp(argv[i],"-f")==0){

      //set the decoy filter
      matrix.setFilter(argv[i+1]);
      i++;

    } else if (strcmp(argv[i], "-g") == 0){

      //set global protXML file name.
      global=argv[i+1];
      i++;

    } else if (strcmp(argv[i], "-l") == 0){

      //read in a list of input files. Note that these are IN ADDITION TO any
      //files listed explicitly in the commandline. Also note the possibility
      //to use multiple list files.
      FILE* f = fopen(argv[i+1], "rt");
      char st[1024];
      char* tok;
      if(f==NULL) {
        cout << "ERROR: Cannot read " << argv[i+1] << " for protXML input." << endl;
        exit(-1);
      }
      while (!feof(f)){
        if (fgets(st, 1024, f) == NULL) continue;
        if (strlen(st)<2) continue;
        tok = strtok(st, " \t\n\r");
        if(tok!=NULL){
          str=tok;
          files.push_back(str);
        }
      }
      fclose(f);
      i++;

    } else { 

      //if we get here, assume an input file (i.e. no indicator flag)
      str=argv[i];
      files.push_back(str);

    }
  }

  //output file is ALWAYS the last word on the commandline.
  matrix.setOutput(argv[argc-1]);
}

void usage(){
  cout << "\nUSAGE: StPeter2Matrix [options] <protXML ...> <output file>" << endl;
  cout << "OPTIONS:" << endl;
  cout << "  -f <text>    :  filter out all protein identifications containing <text>" << endl;
  cout << "  -g <protXML> :  filter proteins using a global protXML analysis." << endl;
  cout << "  -l <file>    :  reads list of protXML input files from <file>." << endl;
  cout << "                  <file> is text-based, with one protXML name per line." << endl;
}
