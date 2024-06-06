#include "Parsers/mzParser/mzParser.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#ifdef _MSC_VER
#include <direct.h>
#include <Windows.h>
#define getcwd _getcwd
#define slashdir '\\'
#else
#include <unistd.h>
#define slashdir '/'
#endif

using namespace std;
using namespace mzParser;

//a sortable index structure
typedef struct sIndex{
  string id;
  int scanNum;
} sIndex;

bool compareID(const sIndex& a, const sIndex& b);
string extractMzML(string s, MzParser& m);
size_t findIndex(string& s);
void fixMe(const char* in, const char* out, MzParser& m);
string fixMPA(string s, MzParser& m);
string fixQuery(string s, MzParser& m);

vector<sIndex> gIndex;
string outFile;

int main(int argc, char* argv[]){
  cout << "FixPepXMLIndex v1.0, March 26 2020, Mike Hoopmann, ISB" << endl;
  if(argc<3) {
    cout << "USAGE:" << endl;
    cout << "  FixPepXMLIndex [old pepXML file] [new pepXML file]" << endl;
    return 1;
  }

  BasicSpectrum s;
  MzParser m(&s);
  outFile=argv[2];
  fixMe(argv[1],argv[2],m);

  return 0;
}

bool compareID(const sIndex& a, const sIndex& b) {
  return (a.id.compare(b.id) < 0);
}

string extractMzML(string s, MzParser& m){

  //parse the run summary
  vector<string> tokens;
  string str;
  char tmp[10000];
  char* tok;
  strcpy(tmp,s.c_str());
  tok=strtok(tmp," =<>\"\n\r\t");
  while(tok!=NULL){
    str=tok;
    tokens.push_back(str);
    tok=strtok(NULL, " =<>\"\n\r\t");
  }

  //extract file name
  string fname;
  for(size_t i=0;i<tokens.size();i++){
    if(tokens[i].compare("base_name")==0) fname=tokens[++i];
    else if(tokens[i].compare("raw_data")==0) fname+=tokens[++i];
  }

  //read the file
  if (!m.load(fname.c_str())) {
    //temporary for testing
    fname=fname.substr(fname.find_last_of('/')+1,fname.size());
    if(!m.load(fname.c_str())){
    cout << "ERROR: Failure reading input file " << fname << endl;
    exit(-1);
    }
  }

  //create our sortable index
  gIndex.clear();
  sIndex x;
  vector<cindex>* index = m.getSpectrumIndex();
  for(size_t i=0;i<index->size();i++){
    x.id=index->at(i).idRef;
    x.scanNum=(int)i+1;
    gIndex.push_back(x);
  }
  sort(gIndex.begin(),gIndex.end(),compareID);

  //return the filename for diagnostics
  return fname;

}

size_t findIndex(string& s){
  size_t sz = gIndex.size();
  size_t lower = 0;
  size_t mid = sz / 2;
  size_t upper = sz;
  int i;

  if (sz > 0) {
    i = gIndex[mid].id.compare(s);
    while (i != 0) {
      if (lower >= upper) break;
      if (i > 0) {
        if (mid == 0) break;
        upper = mid - 1;
        mid = (lower + upper) / 2;
      } else {
        lower = mid + 1;
        mid = (lower + upper) / 2;
      }
      if (mid == sz) break;
      i = gIndex[mid].id.compare(s);
    }

    if (i == 0) return mid;
  }

  //SIZE_MAX indicates index not found
  return SIZE_MAX;
}

void fixMe(const char* in, const char* out, MzParser& m){
  FILE* fin=fopen(in,"rt");
  FILE* fout=fopen(out,"wt");

  char str[10000];
  string s;

  string mzml;
  string newQuery;
  while(!feof(fin)){
    if(fgets(str,10000,fin)==NULL) continue;
    s=str;
    newQuery.clear();
    if(s.find("<msms_run_summary")!=string::npos) mzml=extractMzML(s,m);
    else if (s.find("<msms_pipeline_analysis") != string::npos) newQuery = fixMPA(s, m);
    else if(s.find("<spectrum_query")!=string::npos) newQuery=fixQuery(s,m);
    
    if(newQuery.size()==0) fprintf(fout,"%s",str);
    else fprintf(fout,"%s",newQuery.c_str());
  }

  fclose(fout);
  fclose(fin);

}

string fixMPA(string s, MzParser& m) {
  char* ret;
  char cwd[1024];
  string st;

  //set the new summary_xml string
  ret = getcwd(cwd, 1024);
  st="summary_xml=\"";
  st+=cwd;
  st+=slashdir;
  st+=outFile;
  st+="\"";

  //find location of summary_xml in current string
  size_t start=s.find("summary_xml");
  int count=0;
  size_t i;
  for(i=start;i<s.size();i++){
    if(s[i]=='\"') count++;
    if(count==2) break;
  }
  size_t end=i;

  //swap and return updated string
  s.replace(start,end-start+1,st);
  return s;
}

string fixQuery(string s, MzParser& m) {
  //parse the run spectrum_query
  //count leading spaces
  int spaces=0;
  for(size_t i=0;i<s.size();i++){
    if(s[i]==' ') spaces++;
    else break;
  }
  vector<string> tokens;
  string word;
  bool bInQuotes=false;
  for(size_t i=0;i<s.size();i++){
    if(s[i]==' '){
      if(word.size()>0 && word[0]=='<') {
        tokens.push_back(word);
        word.clear();
      } else if(bInQuotes){
        word+=s[i];
      } else {
        //do nothing
      }
    } else if(s[i]=='='){
      word+=s[i];
      if(!bInQuotes) {
        tokens.push_back(word);
        word.clear();
      }
    } else if(s[i]=='\"'){
      if(!bInQuotes){
        if(word.size()>0){
          tokens.push_back(word);
          word.clear();
        }
        word+=s[i];
        bInQuotes=true;
      } else {
        word+=s[i];
        tokens.push_back(word);
        word.clear();
        bInQuotes=false;
      }
    } else {
      word+=s[i];
    }
  }
  if(word.size()>0) tokens.push_back(word);

  //mark important tokens
  size_t idIndex=0;
  size_t specIndex=0;
  size_t startIndex=0;
  size_t endIndex=0;
  for(size_t i=0;i<tokens.size();i++){
    if(tokens[i].compare("spectrum=")==0){
      specIndex=++i;
      continue;
    } else if(tokens[i].compare("spectrumNativeID=") == 0) {
      idIndex = ++i;
      continue;
    } else if (tokens[i].compare("start_scan=") == 0) {
      startIndex = ++i;
      continue;
    } else if (tokens[i].compare("end_scan=") == 0) {
      endIndex = ++i;
      continue;
    }
  }

  //find the idIndex in the mzML file index
  //remove quotes
  string str=tokens[idIndex].substr(1,tokens[idIndex].size()-2);
  size_t index=findIndex(str);
  if(index==SIZE_MAX) {
    cout << "Error" << endl;
    exit(1);
  }

  //make new strings
  string newSpec;
  int counter=0;
  char cs[64];
  for(size_t i=0;i<tokens[specIndex].size();i++){
    if(tokens[specIndex][i]=='.'){
      sprintf(cs,".%06d.%06d.",gIndex[index].scanNum, gIndex[index].scanNum);
      newSpec+=cs;
      counter++;
      while(counter<3){
        i++;
        if(tokens[specIndex][i] == '.') counter++;
      }
    } else {
      newSpec+=tokens[specIndex][i];
    }
  }

  string newStart;
  sprintf(cs,"\"%d\"",gIndex[index].scanNum);
  newStart=cs;

  //note that the new end_scan is currently assumed the same as the start scan
  string newEnd;
  sprintf(cs, "\"%d\"", gIndex[index].scanNum);
  newEnd = cs;

  string newQuery;
  for(int i=0;i<spaces;i++) newQuery+=' ';
  for(size_t i=0;i<tokens.size()-2;i++){
    if(i==specIndex) newQuery+=newSpec;
    else if(i==startIndex) newQuery+=newStart;
    else if(i==endIndex) newQuery+=newEnd;
    else newQuery+=tokens[i];
    if(tokens[i][tokens[i].size()-1]!='=') newQuery+=" ";
  }
  newQuery+=tokens[tokens.size()-2];
  newQuery+=tokens[tokens.size()-1];

  return newQuery;

}
