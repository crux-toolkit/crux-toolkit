#include "CombineOut.h"
#include "Parsers/mzParser/mzParser.h"

#include "Common/tpp_hashmap.h" // hashmap isn't in std space for some compilers
typedef TPP_HASHMAP_T<char, char> hashmap_char_char_t;

using namespace mzParser;

CombineOut::CombineOut(char* path1, char* path2, char* outpath) : Out2XML () {
#define OPATHLEN 4086
  pepXmlPath_ = new char [OPATHLEN];
  paramsPath_ = new char [OPATHLEN];
  mzXmlPath_ = new char [OPATHLEN];
  baseName_ = new char [OPATHLEN];
  baseDir_ = new char [OPATHLEN];

  numHitsReport_ = 1; 
  if (numHitsReport_ > 20) 
    numHitsReport_ = 20;
  
  //TODO: Have to support multiple enzymes
  enzyme_ = (new ProteolyticEnzymeFactory())->getProteolyticEnzyme("trypsin"); 
  
  strcpy(outpath_, outpath);
  if (outpath_[strlen(outpath_)-1] == '/' || outpath_[strlen(outpath_)-1] == '\\') {
    outpath_[strlen(outpath_)-1] = '\0';
  }

  strcpy(inpath1_, path1);
  strcpy(inpath2_, path2);
  
  tgzfile1_ = False;
  tgzfile2_ = False;

  if (strstr(inpath1_, ".tgz") != NULL) {
     tgzfile1_ = True;
  }
  if (strstr(inpath2_, ".tgz") != NULL) {
     tgzfile2_ = True;
  }
  
  if (tgzfile1_) {
    readTgzFile(path1, outFiles1_, seqParams1_);
  }
  else {
    readOutDir(path1, outFiles1_, seqParams1_);
  }  

  if (tgzfile2_) {
    readTgzFile(path2, outFiles2_, seqParams2_);
  }
  else {
    readOutDir(path2, outFiles2_, seqParams2_);
  }

  sprintf(pepXmlPath_, "%s%s", outpath_,get_pepxml_dot_ext());
  rampConstructInputFileName(mzXmlPath_, OPATHLEN, outpath_);
  sprintf(paramsPath_, "%s/sequest.params", outpath_);
  sprintf(baseName_, "%s", outpath_);
}

CombineOut::~CombineOut() {
  //TODO There's more to clean up
  //delete [] mzXmlPath_;
  //delete [] paramsPath_;
  //delete [] baseName_;
  //delete [] baseDir_;
  //delete enzyme_;
}

void CombineOut::processData() {
  combine();
  if( (pepXmlFile_ = ogzfopen(pepXmlPath_)) == NULL) { // will write as gzip is filename ends in .gz
      cerr << "cannot open pepXML output file " << pepXmlPath_ << endl;
      exit(1);
  }
  writeData();

  //  writeXMLHeader();
  //writeOutData();
  //writeXMLFooter();
}



void CombineOut::combine() {
  SequestOutPtrMap::iterator curr1; 
  SequestOutPtrMap::iterator curr2;
  SequestOutPtrMap::iterator end1 = outFiles1_.end();
  SequestOutPtrMap::iterator end2 = outFiles2_.end();
  combineParams();
  for (curr1 = outFiles1_.begin(); curr1 != end1; curr1++) {
    if ( (curr2 = outFiles2_.find(curr1->first)) != end2 ) {
      //printf("%s\n", curr2->first);
      updateHits(curr2->second);
      SequestOut* combo = mergeOuts(curr1->second, curr2->second);
      outFiles_.insertAtEnd(combo);
    }
    else {
      outFiles_.insertAtEnd(curr1->second);
    }
  }
  for (curr2 = outFiles2_.begin(); curr2 != end2; curr2++) {
    if ( (curr1 = outFiles1_.find(curr2->first)) == end1 ) {
      updateHits(curr2->second);
      outFiles_.insertAtEnd(curr2->second);
    }
  }
}

SequestOut* CombineOut::mergeOuts(SequestOut* out1, SequestOut* out2) {
  SequestOut* tmp = new SequestOut(*out1); //TODO Need to correct the dMasses
  SequestHit* hit1;
  SequestHit* hit2;

  double topXC = -1;

  int totHits1 = out1->sequestHits_.size();
  int totHits2 = out2->sequestHits_.size();

  int rank, i, i1, i2;
  rank = 1; i = 0; i1 = 0; i2 = 0; 
  while (i < totHits1 + totHits2) {   
    hit1 = out1->getHitByIndex(i1);
    hit2 = out2->getHitByIndex(i2);

    if (topXC < 0) {
      double nextXC;
      if (hit1 != NULL && hit2 != NULL) 
	nextXC = hit1->dXC > hit2->dXC ? hit1->dXC : hit2->dXC;
      else if (hit1 != NULL) 
      	nextXC = hit1->dXC;
      else if (hit2 != NULL) 
      	nextXC = hit2->dXC;
      topXC = nextXC;
    }

    if (hit1 != NULL && (i2 >= totHits2 || hit1->dXC > hit2->dXC) ) {
      hit1->dDeltCn = topXC > 0 ? (topXC - hit1->dXC) / topXC : 9999;
      hit1->iRank = rank++;
      i1++; 
      i++;
      tmp->sequestHits_.insertAtEnd(hit1);
    }
    else if (hit2 != NULL && (i1 >= totHits1 || hit1->dXC < hit2->dXC) ) {
      hit2->dDeltCn = topXC > 0 ? (topXC - hit2->dXC) / topXC : 9999;
      hit2->iRank = rank++;
      i2++; 
      i++;
      tmp->sequestHits_.insertAtEnd(hit2);
    }
    else if (hit1 != NULL && hit2 != NULL) { //equal
      hit1->dDeltCn = topXC > 0 ? (topXC - hit1->dXC) / topXC : 9999;
      hit2->dDeltCn = topXC > 0 ? (topXC - hit2->dXC) / topXC : 9999;
      
      i1++; i2++;
      if (!strcmp(hit1->szSubPep, hit2->szSubPep)) {
	tmp->sequestHits_.insertAtEnd(hit1);
	hit1->iRank = rank++; 
	i+=2; 
	
      }
      else {
	tmp->sequestHits_.insertAtEnd(hit1);
	tmp->sequestHits_.insertAtEnd(hit2);
	hit1->iRank = rank++;
	i+=2; 
	hit2->iRank = rank++;
      }   
    } 
  }
  tmp->getDeltaCn();
  return tmp;
}

void CombineOut::updateHits(SequestOut* out) {
  // new modifications require new symbols in the peptide string
  SequestHit* hit;
  int totHits = out->sequestHits_.size();
  char* newPep = new char[SIZE_PEP];
  for (int i = 0; i < totHits ; i++) {
    hit = out->getHitByIndex(i);
    int iLenPep = (int)strlen(hit->szSubPep);
    int pepI = 0;
    Boolean subsym = False;
    for (int j = 0; j < iLenPep; j++) {
      if (!subsym) {
	newPep[pepI] = hit->szSubPep[j];
	pepI++;
      }
      subsym = False;
      if (isalpha(hit->szSubPep[j])) {
	char* sym = NULL;
	if (newModsfromVariable_[hit->szSubPep[j]] != NULL) {
	  sym = &(*newModsfromVariable_[hit->szSubPep[j]])[hit->szSubPep[j+1]];
	  if ( *sym != '\0') {
	    subsym = True;
	    newPep[pepI] = *sym;
	    pepI++;
	  }
	}
	else {
	  sym = &newModsfromStatic_[hit->szSubPep[j]];
	  if ( *sym != '\0') {
	    newPep[pepI] = *sym;
	    pepI++;
	  }
	}

	if (j == 0) {
	  sym = NULL;
	  sym = &newModsfromStatic_['n'];
	  if ( *sym != '\0') {
	    newPep[pepI] = *sym;
	    pepI++;
	  }
	}
	if (j == iLenPep-1) {
	  sym = NULL;
	  sym = &newModsfromStatic_['c'];
	  if ( *sym != '\0') {
	    newPep[pepI] = *sym;
	    pepI++;
	  }
	}
      }
    }
    newPep[pepI] = '\0';
    strcpy(hit->szSubPep, newPep);
  }
  //  printf("%s\n", newPep);
  delete [] newPep;

}

void CombineOut::writeParams() {
  if( (paramsFile_ = fopen(paramsPath_, "w")) == NULL) {
    cerr << "cannot open sequest.params output file " << paramsPath_ << endl;
    exit(1);
  }
  sequestParams_->writeParams(paramsFile_);
  fclose(paramsFile_);

}

void CombineOut::combineParams() {
  //TODO can be made faster by using hashes
  int totMods1 = seqParams1_->getNumModifiedAAs();
  int totMods2 = seqParams2_->getNumModifiedAAs();
  for (int i = 0; i < totMods1; i++) {
    //TODO For now assume that the variable mods are the same and combine only static mods
    if (!(*(seqParams1_->modifications_))[i]->variable) {
      //Found a static mod in the first file, look for a static mod on same aa in other file
      Boolean found = False;
      for (int j=0; j < totMods2; j++) {
	if (seqParams1_->getModifiedAA(i) == seqParams2_->getModifiedAA(j)) { 
	  if ((*(seqParams1_->modifications_))[i]->protein_terminus == (*(seqParams2_->modifications_))[j]->protein_terminus &&
	      DABS((*(seqParams1_->modifications_))[i]->massdiff - (*(seqParams2_->modifications_))[j]->massdiff) > 0.00001) {
	    found = True;
	    Modification* next = new Modification();
	    next->terminal = (*(seqParams1_->modifications_))[i]->terminal;
	    next->aa = (*(seqParams1_->modifications_))[i]->aa;
	    next->massdiff = (*(seqParams2_->modifications_))[j]->massdiff - (*(seqParams1_->modifications_))[i]->massdiff;
	    next->variable = True;
	    next->protein_terminus = (*(seqParams1_->modifications_))[i]->protein_terminus;
	    seqParams1_->setModificationSymbol(next, True);
	    (*(seqParams1_->modifications_)).insertAtEnd(next);
	    if ((*(seqParams2_->modifications_))[j]->variable) {
	      hashmap_char_char_t * tmp = newModsfromVariable_[(*(seqParams2_->modifications_))[j]->aa];
	      if (tmp == NULL) {
		tmp = new hashmap_char_char_t;
	      }
	      newModsfromVariable_[(*(seqParams2_->modifications_))[j]->aa] = tmp;
	      (*newModsfromVariable_[(*(seqParams2_->modifications_))[j]->aa])[(*(seqParams2_->modifications_))[j]->symbol[0]] = *next->symbol;
	    }
	    else {
	      newModsfromStatic_[(*(seqParams2_->modifications_))[j]->aa] = *next->symbol;
	    }
	  }
	  else if ((*(seqParams1_->modifications_))[i]->protein_terminus == (*(seqParams2_->modifications_))[j]->protein_terminus &&
		   DABS((*(seqParams1_->modifications_))[i]->massdiff - (*(seqParams2_->modifications_))[j]->massdiff) <= 0.00001) {
	    found = True;
	  }
	}
      }
      
      if (!found) {
	(*(seqParams1_->modifications_))[i]->variable = True;
	(*(seqParams1_->modifications_))[i]->massdiff *= -1;
	seqParams1_->setModificationSymbol((*(seqParams1_->modifications_))[i], True);
	if ((*(seqParams1_->modifications_))[i]->aa == 'n') {
	  newModsfromStatic_[(*(seqParams1_->modifications_))[i]->aa] = ']';
	}
	else if ((*(seqParams1_->modifications_))[i]->aa == 'c') {
	  newModsfromStatic_[(*(seqParams1_->modifications_))[i]->aa] = '[';
	}
	else {
	  newModsfromStatic_[(*(seqParams1_->modifications_))[i]->aa] = *(*(seqParams1_->modifications_))[i]->symbol;
	}      
      }    
    }
    else {
      // Mod1 is variable
      for (int k=0; k < totMods2; k++) {
	if (seqParams1_->getModifiedAA(i) == seqParams2_->getModifiedAA(k) && 
	    DABS((*(seqParams1_->modifications_))[i]->massdiff - (*(seqParams2_->modifications_))[k]->massdiff) > 0.00001 &&
	    (*(seqParams1_->modifications_))[i]->protein_terminus == (*(seqParams2_->modifications_))[k]->protein_terminus &&
	    (*(seqParams1_->modifications_))[i]->terminal == (*(seqParams2_->modifications_))[k]->terminal) {
	  // check for a static mod in params2
	  double staticmod = 0;
	  
	  if ((*(seqParams2_->modifications_))[k]->variable) {
	    for (int j=0; j < totMods2; j++) {
	      if (j !=  k && seqParams2_->getModifiedAA(j) == seqParams2_->getModifiedAA(k) && !(*(seqParams2_->modifications_))[j]->variable) {
		staticmod = (*(seqParams2_->modifications_))[j]->massdiff;
	      }
	    }
	  }

	  Modification* next = new Modification();
	  next->terminal = (*(seqParams1_->modifications_))[i]->terminal;
	  next->aa = (*(seqParams1_->modifications_))[i]->aa;
	  next->massdiff = (*(seqParams2_->modifications_))[k]->massdiff + staticmod;
	  next->variable = True;
	  next->protein_terminus = (*(seqParams1_->modifications_))[i]->protein_terminus;
	  seqParams1_->setModificationSymbol(next, True);
	  (*(seqParams1_->modifications_)).insertAtEnd(next);

	  if ((*(seqParams2_->modifications_))[k]->variable) {
	      hashmap_char_char_t* tmp =  newModsfromVariable_[(*(seqParams2_->modifications_))[k]->aa];
	      if (tmp == NULL) {
		tmp = new hashmap_char_char_t;
	      }
	      newModsfromVariable_[(*(seqParams2_->modifications_))[k]->aa] = tmp;
	      (*newModsfromVariable_[(*(seqParams2_->modifications_))[k]->aa])[(*(seqParams2_->modifications_))[k]->symbol[0]] = *next->symbol;
	  }
	  else {
	    newModsfromStatic_[(*(seqParams2_->modifications_))[k]->aa] = *next->symbol;
	  }
	}
	else if (seqParams1_->getModifiedAA(i) == seqParams2_->getModifiedAA(k) && 
	    DABS((*(seqParams1_->modifications_))[i]->massdiff - (*(seqParams2_->modifications_))[k]->massdiff) <= 0.00001 &&
	    (*(seqParams1_->modifications_))[i]->protein_terminus == (*(seqParams2_->modifications_))[k]->protein_terminus &&
	    (*(seqParams1_->modifications_))[i]->terminal == (*(seqParams2_->modifications_))[k]->terminal) {
	  // Must update the symbol in params2
	  if ((*(seqParams2_->modifications_))[k]->variable) {
	    hashmap_char_char_t* tmp =  newModsfromVariable_[(*(seqParams2_->modifications_))[k]->aa];
	    if (tmp == NULL) {
	      tmp = new hashmap_char_char_t;
	    }
	    newModsfromVariable_[(*(seqParams2_->modifications_))[k]->aa] = tmp;
	    (*newModsfromVariable_[(*(seqParams2_->modifications_))[k]->aa])[(*(seqParams2_->modifications_))[k]->symbol[0]] = *(*(seqParams1_->modifications_))[i]->symbol;
	  }
	  else {
	    newModsfromStatic_[(*(seqParams2_->modifications_))[k]->aa] = *(*(seqParams1_->modifications_))[i]->symbol;
	  }
	}
      }

    }
  }

 
  for (int i = 0; i < totMods2; i++) {
    //DDS Take care of mods only present in params 2
    totMods1 = seqParams1_->getNumModifiedAAs(); 
    Boolean found = False;
    for (int j = 0; j < totMods1; j++) {
      if ( seqParams2_->getModifiedAA(i) == seqParams1_->getModifiedAA(j) && 
	   (*(seqParams2_->modifications_))[i]->protein_terminus == (*(seqParams1_->modifications_))[j]->protein_terminus &&
	   (*(seqParams2_->modifications_))[i]->terminal == (*(seqParams1_->modifications_))[j]->terminal) {
	found = True;
      }           
    }
    
    if (!found) {
      Modification* next = new Modification();
      next->terminal = (*(seqParams2_->modifications_))[i]->terminal;
      next->aa = (*(seqParams2_->modifications_))[i]->aa;
      next->massdiff = (*(seqParams2_->modifications_))[i]->massdiff;
      next->variable = True;
      next->protein_terminus = (*(seqParams2_->modifications_))[i]->protein_terminus;
      seqParams1_->setModificationSymbol(next, True);
      (*(seqParams1_->modifications_)).insertAtEnd(next);

      if ((*(seqParams2_->modifications_))[i]->variable) {
	hashmap_char_char_t* tmp =  newModsfromVariable_[(*(seqParams2_->modifications_))[i]->aa];
	if (tmp == NULL) {
	  tmp = new hashmap_char_char_t;
	}
	newModsfromVariable_[(*(seqParams2_->modifications_))[i]->aa] = tmp;
	(*newModsfromVariable_[(*(seqParams2_->modifications_))[i]->aa])[(*(seqParams2_->modifications_))[i]->symbol[0]] = *next->symbol;
      }
      else {
	newModsfromStatic_[(*(seqParams2_->modifications_))[i]->aa] = *next->symbol;
      }
    }

  }
  sequestParams_ = seqParams1_;

}

void CombineOut::writeData() {
  //write XML
  writeXMLHeader();
  writeOutData();
  writeXMLFooter();
  
  //write Out files
  writeOutFiles();
}

void CombineOut::writeOutFiles() {

  FILE* outFile;
  char command[OPATHLEN];
  char fileName[OPATHLEN];
  sprintf(command, "mkdir -p %s", outpath_);
  verified_system(command);

  for (int i=0; i < outFiles_.size(); i++) {
    sprintf(fileName, "%s/%s.out", outpath_, outFiles_[i]->szBaseFileName);
    sprintf(command, "%s/%s.dta",outpath_, outFiles_[i]->szBaseFileName);
    
    if( (outFile = fopen(command, "r")) != NULL) {
      sprintf(command, "cp %s/%s.dta %s", inpath1_, outFiles_[i]->szBaseFileName, outpath_);
      verified_system(command);
      fclose(outFile);
    }
    
    if( (outFile = fopen(fileName, "w")) == NULL) {
      cerr << "WARNING: cannot open output file " << fileName << endl;
    }
    else {
      outFiles_[i]->writeOutFile(outFile, sequestParams_);
      fclose(outFile);
    }
  }
  writeParams();

}


void  CombineOut::readOutDir(char* path, SequestOutPtrMap& outFiles, SequestParams*& seqParam) 
{
  struct stat fileStat;
#ifdef _MSC_VER
   struct _finddata_t c_file;
   long hFile;
   char szFileMask[SIZE_BUF];
#else
  struct dirent *dirEntry;
  DIR *dirStream;
#endif
  char * pathName = new char [OPATHLEN];
  char * fileName = new char [OPATHLEN];
  char * paramsFile = new char [OPATHLEN];

  if (!stat(path, &fileStat))
    {
      if (S_ISDIR(fileStat.st_mode))
	{
	  // Since this is a directory name and we
	  // will be appending filenames to it do we
	  // have a separator at the end?
	  if (!isPathSeperator(*(path + strlen(path))))
            {
	      // Append separator
	      sprintf(pathName, "%s/", path);
            }
	  else
            {
	      strcpy(pathName, path);
            }
	  
	  // sprintf(baseName_, "%s", pathName);
	  sprintf(baseDir_, "%s", pathName);
	  
	  // baseName_[strlen(baseName_)-1] = '\0';

	  // sprintf(pepXmlPath_, "%s%s", baseName_,get_pepxml_dot_ext());
	  // sprintf(mzXmlPath_, "%s.mzXML", baseName_);

	  //TODO: need an option for this 
	  sprintf(paramsFile, "%ssequest.params", baseDir_); 
	  seqParam = new SequestParams(paramsFile);
	  
	  //if( (pepXmlFile_ = fopen(pepXmlPath_, "w")) == NULL) {
	  //  cerr << "cannot open pepXML output file " << pepXmlPath_ << endl;
	  //  exit(1);
	  //}
	  const char *fname;

#ifdef _MSC_VER
	  sprintf(szFileMask,"%s*.out",pathName);
	  hFile = _findfirst( szFileMask, &c_file );
	  if (hFile!=-1) {
		  do {
			  fname = c_file.name;
#else
	  dirStream = (DIR *) opendir(pathName);
	  // Process all files in the directory
	  while ((dirEntry = readdir((DIR *) dirStream)))
            {
	      // If it contains ".out" process it
	      if (strstr(fname=dirEntry->d_name, ".out"))
		{
#endif
		  //struct SequestOutStruct* tmp = (SequestOutStruct *)malloc(sizeof(SequestOutStruct));
		  SequestOut* tmp = new SequestOut();
		  struct HeaderStruct * hdr = (HeaderStruct *)malloc(sizeof(HeaderStruct));
                  sprintf(fileName, "%s%s", pathName, fname);
		  this->readOutFile(fileName, tmp, hdr, seqParam);
		  char* key = new char[OPATHLEN];
		  strcpy(key, fname);
		  SequestOutPtrMap::iterator curr; 
		  SequestOutPtrMap::iterator end = outFiles.end();
		  outFiles.insert(pair<const char*, SequestOut*>(key, tmp));
#ifndef _MSC_VER
               }
            }
#else
		  } while( _findnext( hFile, &c_file ) == 0 );
	  }
#endif
	}
      else 
	{
	  //Do nothing right now
	  //TODO: do something meaningful
	}
      
    }
  delete [] pathName;
  delete [] fileName;
  delete [] paramsFile;
}


void  CombineOut::readTgzFile(char* path, SequestOutPtrMap& outFiles, SequestParams*& seqParam) 
{
  struct stat fileStat;
  FILE *tgzpipe, *outpipe;
#ifdef __MINGW__
  char SEP = '\\';
#else
  char SEP = '/';
#endif
  char * pathName = new char [OPATHLEN];
  char * fileName = new char [OPATHLEN];
  char * paramsFile = new char [OPATHLEN];
  char * command = new char [OPATHLEN];
  char * next = new char [OPATHLEN];
  if (!stat(path, &fileStat))
    {
      if (S_ISREG(fileStat.st_mode))
	{
	  sprintf(command, "tar tfz %s *.out", path);
	  if ((tgzpipe = tpplib_popen(command, "r"))!=NULL) {
	    while (fgets(next, OPATHLEN, tgzpipe)!=NULL) {
	      next[strlen(next)-1] = '\0';
	      sprintf(command, "tar Oxfz %s %s", path, next);
	      if ((outpipe = tpplib_popen(command, "r"))!=NULL) {
		SequestOut* tmp = new SequestOut();
		struct HeaderStruct * hdr = (HeaderStruct *)malloc(sizeof(HeaderStruct));
		this->readOutFile(outpipe, tmp, hdr, seqParam);
		char* key = new char[OPATHLEN];
		sprintf(key, "%s", next);
		SequestOutPtrMap::iterator curr; 
		SequestOutPtrMap::iterator end = outFiles.end();
		outFiles.insert(pair<const char*, SequestOut*>(key, tmp));
	      }
	      pclose(outpipe);
	    }
	  }
	  else
	    {
	      //Do nothing right now
	      //TODO: do something meaningful
	    }	
	}    
      else 
	{
	  //Do nothing right now
	  //TODO: do something meaningful
	}
      
    }
  delete [] next;
  delete [] command;
  delete [] pathName;
  delete [] fileName;
  delete [] paramsFile;
}
