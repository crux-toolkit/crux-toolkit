#include "QuanticParser.h"

/*

Program       : Quantic
Author        : David Shteynberg <dshteynb  AT systemsbiology.org>
Date          : 04.20.2018
SVN Info      : $Id$

Copyright (C) 2018 David Shteynberg

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

David Shteynberg
Institute for Systems Biology
401 Terry Avenue North
Seattle, WA  98109  USA

*/
#include <sstream>
#include <string>
#ifndef ALL_AA
#define ALL_AA "GASPVTCLINDQKEMHFRYWOXBZ"
#endif

typedef struct str_Qdata
{
  int thread_no;
  int _index;
  QuanticParser* par;
} Qdata;


#ifdef MSVC
DWORD WINAPI ParserReadThread(LPVOID ptr) 
#else
void* ParserReadThread(void* ptr) 
#endif
{
  int step = 0;
  int tot = 0;
  int offset = -1;
  Qdata *data;            
  int inc = 1;
  int b = 0;



  data = (Qdata *) ptr; 
  inc = data->par->max_threads_;
  offset = data->thread_no;
  b = data->_index;  
 

  // 
  
  data->par->parseRead(data->par->file_, offset);
 
#ifdef MSVC
   ExitThread(0);
#else
   pthread_exit(0); 
#endif
}

#ifdef MSVC
DWORD WINAPI ParserWriteThread(LPVOID ptr) 
#else
void* ParserWriteThread(void* ptr) 
#endif
{
  int step = 0;
  int tot = 0;
  int offset = -1;
  Qdata *data;            
  int inc = 1;
  int b = 0;

  data = (Qdata *) ptr; 
  inc = data->par->max_threads_;
  offset = data->thread_no;
  b = data->_index;  

  // 
  
  data->par->parseWrite(data->par->file_, offset);
 
#ifdef MSVC
   ExitThread(0);
#else
   pthread_exit(0); 
#endif
}




//QuanticParser::QuanticParser(double mztol, double ppmtol, double massdiff_offset, unsigned int verbose) : Parser("quantic") {
QuanticParser::QuanticParser(double mztol, double ppmtol,  unsigned int verbose) : Parser("quantic") {
  update_mod_tags_ = false;
  verbose_ = verbose;
  keepOld_ = false;
  annotate_ = false;
  pep_coupled_ = false;

  diaMode_ = false;

  modstring_ = "";
  //  massOffset_ = massdiff_offset;
  total_mods_ = 0; 
  total_sites_ = 0;

  Peptide::defaultTables();  
  neutlosses_.push_back(new vector<double>());
  mzTol_ = mztol;
  ppmTol_ = ppmtol;
  out_file_ = "";

  string model = "Quantic_Oscore";
  ptm_model_ = new KDModel(model.c_str()); 

  model = "Quantic_Mscore";
  mat_model_ = new KDModel(model.c_str());

  //  SQoffsets_ = new stream_vec();
#ifdef MSVC
  _mutex = CreateMutex( 
			    NULL,              // default security attributes
			    FALSE,             // initially not owned
			    NULL);             // unnamed mutex
  
  if (_mutex == NULL) 
    {
      printf("CreateMutex error: %d\n", (int)GetLastError());
      exit(1);
    }
#else
  pthread_mutex_init(&_mutex, NULL);

#endif

#ifdef MSVC
  run_mutex = CreateMutex( 
			    NULL,              // default security attributes
			    FALSE,             // initially not owned
			    NULL);             // unnamed mutex
  
  if (run_mutex == NULL) 
    {
      printf("CreateMutex error: %d\n", (int)GetLastError());
      exit(1);
    }
#else
  pthread_mutex_init(&run_mutex, NULL);

#endif
  
#ifdef MSVC
  hdr_mutex = CreateMutex( 
			    NULL,              // default security attributes
			    FALSE,             // initially not owned
			    NULL);             // unnamed mutex
  
  if (hdr_mutex == NULL) 
    {
      printf("CreateMutex error: %d\n", (int)GetLastError());
      exit(1);
    }
#else
  pthread_mutex_init(&hdr_mutex, NULL);

#endif

#ifdef MSVC
  fout_mutex = CreateMutex( 
			    NULL,              // default security attributes
			    FALSE,             // initially not owned
			    NULL);             // unnamed mutex
  
  if (fout_mutex == NULL) 
    {
      printf("CreateMutex error: %d\n", (int)GetLastError());
      exit(1);
    }
#else
  pthread_mutex_init(&fout_mutex, NULL);

#endif

#ifdef MSVC
  sq_mutex = CreateMutex( 
			    NULL,              // default security attributes
			    FALSE,             // initially not owned
			    NULL);             // unnamed mutex
  
  if (sq_mutex == NULL) 
    {
      printf("CreateMutex error: %d\n", (int)GetLastError());
      exit(1);
    }
#else
  pthread_mutex_init(&sq_mutex, NULL);

#endif
}
void YIELD()
{
#ifdef MSVC
		       SwitchToThread();
#else
		       sched_yield();
#endif
}
void LOCK(
#ifdef MSVC
	  HANDLE* mutex
#else
	  pthread_mutex_t* mutex
#endif
	  )  {
  
#ifdef MSVC
  WaitForSingleObject( *mutex,    // handle to mutex
		       INFINITE);      //#include "windows.h"
#else	
  pthread_mutex_lock( &*mutex );//#include <pthread.h>
#endif 
}
int TRYLOCK(
#ifdef MSVC
	  HANDLE* mutex
#else
	  pthread_mutex_t* mutex
#endif
	  )  {
  
#ifdef MSVC
  DWORD rtn = WaitForSingleObject( *mutex,    // handle to mutex
		       0);      //#include "windows.h"

  if (rtn == WAIT_OBJECT_0) {
    return 1;
  }

  return 0
  
#else	
  return pthread_mutex_trylock( &*mutex );//#include <pthread.h>
#endif 
}

void UNLOCK(
#ifdef MSVC
	    HANDLE* mutex
#else
	    pthread_mutex_t* mutex
#endif
	    )  {
#ifdef MSVC
  ReleaseMutex( *mutex);      //#include "windows.h"
#else
  pthread_mutex_unlock( &*mutex );//#include <pthread.h>
#endif
}

void QuanticParser::setOutFile(string name) {
  out_file_ = name;
}

void QuanticParser::setModString(string& modstring) {
  int start = 0;
  //  int comma = modstring.find(",", start);
  int sep = modstring.find(":", start);
  string aa = "";
  double shift = 0.;
  double loss = 0.;

 
  int last_start = 0;
 
  //  SQoffsets_ = new stream_vec();
   
  while (sep != string::npos) {

  
    ptmtypes_.push_back( modstring.substr(start) );
     

    aa = modstring.substr(start, sep-start);

    //start = comma + 1;
    

    start = sep + 1;
    

    if (sep != string::npos) {
      shift = atof(modstring.substr(start, sep-start).c_str());
      last_start = start;
      start = sep + 1;

      aminoacids_.push_back(aa);
      massshift_.push_back(shift);
      //neutlosses_.push_back(new vector<double>());


      sep = modstring.find(":", start);    
      if (sep <0) {
	break;
      }
      else if (sep >= 0)  {
	start = sep+1;
	sep = modstring.find_first_of(":", start);    
	while (sep >= 0)  {
	  loss = atof(modstring.substr(start, sep-start).c_str());
	  last_start = start;
	  start = sep + 1;
	  neutlosses_[neutlosses_.size()-1]->push_back(loss);
	  sep = modstring.find(":", start);
	  if (sep == -1)
	    break;
	}
	
	loss = atof(modstring.substr(start, sep-start).c_str());
	neutlosses_[neutlosses_.size()-1]->push_back(loss);
      }
    
  
      break;
      
    }
    
    sep = modstring.find(":", start);
    
    
  }
  for (int x=0; x<neutlosses_.size(); x++) {
    sort(neutlosses_.begin(), neutlosses_.end());
  }
  modstring_ = modstring;
}

void QuanticParser::setUpdate(bool up) { 
  update_mod_tags_ = up; 
}

void QuanticParser::run(const char* c, const char* opts, int max_threads, unsigned int top_peaks) {
  opts_ = string(opts);
  
  time_ = getDateTime(USE_LOCAL_TIME);
  success_ = false;

  top_peaks_ = top_peaks;
  if(time_ == NULL) {
    cout << "error: null analysis time" << endl;
    exit(1);
  }  
  filter_ = False;
  filter_memory_ = False;

  max_threads_ = max_threads < 1 ? 1 : max_threads;
  header_thread_ = max_threads_;
  if(c != NULL) {
    file_ = new char[strlen(c)+1];
    strcpy(file_, c);
  }
  else {
    file_ = NULL;
  }


  if (out_file_.empty()) {
    tmp_file_ = make_tmpfile_name(c);
  }
  else {
    tmp_file_ = out_file_;
  }
 
  fout.open(tmp_file_.c_str());  
  
  if (spec_prob_hash_.size() == 0) {
    parseSpecProbs(c);
  }


  if(! fout) {
    cerr << "cannot write output to file " << tmp_file_ << endl;
    exit(1);
  }
  else {
    cerr << "INFO: Writing file " << tmp_file_ << " ..." << endl;
  }

  std::cerr << "\nCreating " << max_threads_ << " threads " << endl;
  std::cerr << "Wait for threads to finish ..." << endl;
  //  std::cerr << "0------------------------------------------------50"
  //	    << "------------------------------------------------100%" << endl;     

  

  std::cerr << std::flush;

#ifdef MSVC
    DWORD *pId = new DWORD[max_threads_];
    HANDLE *pHandle = new HANDLE[max_threads_];
#else
    int *pId = new int[max_threads_];
    int *pHandle = new int[max_threads_];
    pthread_t pThreads[max_threads_];
#endif

    Qdata data[max_threads_];
    
  int a = 0;
  int b = 0; 
  int live_threads = 0;
  int done = 0;
  
  //DISABLED
  while (0 && !done) { 
    while (a < max_threads_) { // && b < ms1_lines) { 
      
      data[a].thread_no = a;
      data[a].par = this;
      data[a]._index = b;
      
#ifdef MSVC
      pHandle[a] = CreateThread(NULL,0,ParserReadThread,
				(void*) &data[a],0, 
				NULL);
#else
      pthread_create(&pThreads[a],NULL,ParserReadThread, 
		     (void*) &data[a]);
#endif
      live_threads++;
      a++;
      //b++;
    }
    a = 0;
    while(a < live_threads)  {
      void * ignore = 0;
#ifdef MSVC
      WaitForSingleObject(pHandle[a],INFINITE);
#else
      pthread_join(pThreads[a],&ignore);
#endif
      a++;
    }
    std::cerr << endl << "INFO: done ... " << endl;   
    live_threads = 0;
    done = 1;
  }

  std::cerr << "INFO: done ... " << endl;   
  done = 0;
  a = 0;
  current_run_ = 1;
  threads_done_run_ = 0;
  threads_start_run_ = 0;
  header_thread_ = max_threads_;
 std::cerr << "\nCreating " << max_threads_ << " threads " << endl;
  std::cerr << "Wait for threads to finish ..." << endl;
  while (!done) {
    while (a < max_threads_) { // && b < ms1_lines) { 
      
      data[a].thread_no = a;
      data[a].par = this;
      data[a]._index = b;

#ifdef MSVC
      pHandle[a] = CreateThread(NULL,0,ParserWriteThread,
				(void*) &data[a],0, 
				NULL);
#else
      pthread_create(&pThreads[a],NULL,ParserWriteThread, 
		     (void*) &data[a]);
#endif
      live_threads++;
      a++;
      //b++;
    }
    a = 0;
    while(a < live_threads)  {
      void * ignore = 0;
#ifdef MSVC
      WaitForSingleObject(pHandle[a],INFINITE);
#else
      pthread_join(pThreads[a],&ignore);
#endif
      a++;
    }
    live_threads=0;
    done = 1;
  }
  //Closing tags
  Tag* tag;

  // = new Tag("msms_run_summary", false, true);
  // tag->write(fout);
  //delete tag;
  tag = new Tag("msms_pipeline_analysis", false, true);
  tag->write(fout);
  delete tag;
  fout.close();

    
  if (out_file_.empty()) {
    //unlink(c); 
    
    if (copy_file(tmp_file_.c_str(), c)) {
      unlink(tmp_file_.c_str());
    }
    cerr << endl << "INFO: Renaming output file to " << c << " ..." << endl;
  }
  
  std::cerr << "INFO: done ... " << endl;   
  
}






void QuanticParser::parseWrite(const char* c, int t) {
      parseWriteUpdate(c, update_mod_tags_, t);
}

void QuanticParser::parseRead(const char* c, int t) {

  Array<Tag*>* tags = NULL;
  Tag* tag = NULL;

  char *nextline = new char[line_width_];
  char* data = NULL;


  string spectrum_name = "";
  string exp_lbl = "";
  string charge = "";
  double prob = 0;
  Array<double>* allntt_prob=NULL ;
  double calcnmass = -1;
  double rt = -1;
  string pep_seq = "";
  string mod_pep = "";
 
  string mod_tags = "";

  bool in_mod= false;
  
  TPP_HASHMAP_T<int, double> pos_mod_hash;

  bool is_nterm_pep= false;
  bool is_cterm_pep= false;
   
  bool get_pep = false;
  int k=0;

  long scan = -1;
  vector<string> mtvec;
  SpectraSTLib* mtlib= NULL;
  //  SpectraSTCreateParams* mtparams=NULL;
  //SpectraSTPepXMLLibImpo.rter*  specLib = new SpectraSTPepXMLLibImporter(mtvec, mtlib, *mtparams);
  SpectraSTLibEntry* entry = NULL;
  PeptideUser* specPep = NULL;
  Quantic* quant = NULL ;

  string QuanticReport = "";
  cRamp* cramp = NULL;

  string dataFile = "";
  string dataExt = "";
  
  bool has_mod= false;

  unsigned long SQcount = 0;
  unsigned long SQtot = SQtot_;

  int hit_rank = 0;
  double massdiff = 0;
  
  //  tmp_file_ = make_tmpfile_name(c);
  //  ofstream fout(tmp_outfile_.c_str());  
  //  if(! fout) {
  //    cerr << "cannot write output to file " << tmp_file_ << endl;
  //    exit(1);
  //  }

  // TODO  double allntt_prob[3] = {-100, -100, -100};
  // for(k = 0; k < input_files_->length(); k++) {
  RACI fin(c); // read possibly gzipped files
 
  if(! fin) {
    cerr << "fin: error opening " << c << endl;
    exit(1);
  }
  
  cerr << "INFO: Reading file " << c << " ..." << endl;
  bool nextoffset = false;

  for (size_t X=t; X < SQoffsets_.size(); X+=max_threads_) {     
    fin.seekg(SQoffsets_[X]);
    SQcount = X+1;
    nextoffset = false;
    while(fin.getline(nextline, line_width_)) {
      if (nextoffset) break;
   

    data = strstr(nextline, "<");
    while(data != NULL) {
      tag = new Tag(data);
      //tag->write(fout);
      if (in_mod) {
	mod_tags += data;
      }
      if (tag != NULL) {	

	  if (tag->isStart() && 
		   ! strcmp(tag->getName(), "spectrum_query")) {
	    //	    get_pep = true;
	    //	    SQcount++;
	    QuanticReport = "";
	    if ((SQcount-1)%max_threads_==t && SQcount%10 == 0)
	      cerr << "\rINFO: processed " << SQcount << "/" << SQtot << " spectrum_queries\t\t ";

	    pep_seq = "";
	    mod_pep = "";
	    pos_mod_hash.clear();;
	    spectrum_name = "";
	    rt = -1;
	    if (tag->getAttributeValue("retention_time_sec") != NULL) 
	      rt = atof(tag->getAttributeValue("retention_time_sec"));

	    // store the spectrum name for later
	    int len = strlen(tag->getAttributeValue("spectrum"));
	    
	    spectrum_name += tag->getAttributeValue("spectrum");
	    
	    if (tag->getAttributeValue("experiment_label") != NULL) {
	      len = strlen(tag->getAttributeValue("experiment_label"));
	    }
	    else {
	      len = 0;
	    }
	    exp_lbl = "";
	    if (len > 0)
	      exp_lbl += tag->getAttributeValue("experiment_label");
	    //char* tmp = strrchr(spectrum_name, '.');
	    //tmp = '\0';
	    charge = tag->getAttributeValue("assumed_charge");
	    scan = atoi(tag->getAttributeValue("start_scan"));
	    
	  }
	  else if ((SQcount-1) % max_threads_ == t && tag->isStart() && 
		   ! strcmp(tag->getName(), "search_hit") && 
		   atoi(tag->getAttributeValue("hit_rank")) == 1) {
	    hit_rank = 1;
	
	    pep_seq = "";
	    pep_seq += tag->getAttributeValue("peptide");
	    
	    mod_pep = pep_seq;
	    
	    calcnmass = atof(tag->getAttributeValue("calc_neutral_pep_mass"));
	    //get_pep = false;
	    massdiff  = atof(tag->getAttributeValue("massdiff"));

	    is_nterm_pep = false;
	    is_cterm_pep = false;

	    if (! strcmp(tag->getAttributeValue("peptide_prev_aa") , "-")) {
	      is_nterm_pep = true;
	    }
	    if (! strcmp(tag->getAttributeValue("peptide_next_aa") , "-")) {
	      is_cterm_pep = true;
	    }
	    
	  }
	  else if ((SQcount-1) % max_threads_ == t && tag->isStart() && 
		   ! strcmp(tag->getName(), "search_hit") && 
		   atoi(tag->getAttributeValue("hit_rank")) != 1) {
	    hit_rank = atoi(tag->getAttributeValue("hit_rank"));
	  }
	  else if ( (SQcount-1) % max_threads_ == t && tag->isStart() && ! strcmp(tag->getName(), "peptideprophet_result" )  ) {
	    double pepprob = atof(tag->getAttributeValue("probability"));
	 


	  }
	  else if ( (SQcount-1) % max_threads_ == t && hit_rank == 1 && ( (tag->isEnd() && //get_pep &&
									   ! strcmp(tag->getName(), "modification_info") )|| QuanticReport.empty() && ! strcmp(tag->getName(), "search_score")  ) ) {
	    //	    || ( tag->isEnd() &&  ! strcmp(tag->getName(), "search_hit") && unknown_mode_ && !has_mod)) ) {      
	    
	    //DDS: TODO IT HERE!!!!
	    // specLib->SpectraSTPepXMLLibImporter::addModificationsToPeptide(pep_seq, mod_tags);
	   
	    if (tag->isStart()) {
	      mod_pep = "";
	      if (tag->getAttributeValue("modified_peptide")) 
		mod_pep += tag->getAttributeValue("modified_peptide");
	    }
	    if (mod_pep == "") {
	      mod_pep = pep_seq;
	    }
	    
	
	    if (spec_prob_hash_.find(exp_lbl+spectrum_name) != spec_prob_hash_.end() && spec_prob_hash_.find(exp_lbl+spectrum_name)->second  >= minProb_) {


	      if (quant) {
		delete quant;
		quant = NULL;
	      }
	      quant = new Quantic(mod_pep, atoi(charge.c_str()), calcnmass, cramps_[X], scan,
				  aminoacids_, massshift_, neutlosses_, 
					  &stat_mods_hash_, &var_mods_hash_, &stat_prot_termods_hash_, 
				  &var_prot_termods_hash_,&_mutex,is_nterm_pep, is_cterm_pep, pos_mod_hash);
	      // unknown_mode_, labile_mode_, direct_mode_);
	      quant->setMolForms(ptmMolForms_);
	      if (quant->init()) {
		quant->processPeakList();
		quant->setVerbose(verbose_);
		//if (mzTol_ > massdiff) {
		quant->setMZTolerance(mzTol_);
		quant->setTopPeaks(top_peaks_);
		quant->setPeptideCoupled(pep_coupled_);

		QuanticReport = "yes";
		  //}
		  //else {
		  //quant->setMZTolerance(massdiff);
		  //}

		  //quant->setPPMTolerance(ppmTol_);
	      }
	      else {
#ifdef MSVC
	      ReleaseMutex( _mutex);      //#include "windows.h"
#else
	      pthread_mutex_unlock( &_mutex );//#include <pthread.h>
#endif	      
		if (mod_pep != "") {
		  cerr << "ERROR: Cannot initialize for sequence: " << mod_pep.c_str() 
		       << ", unknown mods may exist in spectrum " << spectrum_name.c_str() << endl;
		  cerr << "Please specify modification and rerun Quantic ..." << endl;
		  exit(1);
		}

		delete quant;
		quant = NULL;
	      }
	    }
	    //	    quant->evaluateModSites(mzTol_);

	    mod_tags = "";
	    in_mod = false;
	  }
	  else if ( (SQcount-1) % max_threads_ == t && tag->isStart() && !tag->isEnd() && //get_pep &&
		   ! strcmp(tag->getName(), "modification_info")) {      
	    in_mod = true;
	    has_mod = true;
       	    mod_tags = data; 
	    //	    get_pep = false;
	    mod_pep = "";
	    if (tag->getAttributeValue("modified_peptide"))
	     mod_pep += tag->getAttributeValue("modified_peptide");

	    
	    double mass = 0;
	    if (tag->getAttributeValue("mod_nterm_mass")) {
	      mass = atof(tag->getAttributeValue("mod_nterm_mass"));
	      pos_mod_hash.insert(make_pair(0, mass));
	    }
	    if (tag->getAttributeValue("mod_cterm_mass")) {
	      mass = atof(tag->getAttributeValue("mod_cterm_mass"));
	      pos_mod_hash.insert(make_pair(pep_seq.length()+1, mass));
	    }

	  }
	  else if ( (SQcount-1) % max_threads_ == t && tag->isStart() && 
		   ! strcmp(tag->getName(), "mod_aminoacid_mass")) {      
	    int pos = atol( tag->getAttributeValue("position") );
	    double mass = atof( tag->getAttributeValue("mass") );
			   
	    pos_mod_hash.insert(make_pair(pos, mass));
	  }


	  else if (tag->isEnd() && 
		   ! strcmp(tag->getName(), "spectrum_query")) {
	    has_mod = true;
	    exp_lbl = "";
	    spectrum_name = "";
	    charge = "";
	    //	    get_pep = false;
	    pep_seq = "";
	    mod_pep = "";
	    prob = 0;
	    scan = -1;
	    if ((SQcount-1) % max_threads_ == t && quant) {
	      delete quant;
	      quant = NULL;
	    }
	    nextoffset = true;
	     // TODO  allntt_prob[0] = -100; allntt_prob[1] = -100; allntt_prob[2] = -100;
	    hit_rank = 0;
	  }
		 	
	}
	delete tag;
	tag = NULL;
	data = strstr(data+1, "<");      
      }
    }
  }
  fin.close();
    //    fout.close();

    //unlink(c); 
    //rename(tmp_file_.c_str(), c);

    delete [] nextline;
    //    inter_proph_->computeModels();
}




void QuanticParser::parseSpecProbs(const char* c) {
  Array<Tag*>* tags = NULL;
  Tag* tag = NULL;

  char *nextline = new char[line_width_];
  char* data = NULL;


  string spectrum_name = "";
  string exp_lbl = "";
  string charge = "";
  double spec_prob = 0;
  Array<double>* allntt_prob=NULL ;
  double calcnmass = -1;
  double rt = -1;
  string pep_seq = "";
  string mod_pep = "";
  unsigned int run = 0;
  streampos lastpos;
  // string mod_tags = "";
  //  Array<Tag*>* mod_tags = new Array<Tag*>();//NULL;

  bool in_mod= false;

  TPP_HASHMAP_T<int, double> pos_mod_hash;


  bool get_pep = false;
  int k=0;

  long scan = -1;
  vector<string> mtvec;
  SpectraSTLib* mtlib= NULL;
  //SpectraSTCreateParams* mtparams=new SpectraSTCreateParams();
  // SpectraSTPepXMLLibImporter*  specLib = new SpectraSTPepXMLLibImporter(mtvec, mtlib, *mtparams);
  SpectraSTLibEntry* entry = NULL;
  PeptideUser* specPep = NULL;
  
  cRamp* cramp = NULL;

  string dataFile = "";
  string dataExt = "";

  bool has_mod = false;
  bool reported = false;
  int hit_rank = 0;
  Quantic* quant = NULL;;
  bool in_oldptm_result = false;

  bool skip_run = false;

  bool is_nterm_pep = false;
  
  bool is_cterm_pep = false;
  SQtot_ = 0;
  // TODO 

  // TODO  double allntt_prob[3] = {-100, -100, -100};
  // for(k = 0; k < input_files_->length(); k++) {
    RACI fin(c); // read possibly gzipped files
    if(! fin) {
      cerr << "fin: error opening " << c << endl;
      exit(1);
    }
    //cerr << "INFO: Reading file " << c << " ..." << endl;
    lastpos = fin.tellg();
    while(fin.getline(nextline, line_width_)) {
      data = strstr(nextline, "<");
      while(data != NULL) {
	if (tag) {
	  delete tag;
	  tag = NULL;
	}
	
	tag = new Tag(data);
	  if (tag->isStart() && 
	      ! strcmp(tag->getName(), "msms_pipeline_analysis")) {
	    //(*input_files_)[k]->assign(tag->getAttributeValue("summary_xml"));
	    //inter_proph_->addSearch((*input_files_)[k]);
	    //	    fout << "<analysis_summary analysis=\"quantic\" time=\"" << time_ << "\"/>" << endl;
	  }
	  else if (! strcmp(tag->getName(), "aminoacid_modification") && !strcmp(tag->getAttributeValue("variable"), "N")) {
	    if (stat_mods_hash_.find(*tag->getAttributeValue("aminoacid")) != stat_mods_hash_.end() 
		&& fabs(atof(tag->getAttributeValue("massdiff"))-stat_mods_hash_[*tag->getAttributeValue("aminoacid")])>0.001) {
	      cerr << "ERROR: multiple static mods with different masses found on aminoacid " << tag->getAttributeValue("aminoacid") << ", cannot be processed together with Quantic." << endl;
	      exit(1);
	    }
	    stat_mods_hash_[*tag->getAttributeValue("aminoacid")] = atof(tag->getAttributeValue("massdiff"));
	  }
	  else if (! strcmp(tag->getName(), "aminoacid_modification") && !strcmp(tag->getAttributeValue("variable"), "Y")) {
	    if (var_mods_hash_.find(*tag->getAttributeValue("aminoacid")) != var_mods_hash_.end() ) {
	
	      bool unseen = true;
	      for (int m =0; m < var_mods_hash_[*tag->getAttributeValue("aminoacid")]->size(); m++) {
		if (fabs(atof(tag->getAttributeValue("massdiff"))-(*var_mods_hash_[*tag->getAttributeValue("aminoacid")])[m])<0.00001) {
		  unseen = false;
		}
	      }
	      if(unseen) var_mods_hash_[*tag->getAttributeValue("aminoacid")]->push_back(atof(tag->getAttributeValue("massdiff")));	      

	    }
	    else if ( var_mods_hash_.find(*tag->getAttributeValue("aminoacid")) == var_mods_hash_.end() ) {
	
	      var_mods_hash_[*tag->getAttributeValue("aminoacid")] = new vector<double>();
	      var_mods_hash_[*tag->getAttributeValue("aminoacid")]->push_back(atof(tag->getAttributeValue("massdiff")));	      

	    }

	  }
	  else if (! strcmp(tag->getName(), "terminal_modification") && !strcmp(tag->getAttributeValue("variable"), "N") && !strcmp(tag->getAttributeValue("protein_terminus"), "N")) {
	    if (stat_mods_hash_.find(tolower(*tag->getAttributeValue("terminus"))) != stat_mods_hash_.end()
		&& fabs(atof(tag->getAttributeValue("massdiff"))-stat_mods_hash_[tolower(*tag->getAttributeValue("terminus"))])>0.001) {
	      cerr << "ERROR: multiple static mods with different masses found on terminus " << tag->getAttributeValue("terminus") << ", cannot be processed together with Quantic." << endl;
	      exit(1);
	    }
	    stat_mods_hash_[tolower(tolower(*tag->getAttributeValue("terminus")))] = atof(tag->getAttributeValue("massdiff"));
	  }
	  else if (! strcmp(tag->getName(), "terminal_modification") && !strcmp(tag->getAttributeValue("variable"), "Y") && !strcmp(tag->getAttributeValue("protein_terminus"), "N")) {
	    if (var_mods_hash_.find(tolower(*tag->getAttributeValue("terminus"))) != var_mods_hash_.end() ) {

	      bool unseen = true;
	      if (var_mods_hash_.find(tolower(*tag->getAttributeValue("terminus")))->second) {
		for (int m =0; m < var_mods_hash_[tolower(*tag->getAttributeValue("terminus"))]->size(); m++) {
		  if (fabs(atof(tag->getAttributeValue("massdiff"))-(*var_mods_hash_[tolower(*tag->getAttributeValue("terminus"))])[m])<0.00001) {
		    unseen = false;
		  }
		}
	      }
	      else {
		var_mods_hash_[tolower(*tag->getAttributeValue("terminus"))] = new vector<double>();
	      }
		
		
		
	      if (unseen) var_mods_hash_[tolower(*tag->getAttributeValue("terminus"))]->push_back(atof(tag->getAttributeValue("massdiff")));	      

	    }
	    else if ( var_mods_hash_.find(tolower(*tag->getAttributeValue("terminus"))) == var_mods_hash_.end() ) {
	
	      var_mods_hash_[tolower(*tag->getAttributeValue("terminus"))] = new vector<double>();
	      var_mods_hash_[tolower(*tag->getAttributeValue("terminus"))]->push_back(atof(tag->getAttributeValue("massdiff")));	      

	    }


	  }
	  else if (! strcmp(tag->getName(), "terminal_modification") && !strcmp(tag->getAttributeValue("variable"), "N") && !strcmp(tag->getAttributeValue("protein_terminus"), "Y")) {
	    if (stat_prot_termods_hash_.find(tolower(*tag->getAttributeValue("terminus"))) != stat_prot_termods_hash_.end()
		&& fabs(atof(tag->getAttributeValue("massdiff"))-stat_prot_termods_hash_[tolower(*tag->getAttributeValue("terminus"))])>0.001) {
	      cerr << "ERROR: multiple static mods with different masses found on terminus " << tag->getAttributeValue("terminus") << ", cannot be processed together with Quantic." << endl;
	      exit(1);
	    }
	    stat_prot_termods_hash_[tolower(*tag->getAttributeValue("terminus"))] = atof(tag->getAttributeValue("massdiff"));
	  }
	  else if (! strcmp(tag->getName(), "terminal_modification") && !strcmp(tag->getAttributeValue("variable"), "Y") && !strcmp(tag->getAttributeValue("protein_terminus"), "Y")) {
	    if (var_prot_termods_hash_.find(tolower(*tag->getAttributeValue("terminus"))) != var_prot_termods_hash_.end() ) {
	      
	        bool unseen = true;
	      for (int m =0; var_prot_termods_hash_[tolower(*tag->getAttributeValue("terminus"))] && m < var_prot_termods_hash_[tolower(*tag->getAttributeValue("terminus"))]->size(); m++) {
		if (fabs(atof(tag->getAttributeValue("massdiff"))-(*var_prot_termods_hash_[tolower(*tag->getAttributeValue("terminus"))])[m])<0.00001) {
		    unseen = false;
		  }
	      }
		
	      
	      if (unseen) var_prot_termods_hash_[tolower(*tag->getAttributeValue("terminus"))]->push_back(atof(tag->getAttributeValue("massdiff")));	      

	    }
	    else if ( var_prot_termods_hash_.find(tolower(*tag->getAttributeValue("terminus"))) == var_prot_termods_hash_.end() ) {
	
	      var_prot_termods_hash_[tolower(*tag->getAttributeValue("terminus"))] = new vector<double>();
	      var_prot_termods_hash_[tolower(*tag->getAttributeValue("terminus"))]->push_back(atof(tag->getAttributeValue("massdiff")));	      

	    }

	  }

	  else if (tag->isStart() && 
	      ! strcmp(tag->getName(), "msms_run_summary")) {
	    //(*input_files_)[k]->assign(tag->getAttributeValue("summary_xml"));
	    //inter_proph_->addSearch((*input_files_)[k]);
	    run++;
	    dataFile = tag->getAttributeValue("base_name");
	    dataExt = tag->getAttributeValue("raw_data");
	    if (dataExt[0] != '.') {
	      dataFile += ".";
	    }
	    dataFile += dataExt;
#ifdef MSVC
	      WaitForSingleObject( _mutex,    // handle to mutex
				   INFINITE);      //#include "windows.h"
#else	
	      pthread_mutex_lock( &_mutex );//#include <pthread.h>
#endif   
	      //if (cramp != NULL)
	      //delete cramp;
	      //else
	      //	cramps_.clear();

	      cramp = new cRamp(dataFile.c_str());
	    if (!cramp->OK()) {
	      cerr << "ERROR: cannot read scan in data file " << dataFile << " exiting ..." << endl;                
	      exit(1);
	    }
#ifdef MSVC
		ReleaseMutex( _mutex);      //#include "windows.h"
#else
		pthread_mutex_unlock( &_mutex );//#include <pthread.h>
#endif

	    //stat_mods_hash_.clear();	    
	    //stat_prot_termods_hash_.clear();
	  }
	  
	  else if (tag->isStart() && 
	    ! strcmp(tag->getName(), "spectrum_query")) {
	    SQtot_++;
	    spec_prob = 0.;
	    if (tag->getAttributeValue("experiment_label") != NULL) {
	      exp_lbl = tag->getAttributeValue("experiment_label");
	    }
	    pos_mod_hash.clear();
	    spectrum_name = tag->getAttributeValue("spectrum");
	    SQoffsets_.push_back(lastpos);
	    SQruns_.push_back(run);
	    cramps_.push_back(cramp);
	    charge = tag->getAttributeValue("assumed_charge");
	    scan = atoi(tag->getAttributeValue("start_scan"));
	  }
	  else if ( tag->isStart() && 
		    ! strcmp(tag->getName(), "search_hit") && 
		    atoi(tag->getAttributeValue("hit_rank")) == 1) {
	    
	    hit_rank = 1;
	    reported = false;
	    pep_seq = "";
	    pep_seq += tag->getAttributeValue("peptide");

	    mod_pep = pep_seq;

	    calcnmass = atof(tag->getAttributeValue("calc_neutral_pep_mass"));
	    
	    //	    if (unknown_mode_) massshift_[0] = atof(tag->getAttributeValue("massdiff"));
	    
	    is_nterm_pep = false;
	    is_cterm_pep = false;
	    
	    if (! strcmp(tag->getAttributeValue("peptide_prev_aa") , "-")) {
	      is_nterm_pep = true;
	    }
	    if (! strcmp(tag->getAttributeValue("peptide_next_aa") , "-")) {
	      is_cterm_pep = true;
	    }
	  }
	  else if (tag->isStart() && 
		   ! strcmp(tag->getName(), "search_hit") && 
		   atoi(tag->getAttributeValue("hit_rank")) != 1) {
	    hit_rank = atoi(tag->getAttributeValue("hit_rank"));
	  }
	  else if (  hit_rank == 1 && ( (tag->isEnd() && //get_pep &&
					 ! strcmp(tag->getName(), "modification_info") ) ) ) {
					//|| 
					//( tag->isEnd() && //get_pep &&
					// ! strcmp(tag->getName(), "search_hit") && 
					//  unknown_mode_ && !has_mod)) ) {      
	    
	    //DDS: TODO IT HERE!!!!
	    // specLib->SpectraSTPepXMLLibImporter::addModificationsToPeptide(pep_seq, mod_tags);
	    if (quant) {
	      delete quant;
	      quant = NULL;
	    }
	    if (tag->isStart()) {
	      mod_pep = "";
	      if (tag->getAttributeValue("modified_peptide")) 
		mod_pep += tag->getAttributeValue("modified_peptide");
	    }
	    if (mod_pep == "") {
	      mod_pep = pep_seq;
	    }
	    in_mod = false;
	  }
	  else if ( hit_rank == 1 && tag->isStart() && !tag->isEnd() && //get_pep &&
		     ! strcmp(tag->getName(), "modification_info")) {      
	    in_mod = true;
	    has_mod = true;
	    //       	    mod_tags = data; 
	    //	    get_pep = false;
	    mod_pep = "";
	    if (tag->getAttributeValue("modified_peptide"))
	      mod_pep += tag->getAttributeValue("modified_peptide");

	    double mass = 0;
	    if (tag->getAttributeValue("mod_nterm_mass")) {
	      mass = atof(tag->getAttributeValue("mod_nterm_mass"));
	      pos_mod_hash.insert(make_pair(0, mass));
	    }
	    if (tag->getAttributeValue("mod_cterm_mass")) {
	      mass = atof(tag->getAttributeValue("mod_cterm_mass"));
	      pos_mod_hash.insert(make_pair(pep_seq.length()+1, mass));
	    }
	    
	  }
	  else if ( hit_rank == 1 && tag->isStart() &&
		    ! strcmp(tag->getName(), "mod_aminoacid_mass")) {      
	    
	    int pos = atoi(tag->getAttributeValue("position"));
	    double mass = atof(tag->getAttributeValue("mass"));;
	    pos_mod_hash.insert(make_pair(pos, mass));
	    
	  }

	  else if ( ( tag->isEnd() &&  ! strcmp(tag->getName(), "spectrum_query") ) ||
		    ( tag->isEnd() &&  ! strcmp(tag->getName(), "search_hit") && hit_rank == 1) ) {
	    
	    if (spec_prob  >= 0.9) {
	      if (quant) {
		delete quant;
		quant = NULL;
	      }
	      
	      quant = new Quantic(mod_pep, atoi(charge.c_str()), calcnmass, cramp, scan,
				   aminoacids_, massshift_, neutlosses_, 
					  &stat_mods_hash_, &var_mods_hash_, &stat_prot_termods_hash_, 
				  &var_prot_termods_hash_,&_mutex,is_nterm_pep, is_cterm_pep, pos_mod_hash);// unknown_mode_, labile_mode_, direct_mode_);
	      quant->setMolForms(ptmMolForms_);
	      if (mzTol_ < 0 && quant->init())
		quant->computeMassDiffs(&massdiffs_);
	      
	      delete quant;
	      quant = NULL;
	      
	      hit_rank = 0;
	      //	    quant->evaluateModSites(mzTol_);
	      
	      //mod_tags->clear();
	    }
	    
	      spec_prob_hash_.insert(make_pair(exp_lbl+spectrum_name, spec_prob));
	      
	      has_mod = false;
	      exp_lbl = "";
	      hit_rank = 0;
	      spectrum_name = "";
	      spec_prob = 0.;
	}
	if (tag->isStart() && ! strcmp(tag->getName(), "peptideprophet_result")) { 
	  
	    spec_prob = atof( tag->getAttributeValue("probability") );
	    
	  
	
	}
	if ( tag->isStart() && ! strcmp(tag->getName(), "interprophet_result")) { 

	  
	    spec_prob = atof( tag->getAttributeValue("probability") );
	    

	}
	data = strstr(data+1, "<");      
      }
      delete tag;
      tag = NULL;
      lastpos = fin.tellg();
    }
    fin.close();
    delete [] nextline;
    
    
    if (mzTol_ < 0) {
      mzTol_ = 0.1;
      double currentTol = mzTol_;
      double lastTol = 0;
      double mean = 0;
      double stddev = 0;
      gsl_vector* r;
      size_t size = 0;
      for (int k = 0; k < massdiffs_.size(); k++) {
	if (fabs(massdiffs_[k]) < currentTol)
	  size++;
      }
      while (size && fabs(currentTol-lastTol) > 0.00001) {
	r = gsl_vector_calloc(size);
	int s=0;
	for (int k = 0; k < massdiffs_.size(); k++) {
	  if (fabs(massdiffs_[k]) < currentTol)
	    gsl_vector_set(r, s++, massdiffs_[k]);
	}
	mean = gsl_stats_mean(gsl_vector_ptr(r, 0), 1, size);
	stddev = size > 1 ? gsl_stats_sd(gsl_vector_ptr(r, 0), 1, size) : 0;
	
	lastTol = currentTol;
	currentTol = 3*stddev;
	
	gsl_vector_free(r);
	size = 0;
	for (int k = 0; k < massdiffs_.size(); k++) {
	  if (fabs(massdiffs_[k]) < currentTol) 
	    size++;
	}
      }
      
      cout << "INFO: mean massdiff: " << mean << endl;
      cout << "INFO: stdev massdiff: " << stddev << endl;
      
      if (stddev > 0) {
	cout << "INFO: bootstrapped mzTOL 3*stdev: " << 3*stddev << endl;
	mzTol_ = 3*stddev;
      }
      else {
	cout << "INFO: cannot bootstrap  mzTOL, using: " << mzTol_ << endl;
      }

    }
    

    //if (size)
    //  gsl_vector_free(r);
}

void QuanticParser::parseWriteUpdate(const char* c, bool update, int t) {

  Array<Tag*>* tags = NULL;
  Tag* tag = NULL;

  char *nextline = new char[line_width_];
  char* data = NULL;

  update = false;


  CMercury8* mercury = new CMercury8(NULL);

  string spectrum_name = "";
  string exp_lbl = "";
  string charge = "";
  double spec_prob = 0;
  Array<double>* allntt_prob=NULL ;
  double calcnmass = -1;
  double rt = -1;
  string pep_seq = "";
  string mod_pep = "";
  TPP_HASHMAP_T<int, double> pos_mod_hash;
  // string mod_tags = "";
  Array<Tag*>* mod_tags = new Array<Tag*>();//NULL;

  bool in_mod= false;

  bool has_mod= false;
  
  bool get_pep = false;
  int k=0;

  long scan = -1;
  vector<string> mtvec;
  SpectraSTLib* mtlib= NULL;
  //SpectraSTCreateParams* mtparams=new SpectraSTCreateParams();
  // SpectraSTPepXMLLibImporter*  specLib = new SpectraSTPepXMLLibImporter(mtvec, mtlib, *mtparams);
  SpectraSTLibEntry* entry = NULL;
  PeptideUser* specPep = NULL;
  
  cRamp* cramp = NULL;

  string dataFile = "";
  string dataExt = "";

  bool in_oldptm_result = false;

  bool skip_run = false;
  
  bool reported = false;



  string SQstring = "";

  unsigned long SQtot = SQtot_;

  string QuanticReport = "";
  Quantic* quant = NULL; 

  bool is_nterm_pep = false;
  
  bool is_cterm_pep = false;
  unsigned long SQcount = 0;
  unsigned int hit_rank = 0;
  unsigned int run = 0;

  double massdiff = 0;
  // TODO 
  
  // TODO  double allntt_prob[3] = {-100, -100, -100};
  // for(k = 0; k < input_files_->length(); k++) {
  RACI fin(c); // read possibly gzipped files
  if(! fin) {
    cerr << "fin: error opening " << c << endl;
    exit(1);
  }
  cerr << "INFO: Reading file " << c << " ..." << endl;
  
  unsigned int thread_run_idx = 0;
  bool writeHeaders = false;
  
  while(fin.getline(nextline, line_width_)) {
      data = strstr(nextline, "<");
      
      while(data != NULL) {
	
	if (!in_mod && tag) {
	  delete tag;
	}
	
	tag = new Tag(data);
	
	if (!keepOld_) {
	  //leave out old Qresults
	  if (in_oldptm_result && !strcmp(tag->getName(), "analysis_result") && tag->isEnd()) {
	    in_oldptm_result = false;
	    data = NULL;
	    continue;
	  } 
	  if (!strcmp(tag->getName(), "analysis_result") && tag->isStart() && !strcmp(tag->getAttributeValue("analysis"), "quantic")) {
	    in_oldptm_result = true;
	  }
	  
	  if (in_oldptm_result) {
	    data = NULL;
	    continue;
	  }
	}
	if (!update ||
	    (!in_mod && strcmp(tag->getName(), "modification_info") )) { 
	  /////&& !unknown_mode_) ||
	  //  (!in_mod && !(tag->isEnd() && ! strcmp(tag->getName(), "search_hit") && unknown_mode_ && !has_mod) ) ) {   
	  if (!update || (!in_mod && strcmp(tag->getName(), "modification_info"))) {
	    if (hit_rank <= 1 && !strcmp(tag->getName(), "search_hit") && tag->isEnd() && !QuanticReport.empty() ) {
	      //hit_rank = 0;
	      
	      
	      // fout << QuanticReport;
	      SQstring += QuanticReport;
	      
	    }
	    if ( !( !strcmp(tag->getName(), "spectrum_query") && tag->isStart() ) ) {
	      if (!strcmp(tag->getName(), "msms_run_summary") && tag->isEnd()) {
		writeHeaders = false;
	      }
	    
	      if (writeHeaders) {
		LOCK(&fout_mutex);
		tag->write(fout);
		UNLOCK(&fout_mutex);

	      }
	      else if (!strcmp(tag->getName(), "xml")  ) {
		LOCK(&run_mutex);
		LOCK(&hdr_mutex);
		if (header_thread_ == max_threads_ ) {

		  header_thread_ = t;
		  run=1;
		  threads_start_run_=1;
		  threads_done_run_=0;
		  writeHeaders = true;
		  LOCK(&fout_mutex);
		  tag->write(fout);
		  UNLOCK(&fout_mutex);
		}
		else {
		  UNLOCK(&hdr_mutex);
		}
		UNLOCK(&run_mutex);

	// run++;
	// 	LOCK(&run_mutex);
	// 	threads_start_run_++;
	// 	LOCK(&fout_mutex);
	// 	tag->write(fout);
	// 	UNLOCK(&fout_mutex);
	// 	threads_done_run_=0;
	// 	writeHeaders = true;
	      }
	      else if (!strcmp(tag->getName(), "msms_run_summary") ) {
		if (tag->isEnd()) {

  
		    writeHeaders = false;

		    while (threads_start_run_ != max_threads_) {	      //DDS: DEBUG this infinite loop
		      YIELD();
		    }
		    QuanticReport = "";


		  LOCK(&run_mutex);
		  
		   //cerr << "INFO: thread " << t << " finished with run " << run << ".  Current Run is " << current_run_ << endl;
		
		   threads_done_run_++;
		  
		    
		   if (threads_done_run_ == max_threads_) {
		     LOCK(&fout_mutex);
		     tag->write(fout);//last one done closes tag // unlocks 
		     UNLOCK(&fout_mutex);
		     threads_start_run_=0;
		     current_run_++;	
		     UNLOCK(&run_mutex);
		   }
		   else {
		     UNLOCK(&run_mutex);
		     //Need to wait
		     while (threads_done_run_ != max_threads_ && run == current_run_) {
		       YIELD();
		     }
		   }
		
		}
		else	if (tag->isStart() ) {
		  LOCK(&run_mutex);
		  //		  writeHeaders = false;
		  run++;
		   
		  threads_start_run_++;
		
		  LOCK(&hdr_mutex);
		  if (threads_start_run_ == 1) {
		    LOCK(&fout_mutex);
		    tag->write(fout);//first one to start writes run headers  //keeps the Lock  writing until the header
		    UNLOCK(&fout_mutex);
		    threads_done_run_=0;
		    writeHeaders= true;
		    header_thread_ = t;
		  }
		  else {
		    UNLOCK(&hdr_mutex);

		  }
		  UNLOCK(&run_mutex);		  
		  
		}		  
 

	      }
	      else if (!SQstring.empty() && (SQcount-1) % max_threads_ == t)   {
		ostringstream tagout;
		if (  (!strcmp(tag->getName(), "search_result") && tag->isEnd()) ||  (!strcmp(tag->getName(), "spectrum_query") && tag->isEnd()) || 
		      ( hit_rank <= 1 && (strcmp(tag->getName(), "search_hit") || 
					  (!strcmp(tag->getName(), "search_hit") && tag->isStart() && atoi(tag->getAttributeValue("hit_rank"))==1) ||
					  (!strcmp(tag->getName(), "search_hit") && tag->isEnd())) ) ) {
		  
		  
		  tag->write(tagout);
		  SQstring += tagout.str();
		  
		}
	      }	
		
	    }
	  }
	}	  

	
      
	if (in_mod &&  (SQcount-1) % max_threads_ == t ) {
	  if (!mod_tags->size() || tag != (*mod_tags)[mod_tags->size()-1])
	    mod_tags->insertAtEnd(tag);
	}
      
	if (tag != NULL) {
	  if (tag->isStart() && !strcmp(tag->getName(), "search_summary") && 
	      ( !strcmp(tag->getAttributeValue("search_engine"), "SpectraST") ) ) { //&& !unknown_mode) ) {
	    skip_run = true;
	    
	  }
	  else if (tag->isStart() && !strcmp(tag->getName(), "search_summary")) {
	    skip_run = false;
	    
	  }
	  
	  if ( tag->isStart() && ! strcmp(tag->getName(), "peptideprophet_result")) { 
	    
	    spec_prob = atof( tag->getAttributeValue("probability") );
	    
	    
	  }
	  if ( tag->isStart() && ! strcmp(tag->getName(), "interprophet_result")) { 
	    
	    spec_prob = atof( tag->getAttributeValue("probability") );
	    
	  }
	  
	  if (t==header_thread_ && !skip_run && tag->isEnd() && !strcmp(tag->getName(), "search_summary")) {
	    Tag* ts_tag = new Tag("analysis_timestamp", True, True);
	    ts_tag->setAttributeValue("analysis", "quantic");
	    ts_tag->setAttributeValue("time", time_);
	    ts_tag->setAttributeValue("id", "1");
	    ostringstream tagout;
	    ts_tag->write(tagout);
	    SQstring += tagout.str();
	    delete ts_tag;
	    
	  }
	  
	  if (!skip_run && tag->isStart() && 
	      ! strcmp(tag->getName(), "msms_pipeline_analysis")) {
	    //(*input_files_)[k]->assign(tag->getAttributeValue("summary_xml"));
	    //inter_proph_->addSearch((*input_files_)[k]);
	    if (t==header_thread_) {
		fout << "<analysis_summary analysis=\"quantic\" time=\"" << time_   << "\">" << endl
		     << "<quantic_summary version=\"" << szTPPVersionInfo  << "\""
		     << " options=\"" << opts_.c_str() 
		     << "\">" << endl;
		
		fout << "<inputfile name=\"" << c << "\"/>" << endl;
	    
				
		fout << "</quantic_summary>" << endl;
	    
	    
		fout << "</analysis_summary>" << endl;
	      }
	  }
	  else if (! strcmp(tag->getName(), "aminoacid_modification") && !strcmp(tag->getAttributeValue("variable"), "N")) {
	    if (stat_mods_hash_.find(*tag->getAttributeValue("aminoacid")) != stat_mods_hash_.end() 
		&& fabs(atof(tag->getAttributeValue("massdiff"))-stat_mods_hash_[*tag->getAttributeValue("aminoacid")])>0.001) {
	      cerr << "ERROR: multiple static mods with different masses found on aminoacid " << tag->getAttributeValue("aminoacid") << ", cannot be processed together with Quantic." << endl;
	      exit(1);
	    }
	    stat_mods_hash_[*tag->getAttributeValue("aminoacid")] = atof(tag->getAttributeValue("massdiff"));
	  }
	  else if (! strcmp(tag->getName(), "aminoacid_modification") && !strcmp(tag->getAttributeValue("variable"), "Y")) {
	    if (var_mods_hash_.find(*tag->getAttributeValue("aminoacid")) != var_mods_hash_.end() ) {
	
	      bool unseen = true;
	      for (int m =0; m < var_mods_hash_[*tag->getAttributeValue("aminoacid")]->size(); m++) {
		if (fabs(atof(tag->getAttributeValue("massdiff"))-(*var_mods_hash_[*tag->getAttributeValue("aminoacid")])[m])<0.00001) {
		  unseen = false;
		}
	      }
	      if(unseen) var_mods_hash_[*tag->getAttributeValue("aminoacid")]->push_back(atof(tag->getAttributeValue("massdiff")));	      

	    }
	    else if ( var_mods_hash_.find(*tag->getAttributeValue("aminoacid")) == var_mods_hash_.end() ) {
	
	      var_mods_hash_[*tag->getAttributeValue("aminoacid")] = new vector<double>();
	      var_mods_hash_[*tag->getAttributeValue("aminoacid")]->push_back(atof(tag->getAttributeValue("massdiff")));	      

	    }

	  }
	  else if (! strcmp(tag->getName(), "terminal_modification") && !strcmp(tag->getAttributeValue("variable"), "N") && !strcmp(tag->getAttributeValue("protein_terminus"), "N")) {
	    if (stat_mods_hash_.find(tolower(*tag->getAttributeValue("terminus"))) != stat_mods_hash_.end()
		&& fabs(atof(tag->getAttributeValue("massdiff"))-stat_mods_hash_[tolower(*tag->getAttributeValue("terminus"))])>0.001) {
	      cerr << "ERROR: multiple static mods with different masses found on terminus " << tag->getAttributeValue("terminus") << ", cannot be processed together with Quantic." << endl;
	      exit(1);
	    }
	    stat_mods_hash_[tolower(tolower(*tag->getAttributeValue("terminus")))] = atof(tag->getAttributeValue("massdiff"));
	  }
	  else if (! strcmp(tag->getName(), "terminal_modification") && !strcmp(tag->getAttributeValue("variable"), "Y") && !strcmp(tag->getAttributeValue("protein_terminus"), "N")) {
	    if (var_mods_hash_.find(tolower(*tag->getAttributeValue("terminus"))) != var_mods_hash_.end() ) {
	      
	      bool unseen = true;
	      if (var_mods_hash_.find(tolower(*tag->getAttributeValue("terminus")))->second) {
		for (int m =0; m < var_mods_hash_[tolower(*tag->getAttributeValue("terminus"))]->size(); m++) {
		  if (fabs(atof(tag->getAttributeValue("massdiff"))-(*var_mods_hash_[tolower(*tag->getAttributeValue("terminus"))])[m])<0.00001) {
		    unseen = false;
		  }
		}
	      }
	      else {
		var_mods_hash_[tolower(*tag->getAttributeValue("terminus"))] = new vector<double>();
	      }
	      
	      
	      
	      if (unseen) var_mods_hash_[tolower(*tag->getAttributeValue("terminus"))]->push_back(atof(tag->getAttributeValue("massdiff")));	      
	      
	    }
	    else if ( var_mods_hash_.find(tolower(*tag->getAttributeValue("terminus"))) == var_mods_hash_.end() ) {
	      
	      var_mods_hash_[tolower(*tag->getAttributeValue("terminus"))] = new vector<double>();
	      var_mods_hash_[tolower(*tag->getAttributeValue("terminus"))]->push_back(atof(tag->getAttributeValue("massdiff")));	      

	    }
	    
	    
	  }
	  else if (! strcmp(tag->getName(), "terminal_modification") && !strcmp(tag->getAttributeValue("variable"), "N") && !strcmp(tag->getAttributeValue("protein_terminus"), "Y")) {
	    if (stat_prot_termods_hash_.find(tolower(*tag->getAttributeValue("terminus"))) != stat_prot_termods_hash_.end()
		&& fabs(atof(tag->getAttributeValue("massdiff"))-stat_prot_termods_hash_[tolower(*tag->getAttributeValue("terminus"))])>0.001) {
	      cerr << "ERROR: multiple static mods with different masses found on terminus " << tag->getAttributeValue("terminus") << ", cannot be processed together with Quantic." << endl;
	      exit(1);
	    }
	    stat_prot_termods_hash_[tolower(*tag->getAttributeValue("terminus"))] = atof(tag->getAttributeValue("massdiff"));
	  }
	  else if (! strcmp(tag->getName(), "terminal_modification") && !strcmp(tag->getAttributeValue("variable"), "Y") && !strcmp(tag->getAttributeValue("protein_terminus"), "Y")) {
	    if (var_prot_termods_hash_.find(tolower(*tag->getAttributeValue("terminus"))) != var_prot_termods_hash_.end() ) {
	      
	      bool unseen = true;
	      for (int m =0; var_prot_termods_hash_[tolower(*tag->getAttributeValue("terminus"))] && m < var_prot_termods_hash_[tolower(*tag->getAttributeValue("terminus"))]->size(); m++) {
		if (fabs(atof(tag->getAttributeValue("massdiff"))-(*var_prot_termods_hash_[tolower(*tag->getAttributeValue("terminus"))])[m])<0.00001) {
		  unseen = false;
		}
	      }
	      
	      
	      if (unseen) var_prot_termods_hash_[tolower(*tag->getAttributeValue("terminus"))]->push_back(atof(tag->getAttributeValue("massdiff")));	      
	      
	    }
	    else if ( var_prot_termods_hash_.find(tolower(*tag->getAttributeValue("terminus"))) == var_prot_termods_hash_.end() ) {
	      
	      var_prot_termods_hash_[tolower(*tag->getAttributeValue("terminus"))] = new vector<double>();
	      var_prot_termods_hash_[tolower(*tag->getAttributeValue("terminus"))]->push_back(atof(tag->getAttributeValue("massdiff")));	      
	      
	    }
	    
	  }
	  else if (tag->isStart() && 
	      ! strcmp(tag->getName(), "msms_run_summary")) {
	    //(*input_files_)[k]->assign(tag->getAttributeValue("summary_xml"));
	    //inter_proph_->addSearch((*input_files_)[k]);
	    dataFile = tag->getAttributeValue("base_name");
	    dataExt = tag->getAttributeValue("raw_data");
	    if (dataExt[0] != '.') {
	      dataFile += ".";
	    }
	    dataFile += dataExt;
	    
	    LOCK(&_mutex);

	    if (cramp != NULL)
	      delete cramp;
	    else 
	      cramps_.clear();
	    
	    cramp = new cRamp(dataFile.c_str());
	    if (!cramp->OK()) {
	     cerr << "ERROR: cannot read scan in data file " << dataFile << " exiting ..." << endl;
	     exit(1);
	     }

	    UNLOCK(&_mutex);

	    //stat_mods_hash_.clear();
	    //stat_prot_termods_hash_.clear();
	  }
	  else if (tag->isStart() && 
		   ! strcmp(tag->getName(), "spectrum_query")) {

	    if (header_thread_ != t) {
	      while (!TRYLOCK(&hdr_mutex));
	    }
	    UNLOCK(&hdr_mutex);
	    
	    if (writeHeaders) {
	      UNLOCK(&run_mutex);
	      writeHeaders = false;
	    

	      while (threads_start_run_ != max_threads_) {
		YIELD();
	      }
	    }

	    
	    SQcount++;
	    
	    
	    //	    get_pep = true;

	    QuanticReport = "";
	    SQstring = "";
	    if ((SQcount-1)%max_threads_ == t) {
	      ostringstream tagout;
	      tag->write(tagout);
	      SQstring += tagout.str();
	    }
	    spec_prob = 0.;
	    pep_seq = "";
	    mod_pep = "";
	    spectrum_name = "";
	    pos_mod_hash.clear();
	    rt = -1;
	    if (tag->getAttributeValue("retention_time_sec") != NULL) 
	      rt = atof(tag->getAttributeValue("retention_time_sec"));
	    
	    // store the spectrum name for later
	    int len = strlen(tag->getAttributeValue("spectrum"));
	    
	    spectrum_name += tag->getAttributeValue("spectrum");
	    
	    if (tag->getAttributeValue("experiment_label") != NULL) {
	      len = strlen(tag->getAttributeValue("experiment_label"));
	    }
	    else {
	      len = 0;
	    }
	    exp_lbl = "";
	    if (len > 0)
	      exp_lbl += tag->getAttributeValue("experiment_label");
	    //char* tmp = strrchr(spectrum_name, '.');
	    //tmp = '\0';
	    charge = tag->getAttributeValue("assumed_charge");
	    scan = atoi(tag->getAttributeValue("start_scan"));
	    
	    
	    
	  }
	  else if ((SQcount-1) % max_threads_ == t && tag->isStart() && 
		   ! strcmp(tag->getName(), "search_hit") && 
		   atoi(tag->getAttributeValue("hit_rank")) == 1) {

	    hit_rank = 1;
	    reported = false;
	    pep_seq = "";
	    pep_seq += tag->getAttributeValue("peptide");
	    mod_pep = pep_seq;
	    calcnmass = atof(tag->getAttributeValue("calc_neutral_pep_mass"));
	    massdiff =  atof(tag->getAttributeValue("massdiff"));
	    //if (unknown_mode_) massshift_[0] = atof(tag->getAttributeValue("massdiff"));
	    
	    is_nterm_pep = false;
	    is_cterm_pep = false;
	    
	    if (! strcmp(tag->getAttributeValue("peptide_prev_aa") , "-")) {
	      is_nterm_pep = true;
	    }
	    if (! strcmp(tag->getAttributeValue("peptide_next_aa") , "-")) {
	      is_cterm_pep = true;
	    }
	  }
	  else if ((SQcount-1) % max_threads_ == t && tag->isStart() && 
		   ! strcmp(tag->getName(), "search_hit") && 
		   atoi(tag->getAttributeValue("hit_rank")) != 1) {
	    hit_rank = atoi(tag->getAttributeValue("hit_rank"));
	  }
	  else if ( (SQcount-1) % max_threads_ == t &&
		    hit_rank == 1 && ( (tag->isEnd() && //get_pep &&
					! strcmp(tag->getName(), "modification_info") ) || QuanticReport.empty() && ! strcmp(tag->getName(), "search_score")  ) ) {

	    //||      ( tag->isEnd() &&  ! strcmp(tag->getName(), "search_hit") && unknown_mode_ && !has_mod)) ) {   
	    if (tag->isStart()) {
	      mod_pep = "";	
	      if (tag->getAttributeValue("modified_peptide"))
		mod_pep += tag->getAttributeValue("modified_peptide");
	    }
	    if (mod_pep == "") {
	      mod_pep = pep_seq;
	    }
	    if (spec_prob_hash_.find(exp_lbl+spectrum_name) != spec_prob_hash_.end() && spec_prob_hash_.find(exp_lbl+spectrum_name)->second  >= minProb_) {
	      
	      if ((SQcount-1)%max_threads_==t && SQcount%1000 == 0)
		cerr << "\rINFO: processed " << SQcount << "/" << SQtot << " spectrum_queries\t\t ";
	      
	      //}
	      //else if ( ( tag->isEnd() && //get_pe
	      //	      ! strcmp(tag->getName(), "analysis_result") && unknown_mode_ && !has_mod && !reported) ) { 
	      
	      reported = true;

	      if (quant) {
		delete quant;
		quant = NULL;
		
	      }
	      
	      quant = new Quantic(mod_pep, atoi(charge.c_str()),   calcnmass, cramp,
					   scan, 
				  aminoacids_, massshift_, neutlosses_, 
					   &stat_mods_hash_, &var_mods_hash_, 
				  &stat_prot_termods_hash_, &var_prot_termods_hash_, &_mutex, is_nterm_pep, is_cterm_pep, pos_mod_hash); 
	      // unknown_mode_, labile_mode_, direct_mode_);

	      quant->setMolForms(ptmMolForms_);
	      if (diaMode_) {
		quant->setDiaMode(true);	   
		//string molecform = quant->pepToMolecForm(mod_pep);
		//mercury->GoMercury("C44H68N12O14", 2);
		//mercury->GoMercury((char*)molecform.c_str(), atoi(charge.c_str()));
	      }
	      
	      if (quant->init()) {
		quant->processPeakList();
		quant->setVerbose(verbose_);
		//if (mzTol_ > massdiff) {
		quant->setMZTolerance(mzTol_);
		quant->setTopPeaks(top_peaks_);
		quant->setPeptideCoupled(pep_coupled_);
		//}
		//else {
		//quant->setMZTolerance(massdiff);
		//}
		quant->evaluatePeptideAnnotatedTIC(NULL);

		ostringstream modout;
		
		bool anal_report = false;
		
		if (hit_rank <= 1) {
		  ostringstream rep;
		  
		  rep << "<analysis_result analysis=\"quantic\">\n";		  
		  
		  double antic = quant->getAnnotatedTIC();
		  
		  rep << "<quantic_result antic=\"" << antic << "\"" ;
		  rep << " tic=\"" << quant->getTIC() << "\"" ;
		  if (annotate_)
		    rep << " string=\"" << quant->getAnnotatedString().c_str() << "\"" ;
		  if (!modstring_.empty()) {
		    rep << " quant=\"";
		    //if (diaMode_) {
		    //  quant->correctIsotopeErrors(0, mercury);
		    //}
		    for (int q=0; q<neutlosses_[0]->size(); q++) {
		      rep << std::fixed << std::setprecision(1);
		      rep << quant->getQuant(0, q);
		      if (q<neutlosses_[0]->size()-1) {
			rep << ",";
		      }
		    }
		    rep  << "\"";

		    rep << " norm_quant_mean=\"";
		    for (int q=0; q<neutlosses_[0]->size(); q++) {
		      rep << std::fixed << std::setprecision(3);
		      rep << quant->getQuantMean(0, q);
		      if (q<neutlosses_[0]->size()-1) {
			rep << ",";
		      }
		    }
		    rep  << "\"";

		    rep << " norm_quant_stdev=\"";
		    for (int q=0; q<neutlosses_[0]->size(); q++) {
		      rep << std::fixed << std::setprecision(3);
		      rep << quant->getQuantStdev(0, q);
		      if (q<neutlosses_[0]->size()-1) {
			rep << ",";
		      }
		    }
		    rep  << "\"";
		  }

		  rep << "/>\n";		
		  
		  
		  rep << "</analysis_result>\n";
		  
		  QuanticReport = rep.str();
		}
		
		
		delete quant;
		quant = NULL;
	
		
	      }
	      else if ( (SQcount-1) % max_threads_ == t ){  

		if (mod_pep != "") {
		  cerr << "ERROR: Cannot initialize for sequence: " << mod_pep.c_str() 
		       << ", unknown mods may exist in spectrum " << spectrum_name.c_str() << endl;
		  cerr << "Please specify modification and rerun Quantic ..." << endl;
		  exit(1);
		}
     		delete quant;
		quant = NULL;
		
	
		
	      }
	      
	    }

	    
	    //mod_tags->insertAtEnd(tag);
	    for (int m = 0; m < mod_tags->size(); m++) {
	      //ostringstream modout;
	      //(*mod_tags)[m]->write(modout);
	      //SQstring += modout.str();
	      if (m != mod_tags->size()-1) {
		delete (*mod_tags)[m];
	      }
	    }
	    
	    //tag->write(fout);
	    in_mod = false;
	    mod_tags->clear();
	      
	
	    
	  }
	  else if ( (SQcount-1) % max_threads_ == t  && //get_pep &&
		    ! strcmp(tag->getName(), "modification_info")) {      
	    if (tag->isStart()) {
	      in_mod = true;
	      has_mod = true;
	      
	      if (!mod_tags->size() || tag != (*mod_tags)[mod_tags->size()-1])
		mod_tags->insertAtEnd(tag); 
	      
	      //	    get_pep = false;
	      mod_pep = "";
	      if (tag->getAttributeValue("modified_peptide"))
		mod_pep += tag->getAttributeValue("modified_peptide");
	    
	      double mass = 0;
	      if (tag->getAttributeValue("mod_nterm_mass")) {
		mass = atof(tag->getAttributeValue("mod_nterm_mass"));
		pos_mod_hash.insert(make_pair(0, mass));
	      }
	      if (tag->getAttributeValue("mod_cterm_mass")) {
		mass = atof(tag->getAttributeValue("mod_cterm_mass"));
		pos_mod_hash.insert(make_pair(pep_seq.length()+1, mass));
	      }


	    }
	    
	    if (tag->isEnd()) {
	      in_mod = false;
	      has_mod = false;
	      mod_pep = "";
	    }
	    
	  }

	 else if ( (SQcount-1) % max_threads_ == t  &&  hit_rank == 1 && tag->isStart() &&//get_pep &&
		   ! strcmp(tag->getName(), "mod_aminoacid_mass")) {      
	   int pos = atoi(tag->getAttributeValue("position"));
	   double mass = atof(tag->getAttributeValue("mass"));;
	   pos_mod_hash.insert(make_pair(pos, mass));
	   
	 }     


	  else if (tag->isEnd() && 
		   ! strcmp(tag->getName(), "spectrum_query"))  {
	    //if (t == 0 && SQcount < 2) {

	    // UNLOCK(&fout_mutex);

	    //}

	    if (!SQstring.empty() && SQcount && (SQcount-1)%max_threads_ == t ) {

	      LOCK(&fout_mutex);
	      if (SQcount%1000 == 0)
		cerr << "\rINFO: processed " << SQcount << "/" << SQtot << " spectrum_queries\t\t ";
	      
	      fout << SQstring; 

	      UNLOCK(&fout_mutex);
	      


	    }

	    has_mod = false;
	    exp_lbl = "";
	    spectrum_name = "";
	    charge = "";
	    //	    get_pep = false;
	    pep_seq = "";
	    mod_pep = "";
	    spec_prob = 0;
	    scan = -1;
	    // TODO  allntt_prob[0] = -100; allntt_prob[1] = -100; allntt_prob[2] = -100;
	    hit_rank = 0;
	  }
	}
	
	
	if (!in_mod) {
	  delete tag;
	  tag = NULL;
	}
	
	
	data = strstr(data+1, "<");      
      }
    }

    
    
    fin.close();
    //    fout.close();

    
    mod_tags->clear();
    delete cramp;
    delete mod_tags;
    delete [] nextline;
    //    inter_proph_->computeModels();
}





bool compare_probs (ProbPos* i,ProbPos* j) { 
  return (i->prob_>j->prob_); 
}

