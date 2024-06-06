/*
Program       : PeptideMapper
Author        : Andrew Keller <akeller@systemsbiology.org>, Robert Hubley, and 
open source code
Date          : 11.02.2015
SVN Info      : $Id$

PeptideMapper implementation

Copyright (C) 2003 Andrew Keller
Copyright (C) 2015 David Shteynberg

*/

#ifndef __LGPL__
/*
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#else
/*
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
#endif

/*
Andrew Keller
Institute for Systems Biology
401 Terry Avenue North
Seattle, WA  98109  USA
akeller@systemsbiology.org
*/

#include <time.h>
#include <fstream>
#include "PeptideMapper.h"
#include "Parsers/Parser/TagListComparator.h" // regression test stuff - bpratt Insilicos LLC, Nov 2005
#include "Common/util.h"
#include <errno.h>



#ifdef __LGPL__
// SPARE parts string match stuff
#include "ac/acs.hpp"
#include "cw/cws.hpp"
#include "kmp/pm-kmp.hpp"
#include "bm/bms.hpp"
// dynamically selecting the smaller peptide alphabet when possible
#define STRINGFINDER_std PMACFail_std  
#define STRINGFINDER_pep PMACFail_peptides  
//#define STRINGFINDER PMACOpt_std  
//#define STRINGFINDER PMACKMPFail_std
//#define STRINGFINDER PMCWNaive_std
//#define STRINGFINDER PMCWNLA_std
static kwset_t glob_matches;
static bool spareparts_callback( int a, const kwset_t& O ) {
  kwset_t::const_iterator iter;
  for( iter = O.begin(); iter != O.end(); ++iter) {
    glob_matches.insert(*iter);
  }
  return( true ); // keep seaching for all matches
}
#endif //__LGPL__
//
// do a bit of cacheing of text for help with large redundant databases
//
#ifdef _DEBUG
static unsigned long bytes_unique=0;
static unsigned long bytes_saved=0;
#endif
void consolidate_text_ptr(char * &ptr,stringSet &strset) {
  stringSet::iterator s=strset.find(ptr);
  if (s != strset.end()) {
    delete ptr; // let go of redundant copy
    ptr = *s;
#ifdef _DEBUG
    bytes_saved+=(int)strlen(ptr)+1;
#endif
  } else {
    strset.insert(ptr);
#ifdef _DEBUG
    bytes_unique+=(int)strlen(ptr)+1;
#endif
  }
}
static stringSet aliases;
static stringSet descriptions;
static stringSet prev_aas;
static stringSet next_aas;

PeptideMapper::PeptideMapper(const char* xmlfile, const char* database, int min_num_tol_term,
                             int n_desired_prev_aas, int n_desired_next_aas, 
			     Boolean degen_only, Boolean prot_wt, const char *testMode, const char *restore_prefix) : Parser("database_refresh") {
  database_ = new char[strlen(database)+1];
  strcpy(database_, database);
  unCygwinify(database_); // remove any /cygdrive/ stuff etc

  degen_only_ = degen_only;
  enzymes_ = new StringArray;
  enzyme_digestions_ = new Array<ProteolyticEnzyme*>;
  enzyme_index_ = 0;

  testMode_ = testMode?strdup(testMode):NULL;

  restore_prefix_ = restore_prefix ? strdup(restore_prefix) : NULL;

  calc_prot_wt_ = prot_wt;

  min_num_tol_term_ = min_num_tol_term; // all peptides must have at least this number of tol term to be reported
  use_default_min_ntt_ = min_num_tol_term_ < 0;
  n_desired_prev_aas_ = n_desired_prev_aas; // resize the prev aa lists?
  n_desired_next_aas_ = n_desired_next_aas; // resize the follow aa lists?

  init(xmlfile);
}

PeptideMapper::~PeptideMapper() {
  delete enzymes_;
  if(enzyme_digestions_ != NULL) {
    for(int k = 0; k < enzyme_digestions_->length(); k++)
      if((*enzyme_digestions_)[k] != NULL)
	delete (*enzyme_digestions_)[k];
    delete enzyme_digestions_;
  }
  if(database_ != NULL)
    delete database_;
  //   if (szBuf_ != NULL)
  //     delete szBuf_;
  //   if (szSeq_ != NULL)
  //     delete szSeq_;
  //   if (szHdr_ != NULL)
  //     delete szHdr_;
  //   if (szBuf_ != NULL)
  //     delete szOutputDb_;
  free(testMode_);
}

void PeptideMapper::parse(const char* xmlfile) {
  char* engine = NULL;
  char* enzyme = NULL;
  char* massspec = NULL;
  Tag* tag = NULL;


  char* data = NULL;
  Array<char*>* inputfiles = new Array<char*>;
  
  //
  // regression test stuff - bpratt Insilicos LLC, Nov 2005
  //
  Array<Tag*> test_tags;
  eTagListFilePurpose testType;
  char *testFileName=NULL;
  checkRegressionTestArgs(testMode_,testType);
  if (testType!=NO_TEST) {
    testFileName = constructTagListFilename(xmlfile, // input file
					    testMode_, // program args
					    "PeptideMapper", // program name
					    testType); // user info output
  }
#define RECORD(tag) {(tag)->write(fout);if (testType!=NO_TEST) {test_tags.insertAtEnd(new Tag(*(tag)));}}


  double MIN_PROB =  0.0; //modelOpts_.minprob_; //0.0; //0.05; // for now

  // must do first run to determine whether icat or not, prior to initializing MixtureModel with options...
#ifndef __LGPL__
  kwset_t kwset;
  //  size_t offset;
  struct kwsmatch * match;
#endif

  char * fasta_db_file = NULL;
  char * out_file = NULL;
  char * interact_file = NULL;
  char * seq_string = NULL;
  FILE *fp;
  interact_data_line_t ** interact_data = NULL;
  interact_data_line_t ** uniq_interact_data = NULL;
  interact_data_line_t * interact_line = NULL;
  db_ref_t * new_db_ref = NULL;
  db_ref_t ** ref_list = NULL;
  long num_interact_lines = 0L;
  long num_hits = 0L;
  long uniq_lines = 0L;
  int i = 0;
  int opt = 0;
  int prev_index;
  time_t start_time = 0;
  time_t incr_time = 0;

  Boolean monoisotopic = False; // for calculating protein MW's
  double* prot_wt_list = NULL;

  int* ntt_list = NULL;
  int ntt;
  long line_num = 0L;
  int interact_lines_incr = 2000;
  interact_data_line_t ** interact_file_lines = NULL;


  Boolean enzyme_data = False;

  RACI fin(xmlfile); // can read gzipped xml
  if(! fin) {
    cerr << "PeptideMapper: error opening " << xmlfile << endl;
    exit(1);
  }
  char *nextline_tag = new char[line_width_];
  char *nextline_pep = new char[line_width_];
  
  while(fin.getline(nextline_tag, line_width_)) {
    fin.getline(nextline_pep, line_width_);
    //cout << "next: " << nextline << endl;
    if (line_num >= num_interact_lines ) {
      num_interact_lines += interact_lines_incr;
      interact_file_lines = (interact_data_line_t **)
	realloc(interact_file_lines,
		sizeof(interact_data_line_t *) * 
		(num_interact_lines));
      
    }
    interact_file_lines[line_num++] = getDataLine(nextline_tag, nextline_pep); // add on this one
    
    
  } // next line
  fin.close();



  // now do the db search

  uniq_lines = line_num;
  num_interact_lines = line_num;
  uniq_interact_data = build_uniq_kwlist(interact_file_lines,&uniq_lines);

#ifndef __LGPL__
  if ( !(kwset = kwsalloc((char *) 0)) )
    cerr << "ERROR: initializing the keyword set!\n";
 
  cerr << "INFO:  - Building Commentz-Walter keyword tree..." << endl;
  const char* prev_seq = "";
  const char* prev_alt = "";
  int tot_added = 0;
  Array<int>* alt_index_map = new Array<int>();
  for ( i = 0; i < uniq_lines; i++ ) {
    interact_line = (interact_data_line_t *)uniq_interact_data[i];

   
    if ( interact_line->sequence_I2L != NULL ) {
      if (  strcmp( interact_line->sequence_I2L,prev_seq) != 0 ) {
	
	prev_seq = interact_line->sequence_I2L;
	if ( kwsincr(kwset, 
		     interact_line->sequence_I2L,
		     strlen(interact_line->sequence_I2L)) != 0 ) 
	  fprintf(stderr, "Error adding %s to the keywords structure!\n",
		  interact_line->sequence_I2L);

	alt_index_map->insertAtEnd(i);

	tot_added++;
		
      }
    }
    
    if ( interact_line->alt_sequence != NULL ) {
      if ( kwsincr(kwset, 
                   interact_line->alt_sequence,
                   strlen(interact_line->alt_sequence)) != 0 ) 
	fprintf(stderr, "Error adding %s to the keywords structure!\n",
		interact_line->alt_sequence);
      
      alt_index_map->insertAtEnd(i);
      tot_added++;
	             
    }

  }

  /* Prep the keyword trie for searching */
  if ( kwsprep (kwset) != 0 )
    cerr << "Error prepping the structure for a search!\n";
#else
  /* SPARE Parts Section */
  kwset_t kwset;
  kwmap_t kwmap; // preserves order of uniq_lines
  int n_kwsets;
  int max_kwset_size = min(5000,uniq_lines); //limit to avoid trie width explosion, also may adjust downward dynamically
  std::vector<STRINGFINDER_pep *> finder_pep;
  std::vector<STRINGFINDER_std *> finder_std;
  cerr << "INFO:  - Building Commentz-Walter keyword tree for " << uniq_lines << " unique peptides...";
  bool bUsePeptideAlphabet = true; // use smaller alphabet if posssible
  bool firstpass = true; // we may have to retry on trie width exception
  while (true) { // dynamic sizing of tries
    std::vector<kwset_t> kwsets; // break large jobs into multiple sets
    n_kwsets = (uniq_lines/max_kwset_size)+((uniq_lines%max_kwset_size)!=0);
    kwsets.resize(n_kwsets);
    try {
      peptideAlphabet::init(); // initialize lookup tables, etc
      const char *badStr=NULL;
      for ( i = 0; i < uniq_lines; i++ ) {
	interact_line = (interact_data_line_t *)uniq_interact_data[i];
    
	if ( interact_line->sequence != NULL ) {
	  kwsets[i/max_kwset_size].insert(interact_line->sequence_I2L);
	  if (bUsePeptideAlphabet && !peptideAlphabet::legal_string(interact_line->sequence_I2L)) {
	    badStr = interact_line->sequence_I2L;
	    bUsePeptideAlphabet = false;
	  }
	  if (firstpass && (i!=kwmap[interact_line->sequence_I2L])) { // will insert a record
	    kwmap.find(interact_line->sequence_I2L)->second = i;
	  }
	}
    
	if ( interact_line->alt_sequence != NULL ) {
	  kwsets[i/max_kwset_size].insert(interact_line->alt_sequence);
	  if (bUsePeptideAlphabet && !peptideAlphabet::legal_string(interact_line->alt_sequence)) {
	    badStr = interact_line->alt_sequence;
	    bUsePeptideAlphabet = false;
	  }
	  if (firstpass && (i!=kwmap[interact_line->alt_sequence])) { // will insert a record
	    kwmap.find(interact_line->alt_sequence)->second = i;
	  }
	}
      }
      /* Prep the keyword trie for searching */
      if (bUsePeptideAlphabet) {
	firstpass = false; 
	for (int n=n_kwsets;n--;) {
	  finder_pep.push_back(new STRINGFINDER_pep(kwsets[n])); 
	}
      } else {
	if (firstpass) {
	  cerr << "WARNING: unexpected character in peptide sequence \"" << badStr << "\", using more general but less efficient search " << endl;
	}
	firstpass = false; 
	for (int n=n_kwsets;n--;) {
	  finder_std.push_back(new STRINGFINDER_std(kwsets[n])); 
	}
      }
      break; // success
    } // end try
    catch (...) {
      // go back and try it again with smaller trie width
      max_kwset_size /= 2; 
      if (bUsePeptideAlphabet) {
	for (size_t n=finder_pep.size();n--;) {
	  delete finder_pep[n];
	} 
	finder_pep.clear();
      } else {
	for (size_t n=finder_std.size();n--;) {
	  delete finder_std[n];
	}
	finder_std.clear();
      }
    }
  }
#endif //__LGPL__
  /* Perform the search */
  cerr << "INFO:  - Searching the tree..." << endl;
  if ( ! (fp = fopen(database_,"r"))){
    std::string rdb(resolve_root(database_));
    const char *trythis = rdb.c_str();
    if (strcmp(trythis, database_)) {
      int err = errno;
      if (! (fp = fopen(trythis,"r"))) {
	fprintf(stderr, "error: could not open database file \"%s\" (%s)\n", database_, strerror(err));
	fprintf(stderr, "also tried \"%s\" (%s)\n", trythis, strerror(errno));
      } else {
	fprintf(stderr, "opening \"%s\" as \"%s\"\n",database_,trythis);
      }
    } else {
      fprintf(stderr, "error: could not open database file \"%s\" (%s)\n", database_, strerror(errno));
    }
    if (!fp) {
      exit(-1);
    }
  }

  std::string outfilerep = xmlfile;
  std::string outfilemap = xmlfile;
  stringSet* replSet = new stringSet();
  outfilerep += ".replace";
  outfilemap += ".map";

  //cerr << "writing data to " << outfile << endl;

  ofstream foutrep(outfilerep.c_str());

  ofstream foutmap(outfilemap.c_str());
  



  int iLenSeq;
   
  char *szBuf_ = (char*) malloc(sizeof(char)*line_width_);
  char *szSeq_ = (char*) malloc(sizeof(char)*MAX_SEQ);
  char *szOrigSeq_ = (char*) malloc(sizeof(char)*MAX_SEQ);
  char *szHdr_ = (char*) malloc(sizeof(char)*MAX_HEADER_LEN); 
  char *szOutputDb_ = (char*) malloc(sizeof(char)*SIZE_FILE);
  int* sort_idx = NULL;
  int max_num_hits = 0;

  unsigned long maxSeqLen = MAX_SEQ;
   
  szBuf_[0] = '\0';
  szSeq_[0] = '\0';
  szOrigSeq_[0] = '\0';
  szHdr_[0] = '\0';
  szOutputDb_[0] = '\0';



  while (fgets(szBuf_, line_width_, fp))

    {
      if (szBuf_[0]=='>')
	{
	  int cAA;
	 
	  strncpy(szHdr_, szBuf_, MAX_HEADER_LEN);
	  szHdr_[MAX_HEADER_LEN-1] = '\0';
	  if(strlen(szBuf_) < MAX_HEADER_LEN)
	    szHdr_[strlen(szBuf_)-1] = '\0';
	 
	  // now get rid of quotations
	  char* quotematch = strchr(szHdr_, '"');
	  while(quotematch != NULL) {
	    memmove(quotematch,quotematch+1,strlen(quotematch));
	    quotematch = strchr(quotematch, '"');
	  } // while

	  iLenSeq=0;

	  while ( (cAA=fgetc(fp)) )
	    {
	      if (isalpha(cAA) || cAA=='*')
		{
		  szOrigSeq_[iLenSeq]=cAA;
		  if (cAA == 'I') {
		    cAA = 'L';
		  }
		  if (cAA == '*') {
		    cAA = '-';
		  }
		  szSeq_[iLenSeq++]=cAA;
		  if (iLenSeq > maxSeqLen)
		    {
                 
		      szSeq_ =  (char*) realloc(szSeq_,
						sizeof(char)*iLenSeq*2);
		      szOrigSeq_ =  (char*) realloc(szOrigSeq_,
						    sizeof(char)*iLenSeq*2);
		      maxSeqLen = iLenSeq*2;
		      //cerr << "ERROR: - sequence larger than MAX_SEQ\n\n";
		      //              fclose(fp);
		      //exit(1);
		    }
		}
	      else if (feof(fp) || cAA=='>')
		{
		  int iReturn;

		  iReturn=ungetc(cAA, fp);

		  if (iReturn!=cAA)
		    {
		      fprintf(stderr, "Error with ungetc.\n\n");
		      fclose(fp);
		      exit(1);
		    }
		  break;
		}
	    }

	  szSeq_[iLenSeq]='\0';
	  szOrigSeq_[iLenSeq]='\0';
	 
#ifdef __LGPL__	 
	  /* SPARE Parts */
	  glob_matches.clear(); // nothing found yet
	  if (bUsePeptideAlphabet) {
	    for (int n=n_kwsets;n--;) {
	      finder_pep[n]->match(szSeq_,spareparts_callback); // any proteins match?
	    }
	  } else {
	    for (int n=n_kwsets;n--;) {
	      finder_std[n]->match(szSeq_,spareparts_callback); // any proteins match?
	    }
	  }
	  num_hits = (int)glob_matches.size();
	  kwset_t::const_iterator iter = glob_matches.begin();
#else	  
	  num_hits = kwsexec_multiple(kwset, szSeq_, strlen(szSeq_), &match);
#endif //__LGPL__
	 
	  //DDS: sort the matched indices
	  if (num_hits > 0) {
	   
	    if (num_hits>max_num_hits) {
	      sort_idx = (int *)realloc(sort_idx,(max_num_hits=num_hits)*sizeof(int));
	    }
	    for (i=0; i<num_hits; i++){
#ifdef __LGPL__	 
	      sort_idx[i] = kwmap.find(iter->c_str())->second;
	      iter++;
#else	  
	      sort_idx[i] = match[i].index;
#endif 
	    }
	    qsort(sort_idx, num_hits, sizeof(int), compare_indices);
       
	    // Process the hits for this protein sequence
	    prev_index = -1;
	   
	    for ( i=0; i<num_hits; i++){
	     
	      int this_index = sort_idx[i];
	      this_index = (*alt_index_map)[this_index];

	      const char* result = NULL;

	      if (uniq_interact_data[this_index]) {
		result = strstr(szSeq_, ((interact_data_line_t *)
					 uniq_interact_data[this_index])->sequence_I2L);
	      }
	      bool alt = false;
	      if (!result) {
		alt = true;
		result = strstr(szSeq_, ((interact_data_line_t *)
					 uniq_interact_data[this_index])->alt_sequence);
	      }
	     

	      //DDS: For this to work the indices must be sorted
	      if ( this_index != prev_index ) {
	       
		if ( (new_db_ref = (db_ref_t *)malloc(sizeof(db_ref_t))) == NULL )
		  fprintf(stderr, "Error: Could not allocate memory for db_ref_t!\n");
	       
		/*
		 * use space or | as delimiter of accession if IPI database
		 */
		char cDelimiter;

		if (!strncmp(szHdr_, ">IPI", 4) && strchr(szHdr_, '|') != NULL && strchr(szHdr_, '|')<strchr(szHdr_,' ') )
		  cDelimiter = '|';
		else
		  cDelimiter = ' ';
                
		if ( strchr(szHdr_,cDelimiter) != NULL ) {
		  new_db_ref->alias = make_substr(szHdr_,
						  0,
						  strchr(szHdr_,cDelimiter) - 
						  szHdr_ - 1);
		 
		  // get rid of windows control characters (^M, etc) by truncating at the first one
		  size_t len = strlen(new_db_ref->alias);
		  for (size_t i = 0; i < len; i++) {
		    if (new_db_ref->alias[i] < 20 || new_db_ref->alias[i] > 126) {
		      new_db_ref->alias[i] = 0;
		      break;
		    }
		  }
		 
		  new_db_ref->description = make_substr(szHdr_,
							strchr(szHdr_,cDelimiter) - 
							szHdr_,
							(int)strlen(szHdr_));

		  // get rid of windows control characters (^M, etc) by truncating at the first one
		  len = strlen(new_db_ref->description);
		  for (size_t i = 0; i < len; i++) {
		    if (new_db_ref->description[i] < 20 || new_db_ref->description[i] > 126) {
		      new_db_ref->description[i] = 0;
		      break;
		    }
		  }
		 
		 
		  // make sure no xml problem characters in prot descr
		  for(int z = 0; z < (int) strlen(new_db_ref->description); z++)
		    if(new_db_ref->description[z] == '&')
		      new_db_ref->description[z] = '+';
		    else if(new_db_ref->description[z] == '"')
		      new_db_ref->description[z] = '\'';
		    else if(new_db_ref->description[z] == '>') {
		      //if(verbose)
		      // cout << "got one!" << endl;
		      new_db_ref->description[z] = ' ';
		    }
		    else if(new_db_ref->description[z] == '<')
		      new_db_ref->description[z] = ' ';
		    else if(new_db_ref->description[z] == '/')
		      new_db_ref->description[z] = '\\';
	 
		}else {
		  new_db_ref->alias = strdup(szHdr_);
		  // get rid of windows control characters (^M, etc) by truncating at the first one
		  size_t len = strlen(new_db_ref->alias);
		  for (size_t i = 0; i < len; i++) {
		    if (new_db_ref->alias[i] < 20 || new_db_ref->alias[i] > 126) {
		      new_db_ref->alias[i] = 0;
		      break;
		    }
		  }
		  new_db_ref->description = new char[2];
		  new_db_ref->description[0] = 0;
		}
	 
	       
		ref_list = ((interact_data_line_t *)
			    uniq_interact_data[this_index])->updated_refs;
		ntt_list = ((interact_data_line_t *)uniq_interact_data[this_index])->updated_ntts;
	       
		if(calc_prot_wt_)
		  prot_wt_list = ((interact_data_line_t *)uniq_interact_data[this_index])->updated_prot_wts;
	       
	       
		if ( uniq_interact_data[this_index]->updated_hits == 0 ){
		  if ( (ref_list = (db_ref_t **)malloc( sizeof(db_ref_t *))) == NULL ) 
		    fprintf(stderr, "Error: Could not allocate memory for a new ref_list!\n");
		  if ( (ntt_list = (int*)malloc( sizeof(int *))) == NULL ) 
		    fprintf(stderr, "Error: Could not allocate memory for a new ntt_list!\n");
		 
		  if (calc_prot_wt_ && (prot_wt_list = (double*)malloc( sizeof(double *))) == NULL ) 
		    fprintf(stderr, "Error: Could not allocate memory for a new prot_wt_list!\n");
		 
		 
		}else{
		  if ( (ref_list = (db_ref_t **)realloc(ref_list,
							(((interact_data_line_t *)
							  uniq_interact_data[this_index])->updated_hits + 1) * 
							sizeof(db_ref_t *))) == NULL ) 
		    fprintf(stderr, "Error: Could not reallocate memory for ref_list!\n");
		  if ( (ntt_list = (int*)realloc(ntt_list,
						 (((interact_data_line_t *)
						   uniq_interact_data[this_index])->updated_hits + 1) * 
						 sizeof(int *))) == NULL ) 
		    fprintf(stderr, "Error: Could not reallocate memory for ntt_list!\n");
		 
		  if (calc_prot_wt_ && (prot_wt_list = (double*)realloc(prot_wt_list,
									(((interact_data_line_t *)
									  uniq_interact_data[this_index])->updated_hits + 1) * 
									sizeof(double *))) == NULL ) 
		    fprintf(stderr, "Error: Could not reallocate memory for prot_wt_list!\n");
		 
		 
		 
		}
		ref_list[uniq_interact_data[this_index]->updated_hits] = new_db_ref;
		ntt = -1;
	       

		char prev, next;
		ntt = 0;
		bool firstpass = true;
		if(result != NULL) {
		  size_t seqlen=strlen(szSeq_);
		  //DDS: Find all occurrences of the peptide in the protein
		  while (result != NULL) { 
		   
		    size_t reslen = strlen(result);
		    if(reslen == seqlen)
		      prev = '-';
		    else
		      prev = szOrigSeq_[result-szSeq_-1];//prev = (result-1)[0];
		   
		    // change M to - when it's first character in the Protein
		    if(prev == 'M' && reslen >= seqlen - 1) {
		      prev = '-';
		    } 
	
		    size_t nextlen=strlen(((interact_data_line_t *)
					   uniq_interact_data[this_index])->sequence_I2L);
							      
							      
		  
		    if(reslen == nextlen)
		      next = '-';
		    else
		      next = szOrigSeq_[result-szSeq_+strlen(((interact_data_line_t *)
							      uniq_interact_data[this_index])->sequence)];//		     next = result[nextlen]; // next aa

		    char* origPep = ((interact_data_line_t *)uniq_interact_data[this_index])->sequence;
		    char* replPep = make_substr(szOrigSeq_, 
						result-szSeq_, 
						result-szSeq_+strlen(((interact_data_line_t *)
								      uniq_interact_data[this_index])->sequence)-1);

		    char* pepAcs  = new char[strlen("PAp00000000")+1]; //make_substr(((interact_data_line_t *)uniq_interact_data[this_index])->dbname, 1,strlen(((interact_data_line_t *)uniq_interact_data[this_index])->dbname)-1);
		    strcpy(pepAcs,"PAp00000000");
		    if (pepAcsHash_.find(replPep) != pepAcsHash_.end()) {
		      pepAcs = make_substr((char*) pepAcsHash_[replPep],1,strlen(pepAcsHash_[replPep])-1);
		    }

		    string key = origPep;
		    key += replPep;
		    key += pepAcs;

		    char * ckey = new char[key.length()+1];
		    strcpy(ckey, key.c_str());

		    if ( replSet->find(ckey) == replSet->end() ) {
		     
		      foutrep << origPep
			      << "\t"
			      << replPep
			      << "\t" 
			      << pepAcs
			      << "\t" 
			      << make_substr(new_db_ref->alias, 1, strlen(new_db_ref->alias)-1)
			      << endl;
		      replSet->insert(ckey);
		    }
		 
		    foutmap << pepAcs
		   
			    << "\t"
		   
			    << make_substr(szOrigSeq_, result-szSeq_, result-szSeq_+strlen(((interact_data_line_t *)uniq_interact_data[this_index])->sequence)-1)
			    << "\t"
		   
			    << make_substr(new_db_ref->alias, 1, strlen(new_db_ref->alias)-1)
		      //<< result-szSeq_ 
			    << "\t"
			    << result-szSeq_+1
			    << "\t"
			    << result-szSeq_+strlen(((interact_data_line_t *)uniq_interact_data[this_index])->sequence)
			    << "\t"
			    << prev
			    << "\t"		    
			    << next << endl;
		 
	
		   

		    if (firstpass) { // do this just once
		      new_db_ref->prev_aa = (char *)malloc(n_desired_prev_aas_+1);
		      new_db_ref->next_aa = (char *)malloc(n_desired_next_aas_+1);
		     
		      // resize the prev AA string

		      if (NULL != new_db_ref->prev_aa ) {
			const char *aa=result;
			int count = 0;
			while ((count<n_desired_prev_aas_) && (aa>szSeq_)) {
			  aa--;
			  count++;
			}
			char *w=new_db_ref->prev_aa;
			while (count--) {
			  *w++ = *aa++;
			}
			*w = 0;
			if (!new_db_ref->prev_aa[0]) {
			  strcpy(new_db_ref->prev_aa,"-");
			}
		      }
		      // resize the follow AA string
		      if (NULL != new_db_ref->next_aa) {
			const char *aa=result+nextlen;
			int count = 0;
			while ((count<n_desired_next_aas_) && *aa) {
			  aa++;
			  count++;
			}
			char *w=new_db_ref->next_aa;
			aa=result+nextlen;
			while (count--) {
			  *w++ = *aa++;
			}
			*w = 0;
			if (!new_db_ref->next_aa[0]) {
			  strcpy(new_db_ref->next_aa,"-");
			}
		      }
		    }
		    result = strstr(result+strlen(((interact_data_line_t *)uniq_interact_data[this_index])->sequence_I2L), 
				    ((interact_data_line_t *)uniq_interact_data[this_index])->sequence_I2L);
		    firstpass = false;
		  }

 
		} // if not null
		else {
		 
		  cerr << "ERROR: " << ((interact_data_line_t *)uniq_interact_data[this_index])->sequence_I2L << " not found in " <<  szSeq_ << endl; 
		  exit(1);
		}
	       
		// now see if we can reuse existing text
		consolidate_text_ptr(new_db_ref->alias,aliases);
		consolidate_text_ptr(new_db_ref->description,descriptions);
		consolidate_text_ptr(new_db_ref->prev_aa,prev_aas);
		consolidate_text_ptr(new_db_ref->next_aa,next_aas);
	       
	       
		if(calc_prot_wt_)
		  prot_wt_list[uniq_interact_data[this_index]->updated_hits] = ResidueMass::getProteinMass(szSeq_, monoisotopic);
	       
	       
		uniq_interact_data[this_index]->updated_hits += 1;
		((interact_data_line_t *)uniq_interact_data[this_index])->updated_refs = ref_list;

	       
		if(calc_prot_wt_)
		  ((interact_data_line_t *)uniq_interact_data[this_index])->updated_prot_wts = prot_wt_list;
	       
		prev_index = this_index;
	      } // end if this_index != prev_index
	    } // end for num_hits
	  } // end if num_hits
	       
	} // while fgets database
    } // while fin.getline
  free(sort_idx);
#ifdef __LGPL__
  /*SPARE Parts Section*/
  for (int n=n_kwsets;n--;) {
    if (bUsePeptideAlphabet) {
      delete finder_pep[n];
    } else {
      delete finder_std[n];
    }
  }
#endif
  fclose( fp );

  //fprintf(stdout,"\n  - Linking duplicate entries...");
  //link_duplicates(interact_file_lines,num_interact_lines);

  incr_time = time(NULL);
 
  free(szBuf_);
  free(szSeq_);
  free(szHdr_);
  free(szOutputDb_);
  free(uniq_interact_data);
}
  


void PeptideMapper::setFilter(Tag* tag) {
  if(tag == NULL)
    return;

  if(filter_memory_) {
    filter_memory_ = False;
    filter_ = False;
  }

  if(! strcmp(tag->getName(), "spectrum_query")){ 
    if(tag->isStart()) {
      filter_ = True;
    }else{
      filter_memory_ = True;
    }
  }

}

Boolean PeptideMapper::refresh(Tag* tag) {

  return (!degen_only_ || strcmp(tag->getAttributeValue("num_tot_proteins"), "1"));
  
  //  && (!ignore_prefix_ || strncmp(tag->getAttributeValue("protein"), ignore_prefix_, strlen(ignore_prefix_)) != 0));
}

int PeptideMapper::getEnzymeIndex(char* enz) {
  return enzymes_->findByStringValue(enz);
}


interact_data_line_t* PeptideMapper::getDataLine(char* line_tag, char* line_pep ) {
  
  interact_data_line_t* interact_data_line = (interact_data_line_t*)calloc(sizeof(interact_data_line_t),1);

  interact_data_line->hits = 1; // default to 1

  interact_data_line->sequence = new char[strlen(line_pep)+1];
  strcpy(interact_data_line->sequence, line_pep);
  interact_data_line->sequence_I2L = new char[strlen(line_pep)+1];
  strcpy(interact_data_line->sequence_I2L, line_pep);

  interact_data_line->dbname = new char[strlen(line_tag)+1];
  strcpy(interact_data_line->dbname, line_tag);

  if (pepAcsHash_.find(interact_data_line->sequence) == pepAcsHash_.end()) {
    pepAcsHash_.insert(make_pair(interact_data_line->sequence, interact_data_line->dbname));
  }

  interact_data_line->enzyme_ind = enzyme_index_;

  for (int i =0; i < (int)strlen(interact_data_line->sequence_I2L); i++) {
    if (interact_data_line->sequence_I2L[i] == 'I') {
      interact_data_line->sequence_I2L[i] = 'L';
    }
  }

  interact_data_line->alt_sequence = NULL;


  

  interact_data_line->min_num_enz_term = 0; //min_num_tol_term_;


  return interact_data_line;
}

// char * make_substr(char * string, long start, long end)
//
// Create a new string from the substr of an existing
// string.  This routine allocates memory for the new
// string which can be released using free().
char* PeptideMapper::make_substr(const char * string, long start, long end) 
{
  char * substr;
  /*
    int verbose = 0;

    if ( verbose > 1 ) 
    fprintf(stderr, "make_substr: called with s=%s start=%ld end=%ld\n",
    string, start, end);
  */
  if(end > start && string[start] == ' ')
    start++;
  substr = (char *)malloc((end-start+2) * sizeof( char ));
  strncpy(substr,string+start,(end-start+1));
  // Ask paul about this why is this not equivalent to:
  //   *(substr+end-start+1) = '\0';
  substr[end-start+1] = '\0';
  return substr;
}




char * PeptideMapper::replace_substr(char * string, char * orig_substr, char * repl_substr ) 
{
  char * substr_ptr;
  char * suffix_ptr;
  char * search_ptr;
  size_t substr_len_diff = 0L;
  long offset = 0L;
  int verbose = 0;

  if ( verbose > 1 ) 
    //    fprintf(stderr, "replace_substr: called with s=%s o=%s r=%s\n",
    //      string, orig_substr, repl_substr);

    substr_len_diff =  strlen(repl_substr) - strlen(orig_substr);
  search_ptr = string;

  while ( (substr_ptr = strstr(search_ptr,orig_substr)) != NULL ) {

    if ( substr_len_diff > 0 ) {
      offset = substr_ptr - string;
      char* tmpString = string;
      //      string = (char*)realloc(string,
      //                 (strlen(string) + substr_len_diff + 1)*sizeof(char));
      string = (char*)malloc((strlen(string) + substr_len_diff + 1)*sizeof(char));
      memcpy(string, tmpString, (strlen(string) + substr_len_diff + 1)*sizeof(char));
      free(tmpString);
      substr_ptr = string + offset;
    }
    suffix_ptr = substr_ptr + strlen(orig_substr);
    memmove(suffix_ptr+substr_len_diff,
	    suffix_ptr,
	    strlen(suffix_ptr)+1);
    strncpy(substr_ptr,repl_substr,strlen(repl_substr));

    search_ptr = substr_ptr + strlen(repl_substr);
  }
  // TODO: could free some memory if we shrunk the string
  return string;
}



interact_data_line_t ** PeptideMapper::build_uniq_kwlist(interact_data_line_t ** interact_data,long * interact_lines)
{
  int i = 0;
  interact_data_line_t ** sorted_lines;
  interact_data_line_t ** uniq_lines;
  long u_index;
  const char * prev_sequence;
  const char * prev_alt;
  int prev_index = -1;
  int verbose = 0;

  if ( verbose > 1 ) 
    cerr << "INFO: build_uniq_kwlist: called\n";


  //cout << "interact lines: " << *interact_lines << endl;

  sorted_lines = (interact_data_line_t **)malloc(sizeof(interact_data_line_t *)* *interact_lines);
  uniq_lines = (interact_data_line_t **)malloc(sizeof(interact_data_line_t *)* *interact_lines);
  //cout << "0" << endl;

  interact_data_line_t **sp = sorted_lines;
  interact_data_line_t **idp = interact_data;
  for ( i = *interact_lines; i--; ) { 
    *sp++ = *idp++;
  }
  qsort(sorted_lines,*interact_lines,sizeof(interact_data_line_t *),
	compare_interact_lines);

  //for(int k = 0; k < 10; k++)
  // cout << "pep: " << sorted_lines[k]->sequence << ": " << sorted_lines[k]->enzyme_ind << endl;
  //cout << "2" << endl;
  u_index = 0;
  prev_sequence = "";
  sp = sorted_lines;
  interact_data_line_t *id;
  prev_sequence = "";
  prev_alt = NULL;
  for ( i = *interact_lines; i--; ) { 
  
    if ( (id = *sp++)->sequence != NULL ) {
      if ( //prev_index != id->enzyme_ind ||
	  strcmp(id->sequence_I2L,prev_sequence) != 0 ) {
	prev_sequence = id->sequence_I2L; 
	prev_alt = id->alt_sequence; 
	//prev_index = id->enzyme_ind;
	uniq_lines[u_index++] = id;
      }
      else if (id->alt_sequence || prev_alt) {
	if ( (id->alt_sequence && !prev_alt) ||
	     (!id->alt_sequence && prev_alt) ||
	     (id->alt_sequence && prev_alt &&  strcmp(id->alt_sequence,prev_alt) != 0)
	     )
	  {
	    prev_sequence = id->sequence_I2L; 
	    prev_alt = id->alt_sequence; 
	    //prev_index = id->enzyme_ind;
	    uniq_lines[u_index++] = id;
	  }
      }


    }
  }
  free(sorted_lines);
  *interact_lines = u_index;
  return(uniq_lines);

}



int PeptideMapper::link_duplicates(interact_data_line_t ** interact_data,long interact_lines)
{
  int i = 0;
  interact_data_line_t ** sorted_lines;
  db_ref_t ** ref_ptr=NULL;
  char *prev_aa=NULL;
  char *next_aa=NULL;
  const char * sequence;
  const char * alt_sequence = NULL;
  long ref_hits=-1;
  int enzyme_ind;
  int verbose = 0;

  double* prot_wt_ptr = NULL;

  if ( verbose > 1 ) 
    cerr << "link_duplicates: called\n";

  //fprintf(stderr, "Allocating memory for sort\n");
  sorted_lines = (interact_data_line_t **)malloc(sizeof(interact_data_line_t *)*interact_lines);

  //fprintf(stderr, "Creating array for sort\n");
  for ( i = 0; i < interact_lines; i++ )
    sorted_lines[i] = (interact_data_line_t *)interact_data[i];

  //fprintf(stderr, "sorting...\n");
  qsort(sorted_lines,interact_lines,sizeof(interact_data_line_t *),
        compare_interact_lines);

  //fprintf(stderr, "Printing sorted data\n");
  //for ( i = 0; i < interact_lines; i++ ) 
  // fprintf(stderr, "sequence = %s  pointer = %d\n", sorted_lines[i]->sequence,
  //        sorted_lines[i]->updated_refs);


  sequence = "start";
  
  for ( i = 0; i < interact_lines; i++ ) {
    if ( sorted_lines[i]->sequence != NULL ) {   // must be a data line
      if ( sorted_lines[i]->updated_refs != NULL ) {
        ref_ptr = sorted_lines[i]->updated_refs;

	if(calc_prot_wt_)
	  prot_wt_ptr = sorted_lines[i]->updated_prot_wts;

        ref_hits = sorted_lines[i]->updated_hits;

        sequence = sorted_lines[i]->sequence_I2L;	

	alt_sequence = sorted_lines[i]->alt_sequence;

	enzyme_ind = sorted_lines[i]->enzyme_ind;

      }else {

        if ( strcmp(sequence, sorted_lines[i]->sequence_I2L) == 0 && 
	     ((!alt_sequence && !sorted_lines[i]->alt_sequence) ||
	      (alt_sequence && sorted_lines[i]->alt_sequence && strcmp(alt_sequence, sorted_lines[i]->alt_sequence))==0 )) {

          //fprintf(stderr, "updating:\n%s\n",sorted_lines[i]->line);
          sorted_lines[i]->updated_refs = ref_ptr;
	  
	  if(calc_prot_wt_)
	    sorted_lines[i]->updated_prot_wts = prot_wt_ptr;

	  sorted_lines[i]->updated_hits = ref_hits;
          //fprintf(stderr, "now:\n%s\n",sorted_lines[i]->line);
        }
      }
    }
  }
  return 0;
}

int compare_indices(const void * int1, const void * int2) {
  if (*(int*)int1 > *(int*)int2) {
    return 1;
  }
  else {
    return 0;
  }
  
}

// int compare_interact_lines(const void * line1, const void * line2)
//
int compare_interact_lines(const void * line1, const void * line2)
{
  interact_data_line_t * linet1 = *((interact_data_line_t **)line1);
  interact_data_line_t * linet2 = *((interact_data_line_t **)line2);
 
 
  const char * s1, * s2;

  
  if ( (s1 = linet1->sequence_I2L) == NULL || 
       (s2 = linet2->sequence_I2L) == NULL  ) {
    if ( linet1->sequence_I2L != NULL ) 
      return 1;
    if ( linet2->sequence_I2L != NULL  ) 
      return -1;
    return 0;
  }
  else if (0) {
    if ( (s1 = linet1->alt_sequence) == NULL || 
	 (s2 = linet2->alt_sequence) == NULL  ) {
      if ( linet1->alt_sequence != NULL ) 
	return 1;
      if ( linet2->alt_sequence != NULL  ) 
	return -1;
      return 0;
    }
  }
  int return_val;
  
  if (!(return_val = strcmp(s1,s2))) {
    if ( linet1->updated_refs == NULL || 
         linet2->updated_refs == NULL  ) {
      if ( linet1->updated_refs != NULL ) 
	return -1;
      if ( linet2->updated_refs != NULL  ) 
	return 1;
      return 0;
    }
    
  }

  return return_val;

}


// char * long_to_string(long value)
//
char * PeptideMapper::long_to_string(long value) 
{
  long size = 100;
  long conv_size = 0;
  char *string, *tmpstring;
  int verbose = 0;

  if ( verbose > 8 ) 
    cerr <<"long_to_string: called\n";
  
  tmpstring = (char*)malloc(sizeof(char) * size );
  conv_size = snprintf(tmpstring,size,"%ld",value);
  if ( conv_size > 0 ) {
    string = (char*)malloc(sizeof(char) * (conv_size+2));
    memcpy(string, tmpstring, sizeof(char) * (conv_size+2));
    //   string = (char*)realloc(string,conv_size + 2);
  }
  else
    {
      free(tmpstring);
      return NULL;
    }
  free(tmpstring);
  return string;
}


