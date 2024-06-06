/*
  Program       : RefreshParser
  Author        : Andrew Keller <akeller@systemsbiology.org>, Robert Hubley, and open source code
  Date          : 11.27.02
  SVN Info      : $Id: RefreshParser.cpp 9086 2023-12-19 00:14:09Z dshteyn $

  Copyright (C) 2003 Andrew Keller
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
#include "RefreshParser.h"
#include "Parsers/Parser/TagListComparator.h" // regression test stuff - bpratt Insilicos LLC, Nov 2005
#include "Common/util.h"
#include <errno.h>


static int compare_indices(const void * int1, const void * int2);
struct equalstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) == 0;
  }
};

#ifdef _MSC_VER
#include <hash_map> // if you're using VC6 you need STLPort 
#include <hash_set>
#ifdef _STLP_HASH_MAP 
// a bit of ugly cut and paste to deal with VC6 template borkeness
inline size_t __stl_hash_string(const char* __s)
{
  _STLP_FIX_LITERAL_BUG(__s)
    unsigned long __h = 0; 
  for ( ; *__s; ++__s)
    __h = 5*__h + *__s;
  
  return size_t(__h);
}
struct cchash
{
  size_t operator()(const char* __s) const { _STLP_FIX_LITERAL_BUG(__s) return __stl_hash_string(__s); }
};
typedef std::hash<const char *> cchash;
typedef std::hash_map<const char *,int,cchash,equalstr> kwmap_t;
typedef std::hash_set<char *,cchash,equalstr> stringSet;
#else // VC 7 or 8
struct strltpred
{
  bool operator()(const char* a, const char* b) const
  {
    return strcmp(a, b) < 0;
  }
};
typedef stdext::hash_map<const char *,int,stdext::hash_compare<const char*, strltpred>> kwmap_t;
typedef stdext::hash_set<char *,stdext::hash_compare<const char*, strltpred>> stringSet;
#endif
#else // GCC
#if (__GNUC__ < 4)
#include <hash_map.h>
typedef hash<const char *> cchash;
typedef hash_map<const char *,int,cchash,equalstr> kwmap_t;
#include <hash_set.h>
typedef hash_set<char *,cchash,equalstr> stringSet;
#else
#include <ext/hash_map>
typedef __gnu_cxx::hash<const char *> cchash;
typedef __gnu_cxx::hash_map<const char *,int,cchash,equalstr> kwmap_t;
#include <ext/hash_set>
typedef __gnu_cxx::hash_set<char *,cchash,equalstr> stringSet;
#endif
#endif

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
    free(ptr); // let go of redundant copy
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


void save_ptrs_to_free(interact_data_line_t * line, map <db_ref_t**, int> *freed, map <db_ref_t*, int> *freedp, map <char*, int> *freedc) {
  for (int h=0; h < line->updated_hits; h++) {    
    (*freedc)[line->updated_refs[h]->alias] = 1;
    (*freedc)[line->updated_refs[h]->description] = 1;
    (*freedc)[line->updated_refs[h]->prev_aa] = 1;
    (*freedc)[line->updated_refs[h]->next_aa] = 1;
    (*freedp)[line->updated_refs[h]] = 1;
  }
  (*freed)[line->updated_refs] = 1;
}

static stringSet aliases;
static stringSet descriptions;
static stringSet prev_aas;
static stringSet next_aas;


RefreshParser::RefreshParser(const char* xmlfile, const char* database, int min_num_tol_term,
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
  num_mapped_ = 0;
  num_not_mapped_ = 0;

  calc_prot_wt_ = prot_wt;

  min_num_tol_term_ = min_num_tol_term; // all peptides must have at least this number of tol term to be reported
  use_default_min_ntt_ = min_num_tol_term_ < 0;
  n_desired_prev_aas_ = n_desired_prev_aas; // resize the prev aa lists?
  n_desired_next_aas_ = n_desired_next_aas; // resize the follow aa lists?

  init(xmlfile);
}

RefreshParser::~RefreshParser() {
  delete enzymes_;
  if(enzyme_digestions_ != NULL) {
    for(int k = 0; k < enzyme_digestions_->length(); k++)
      if((*enzyme_digestions_)[k] != NULL)
	delete (*enzyme_digestions_)[k];
    delete enzyme_digestions_;
  }
  if(database_ != NULL)
    delete[] database_;

  free(testMode_);
}

void RefreshParser::parse(const char* xmlfile) {
  char* engine = NULL;
  char* enzyme = NULL;
  char* massspec = NULL;
  char* data = NULL;
  Tag*  tag = NULL;
  
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
					    "RefreshParser", // program name
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

  interact_data_line_t ** interact_file_lines = NULL;
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

  Boolean enzyme_data = False;

  RACI fin(xmlfile); // can read gzipped xml
  if(! fin) {
    cerr << "RefreshParser: error opening " << xmlfile << endl;
    exit(1);
  }
  char *nextline = new char[line_width_];
  while(fin.getline(nextline, line_width_)) {
    //cout << "next: " << nextline << endl;

    data = strchr(nextline, '<');
    while(data != NULL) {
      tag = new Tag(data);
      //setFilter(tag);
      //tag->write(cout);
      if(tag != NULL) {
	if ( (tag->isStart() && ! strcmp(tag->getName(), "linked_peptide")) || 
	     (tag->isStart() && ! strcmp(tag->getName(), "search_hit") && 
	      ! strcmp(tag->getAttributeValue("hit_rank"), "1") &&  
	      (!tag->getAttributeValue("is_rejected") || 
	       strcmp(tag->getAttributeValue("is_rejected"), "1")) &&
	      refresh(tag)) ) {
	  if (line_num >= num_interact_lines ) {
	    num_interact_lines += interact_lines_incr;
	    interact_file_lines = (interact_data_line_t **)
	      realloc(interact_file_lines, sizeof(interact_data_line_t *) * (num_interact_lines));
	  }
	  interact_file_lines[line_num++] = getDataLine(tag); // add on this one
	  // check
#ifdef USE_STD_MODS
	  if(strcmp(interact_file_lines[line_num-1]->sequence, tag->getAttributeValue("peptide"))) {
	    cerr << "Error: " << interact_file_lines[line_num-1]->sequence << " does not match " << tag->getAttributeValue("peptide") << endl;
#else
	  if(strcmp(interact_file_lines[line_num-1]->sequence, tag->getAttributeValue("stripped_peptide"))) {
	    cerr << "Error: " << interact_file_lines[line_num-1]->sequence << " does not match " << tag->getAttributeValue("stripped_peptide") << endl;
#endif

	    tag->write(cerr);
	    exit(1);
	  }

	}
	else if(use_default_min_ntt_ && tag->isStart() && ! strcmp(tag->getName(), "enzymatic_search_constraint") &&
		! strcasecmp(tag->getAttributeValue("enzyme"), (*enzymes_)[enzyme_index_])) {
	  min_num_tol_term_ = atoi(tag->getAttributeValue("min_number_termini"));

	}
	else if ( ! strcmp(tag->getName(), "mod_aminoacid_mass") && tag->getAttributeValue("alt_aa")) {
	  char* alt_aa = new char[strlen(tag->getAttributeValue("alt_aa"))+1];
	  strcpy(alt_aa, tag->getAttributeValue("alt_aa"));
	  int pos = atoi(tag->getAttributeValue("position"));
	    
	  if (!interact_file_lines[line_num-1]->alt_sequence) {
	    interact_file_lines[line_num-1]->alt_sequence  = new char[strlen(interact_file_lines[line_num-1]->sequence)+1];
	    strcpy(interact_file_lines[line_num-1]->alt_sequence, interact_file_lines[line_num-1]->sequence_I2L);
	    interact_file_lines[line_num-1]->alt_sequence[pos-1] = *alt_aa;
	  }
	  delete[] alt_aa;
	}
	else if(tag->isStart() && ! strcmp(tag->getName(), "sample_enzyme")) {
	  char* nextenz = new char[strlen(tag->getAttributeValue("name"))+1];
	  strcpy(nextenz, tag->getAttributeValue("name"));
	  int enzIdx = enzymes_->findByStringValue(nextenz);
	  if (enzIdx == -1) {
	    //not an enzyme already stored
	    enzymes_->insertAtEnd(nextenz);
	    enzyme_digestions_->insertAtEnd(new ProteolyticEnzyme(tag));
	    enzyme_index_ = enzyme_digestions_->length()-1;
	    enzyme_data = True;
	    if(use_default_min_ntt_ && min_num_tol_term_ < 0)
	      min_num_tol_term_ = 0;
	  }
	  else {
	    //already listed
	    delete [] nextenz;
	    enzyme_index_ = enzIdx;
	    if(use_default_min_ntt_ && min_num_tol_term_ < 0)
	      min_num_tol_term_ = 0;
	  }
	}
	  
	if(enzyme_data) {
	  if(tag->isEnd() && ! strcmp(tag->getName(), "sample_enzyme")) {
	    (*enzyme_digestions_)[enzyme_index_]->fixSpecificity();
	    enzyme_data = False;
	  }
	  else if(strcmp(tag->getName(), "sample_enzyme")) {
	    (*enzyme_digestions_)[enzyme_index_]->enterSpecificity(tag);
	  }
	    
	}

	delete tag;
      } // if not null

      data = strchr(data+1, '<');
    } // next tag

  } // next line
  fin.close();

  num_interact_lines = line_num; 
  interact_file_lines = (interact_data_line_t**)realloc(interact_file_lines, sizeof(interact_data_line_t*)* num_interact_lines);

  // now do the db search

  uniq_lines = num_interact_lines;
  uniq_interact_data = build_uniq_kwlist(interact_file_lines,&uniq_lines);

#ifndef __LGPL__
  if ( !(kwset = kwsalloc((char *) 0)) )
    printf("Error initializing the keyword set!\n");
 
  fprintf(stderr,"  - Building Commentz-Walter keyword tree...");
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
	  printf("Error adding %s to the keyword structure!\n",
		 interact_line->sequence_I2L);

	alt_index_map->insertAtEnd(i);

	tot_added++;
      }
    }

    if ( interact_line->alt_sequence != NULL ) {
      if ( kwsincr(kwset, 
		   interact_line->alt_sequence,
		   strlen(interact_line->alt_sequence)) != 0 ) 
	printf("Error adding %s to the keyword structure!\n",
	       interact_line->alt_sequence);
      
      alt_index_map->insertAtEnd(i);
      tot_added++;
    }
  }

  /* Prep the keyword trie for searching */
  if ( kwsprep (kwset) != 0 )
    printf("Error prepping the structure for a search!\n");
#else
  /* SPARE Parts Section */
  kwset_t kwset;
  kwmap_t kwmap; // preserves order of uniq_lines
  int n_kwsets;
  int max_kwset_size = min(5000,uniq_lines); //limit to avoid trie width explosion, also may adjust downward dynamically
  std::vector<STRINGFINDER_pep *> finder_pep;
  std::vector<STRINGFINDER_std *> finder_std;
  fprintf(stdout,"  - Building Commentz-Walter keyword tree for %ld unique peptides...",uniq_lines);
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
	  cout << "warning: unexpected character in peptide sequence \"" << badStr << "\", using more general but less efficient search " << endl;
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
  fprintf(stdout,"\n  - Searching the tree...");
  if ( ! (fp = fopen(database_,"r"))){
    std::string rdb(resolve_root(database_));
    const char *trythis = rdb.c_str();
    if (strcmp(trythis, database_)) {
      int err = errno;
      if (! (fp = fopen(trythis,"r"))) {
	printf("error: could not open database file \"%s\" (%s)\n", database_, strerror(err));
	printf("also tried \"%s\" (%s)\n", trythis, strerror(errno));
      } else {
	printf("opening \"%s\" as \"%s\"\n",database_,trythis);
      }
    } else {
      printf("error: could not open database file \"%s\" (%s)\n", database_, strerror(errno));
    }
    if (!fp) {
      exit(-1);
    }
  }

  int iLenSeq;

  char *szBuf_ = (char*) malloc(sizeof(char)*line_width_);
  char *szSeq_ = (char*) malloc(sizeof(char)*MAX_SEQ);
  char *szHdr_ = (char*) malloc(sizeof(char)*MAX_HEADER_LEN); 
  char *szOutputDb_ = (char*) malloc(sizeof(char)*SIZE_FILE);
  int* sort_idx = NULL;
  int max_num_hits = 0;

  szBuf_[0] = '\0';
  szSeq_[0] = '\0';
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

	  while ( (cAA=fgetc(fp)) ) {
	    if (isalpha(cAA) || cAA=='*') {
	      if (cAA == 'I')
		cAA = 'L';
	      if (cAA == '*')
		cAA = '-';

	      szSeq_[iLenSeq++]=cAA;
	      if (iLenSeq > MAX_SEQ) {
		printf(" Error - sequence larger than MAX_SEQ\n\n");
		printf("%s\n\n", szHdr_);
		fclose(fp);
		exit(1);
	      }
	    }
	    else if (feof(fp) || cAA=='>') {
	      int iReturn;
	      iReturn=ungetc(cAA, fp);

	      if (iReturn!=cAA) {
		printf("Error with ungetc.\n\n");
		fclose(fp);
		exit(1);
	      }
	      break;
	    }
	  }

	  szSeq_[iLenSeq]='\0';
	 
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
#ifndef __LGPL__
	    free(match);
#endif
	    qsort(sort_idx, num_hits, sizeof(int), compare_indices);
       
	    // Process the hits for this protein sequence
	    prev_index = -1;
	   
	    for ( i=0; i<num_hits; i++) {
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
		  printf("Error: Could not allocate memory for db_ref_t!\n");

		// use space or | as delimiter of accession if IPI database
		char cDelimiter;

		if (!strncmp(szHdr_, ">IPI", 4) && strchr(szHdr_, '|') != NULL && strchr(szHdr_, '|')<strchr(szHdr_,' ') )
		  cDelimiter = '|';
		else
		  cDelimiter = ' ';

		if ( strchr(szHdr_,cDelimiter) != NULL ) {
		  new_db_ref->alias = make_substr(szHdr_,
						  0,
						  strchr(szHdr_,cDelimiter) - szHdr_ - 1);

		  // get rid of windows control characters (^M, etc) by truncating at the first one
		  size_t len = strlen(new_db_ref->alias);
		  for (size_t i = 0; i < len; i++) {
		    if (new_db_ref->alias[i] < 20 || new_db_ref->alias[i] > 126) {
		      new_db_ref->alias[i] = 0;
		      break;
		    }
		  }

		  new_db_ref->description = make_substr(szHdr_,
							strchr(szHdr_,cDelimiter) - szHdr_,
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
		    else if(new_db_ref->description[z] == '>')
		      new_db_ref->description[z] = ' ';
		    else if(new_db_ref->description[z] == '<')
		      new_db_ref->description[z] = ' ';
		    else if(new_db_ref->description[z] == '/')
		      new_db_ref->description[z] = '\\';
		}
		else {
		  new_db_ref->alias = strdup(szHdr_);

		  // get rid of windows control characters (^M, etc) by truncating at the first one
		  size_t len = strlen(new_db_ref->alias);
		  for (size_t i = 0; i < len; i++) {
		    if (new_db_ref->alias[i] < 20 || new_db_ref->alias[i] > 126) {
		      new_db_ref->alias[i] = 0;
		      break;
		    }
		  }
		  new_db_ref->description = (char *)malloc((2) * sizeof( char )); // allocate as done elsewhere
		  new_db_ref->description[0] = 0;
		}

		ref_list = ((interact_data_line_t *)uniq_interact_data[this_index])->updated_refs;
		ntt_list = ((interact_data_line_t *)uniq_interact_data[this_index])->updated_ntts;
		if(calc_prot_wt_)
		  prot_wt_list = ((interact_data_line_t *)uniq_interact_data[this_index])->updated_prot_wts;

		if ( uniq_interact_data[this_index]->updated_hits == 0 ){
		  if ( (ref_list = (db_ref_t **)malloc( sizeof(db_ref_t *))) == NULL )
		    printf("Error: Could not allocate memory for a new ref_list!\n");
		  if ( (ntt_list = (int*)malloc( sizeof(int *))) == NULL )
		    printf("Error: Could not allocate memory for a new ntt_list!\n");
		  if (calc_prot_wt_ && (prot_wt_list = (double*)malloc( sizeof(double *))) == NULL )
		    printf("Error: Could not allocate memory for a new prot_wt_list!\n");
		}
		else {
		  if ( (ref_list = (db_ref_t **)realloc(ref_list,
							(((interact_data_line_t *)
							  uniq_interact_data[this_index])->updated_hits + 1) * 
							sizeof(db_ref_t *))) == NULL ) 
		    printf("Error: Could not reallocate memory for ref_list!\n");
		  if ( (ntt_list = (int*)realloc(ntt_list,
						 (((interact_data_line_t *)
						   uniq_interact_data[this_index])->updated_hits + 1) * 
						 sizeof(int *))) == NULL ) 
		    printf("Error: Could not reallocate memory for ntt_list!\n");
		  if (calc_prot_wt_ && (prot_wt_list = (double*)realloc(prot_wt_list,
									(((interact_data_line_t *)
									  uniq_interact_data[this_index])->updated_hits + 1) * 
									sizeof(double *))) == NULL ) 
		    printf("Error: Could not reallocate memory for prot_wt_list!\n");
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
		      prev = (result-1)[0];

		    // change M to - when it's first character in the Protein
		    if(prev == 'M' && reslen >= seqlen - 1) {
		      prev = '-';
		    }

		    size_t nextlen=strlen(((interact_data_line_t *)
					   uniq_interact_data[this_index])->sequence_I2L);

		    if(reslen == nextlen)
		      next = '-';
		    else
		      next = result[nextlen]; // next aa

		    if (enzyme_digestions_->size() == 0) {
		      cerr << "ERROR: <sample_enzyme> tag was not found in the pep.xml file." << endl;
		      exit(1);
		    }

		    int new_ntt = (*enzyme_digestions_)[((interact_data_line_t *)
							 uniq_interact_data[this_index])->enzyme_ind]->
		      getNumTolTerm(prev, ((interact_data_line_t *) uniq_interact_data[this_index])->sequence_I2L, next);

		    if (firstpass || new_ntt > ntt) { // do this just once
		      if (firstpass) {
			new_db_ref->prev_aa = (char *)malloc(n_desired_prev_aas_+1);
			new_db_ref->next_aa = (char *)malloc(n_desired_next_aas_+1);
		      }
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
		    ntt = new_ntt > ntt ? new_ntt : ntt;		   
		    result = strstr(result+strlen(((interact_data_line_t *)uniq_interact_data[this_index])->sequence_I2L), 
				    ((interact_data_line_t *)uniq_interact_data[this_index])->sequence_I2L);
		    firstpass = false;
		  }
		} // if not null
		else {
		  cout << "error: " << ((interact_data_line_t *)uniq_interact_data[this_index])->sequence_I2L << " not found in " <<  szSeq_ << endl; 
		  exit(1);
		}

		// now see if we can reuse existing text
		consolidate_text_ptr(new_db_ref->alias,aliases);
		consolidate_text_ptr(new_db_ref->description,descriptions);
		consolidate_text_ptr(new_db_ref->prev_aa,prev_aas);
		consolidate_text_ptr(new_db_ref->next_aa,next_aas);

		ntt_list[uniq_interact_data[this_index]->updated_hits] = ntt;

		if(calc_prot_wt_)
		  prot_wt_list[uniq_interact_data[this_index]->updated_hits] = ResidueMass::getProteinMass(szSeq_, monoisotopic);

		uniq_interact_data[this_index]->updated_hits += 1;
		((interact_data_line_t *)uniq_interact_data[this_index])->updated_refs = ref_list;
		((interact_data_line_t *)uniq_interact_data[this_index])->updated_ntts = ntt_list;

		if(calc_prot_wt_)
		  ((interact_data_line_t *)uniq_interact_data[this_index])->updated_prot_wts = prot_wt_list;

		prev_index = this_index;
	      } // end if this_index != prev_index
	    } // end for num_hits
	  } // end if num_hits
	} // while fgets database
    } // while fin.getline
#ifdef __LGPL__
  /*SPARE Parts Section*/
  for (int n=n_kwsets;n--;) {
    if (bUsePeptideAlphabet)
      delete finder_pep[n];
    else
      delete finder_std[n];
  }
#else
  kwsfree(kwset);
#endif
  fclose( fp );

  fprintf(stdout,"\n  - Linking duplicate entries...");
  link_duplicates(interact_file_lines,num_interact_lines);

  incr_time = time(NULL);
  fprintf(stdout,"\n  - Printing results...");

  // second time through data
  Tag* analysis_summary = new Tag("analysis_summary", True, True);
  analysis_summary->setAttributeValue("analysis", getName());
  analysis_summary->setAttributeValue("time", time_);

  Tag* analysis_timestamp_start = new Tag("analysis_timestamp", True, False);
  analysis_timestamp_start->setAttributeValue("analysis", getName());
  analysis_timestamp_start->setAttributeValue("time", time_);
  analysis_timestamp_start->setAttributeValue("id", "1");

  Tag* analysis_timestamp_stop = new Tag("analysis_timestamp", False, True);

  Tag* timestamp = new Tag("database_refresh_timestamp", True, True);
  // get time info
  timestamp->setAttributeValue("database", database_);
  //  if(! use_default_min_ntt_) {
  
  char next[20];
  sprintf(next, "%d", min_num_tol_term_);
  timestamp->setAttributeValue("min_num_enz_term", next);
  
  //}
  // want to overwrite timestamp to summary.xml file......

  int index = 0;
  double nextprob = -1.0;

  // construct a tmpfile name based on xmlfile
  std::string outfile = make_tmpfile_name(xmlfile);
  //cerr << "writing data to " << outfile << endl;
  ofstream fout(outfile.c_str());
  if(! fout) {
    cerr << "cannot write output to file " << outfile << endl;
    exit(1);
  }

  TagFilter* refresh_filter = new TagFilter("alternative_protein");

  TagFilter* timestamp_filter = new TagFilter("analysis_timestamp");
  timestamp_filter->enterRequiredAttributeVal("analysis", getName());

  TagFilter* summary_filter = new TagFilter("analysis_summary");
  summary_filter->enterRequiredAttributeVal("analysis", getName());

  int result_index = 1;
  char search_result[] = "spectrum_query";
  char attr_name[] = "index";

  RACI fin2(xmlfile); // can read gzipped xml
  if(! fin2) {
    cerr << "RefreshParser(2): error opening " << xmlfile << endl;
    exit(1);
  }

  char current_database[500];
  current_database[0] = 0;
  Array<Tag*>* tags = NULL;
  long line_index = 0L;
  Boolean refreshed = False;

  while(fin2.getline(nextline, line_width_)) {   
    data = strchr(nextline, '<');
    while(data != NULL) {
      tag = new Tag(data);
      refreshed = False;

      //tag->write(cout);
      if(line_index > num_interact_lines) {
	cerr << "error1" << endl;
	exit(1);
      }
      if(tag != NULL && ! timestamp_filter->filter(tag) && ! refresh_filter->filter(tag) &&
	 ! summary_filter->filter(tag)) {
	if ( (tag->isStart() && ! strcmp(tag->getName(), "linked_peptide"))  || 
	     (tag->isStart() && ! strcmp(tag->getName(), "search_hit") && 
	      ! strcmp(tag->getAttributeValue("hit_rank"), "1") &&
	      (!tag->getAttributeValue("is_rejected") || 
	       strcmp(tag->getAttributeValue("is_rejected"), "1")) &&
	      refresh(tag)) ) {
#ifdef USE_STD_MODS
	  if(strcmp(interact_file_lines[line_index]->sequence, tag->getAttributeValue("peptide"))) {
	    cerr << "error2 " << interact_file_lines[line_index]->sequence << " vs " << tag->getAttributeValue("peptide") << endl;
#else
	  if(strcmp(interact_file_lines[line_index]->sequence, tag->getAttributeValue("stripped_peptide"))) {
	    cerr << "error2 " << interact_file_lines[line_index]->sequence << " vs " << tag->getAttributeValue("stripped_peptide") << endl;
#endif
	    exit(1);
	  }
	  
	  tags = getRefreshTags(tag,
				interact_file_lines[line_index++],
				current_database);
	  
	  if(tags != NULL) {
	    for(int k = 0; k < tags->length(); k++)
	      if((*tags)[k] != NULL) {
		RECORD((*tags)[k]);
		delete (*tags)[k];
	      }
	    delete tags;
	    tags = NULL;
	    refreshed = True;
	  }
	}
	else if(tag->isStart() && ! strcmp(tag->getName(), "msms_pipeline_analysis")) {
	  tag->write(fout);
	  // here write the timestamp
	  analysis_summary->write(fout);
	}
	else if(tag->isEnd() && ! strcmp(tag->getName(), "search_summary")) {
	  tag->write(fout);
	  // here write the timestamp
	  analysis_timestamp_start->write(fout);
	  timestamp->write(fout);
	  analysis_timestamp_stop->write(fout);
	}
	else if(tag->isStart() && ! strcmp(tag->getName(), "search_database")) {
	  memcpy(current_database, 
		 tag->getAttributeValue("local_path"), 
		 strlen(tag->getAttributeValue("local_path"))+1);
	  RECORD(tag);
	}
	else {
	  tag->write(fout);
	}
      }
      if(! refreshed && tag != NULL)
	delete tag;
      data = strchr(data+1, '<');
    } // next tag
  } // next line
  fin2.close();
  fout.close();

  delete [] nextline;

  if (num_mapped_ == 0)
    fprintf(stdout,"\nERROR: no entries mapped! Please check input files.");
  else
    fprintf(stdout,"\n  - Mapped %d entries", num_mapped_);
  if (num_not_mapped_ > 0)
    fprintf(stdout,"\nWarning: could not map %d entries", num_not_mapped_);

  fprintf(stdout,"\n\n");

  if (testType!=NO_TEST) {
    //
    // regression test stuff - bpratt Insilicos LLC, Nov 2005
    //
    TagListComparator("RefreshParser",testType,test_tags,testFileName);
    delete[] testFileName;
    for(int k = test_tags.length(); k--;) {
      delete test_tags[k];
    }
  }

  if(! overwrite(xmlfile, outfile.c_str(), "</msms_pipeline_analysis>")) {
    cerr << "error: no RefreshParser data written to file " << xmlfile << endl;
  }


  // clean up...
  map <int*, int> freed_iptrs_map;
  map <char*, int> *freed_cptrs_map = new map <char*, int>;
  map <db_ref_t*, int> *freed_pptrs_map = new map <db_ref_t*, int>;
  map <db_ref_t**, int> *freed_rptrs_map = new map <db_ref_t**, int>;

  if (uniq_interact_data != NULL) {
    for(int k = 0; k < uniq_lines; k++) {
      if (uniq_interact_data[k] != NULL) {
        if (uniq_interact_data[k]->updated_refs != NULL)
          save_ptrs_to_free(uniq_interact_data[k],freed_rptrs_map,freed_pptrs_map,freed_cptrs_map);
        uniq_interact_data[k]->updated_refs = NULL;
      }
    }
    free(uniq_interact_data);
    uniq_interact_data = NULL;
  }

  if (interact_file_lines != NULL) {
    for(int k = 0; k < num_interact_lines; k++) {
      if (interact_file_lines[k] != NULL) {
	delete[] interact_file_lines[k]->sequence;
	delete[] interact_file_lines[k]->sequence_I2L;
	delete[] interact_file_lines[k]->db_ref;
	if ( interact_file_lines[k]->alt_sequence != NULL )
	  delete[] interact_file_lines[k]->alt_sequence;

	if (interact_file_lines[k]->updated_ntts != NULL)
	  freed_iptrs_map[interact_file_lines[k]->updated_ntts] = 1;
	interact_file_lines[k]->updated_ntts = NULL;

	if (interact_file_lines[k]->updated_refs != NULL)
	  save_ptrs_to_free(interact_file_lines[k],freed_rptrs_map,freed_pptrs_map,freed_cptrs_map);
	interact_file_lines[k]->updated_refs = NULL;

	free(interact_file_lines[k]);
	interact_file_lines[k] = NULL;
      }
    }
    free(interact_file_lines);
    interact_file_lines = NULL;
  }

  for (map <int*, int>::const_iterator it = freed_iptrs_map.begin(); it != freed_iptrs_map.end(); it++)
    free(it->first);
  for (map <char*, int>::const_iterator it = (*freed_cptrs_map).begin(); it != (*freed_cptrs_map).end(); it++)
    free(it->first);
  for (map <db_ref_t*, int>::const_iterator it = (*freed_pptrs_map).begin(); it != (*freed_pptrs_map).end(); it++)
    free(it->first);
  for (map <db_ref_t**, int>::const_iterator it = (*freed_rptrs_map).begin(); it != (*freed_rptrs_map).end(); it++)
    free(it->first);
  delete freed_cptrs_map;
  delete freed_pptrs_map;
  delete freed_rptrs_map;

  delete alt_index_map;
  delete analysis_summary;
  delete analysis_timestamp_start;
  delete analysis_timestamp_stop;
  delete timestamp;
  delete refresh_filter;
  delete timestamp_filter;
  delete summary_filter;

  free(sort_idx);
  free(szBuf_);
  free(szSeq_);
  free(szHdr_);
  free(szOutputDb_);

  aliases.clear();
  descriptions.clear();
  prev_aas.clear();
  next_aas.clear();
}


void RefreshParser::setFilter(Tag* tag) {
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

Boolean RefreshParser::refresh(Tag* tag) {
  return (!degen_only_ || strcmp(tag->getAttributeValue("num_tot_proteins"), "1"));
  //  && (!ignore_prefix_ || strncmp(tag->getAttributeValue("protein"), ignore_prefix_, strlen(ignore_prefix_)) != 0));
}

int RefreshParser::getEnzymeIndex(char* enz) {
  return enzymes_->findByStringValue(enz);
}


interact_data_line_t* RefreshParser::getDataLine(Tag* hit_tag) {

  interact_data_line_t* interact_data_line = (interact_data_line_t*)calloc(sizeof(interact_data_line_t),1);

  interact_data_line->hits = 1; // default to 1
  if (hit_tag->getAttributeValue("num_tot_proteins")) {
    interact_data_line->hits = atoi(hit_tag->getAttributeValue("num_tot_proteins"));
  }

#ifdef USE_STD_MODS
  interact_data_line->sequence = new char[strlen(hit_tag->getAttributeValue("peptide"))+1];
  strcpy(interact_data_line->sequence, hit_tag->getAttributeValue("peptide"));
  interact_data_line->sequence_I2L = new char[strlen(hit_tag->getAttributeValue("peptide"))+1];
  strcpy(interact_data_line->sequence_I2L, hit_tag->getAttributeValue("peptide"));
#else
  interact_data_line->sequence = new char[strlen(hit_tag->getAttributeValue("stripped_peptide"))+1];
  strcpy(interact_data_line->sequence, hit_tag->getAttributeValue("stripped_peptide"));
  interact_data_line->sequence_I2L = new char[strlen(hit_tag->getAttributeValue("stripped_peptide"))+1];
  strcpy(interact_data_line->sequence_I2L, hit_tag->getAttributeValue("stripped_peptide"));
#endif
  if (!hit_tag->getAttributeValue("protein")) {
    fprintf(stderr,"  - ERROR: No 'protein' attribute found for search_hit with peptide %s. Exiting...\n",hit_tag->getAttributeValue("peptide"));
    exit(1);
  }
  interact_data_line->db_ref = new char[strlen(hit_tag->getAttributeValue("protein"))+1];
  strcpy(interact_data_line->db_ref, hit_tag->getAttributeValue("protein"));
  interact_data_line->enzyme_ind = enzyme_index_;

  for (int i =0; i < (int)strlen(interact_data_line->sequence_I2L); i++) {
    if (interact_data_line->sequence_I2L[i] == 'I') {
      interact_data_line->sequence_I2L[i] = 'L';
    }
  }

  interact_data_line->alt_sequence = NULL;

  //DDS: Changed
  //  interact_data_line->min_num_enz_term = min_num_tol_term_;
  //HENRY: Check the presence of num_tol_term before doing this -- otherwise
  // it seg-faults for some OMSSA files.
  //DDS: Commented out to make consistent with search constraint parameters
  //  if (hit_tag->getAttributeValue("num_tol_term") && use_default_min_ntt_) {
  //  interact_data_line->min_num_enz_term = atoi( hit_tag->getAttributeValue("num_tol_term") );
  //} else {
  interact_data_line->min_num_enz_term = min_num_tol_term_;
  //}

  //interact_data_line->updated_hits = 0L;
  //interact_data_line->updated_refs = 0L;
  //interact_data_line->updated_ntts = NULL;

  return interact_data_line;
}

// char * make_substr(char * string, long start, long end)
//
// Create a new string from the substr of an existing
// string.  This routine allocates memory for the new
// string which can be released using free().
char* RefreshParser::make_substr(char * string, long start, long end) 
{
  char * substr;
  /*
    int verbose = 0;
    if ( verbose > 1 ) 
    printf("make_substr: called with s=%s start=%ld end=%ld\n",
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


char * RefreshParser::replace_substr(char * string, char * orig_substr, char * repl_substr ) 
{
  char * substr_ptr;
  char * suffix_ptr;
  char * search_ptr;
  size_t substr_len_diff = 0L;
  long offset = 0L;
  int verbose = 0;

  if ( verbose > 1 ) 
    printf("replace_substr: called with s=%s o=%s r=%s\n",
	   string, orig_substr, repl_substr);

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


interact_data_line_t ** RefreshParser::build_uniq_kwlist(interact_data_line_t ** interact_data,long * interact_lines)
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
    printf("build_uniq_kwlist: called\n");

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


int RefreshParser::link_duplicates(interact_data_line_t ** interact_data,long interact_lines)
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
  int* ntt_ptr=NULL;

  double* prot_wt_ptr = NULL;

  if ( verbose > 1 ) 
    printf("link_duplicates: called\n");

  //printf("Allocating memory for sort\n");
  sorted_lines = (interact_data_line_t **)malloc(sizeof(interact_data_line_t *)*interact_lines);

  //printf("Creating array for sort\n");
  for ( i = 0; i < interact_lines; i++ )
    sorted_lines[i] = (interact_data_line_t *)interact_data[i];

  //printf("sorting...\n");
  qsort(sorted_lines,interact_lines,sizeof(interact_data_line_t *),
	compare_interact_lines);

  //printf("Printing sorted data\n");
  //for ( i = 0; i < interact_lines; i++ ) 
  // printf("sequence = %s  pointer = %d\n", sorted_lines[i]->sequence,
  //        sorted_lines[i]->updated_refs);

  sequence = "start";

  for ( i = 0; i < interact_lines; i++ ) {
    if ( sorted_lines[i]->sequence != NULL ) {   // must be a data line
      if ( sorted_lines[i]->updated_refs != NULL ) {
	ref_ptr = sorted_lines[i]->updated_refs;
	ntt_ptr = sorted_lines[i]->updated_ntts;

	if(calc_prot_wt_)
	  prot_wt_ptr = sorted_lines[i]->updated_prot_wts;

	ref_hits = sorted_lines[i]->updated_hits;

	sequence = sorted_lines[i]->sequence_I2L;

	alt_sequence = sorted_lines[i]->alt_sequence;

	enzyme_ind = sorted_lines[i]->enzyme_ind;

      }
      else {
	if ( strcmp(sequence, sorted_lines[i]->sequence_I2L) == 0 && 
	     ((!alt_sequence && !sorted_lines[i]->alt_sequence) ||
	      (alt_sequence && sorted_lines[i]->alt_sequence && strcmp(alt_sequence, sorted_lines[i]->alt_sequence))==0 )) {

	  //printf("updating:\n%s\n",sorted_lines[i]->line);
	  sorted_lines[i]->updated_refs = ref_ptr;
	  sorted_lines[i]->updated_ntts = ntt_ptr;

	  if (sorted_lines[i]->enzyme_ind != enzyme_ind) {
	    for (int t=0; t < ((interact_data_line_t *) sorted_lines[i])->updated_hits; t++) {

	      ((interact_data_line_t *) sorted_lines[i])->updated_ntts[t] = 
		(*enzyme_digestions_)[((interact_data_line_t *)sorted_lines[i])->enzyme_ind]->
		getNumTolTerm(*((interact_data_line_t *) sorted_lines[i])->updated_refs[t]->prev_aa, 
			      ((interact_data_line_t *) sorted_lines[i])->sequence, 
			      *((interact_data_line_t *) sorted_lines[i])->updated_refs[t]->next_aa);

	    }
	  }
	  if(calc_prot_wt_)
	    sorted_lines[i]->updated_prot_wts = prot_wt_ptr;

	  sorted_lines[i]->updated_hits = ref_hits;
	  //printf("now:\n%s\n",sorted_lines[i]->line);
	}
      }
    }
  }

  free(sorted_lines);
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


Array<Tag*>* RefreshParser::getRefreshTags(Tag* search_hit, interact_data_line_t* interact_line, char* database) {
  int i = 0;
  //int j = 0;
  db_ref_t ** db_ref;
  int found_hit = 0;
  const int text_len = 500;
  char text[text_len];
  int verbose = 0;

  Array<Tag*>* tags = new Array<Tag*>;
  Tag* tag = NULL;
 
  int known_proteins = 0;

  // HENRY -- re-calculate num_missed_cleavages regardless
  int new_nmc = (*enzyme_digestions_)[interact_line->enzyme_ind]->getNumMissedCleavages(interact_line->sequence, NULL);
  sprintf(text, "%d", new_nmc);
  search_hit->setAttributeValue("num_missed_cleavages", text);

  if(interact_line->sequence != NULL) {
    //cout << "seq: " << interact_line->sequence << endl;
    if(interact_line->updated_refs != NULL) {
      //cout << "refs not null, with " << interact_line->hits << " hits and " << interact_line->updated_hits << " updated hits" << endl;
      db_ref = (db_ref_t**)interact_line->updated_refs;
      if(interact_line->hits != (interact_line->updated_hits-1)) {
	if(interact_line->hits == 0) {
	  // ???
	}
      }
      sprintf(text, "%ld", interact_line->updated_hits);
      search_hit->setAttributeValue("num_tot_proteins", text);
      found_hit = -1;
      //cout << interact_line->db_ref << ": ";
      int j;
      for(j = 0; j < interact_line->updated_hits; j++ ){
	if ( verbose > 4 ) 
	  fprintf(stdout,"print_interact_file: dbref=%s updrefs[%d]=%s\n",
		  interact_line->db_ref, j, ((db_ref_t *)db_ref[j])->alias);

	if ( strcmp(interact_line->db_ref,((db_ref_t *)db_ref[j])->alias) == 0 )
	  found_hit = j;
	else if(interact_line->db_ref[0] && interact_line->db_ref[0] != '>' &&
		((db_ref_t *)db_ref[j])->alias[0] == '>' &&
		! strcmp(interact_line->db_ref,((db_ref_t *)db_ref[j])->alias+1))
	  found_hit = j;

	//cout << ((db_ref_t *)db_ref[j])->alias << " ";
      }
      
      // cout << "found hit: " << found_hit << endl;

      if(found_hit == -1) {
	// find first hit with acceptable ntt
	found_hit = 0;
	while(found_hit < interact_line->updated_hits && interact_line->updated_ntts[found_hit] < interact_line->min_num_enz_term)
	  found_hit++;
	if(found_hit < interact_line->updated_hits) {
	  // set the parent tag protein to the top hit
	  //found_hit = 0;
	  if(((db_ref_t*)db_ref[found_hit])->alias != NULL) {
	    int offset = 0;
	    if(((db_ref_t*)db_ref[found_hit])->alias[0] == '>')
	      offset++;
	    search_hit->setAttributeValue("protein", ((db_ref_t*)db_ref[found_hit])->alias+offset);
	  }
	  //search_hit->setAttributeValue("protein", ((db_ref_t*)db_ref[found_hit])->alias);
	  if(((db_ref_t*)db_ref[found_hit])->description != NULL && ((db_ref_t*)db_ref[found_hit])->description[0]) {
	    search_hit->setAttributeValue("protein_descr", ((db_ref_t*)db_ref[found_hit])->description);
	  }
	  sprintf(text, "%d", interact_line->updated_ntts[found_hit]);
	  search_hit->setAttributeValue("num_tol_term", text);

	  if(calc_prot_wt_) {
	    sprintf(text, "%0.4f", interact_line->updated_prot_wts[found_hit]);
	    search_hit->setAttributeValue("protein_mw", text);
	  }

	  // adjust prev and follow aa list lengths if desired
	  if (db_ref[found_hit]->prev_aa) {
	    search_hit->setAttributeValue("peptide_prev_aa", db_ref[found_hit]->prev_aa);
	  }
	  if (db_ref[found_hit]->next_aa) {
	    search_hit->setAttributeValue("peptide_next_aa", db_ref[found_hit]->next_aa);
	  }
	
	  known_proteins++;
	  tags->insertAtEnd(search_hit);
	}
      }
      else { // set descr
	//cout << "here..." << endl;
	//cout << "descr: " << ((db_ref_t*)db_ref[found_hit])->description << endl;
	int offset = 0;
	if(((db_ref_t*)db_ref[found_hit])->alias[0] && ((db_ref_t*)db_ref[found_hit])->alias[0] == '>')
	  offset++;
	search_hit->setAttributeValue("protein", ((db_ref_t*)db_ref[found_hit])->alias+offset);
	if(((db_ref_t*)db_ref[found_hit])->description != NULL && ((db_ref_t*)db_ref[found_hit])->description[0]) {
	  //cout << "null description..." << endl; exit(1);
	  //printf("des len: %d\n", strlen(((db_ref_t*)db_ref[found_hit])->description));
	  //cout << ((db_ref_t*)db_ref[found_hit])->description << endl;
	  //int offset = 0;
	  //if(strlen(((db_ref_t*)db_ref[found_hit])->alias) > 0 && ((db_ref_t*)db_ref[found_hit])->alias[0] == '>')
	  //    offset++;
	  //search_hit->setAttributeValue("protein", ((db_ref_t*)db_ref[found_hit])->alias+offset);
	  search_hit->setAttributeValue("protein_descr", ((db_ref_t*)db_ref[found_hit])->description);
	  //sprintf(text, "%d", interact_line->updated_ntts[found_hit]);
	  //search_hit->setAttributeValue("num_tol_term", text);
	  //known_proteins++;
	}
	//else 
	// ; //known_proteins++; //cout << "null description..." << endl;
	//cout << "done!" << endl;
	sprintf(text, "%d", interact_line->updated_ntts[found_hit]);
	search_hit->setAttributeValue("num_tol_term", text);

	if(calc_prot_wt_) {
	  sprintf(text, "%0.4f", interact_line->updated_prot_wts[found_hit]);
	  search_hit->setAttributeValue("protein_mw", text);
	}

	// adjust prev and follow aa list lengths if desired
	if (db_ref[found_hit]->prev_aa) {
	  search_hit->setAttributeValue("peptide_prev_aa", db_ref[found_hit]->prev_aa);
	}
	if (db_ref[found_hit]->next_aa) {
	  search_hit->setAttributeValue("peptide_next_aa", db_ref[found_hit]->next_aa);
	}

	known_proteins++;
	tags->insertAtEnd(search_hit);
      }

      // now an alternative protein for each extra
      for(j = 0; j < interact_line->updated_hits; j++ ){
	if(j != found_hit && interact_line->updated_ntts[j] >= interact_line->min_num_enz_term) {
	  tag = new Tag("alternative_protein", True, True);
	  if(((db_ref_t*)db_ref[j])->alias != NULL) {
	    int offset = 0;
	    if(((db_ref_t*)db_ref[j])->alias[0] == '>')
	      offset++;
	    tag->setAttributeValue("protein", ((db_ref_t*)db_ref[j])->alias+offset);
	  }
	  if(((db_ref_t*)db_ref[j])->description != NULL && ((db_ref_t*)db_ref[j])->description[0]) {
	    tag->setAttributeValue("protein_descr", ((db_ref_t*)db_ref[j])->description);
	  }

	  sprintf(text, "%d", interact_line->updated_ntts[j]);
	  tag->setAttributeValue("num_tol_term", text);

	  if(calc_prot_wt_) {
	    sprintf(text, "%0.4f", interact_line->updated_prot_wts[j]);
	    tag->setAttributeValue("protein_mw", text);
	  }
	  // adjust prev and follow aa list lengths if desired
	  if (db_ref[j]->prev_aa) {
	    tag->setAttributeValue("peptide_prev_aa", db_ref[j]->prev_aa);
	  }
	  if (db_ref[j]->next_aa) {
	    tag->setAttributeValue("peptide_next_aa", db_ref[j]->next_aa);
	  }
	  
	  tags->insertAtEnd(tag);
	  known_proteins++;
	} // if not main result
      } // next

    } // if found some matches in the database
    sprintf(text, "%d", known_proteins);
    search_hit->setAttributeValue("num_tot_proteins", text);

    if(tags->length() == 0 && (interact_line->updated_refs == NULL || ! known_proteins)) {      
      // not found in current database: append "_UNMAPPED" to the protein name, unless the original protein has a special prefix (specified by user), in which case nothing will change. 
      if ( (  !search_hit->getAttributeValue("xl_type") || strcmp(search_hit->getAttributeValue("xl_type"), "xl")  )
	  && ( !restore_prefix_ || strncmp(search_hit->getAttributeValue("protein"), restore_prefix_, strlen(restore_prefix_)) != 0) ) {
	std::string newprotein(search_hit->getAttributeValue("protein"));
	if (newprotein.length() < 9 || newprotein.substr(newprotein.length() - 9) != "_UNMAPPED") {
	  newprotein += "_UNMAPPED";
	}
	search_hit->setAttributeValue("protein", newprotein.c_str());
	char* next_descr = new char[strlen(interact_line->db_ref)+strlen(interact_line->db_ref) + strlen(database)+100];
	sprintf(next_descr, "originally identified as %s in database %s", interact_line->db_ref, database);
	search_hit->setAttributeValue("protein_descr", (XMLEscape(next_descr)).c_str());

	++num_not_mapped_;

	if (!(search_hit->getAttributeValue("num_tol_term"))) {
	  int default_ntt = 0;
	  sprintf(text, "%d", default_ntt);  
	  search_hit->setAttributeValue("num_tol_term", text);
	}

	delete [] next_descr;
      }
      tags->insertAtEnd(search_hit);
    }
    else {
      ++num_mapped_;
    }
  }
  else {
    //cout << "returning with orig" << endl;
    tags->insertAtEnd(search_hit);
  }

  return tags;
}
