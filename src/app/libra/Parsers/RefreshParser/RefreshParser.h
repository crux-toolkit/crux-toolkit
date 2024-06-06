#ifndef REFRESH_PARSER_H
#define REFRESH_PARSER_H

/*

Program       : RefreshParser
Author        : Andrew Keller <akeller@systemsbiology.org>, Robert Hubley, and open source code
Date          : 11.27.02
SVN Info      : $Id: RefreshParser.h 8223 2020-07-18 08:10:23Z real_procopio $

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


#include <stdio.h>
#include <math.h>
#include <time.h>

#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <fstream>
#ifndef __LGPL__
#include "Common/kwset.h"
#endif

#include "Parsers/Parser/Parser.h"
#include "Parsers/Parser/TagFilter.h"
#include "Common/Enzyme/ProteolyticEnzyme/ProteolyticEnzyme.h"
#include "Common/ResidueMass/ResidueMass.h"
#define MAX_SEQ              50000000    /* should really use realloc but I'm lazy */
#define MAX_HEADER_LEN       1024
//#define SIZE_BUF             4096
#define FIND_PROT_MW 1  // undefine to store protein molecular wts in xml




typedef struct
{
  char * alias;
  char * description;
  char * prev_aa; // for tracking prev_aa for alt proteins
  char * next_aa; // for tracking next_aa for alt proteins
} db_ref_t;

typedef struct
{
  char      * line;
  char      * db_ref;
  db_ref_t ** updated_refs;
  int       * updated_ntts; // for ntt purposes
  char      * sequence;
  char      * sequence_I2L;
  char      * alt_sequence;
  
  long      hits;
  long      updated_hits;
  char      * mtype;
  char      * dbname;
  int       enzyme_ind; // use this to sort by for finding unique peptide/enzyme combos
  int       min_num_enz_term;
  double    * updated_prot_wts; // for ntt purposes

} interact_data_line_t;


int compare_match_records(const void * rec1, const void * rec2);
int compare_interact_lines(const void * line1, const void * line2);
#include "Parsers/Parser/TagListComparator.h" // regression test stuff - bpratt Insilicos LLC, Nov 2005
#include "Common/util.h"

class RefreshParser : public Parser {

 public:

  RefreshParser(const char* xmlfile, const char* database, int min_num_tol_term, 
      int n_desired_prev_aas, int n_desired_next_aas, 
		Boolean degen_only, Boolean prot_wt, const char *testMode, const char *restore_prefix);
  ~RefreshParser();
  void setFilter(Tag* tag);

 protected:

  void parse(const char* xmlfile);
  int getEnzymeIndex(char* enz);
  Boolean refresh(Tag* tag);
  interact_data_line_t* getDataLine(Tag* hit_tag);
  interact_data_line_t ** build_uniq_kwlist(interact_data_line_t ** interact_data,long * interact_lines);
  char * replace_substr(char * string, char * orig_substr, char * repl_substr );
  char* make_substr(char * string, long start, long end);
  int link_duplicates(interact_data_line_t ** interact_data,long interact_lines);
  Array<Tag*>* getRefreshTags(Tag* search_hit, interact_data_line_t* interact_line, char* database);

  char* database_;

  char *testMode_; // regression test stuff - bpratt Insilicos LLC, Nov 2005

  char *restore_prefix_;
  int num_mapped_;
  int num_not_mapped_;

  StringArray* enzymes_;
  Array<ProteolyticEnzyme*>* enzyme_digestions_;
  int enzyme_index_;
  Boolean degen_only_;
  int min_num_tol_term_;
  Boolean use_default_min_ntt_;
  Boolean calc_prot_wt_;
  int n_desired_prev_aas_;
  int n_desired_next_aas_;
};

#endif
