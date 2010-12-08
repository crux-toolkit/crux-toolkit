/**
 * \file CruxApplicationList.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE: December 6th, 2010
 * \brief Maintains a list of executable applications
 *****************************************************************************/
#include "CruxApplicationList.h"
#include "DelimitedFile.h"
#include "carp.h"

#include <iostream>

using namespace std;

/**
 * Creates an application list with a listname
 */
CruxApplicationList::CruxApplicationList(const char* list_name) {
  list_name_ = string(list_name);
}


/**
 * Destructor for CruxApplicationList
 */
CruxApplicationList::~CruxApplicationList() {
  
  for (unsigned int idx=0;idx < size();idx++) {
    delete at(idx);
  }
  clear();

}

/**
 * Adds an application pointer to the list of applications
 */
void CruxApplicationList::add(CruxApplication* application) {

  if (find(application->getName()) != NULL) {
    carp(CARP_FATAL, "Name clash! %s",application->getName().c_str());
  }

  push_back(application);
}

/**
 * \returns an application by a name,
 * returns NULL if not found
 */
CruxApplication* CruxApplicationList::find(const string& appname) {
  return find(appname.c_str());
}

/**
 * \returns an application by a name,
 * returns NULL if not found
 */
CruxApplication* CruxApplicationList::find(const char* appname) {

  CruxApplication* crux_application = NULL;

  for (CruxApplicationList::iterator app_iter = this-> begin();
    app_iter != this->end();
    ++app_iter) {
    if (appname == (*app_iter) -> getName()) {
      crux_application = *app_iter;
      break;
    }
  }
  return crux_application;
}

/**
 * prints out the usage statement for this application list.
 * Each applications name is printed along with its description
 */
void CruxApplicationList::usage() {

  CruxApplicationList::iterator iter;

  size_t max_name_length=0;

  for (iter = this->begin();
    iter != this->end();
    ++iter) {

    max_name_length = max(max_name_length, (*iter) -> getName().length());

  }

  cerr <<"max:"<<max_name_length<<endl;


  cerr <<" Usage: " << list_name_ << " <command> [options] <argument>" << endl;
  cerr << endl;
  cerr << list_name_ << " supports the following commands:"<<endl;

  for (CruxApplicationList::iterator iter = this->begin();
    iter != this->end();
    ++iter) {
  
    string name = (*iter)->getName();
    string description = (*iter)->getDescription();
    int name_length = name.length();

    int padding = max_name_length - name_length;

    cerr<<"  "<<name<<"  ";
    for (int idx=0;idx<padding;idx++) {
      cerr<<" ";
    }

    unsigned int max_descr_line = 80-max_name_length-4;

    if (description.length() < max_descr_line) {
      cerr<<description<<endl;
    } else {
      
      vector<string> words;
      DelimitedFile::tokenize(description, words, ' ');

      unsigned int word_index = 0;
      unsigned int line_length = 0;

      //print the first line.
      while(line_length < max_descr_line) {
        cerr << words[word_index];
        cerr << " ";
        line_length += words[word_index].length() + 1;
        word_index++;
      }

      while (word_index < words.size()) {
        line_length = 0;
        cerr << endl;
        for (unsigned int idx =0;idx < max_name_length + 4; idx++) {
          cerr<<" ";
        }
        while ((word_index < words.size()) && (line_length < max_descr_line)) {
          cerr << words[word_index];
          cerr << " ";
          line_length += words[word_index].length() + 1;
          word_index++;
        }
      }
      cerr<<endl;
    }
  }
  cerr << endl;
  cerr << "Options and arguments are specific to each command.";
  cerr << "Type '"<< list_name_ <<" <command>' for details."<<endl;

}

/**
 * the main method for CruxApplicationList.  Attempts to find
 * an application by name from the first argument.  If successful,
 * calls that applications main method with the rest of the parameters.
 */
int CruxApplicationList::main(int argc, char** argv) {

  if (argc < 2) {
    usage();
    return -1;
  }

  string appname = string(argv[1]);
  CruxApplication* crux_application = find(appname);

  if (crux_application == NULL) {
    cerr<< "Cannot find "<<appname<<" in available applications"<<endl;
    usage();
    return -1;
  }

  return (crux_application->main(argc-1, argv+1));
}
