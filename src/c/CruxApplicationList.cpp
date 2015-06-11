/**
 * \file CruxApplicationList.cpp
 * AUTHOR: Sean McIlwain
 * CREATE DATE: December 6th, 2010
 * \brief Maintains a list of executable applications
 *****************************************************************************/
#include "CruxApplicationList.h"
#include "carp.h"
#include "crux-utils.h"

#include <iostream>

using namespace std;

/**
 * Creates an application list with a listname
 */
CruxApplicationList::CruxApplicationList(
  const char* list_name ///<name of the list
) {
  list_name_ = string(list_name);
}


/**
 * Destructor for CruxApplicationList
 */
CruxApplicationList::~CruxApplicationList() {
  
  for (unsigned int idx=0;idx < applications_.size();idx++) {
    delete applications_.at(idx);
  }
  applications_.clear();

}

/**
 * Adds an application pointer to the list of applications
 */
void CruxApplicationList::add(
  CruxApplication* application ///< application to add
  ) {

  if (find(application->getName()) != NULL) {
    carp(CARP_FATAL, "Name clash! %s",application->getName().c_str());
  }

  applications_.push_back(application);
}

/**
 * Adds a message to be printed in the usage statement
 */
void CruxApplicationList::addMessage(
  const string& message ///< message to add
) {
  messages_[applications_.size()] = message;
}

/**
 * \returns an application by a name,
 * returns NULL if not found
 */
CruxApplication* CruxApplicationList::find(
  const string& appname ///<name of the application to find
  ) {
  return find(appname.c_str());
}

/**
 * \returns an application by a name,
 * returns NULL if not found
 */
CruxApplication* CruxApplicationList::find(
  const char* appname ///<name of the application to find
  ) {

  CruxApplication* crux_application = NULL;

  for (vector<CruxApplication*>::iterator app_iter = 
    applications_.begin();
    app_iter != applications_.end();
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

  vector<CruxApplication*>::iterator iter;

  size_t max_name_length=0;

  for (iter = applications_.begin();
    iter != applications_.end();
    ++iter) {

    max_name_length = max(max_name_length, (*iter) -> getName().length());

  }

  cerr << " Usage: " << list_name_ << " <command> [options] <argument>" << endl;

  for (unsigned int i = 0; i < applications_.size(); ++i) {
    map<int, string>::const_iterator msg = messages_.find(i);
    if (msg != messages_.end()) {
      cerr << endl << msg->second << endl << endl;
    }
    if ((applications_[i])->hidden()) { // skip deprecated commands
      continue;
    }
    string name = applications_[i]->getName();
    string description = applications_[i]->getDescription();
    int padding = max_name_length - name.length();

    cerr<<"  "<<name<<" ";
    for (int idx=0;idx<padding;idx++) {
      cerr<<" ";
    }

    const int LINE_WIDTH = 80;
    unsigned int max_descr_line = LINE_WIDTH - (max_name_length+4);


    // If the description is short enough, just print.
    if (description.length() < max_descr_line) {
      cerr<<" "<<description<<endl;
    }

    // Otherwise, insert EOLs.
    else {
      
      vector<string> words;
      tokenize(description, words, ' ');

      unsigned int word_index = 0;
      unsigned int line_length = 0;

      // Print the first line.
      while(line_length + words[word_index].length() + 1 < max_descr_line) {
        cerr << " ";
        cerr << words[word_index];
        line_length += words[word_index].length() + 1;
        word_index++;
      }

      // Print subsequent lines.
      while (word_index < words.size()) {
        line_length = 0;
        cerr << endl;
        for (unsigned int idx =0;idx < max_name_length + 3; idx++) {
          cerr<<" ";
        }
        while ((word_index < words.size()) && 
               (line_length + words[word_index].length() + 1 < max_descr_line)) {
          cerr << " ";
          cerr << words[word_index];
          line_length += words[word_index].length() + 1;
          word_index++;
        }
      }
      cerr<<endl;
    }
  }
  cerr << endl << endl;
  cerr << "Options and arguments are specific to each command."<< endl;
  cerr << "Type '"<< list_name_ <<" <command>' for details."<<endl;

}

/**
 * Gets the name of the list
 */
string CruxApplicationList::getListName() const {
  return list_name_;
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

  int ret = crux_application->main(argc-1, argv+1);

  carp(CARP_INFO, "Elapsed time: %.3g s", wall_clock() / 1e6);
  carp(CARP_INFO, "Finished crux %s.", appname.c_str());
  carp(CARP_INFO, "Return Code:%i", ret);

  return ret;
}
