/** 
 * \file CarpStreamBuf.cpp 
 * AUTHOR: Sean McIlwain
 * CREATE DATE: 10/30/2012
 * $Revision: 1.00 $ 
 * \brief streambuffer that routes c++ stream to the carp logging system 
 *************************************************************************/ 
#include "CarpStreamBuf.h"
#include "carp.h"
#include <iostream>
#include <sstream>

using namespace std;

/** 
 * Handles data being put into the stream.  Fills up the buffer_ string 
 * until an endline is encountered, then passes the buffer to carp. 
 */  
int CarpStreamBuf::overflow( 
  int over_char ///< character to put into the stream 
) {

  if (over_char == EOF) {
    return EOF;
  } else { 

    buffer_.push_back((char)over_char);

    if (buffer_.find(endl_str_) != string::npos) {
      buffer_.erase(buffer_.length()-endl_str_.length(),endl_str_.length());
      carp(CARP_INFO, "%s", buffer_.c_str());
      buffer_.clear();
    }
    return over_char; 
  }
}

/** 
 * Not used, just returns EOF 
 */ 
int CarpStreamBuf::underflow() {

  return EOF;
}

/** 
 * Not used, just returns EOF 
 */ 
int CarpStreamBuf::uflow() {

  return EOF;
}
 
/** 
 * Not used, just returns EOF 
 */  
int CarpStreamBuf::pbackfail(
  int c
  ) {

  return EOF;
}

/** 
 * Not used, just returns EOF 
 */ 
int CarpStreamBuf::sync() {

  return 0;
}

/**
 * Constructor
 */
CarpStreamBuf::CarpStreamBuf() : streambuf() {

  //An attempt to make endline detection cross-platform.
  ostringstream oss;
  oss << endl;
  endl_str_ = oss.str();
}

/**
 * Destructor
 */
CarpStreamBuf::~CarpStreamBuf() { 

  //If there is any buffer left, then just print it out.
  if (buffer_.length() != 0) {
    carp(CARP_INFO, "%s", buffer_.c_str());
  }
}

/*  
 * Local Variables:  
 * mode: c  
 * c-basic-offset: 2  
 * End:  
 */  


