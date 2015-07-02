/**
 * \file CarpStreamBuf.h
 * $Revision: 1.00 $
 * \brief streambuffer that routes c++ stream to the carp logging system
 *************************************************************************/
#ifndef CARPSTREAMBUF_H_
#define CARPSTREAMBUF_H_

#include <iostream>
#include <cstdio>
#include <string>

class CarpStreamBuf : public std::streambuf {
 private:
  std::string endl_str_; ///< Definition of endline
  std::string buffer_; ///< buffer of characters put into the stream

 protected:
  /*
   * Overrides for std::streambuf
   */

  /**
   * Handles data being put into the stream.  Fills up the buffer_ string
   * until an endline is encountered, then passes the buffer to carp.
   */ 
  virtual int overflow(
    int over_char = EOF ///< character to put into the stream
  );

  /**
   * Not used, just returns EOF
   */
  virtual int underflow();

  /**
   * Not used, just returns EOF
   */
  virtual int uflow();

  /**
   * Not used, just returns EOF
   */
  virtual int pbackfail(
    int pback_char = EOF ///< character to popback
  );

  /**
   * Not used, just returns EOF
   */
  virtual int sync();
  
 public:
  /**
   * Constructor
   */
  CarpStreamBuf();


  virtual ~CarpStreamBuf();
};


#endif

/* 
 * Local Variables: 
 * mode: c 
 * c-basic-offset: 2 
 * End: 
 */ 
 
