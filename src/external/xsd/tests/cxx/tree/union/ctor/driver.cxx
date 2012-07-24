// file      : tests/cxx/tree/union/ctor/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test union constructors.
//
#include <string>

#include "test.hxx"

using namespace std;
using namespace test;

int
main ()
{
  // Test ctor(const std::string&)
  //
  {
    string const s ("123");
    int_string_union u (s);
    type t (s);
  }

  // Test ctor(const char*).
  //
  {
    int_string_union u ("123");
    type t ("123");
  }
}
