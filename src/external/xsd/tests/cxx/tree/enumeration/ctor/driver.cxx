// file      : tests/cxx/tree/enumeration/ctor/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test enumeration constructors.
//
#include <string>

#include "test.hxx"

using namespace std;
using namespace test;

int
main ()
{
  // Test ctor(const char*).
  //
  {
    string_enum se ("a");
    type t ("a", 1);
  }

  // Test ctor(const std::string&)
  //
  {
    string const s ("c");
    string_enum se (s);
    type t (s, 3);
  }
}
