// file      : tests/cxx/tree/types-only/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test that code generated without parsing and serialization functions
// still compiles.
//

#include <memory> // std::auto_ptr
#include <iostream>

#include "test.hxx"

using namespace std;
using namespace test;

int
main ()
{
  color_enum red (color_enum::red);
  red_blue_enum blue (red_blue_enum::blue);

  long_string_union num ("123");

  string_list list;
  list.push_back ("Hello");
  list.push_back ("World");

  complex_type t ("Hello, World!", "foo", color_enum::red);

  anon a ("Hello, World!");
}
