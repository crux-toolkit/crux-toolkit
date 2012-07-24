// file      : tests/cxx/tree/encoding/char/lcp/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test local code page encoding (--char-encoding lcp).
// The test just makes sure it still compiles and works.
//

#include <memory> // std::auto_ptr
#include <fstream>
#include <iostream>

#include "test.hxx"

using namespace std;
using namespace test;

int
main (int argc, char* argv[])
{
  if (argc != 2)
  {
    cerr << "usage: " << argv[0] << " test.xml" << endl;
    return 1;
  }

  try
  {
    auto_ptr<type> r (root (argv[1]));

    xml_schema::namespace_infomap map;
    map["t"].name = "test";

    root (std::cout, *r, map);
  }
  catch (xml_schema::exception const& e)
  {
    cerr << "xml_schema::exception: " << e.what () << endl;
    return 1;
  }
}
