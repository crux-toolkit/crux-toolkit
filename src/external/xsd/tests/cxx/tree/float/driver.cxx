// file      : tests/cxx/tree/float/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test floating point (xsd:{float, double, decimal}) type parsing
// and serialization.
//

#include <memory> // std::auto_ptr
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

    r->simple ().push_back (12.129456);
    r->simple ().push_back (123.129456);
    r->simple ().push_back (1234.129456);

    r->s (12.129456);

    r->complex ().push_back (12.129456);
    r->complex ().push_back (123.129456);
    r->complex ().push_back (1234.129456);
    r->complex ().push_back (-12.12);
    r->complex ().push_back (-123.12);

    r->s (12.129456);

    xml_schema::namespace_infomap map;

    map["t"].name = "test";
    root (cout, *r, map);
  }
  catch (xml_schema::exception const& e)
  {
    cerr << e << endl;
    return 1;
  }
}
