// file      : tests/cxx/parser/name-clash/inheritance/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test for name clashes across inheritance hierarchy.
//

#include <string>
#include <iostream>

#include "test-pskel.hxx"

using namespace std;
using namespace test;

struct derived_pimpl: derived_pskel
{
  virtual void
  e (string const& v)
  {
    cout << "e: " << v << endl;
  }

  virtual void
  e1 (string const& v)
  {
    cout << "e1: " << v << endl;
  }
};

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
    xml_schema::string_pimpl string_p;
    derived_pimpl derived_p;

    derived_p.parsers (string_p, string_p);

    xml_schema::document doc_p (derived_p, "test", "root");

    derived_p.pre ();
    doc_p.parse (argv[1]);
    derived_p.post_derived ();
  }
  catch (xml_schema::exception const& e)
  {
    cerr << e << endl;
    return 1;
  }
  catch (std::ios_base::failure const&)
  {
    cerr << "io failure" << endl;
    return 1;
  }
}
