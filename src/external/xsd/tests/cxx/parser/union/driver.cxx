// file      : tests/cxx/parser/union/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test xsd:union parsing.
//

#include <string>
#include <iostream>

#include "test-pskel.hxx"

using namespace std;
using namespace test;

struct int_string_union_pimpl: int_string_union_pskel
{
  virtual void
  _characters (const xml_schema::ro_string& s)
  {
    cout << "'" << s << "'" << endl;
  }
};

struct type_pimpl: type_pskel
{
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
    int_string_union_pimpl int_string_union_p;
    type_pimpl type_p;

    type_p.parsers (int_string_union_p);

    xml_schema::document doc_p (type_p, "test", "root");

    type_p.pre ();
    doc_p.parse (argv[1]);
    type_p.post_type ();
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
