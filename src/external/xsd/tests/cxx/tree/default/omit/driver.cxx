// file      : tests/cxx/tree/default/omit/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test default attribute omission from the output.
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
    auto_ptr<type> r (root (argv[1], xml_schema::flags::dont_validate));

    cout << *r << endl
         << "default x: " << derived::x_default_value () << endl
         << "default y: " << derived::y_default_value () << endl
	 << "fixed p: " << derived::p_default_value () << endl
	 << "fixed q1: " << derived::q1_default_value () << endl
	 << "fixed q2: " << derived::q2_default_value () << endl;

    // Serialize.
    //
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
