// file      : tests/cxx/tree/encoding/char/utf-8/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test UTF-8 encoding.
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

    {
      type::a_sequence const& s (r->a ());

      if (s[0] != "abc" ||
          s[1] != "\xD5\x95" ||
          s[2] != "\xEA\xAA\xAA" ||
          s[3] != "\xF2\xAA\xAA\xAA")
      {
        cerr << "invalid encoding" << endl;
        return 1;
      }
    }

    {
      type::b_sequence const& s (r->b ());

      if (s[0] != strenum::abc ||
          s[1] != strenum::a_c ||
          s[2] != strenum::cxx_bc ||
          s[3] != strenum::ab_)
      {
        cerr << "invalid encoding" << endl;
        return 1;
      }
    }

    xml_schema::namespace_infomap map;
    map["t"].name = "test";

    root (std::cout, *r, map, "UCS-4LE");
  }
  catch (xml_schema::exception const& e)
  {
    cerr << "xml_schema::exception: " << e.what () << endl;
    return 1;
  }
}
