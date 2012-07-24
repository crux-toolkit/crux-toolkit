// file      : tests/cxx/tree/encoding/char/iso-8859-1/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test ISO-8859-1 encoding.
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
    try
    {
      root (argv[1]);
      return 1;
    }
    catch (xsd::cxx::xml::iso8859_1_unrepresentable const&)
    {
    }

    xsd::cxx::xml::char_transcoder::unrep_char ('?');
    auto_ptr<type> r (root (argv[1]));

    {
      type::a_sequence const& s (r->a ());

      if (s[0] != "abc" ||
          s[1] != "\xE6" ||
          s[2] != "\xA2\xA3\xA4\xA5" ||
          s[3] != "???")
      {
        cerr << "invalid encoding" << endl;
        return 1;
      }
    }

    {
      type::b_sequence const& s (r->b ());

      if (s[0] != strenum::abc ||
          s[1] != strenum::a_c ||
          s[2] != strenum::cxx__bc)
      {
        cerr << "invalid encoding" << endl;
        return 1;
      }
    }

    xml_schema::namespace_infomap map;
    map["t"].name = "test";

    root (std::cout, *r, map, "ISO-8859-1");
  }
  catch (xml_schema::exception const& e)
  {
    cerr << "xml_schema::exception: " << e.what () << endl;
    return 1;
  }
}
