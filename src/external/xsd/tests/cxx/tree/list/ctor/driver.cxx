// file      : tests/cxx/tree/list/ctor/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test list constructors.
//
#include <string>

#include "test.hxx"

using namespace std;
using namespace test;

int
main ()
{
  // Test ctor()
  //
  {
    string_list sl;

    xml_schema::nmtokens nt;
    xml_schema::idrefs id;
  }

  // Test ctor(size_type, const X&)
  //
  {
    string_list sl (10, "abc");
    size_type st (10, 123);

    xml_schema::nmtokens nt (10, "abc");
    xml_schema::idrefs id (10, "abc");
  }

  // Test ctor(const I& begin, const I& end)
  //
  {
    string_list sl1 (10, "abc");
    string_list sl2 (sl1.begin (), sl1.end ());

    I i1 (10, 123);
    I i2 (i1.begin (), i1.end ());

    xml_schema::nmtokens nt1 (10, "abc");
    xml_schema::nmtokens nt2 (nt1.begin (), nt1.end ());

    xml_schema::idrefs id1 (10, "abc");
    xml_schema::idrefs id2 (id1.begin (), id1.end ());
  }
}
