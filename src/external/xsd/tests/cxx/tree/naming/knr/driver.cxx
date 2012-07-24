// file      : tests/cxx/tree/naming/knr/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test K&R naming style.
//

#include <memory> // std::auto_ptr
#include <sstream>
#include <iostream>

#include <xercesc/util/PlatformUtils.hpp>

#include "test.hxx"

using namespace std;
using namespace test;

int
main ()
{
  xercesc::XMLPlatformUtils::Initialize ();

  try
  {
    // Enum 'value' type.
    //
    {
      gender::value v;
      v = gender::female;
    }

    // Anonymous type.
    //
    {
      foo f ("a", "b");

      if (f.a () != "a" || f.b () != "b")
        return 1;
    }

    // Type name and accessors/modifiers.
    //
    {
      type t ("bar");

      // foo
      //
      {
        type::foo_type* p = 0;
        type::foo_optional o;

        if (t.foo ().present ())
          return 1;

        t.foo (o);
      }

      // bar
      //
      {
        type::bar_type* p = 0;

        if (t.bar () != "bar")
          return 1;

        t.bar ("barbar");
      }

      // baz
      //
      {
        type::baz_type* p = 0;
        type::baz_sequence s;
        type::baz_iterator i (s.begin ());
        type::baz_const_iterator ci (s.begin ());

        if (t.baz () != s)
          return 1;

        t.baz (s);
      }

      // any
      //
      {
        type::any_sequence s (t.dom_document ());
        type::any_iterator i (s.begin ());
        type::any_const_iterator ci (s.begin ());

        if (t.any () != s)
          return 1;

        t.any (s);
      }

      // foo
      //
      {
        type::fox_type x = type::fox_default_value ();

        if (t.fox () != x)
          return 1;

        t.fox ("fox");
      }

      // any_attribute
      //
      {
        type::any_attribute_set s (t.dom_document ());
        type::any_attribute_iterator i (s.begin ());
        type::any_attribute_const_iterator ci (s.begin ());

        if (t.any_attribute () != s)
          return 1;

        t.any_attribute (s);
      }
    }

    // Parsing/serialization functions.
    //
    {
      istringstream is ("<t:root xmlns:t='test'>foo</t:root>");
      root (is, xml_schema::flags::dont_validate);
    }

    {
      ostringstream os;
      xml_schema::namespace_infomap m;
      m["t"].name = "test";

      root (os, "foo", m);
    }
  }
  catch (xml_schema::exception const& e)
  {
    cerr << e << endl;
    return 1;
  }

  xercesc::XMLPlatformUtils::Terminate ();
}
