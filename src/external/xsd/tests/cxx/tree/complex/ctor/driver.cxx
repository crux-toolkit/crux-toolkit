// file      : tests/cxx/tree/complex/ctor/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test generation of varous complex type constructors.
//

#include <cassert>
#include <memory>

#include "test.hxx"

using namespace std;
using namespace test;

int
main ()
{
  // Test case A.
  //
  {
    a_base b1;
    a_base b2 ("abc");        // empty ultimate base + required
    a_base b3 ("abc", "foo"); // ultimate base + required

    a_derived a1;
    a_derived a2 ("foo", "bar");        // empty ultimate base + required
    a_derived a3 (b3, "bar");           // base + required
    a_derived a4 ("abc", "foo", "bar"); // ultimate base + required
  }

  // Test case B.
  //
  {
    b_simple s ("base");
    b_base b ("base", "foo");
    b_derived d ("base", "foo", "bar");
    b_type t ("base");
  }

  // Test case C.
  //
  {
    c_simple s (c_enum::a);
    c_base b (c_enum::a, "foo");
    c_derived d (c_enum::a, "foo", "bar");
    c_type t (c_enum::a);
  }

  // Test case D.
  //
  {
    d_simple s (1);
    d_base b (1, "foo");
    d_derived d (1, "foo", "bar");
    d_type t (1);
  }

  // Test case E.
  //
  {
    // e_base
    //
    e_base b1 (1, "foo", e_complex_type ("bar"));

    auto_ptr<e_complex_type> c2 (new e_complex_type ("bar"));
    e_base b2 (1, "foo", c2);

    auto_ptr<e_simple_type> s3 (new e_simple_type ("foo"));
    auto_ptr<e_complex_type> c3 (new e_complex_type ("bar"));
    e_base b3 (1, s3, c3);

    assert (b1 == b2);
    assert (b1 == b3);

    // e_derived
    //
    e_derived d1 (1, "foo", e_complex_type ("bar"),
                  true, "baz", e_complex_type ("biz"));

    auto_ptr<e_complex_type> c2a (new e_complex_type ("bar"));
    auto_ptr<e_complex_type> c2b (new e_complex_type ("biz"));
    e_derived d2 (1, "foo", c2a, true, "baz", c2b);

    auto_ptr<e_simple_type> s3a (new e_simple_type ("foo"));
    auto_ptr<xml_schema::string> s3b (new xml_schema::string ("baz"));
    auto_ptr<e_complex_type> c3a (new e_complex_type ("bar"));
    auto_ptr<e_complex_type> c3b (new e_complex_type ("biz"));
    e_derived d3 (1, s3a, c3a, true, s3b, c3b);

    assert (d1 == d2);
    assert (d1 == d3);

  }

  // Test case F.
  //
  {
    f_type f1 (xml_schema::type (), 1, "foo", f_complex_type ("bar"));

    auto_ptr<f_complex_type> c2 (new f_complex_type ("bar"));
    f_type f2 (1, "foo", c2);

    auto_ptr<f_simple_type> s3 (new f_simple_type ("foo"));
    auto_ptr<f_complex_type> c3 (new f_complex_type ("bar"));
    f_type f3 (1, s3, c3);

    assert (f1 == f2);
    assert (f1 == f3);
  }
}
