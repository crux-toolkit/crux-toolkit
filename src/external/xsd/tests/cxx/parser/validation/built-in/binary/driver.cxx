// file      : tests/cxx/parser/validation/built-in/binary/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test the built-in base64Binary and hexBinary types validation.
//
#include <cassert>

#include <xsd/cxx/parser/validating/exceptions.hxx>
#include <xsd/cxx/parser/validating/xml-schema-pimpl.hxx>

using namespace xsd::cxx::parser::validating;

template <typename T>
bool
test_post_fail (T& p)
{
  try
  {
    p._post ();
  }
  catch (invalid_value<char> const&)
  {
    return true;
  }

  return false;
}

int
main ()
{
  typedef xsd::cxx::parser::buffer buffer;

  // Good.
  //

  // hexBinary
  //
  {
    hex_binary_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \t\n ");
    p._characters ("  ");
    p._post ();
    assert (*p.post_hex_binary () == buffer ());
  }

  {
    hex_binary_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \t\n313");
    p._characters ("23334356162636a6b  ");
    p._post ();
    assert (*p.post_hex_binary () == buffer ("12345abcjk", 10));
  }

  // base64Binary
  //
  {
    base64_binary_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \t\n ");
    p._characters ("MTIzND  ");
    p._characters ("VhYmNqaw = = ");
    p._post ();
    assert (*p.post_base64_binary () == buffer ("12345abcjk", 10));
  }

  {
    base64_binary_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("YQ==");
    p._post ();
    assert (*p.post_base64_binary () == buffer ("a", 1));
  }

  {
    base64_binary_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("YWI=");
    p._post ();
    assert (*p.post_base64_binary () == buffer ("ab", 2));
  }

  {
    base64_binary_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("YWJj");
    p._post ();
    assert (*p.post_base64_binary () == buffer ("abc", 3));
  }

  // Bad
  //

  // hexBinary
  //
  {
    hex_binary_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("313");
    assert (test_post_fail (p));
  }

  {
    hex_binary_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("313233343X6162636a6b");
    assert (test_post_fail (p));
  }

  // base64Binary
  //
  {
    base64_binary_pimpl<char> p;
    p.pre ();
    p._pre ();
    // p._characters ("");
    assert (test_post_fail (p));
  }

  {
    base64_binary_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("YQ");
    assert (test_post_fail (p));
  }

  {
    base64_binary_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("==");
    assert (test_post_fail (p));
  }

  {
    base64_binary_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("MTIzNDVhYmNqaw=A");
    assert (test_post_fail (p));
  }
}
