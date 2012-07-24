// file      : tests/cxx/parser/validation/built-in/int/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test the built-in int and unsigned int types validation.
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
  // Good.
  //
  {
    int_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-2147483648");
    p._post ();
    assert (p.post_int () == -2147483648);
  }

  {
    int_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("0");
    p._post ();
    assert (p.post_int () == 0);
  }

  {
    int_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2147483647");
    p._post ();
    assert (p.post_int () == 2147483647);
  }

  {
    unsigned_int_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("0");
    p._post ();
    assert (p.post_unsigned_int () == 0);
  }

  {
    unsigned_int_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("4294967295");
    p._post ();
    assert (p.post_unsigned_int () == 4294967295);
  }

  // Bad
  //

  {
    unsigned_int_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-123");
    assert (test_post_fail (p));
  }


  // Ranges
  //
  {
    int_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-2147483649");
    assert (test_post_fail (p));
  }

  {
    int_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2147483648");
    assert (test_post_fail (p));
  }

  {
    unsigned_int_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("4294967296");
    assert (test_post_fail (p));
  }
}
