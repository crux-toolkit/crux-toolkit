// file      : tests/cxx/parser/validation/built-in/long/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test the built-in long and unsigned long types validation.
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
    long_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-9223372036854775808");
    p._post ();
    assert (p.post_long () == (-9223372036854775807LL - 1));
  }

  {
    long_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("0");
    p._post ();
    assert (p.post_long () == 0);
  }

  {
    long_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("9223372036854775807");
    p._post ();
    assert (p.post_long () == 9223372036854775807LL);
  }

  {
    unsigned_long_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("0");
    p._post ();
    assert (p.post_unsigned_long () == 0);
  }

  {
    unsigned_long_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("18446744073709551615");
    p._post ();
    assert (p.post_unsigned_long () == 18446744073709551615ULL);
  }

  // Bad
  //

  {
    unsigned_long_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-123");
    assert (test_post_fail (p));
  }


  // Ranges
  //
  {
    long_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-9223372036854775809");
    assert (test_post_fail (p));
  }

  {
    long_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("9223372036854775808");
    assert (test_post_fail (p));
  }

  {
    unsigned_long_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("18446744073709551616");
    assert (test_post_fail (p));
  }
}
