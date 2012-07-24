// file      : tests/cxx/parser/validation/built-in/int/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test the built-in integer & friends types validation.
//
#include <limits.h>

#include <string>
#include <sstream>
#include <cassert>

#include <xsd/cxx/parser/validating/exceptions.hxx>
#include <xsd/cxx/parser/validating/xml-schema-pimpl.hxx>

using namespace std;
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

  std::string min;
  std::string max;
  std::string umax;

  {
    ostringstream ostr;
    ostr << LLONG_MIN;
    min = ostr.str ();
  }

  {
    ostringstream ostr;
    ostr << LLONG_MAX;
    max = ostr.str ();
  }

  {
    ostringstream ostr;
    ostr << ULLONG_MAX;
    umax = ostr.str ();
  }

  // integer
  //
  {
    integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (min.c_str ());
    p._post ();
    assert (p.post_integer () == LLONG_MIN);
  }

  {
    integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("0");
    p._post ();
    assert (p.post_integer () == 0);
  }

  {
    integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (max.c_str ());
    p._post ();
    assert (p.post_integer () == LLONG_MAX);
  }

  // negative_integer
  //
  {
    negative_integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (min.c_str ());
    p._post ();
    assert (p.post_negative_integer () == LLONG_MIN);
  }

  {
    negative_integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-1");
    p._post ();
    assert (p.post_negative_integer () == -1);
  }

  // non_positive_integer
  //
  {
    non_positive_integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (min.c_str ());
    p._post ();
    assert (p.post_non_positive_integer () == LLONG_MIN);
  }

  {
    non_positive_integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("+0");
    p._post ();
    assert (p.post_non_positive_integer () == 0);
  }

  // positive_integer
  //
  {
    positive_integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("1");
    p._post ();
    assert (p.post_positive_integer () == 1);
  }

  {
    positive_integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (umax.c_str ());
    p._post ();
    assert (p.post_positive_integer () == ULLONG_MAX);
  }

  // non_negative_integer
  //
  /*
  {
    non_negative_integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-0");
    p._post ();
    assert (p.post_non_negative_integer () == 0);
  }
  */

  {
    non_negative_integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("0");
    p._post ();
    assert (p.post_non_negative_integer () == 0);
  }

  {
    non_negative_integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (umax.c_str ());
    p._post ();
    assert (p.post_non_negative_integer () == ULLONG_MAX);
  }


  // Bad
  //

  std::string past_min (min);
  std::string past_max (max);
  std::string past_umax (umax);

  assert (*past_min.rbegin () != '9');
  assert (*past_max.rbegin () != '9');
  assert (*past_umax.rbegin () != '9');

  (*past_min.rbegin ())++;
  (*past_max.rbegin ())++;
  (*past_umax.rbegin ())++;

  // integer
  //
  {
    integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (past_min.c_str ());
    assert (test_post_fail (p));
  }

  {
    integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (past_max.c_str ());
    assert (test_post_fail (p));
  }

  // negative_integer
  //
  {
    negative_integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (past_min.c_str ());
    assert (test_post_fail (p));
  }

  {
    negative_integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-0");
    assert (test_post_fail (p));
  }

  {
    negative_integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("1");
    assert (test_post_fail (p));
  }

  // non_positive_integer
  //
  {
    non_positive_integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (past_min.c_str ());
    assert (test_post_fail (p));
  }

  {
    non_positive_integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("1");
    assert (test_post_fail (p));
  }

  // positive_integer
  //
  {
    positive_integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-1");
    assert (test_post_fail (p));
  }

  {
    positive_integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("+0");
    assert (test_post_fail (p));
  }

  {
    positive_integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (past_umax.c_str ());
    assert (test_post_fail (p));
  }

  // non_negative_integer
  //
  {
    non_negative_integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-1");
    assert (test_post_fail (p));
  }

  {
    non_negative_integer_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (past_umax.c_str ());
    assert (test_post_fail (p));
  }
}
