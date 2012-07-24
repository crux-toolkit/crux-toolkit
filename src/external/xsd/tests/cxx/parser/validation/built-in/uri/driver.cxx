// file      : tests/cxx/parser/validation/built-in/uri/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test the built-in anyURI type validation.
//
#include <cassert>

#include <xsd/cxx/parser/validating/exceptions.hxx>
#include <xsd/cxx/parser/validating/xml-schema-pimpl.hxx>

using namespace xsd::cxx::parser::validating;

int
main ()
{
  // Good.
  //
  {
    uri_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("  ");
    p._post ();
    assert (p.post_uri () == "");
  }

  {
    uri_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("relative");
    p._post ();
    assert (p.post_uri () == "relative");
  }

  {
    uri_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("#id");
    p._post ();
    assert (p.post_uri () == "#id");
  }

  {
    uri_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("http://www.example.com/foo#bar");
    p._post ();
    assert (p.post_uri () == "http://www.example.com/foo#bar");
  }
}
