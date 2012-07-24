// file      : tests/cxx/parser/validation/built-in/string/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test the built-in string & friends types validation.
//
#include <string>
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
    p._post_impl (); // List implementation needs this to be post_impl.
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
  typedef xsd::cxx::parser::string_sequence<char> strings;

  // Good.
  //

  // string
  //
  {
    string_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \n\t");
    p._characters (" aaa ");
    p._characters ("bbb");
    p._characters (" ");
    p._post ();
    assert (p.post_string () == " \n\t aaa bbb ");
  }

  // normalized_string
  //
  {
    normalized_string_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \n\t");
    p._characters (" aaa \n\t ");
    p._characters (" bbb");
    p._characters ("  ");
    p._post ();
    assert (p.post_normalized_string () == "    aaa     bbb  ");
  }

  // token
  //
  {
    token_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \n\t");
    p._characters (" aaa \n\t ");
    p._characters (" bbb \n\t");
    p._characters ("  ");
    p._post ();
    assert (p.post_token () == "aaa bbb");
  }

  // name
  //
  {
    name_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \n\t");
    p._characters (" a:b-c_d123 ");
    p._characters ("  ");
    p._post ();
    assert (p.post_name () == "a:b-c_d123");
  }

  {
    name_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \n\t");
    p._characters (" _12 ");
    p._characters ("  ");
    p._post ();
    assert (p.post_name () == "_12");
  }

  {
    name_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \n\t");
    p._characters (" :12 ");
    p._characters ("  ");
    p._post ();
    assert (p.post_name () == ":12");
  }

  // nmtoken
  //
  {
    nmtoken_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \n\t");
    p._characters (" 123a:b-c_d123 ");
    p._characters (" \n\t");
    p._characters ("  ");
    p._post ();
    assert (p.post_nmtoken () == "123a:b-c_d123");
  }

  // nmtokens
  //
  {
    strings s;
    s.push_back ("123");
    s.push_back ("abc");

    nmtokens_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \n\t");
    p._characters (" 123 ");
    p._characters (" \n\t abc ");
    p._characters ("  ");
    p._post ();
    assert (p.post_nmtokens () == s);
  }

  // ncname
  //
  {
    ncname_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \n\t");
    p._characters (" a.b-c_d123 ");
    p._characters ("  ");
    p._post ();
    assert (p.post_ncname () == "a.b-c_d123");
  }

  // id
  //
  {
    id_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \n\t");
    p._characters (" a.b-c_d123 ");
    p._characters ("  ");
    p._post ();
    assert (p.post_id () == "a.b-c_d123");
  }

  // idref
  //
  {
    idref_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \n\t");
    p._characters (" a.b-c_d123 ");
    p._characters ("  ");
    p._post ();
    assert (p.post_idref () == "a.b-c_d123");
  }

  // idrefs
  //
  {
    strings s;
    s.push_back ("a123");
    s.push_back ("abc");

    idrefs_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \n\t");
    p._characters (" a123 ");
    p._characters (" \n\t abc ");
    p._characters ("  ");
    p._post ();
    assert (p.post_idrefs () == s);
  }

  // language
  //
  {
    language_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" x ");
    p._post ();
    assert (p.post_language () == "x");
  }

  {
    language_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" en ");
    p._post ();
    assert (p.post_language () == "en");
  }

  {
    language_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" en");
    p._characters ("-us ");
    p._post ();
    assert (p.post_language () == "en-us");
  }

  {
    language_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("one-two-three-four44-seven77-eight888");
    p._post ();
    assert (p.post_language () == "one-two-three-four44-seven77-eight888");
  }


  // Bad
  //

  // name
  //
  {
    name_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("");
    assert (test_post_fail (p));
  }

  {
    name_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (".a");
    assert (test_post_fail (p));
  }

  {
    name_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-a");
    assert (test_post_fail (p));
  }

  {
    name_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("1a");
    assert (test_post_fail (p));
  }

  {
    name_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("a,b");
    assert (test_post_fail (p));
  }

  {
    name_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("a b");
    assert (test_post_fail (p));
  }

  {
    name_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("a<b");
    assert (test_post_fail (p));
  }

  // nmtoken
  //
  {
    nmtoken_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("");
    assert (test_post_fail (p));
  }

  {
    nmtoken_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("a,b");
    assert (test_post_fail (p));
  }

  {
    nmtoken_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("a b");
    assert (test_post_fail (p));
  }

  {
    nmtoken_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("a<b");
    assert (test_post_fail (p));
  }

  // nmtokens
  //
  {
    nmtokens_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" ");
    p._characters (" \t\n  ");
    assert (test_post_fail (p));
  }

  {
    nmtokens_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("ab a,b");
    assert (test_post_fail (p));
  }

  // ncname
  //
  {
    ncname_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("");
    assert (test_post_fail (p));
  }

  {
    ncname_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (".a");
    assert (test_post_fail (p));
  }

  {
    ncname_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-a");
    assert (test_post_fail (p));
  }

  {
    ncname_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (":a");
    assert (test_post_fail (p));
  }

  {
    ncname_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("1a");
    assert (test_post_fail (p));
  }

  {
    ncname_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("1:a");
    assert (test_post_fail (p));
  }

  {
    ncname_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("a,b");
    assert (test_post_fail (p));
  }

  {
    ncname_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("a b");
    assert (test_post_fail (p));
  }

  {
    ncname_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("a<b");
    assert (test_post_fail (p));
  }

  // id
  //
  {
    id_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("a b");
    assert (test_post_fail (p));
  }

  // idref
  //
  {
    idref_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("a b");
    assert (test_post_fail (p));
  }

  // idrefs
  //
  {
    idrefs_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("  ");
    p._characters (" \t\n ");
    assert (test_post_fail (p));
  }

  {
    idrefs_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("ab a<b");
    assert (test_post_fail (p));
  }

  // language
  //
  {
    language_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" ");
    assert (test_post_fail (p));
  }

  {
    language_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("en-");
    assert (test_post_fail (p));
  }

  {
    language_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("a1");
    assert (test_post_fail (p));
  }

  {
    language_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("en+us");
    assert (test_post_fail (p));
  }

  {
    language_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("en-nine99999");
    assert (test_post_fail (p));
  }
}
