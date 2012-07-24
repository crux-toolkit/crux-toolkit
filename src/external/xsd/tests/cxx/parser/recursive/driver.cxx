// file      : tests/cxx/parser/recursive/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test recursive parser invocation.
//

#include <iostream>
#include <string>

#include "test-pskel.hxx"

using namespace std;

struct sub_pimpl: sub_type_pskel
{
  virtual void
  pre ()
  {
    cout << "sub::pre" << endl;
  }

  virtual void
  sub ()
  {
    cout << "sub::sub" << endl;
  }

  virtual void
  sub2 ()
  {
    cout << "sub::sub2" << endl;
  }

  virtual void
  name (string const& n)
  {
    cout << "sub::name: " << n << endl;
  }

  virtual void
  post_sub_type ()
  {
    cout << "sub::post" << endl;
  }
};

struct indir_pimpl: indir_type_pskel
{
  virtual void
  pre ()
  {
    cout << "indir::pre" << endl;
  }

  virtual void
  sub ()
  {
    cout << "indir::sub" << endl;
  }

  virtual void
  name (string const& n)
  {
    cout << "indir::name: " << n << endl;
  }

  virtual void
  post_indir_type ()
  {
    cout << "indir::post" << endl;
  }
};

struct test_pimpl: test_type_pskel
{
  virtual void
  pre ()
  {
    cout << "test::pre" << endl;
  }

  virtual void
  sub ()
  {
    cout << "test::sub" << endl;
  }

  virtual void
  name (string const& n)
  {
    cout << "test::name: " << n << endl;
  }

  virtual void
  post_test_type ()
  {
    cout << "test::post" << endl;
  }
};


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
    xml_schema::string_pimpl string_p;

    sub_pimpl sub_p;
    indir_pimpl indir_p;
    test_pimpl test_p;

    sub_p.parsers (sub_p, indir_p, sub_p, string_p);
    indir_p.parsers (sub_p, string_p);
    test_p.parsers (sub_p, string_p);

    xml_schema::document doc_p (test_p, "test");

    test_p.pre ();
    doc_p.parse (argv[1]);
    test_p.post_test_type ();
  }
  catch (xml_schema::exception const& e)
  {
    cerr << e << endl;
    return 1;
  }
  catch (ios_base::failure const&)
  {
    cerr << "io failure" << endl;
    return 1;
  }
}
