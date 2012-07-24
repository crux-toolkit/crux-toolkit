// file      : tests/cxx/parser/enumeration/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test xsd:enumeration parsing.
//

#include <string>
#include <iostream>

#include "test-pskel.hxx"

using namespace std;
using namespace xml_schema;

struct digit_pimpl: test::digit_pskel, int_pimpl
{
};

struct gender_pimpl: test::gender_pskel, string_pimpl
{
  virtual ::gender
  post_gender ()
  {
    std::string str (post_string ());

    if (str == "male")
      return male;
    else
      return female;
  }
};

struct type_pimpl: test::type_pskel
{
  virtual void
  digit (int i)
  {
    cout << i << endl;
  }

  virtual void
  gender (::gender g)
  {
    cout << g << endl;
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
    digit_pimpl digit_p;
    gender_pimpl gender_p;
    type_pimpl type_p;

    type_p.parsers (digit_p, gender_p);

    document doc_p (type_p, "test", "root");

    type_p.pre ();
    doc_p.parse (argv[1]);
    type_p.post_type ();
  }
  catch (xml_schema::exception const& e)
  {
    cerr << e << endl;
    return 1;
  }
  catch (std::ios_base::failure const&)
  {
    cerr << "io failure" << endl;
    return 1;
  }
}
