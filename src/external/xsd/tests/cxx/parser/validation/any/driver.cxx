// file      : tests/cxx/parser/validation/any/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test the any particle validation.
//

#include <string>
#include <fstream>
#include <iostream>

#include "test-pskel.hxx"

using namespace std;
using namespace test;
using xml_schema::ro_string;

struct any_a_pimpl: any_a_pskel
{
  virtual void
  pre ()
  {
    cout << "{" << endl;
  }

  virtual void
  a (string const& v)
  {
    cout << "  a = " << v << endl;
  }

  virtual void
  x (string const& v)
  {
    cout << "  x = " << v << endl;
  }

  virtual void
  _start_any_element (ro_string const&, 
                      ro_string const& n,
		      ro_string const*)
  {
    cout << "  start any element '" << n << "'" << endl;
  }

  virtual void
  _end_any_element (ro_string const&, ro_string const& n)
  {
    cout << "  end any element '" << n << "'" << endl;
  }

  virtual void
  _any_attribute (ro_string const&,
                  ro_string const& n,
		  ro_string const& v)
  {
    cout << "  any attribute " << n << " = '" << v << "'" << endl;
  }

  virtual void
  _any_characters (ro_string const& s)
  {
    cout << "  any text: '" << s << "'" << endl;
  }

  virtual void
  post_any_a ()
  {
    cout << "}" << endl
         << endl;
  }
};

struct type_pimpl: type_pskel
{
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
    any_a_pimpl any_a_p;
    type_pimpl type_p;

    any_a_p.parsers (string_p, string_p);
    type_p.parsers (any_a_p);

    xml_schema::document doc_p (type_p, "test", "root");

    try
    {
      ifstream ifs (argv[1]);
      type_p.pre ();
      doc_p.parse (ifs, argv[1], "", xml_schema::flags::dont_validate);
      type_p.post_type ();
    }
    catch (xml_schema::exception const& e)
    {
      cout << "  " << e << endl
           << "}" << endl
           << endl;
    }
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
