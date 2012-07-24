// file      : tests/cxx/parser/validation/built-in/any-type/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test the anyType and anySimpleType validation.
//

#include <string>
#include <fstream>
#include <iostream>

#include "test-pskel.hxx"

using namespace std;
using namespace test;
using xml_schema::ro_string;

struct any_type_pimpl: xml_schema::any_type_pimpl
{
  virtual void
  pre ()
  {
    cout << "{" << endl;
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
  post_any_type ()
  {
    cout << "}" << endl
         << endl;
  }
};

struct any_simple_type_pimpl: xml_schema::any_simple_type_pimpl
{
  virtual void
  pre ()
  {
    cout << "{" << endl;
  }

  virtual void
  _any_characters (ro_string const& s)
  {
    cout << "  any text: '" << s << "'" << endl;
  }

  virtual void
  post_any_simple_type ()
  {
    cout << "}" << endl
         << endl;
  }
};

struct any_extension_pimpl: virtual any_extension_pskel,
                            any_type_pimpl

{
  virtual void
  x (const string& v)
  {
    cout << "  x = " << v << endl;
  }
};

struct any_simple_extension_pimpl: virtual any_simple_extension_pskel,
                                   any_simple_type_pimpl
{
  virtual void
  x (const string& v)
  {
    cout << "  x = " << v << endl;
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

    any_type_pimpl any_type_p;
    any_simple_type_pimpl any_simple_type_p;

    any_extension_pimpl any_extension_p;
    any_simple_extension_pimpl any_simple_extension_p;

    type_pimpl type_p;

    any_extension_p.parsers (string_p);
    any_simple_extension_p.parsers (string_p);

    type_p.parsers (any_type_p,
                    any_extension_p,
                    any_simple_extension_p,
                    any_simple_type_p);

    xml_schema::document doc_p (type_p, "test", "root");

    ifstream ifs (argv[1]);
    type_p.pre ();
    doc_p.parse (ifs, argv[1], "", xml_schema::flags::dont_validate);
    type_p.post_type ();
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
