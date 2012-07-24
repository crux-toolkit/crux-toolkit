// file      : tests/cxx/parser/validation/choice/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test the choice compositor validation.
//

#include <string>
#include <fstream>
#include <iostream>

#include "test-pskel.hxx"

using namespace std;
using namespace test;
using xml_schema::ro_string;

struct choice_pimpl: choice_pskel
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
  b (string const& v)
  {
    cout << "  b = " << v << endl;
  }

  virtual void
  c (string const& v)
  {
    cout << "  c = " << v << endl;
  }

  virtual void
  d (string const& v)
  {
    cout << "  d = " << v << endl;
  }

  virtual void
  _start_any_element (ro_string const& ns, 
                      ro_string const& name,
		      ro_string const*)
  {
    cout << "  any: " << ns << "#" << name << endl
         << "  {" << endl;
  }

  virtual void
  _any_characters (ro_string const& v)
  {
    cout << "    chars = " << v << endl;
  }

  virtual void
  _end_any_element (ro_string const&, ro_string const&)
  {
    cout << "  }" << endl;
  }

  virtual void
  post_choice ()
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
    choice_pimpl choice_p;
    type_pimpl type_p;

    choice_p.parsers (string_p, string_p, string_p, string_p);
    type_p.parsers (choice_p);

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
