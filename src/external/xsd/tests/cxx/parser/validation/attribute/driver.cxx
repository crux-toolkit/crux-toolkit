// file      : tests/cxx/parser/validation/attribute/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test attribute and attribute wildcard (anyAttribute) validation.
//

#include <string>
#include <fstream>
#include <iostream>

#include "test-pskel.hxx"

using namespace std;
using namespace test;
using xml_schema::ro_string;

struct pass_a_pimpl: pass_a_pskel
{
  virtual void
  pre ()
  {
    cout << "pass-a" << endl
         << "{" << endl;
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
  _any_attribute (ro_string const& ns,
                  ro_string const& name,
                  ro_string const& value)
  {
    cout << "  any: " << ns << "#" << name << " = " << value << endl;
  }

  virtual void
  post_pass_a ()
  {
    cout << "}" << endl
         << endl;
  }
};

struct pass_b_pimpl: pass_b_pskel
{
  virtual void
  pre ()
  {
    cout << "pass-b" << endl
         << "{" << endl;
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
  _any_attribute (ro_string const& ns,
                  ro_string const& name,
                  ro_string const& value)
  {
    cout << "  any: " << ns << "#" << name << " = " << value << endl;
  }

  virtual void
  post_pass_b ()
  {
    cout << "}" << endl
         << endl;
  }
};

struct pass_c_pimpl: pass_c_pskel
{
  virtual void
  pre ()
  {
    cout << "pass-c" << endl
         << "{" << endl;
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
  post_pass_c ()
  {
    cout << "}" << endl
         << endl;
  }
};

struct fail_pimpl: fail_pskel
{
  virtual void
  pre ()
  {
    cout << "fail" << endl
         << "{" << endl;
  }

  virtual void
  a (string const& v)
  {
    cout << "  a = " << v << endl;
  }

  virtual void
  post_fail ()
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
    pass_a_pimpl pass_a_p;
    pass_b_pimpl pass_b_p;
    pass_c_pimpl pass_c_p;
    fail_pimpl fail_p;
    type_pimpl type_p;

    pass_a_p.parsers (string_p, string_p);
    pass_b_p.parsers (string_p, string_p);
    pass_c_p.parsers (string_p, string_p);
    fail_p.parsers (string_p);
    type_p.parsers (pass_a_p, pass_b_p, pass_c_p, fail_p);

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
