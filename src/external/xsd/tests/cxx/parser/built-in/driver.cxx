// file      : tests/cxx/parser/built-in/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test built-in type parsing.
//

#include <string>
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

struct type_pimpl: type_pskel
{
  virtual void
  boolean (bool v)
  {
    cout << v << endl;
  }

  virtual void
  byte (signed char v)
  {
    cout << short (v) << endl;
  }

  virtual void
  unsigned_byte (unsigned char v)
  {
    cout << (unsigned short) (v) << endl;
  }

  virtual void
  short_ (short v)
  {
    cout << v << endl;
  }

  virtual void
  unsigned_short (unsigned short v)
  {
    cout << v << endl;
  }

  virtual void
  int_ (int v)
  {
    cout << v << endl;
  }

  virtual void
  unsigned_int (unsigned int v)
  {
    cout << v << endl;
  }

  virtual void
  long_ (long long v)
  {
    cout << v << endl;
  }

  virtual void
  unsigned_long (unsigned long long v)
  {
    cout << v << endl;
  }

  virtual void
  integer (long long v)
  {
    cout << v << endl;
  }

  virtual void
  negative_integer (long long v)
  {
    cout << v << endl;
  }


  virtual void
  non_positive_integer (long long v)
  {
    cout << v << endl;
  }


  virtual void
  positive_integer (unsigned long long v)
  {
    cout << v << endl;
  }

  virtual void
  non_negative_integer (unsigned long long v)
  {
    cout << v << endl;
  }

  virtual void
  float_ (float v)
  {
    cout << v << endl;
  }

  virtual void
  double_ (double v)
  {
    cout << v << endl;
  }

  virtual void
  decimal (double v)
  {
    cout << v << endl;
  }

  virtual void
  string (std::string const& v)
  {
    cout << "'" << v << "'" << endl;
  }

  virtual void
  normalized_string (std::string const& v)
  {
    cout << "'" << v << "'" << endl;
  }

  virtual void
  token (std::string const& v)
  {
    cout << "'" << v << "'" << endl;
  }

  virtual void
  name (std::string const& v)
  {
    cout << "'" << v << "'" << endl;
  }

  virtual void
  nmtoken (std::string const& v)
  {
    cout << "'" << v << "'" << endl;
  }

  virtual void
  nmtokens (xml_schema::string_sequence const& s)
  {
    cout << "'";

    for (xml_schema::string_sequence::const_iterator i (s.begin ());
         i != s.end (); ++i)
      cout << *i << " ";

    cout << "'" << endl;
  }

  virtual void
  ncname (std::string const& v)
  {
    cout << "'" << v << "'" << endl;
  }

  virtual void
  id (std::string const& v)
  {
    cout << "'" << v << "'" << endl;
  }

  virtual void
  idref (std::string const& v)
  {
    cout << "'" << v << "'" << endl;
  }

  virtual void
  idrefs (xml_schema::string_sequence const& s)
  {
    cout << "'";

    for (xml_schema::string_sequence::const_iterator i (s.begin ());
         i != s.end (); ++i)
      cout << *i << " ";

    cout << "'" << endl;
  }

  virtual void
  language (std::string const& v)
  {
    cout << "'" << v << "'" << endl;
  }

  virtual void
  uri (std::string const& v)
  {
    cout << "'" << v << "'" << endl;
  }

  virtual void
  qname (xml_schema::qname const& v)
  {
    cout << "'" << v.prefix () << ":" << v.name () << "'" << endl;
  }

  virtual void
  base64_binary (std::auto_ptr<xml_schema::buffer> v)
  {
    std::string tmp (v->data (), v->size ());
    cout << "'" << tmp << "'" << endl;
  }

  virtual void
  hex_binary (std::auto_ptr<xml_schema::buffer> v)
  {
    std::string tmp (v->data (), v->size ());
    cout << "'" << tmp << "'" << endl;
  }

  virtual void
  gday (xml_schema::gday const& v)
  {
    cout << v.day ();

    if (v.zone_present ())
      cout << (v.zone_hours () < 0 ? "" : "+") << v.zone_hours ()
           << ':' << v.zone_minutes ();

    cout << endl;
  }

  virtual void
  gmonth (xml_schema::gmonth const& v)
  {
    cout << v.month ();

    if (v.zone_present ())
      cout << (v.zone_hours () < 0 ? "" : "+") << v.zone_hours ()
           << ':' << v.zone_minutes ();

    cout << endl;
  }

  virtual void
  gyear (xml_schema::gyear const& v)
  {
    cout << v.year ();

    if (v.zone_present ())
      cout << (v.zone_hours () < 0 ? "" : "+") << v.zone_hours ()
           << ':' << v.zone_minutes ();

    cout << endl;
  }

  virtual void
  gmonth_day (xml_schema::gmonth_day const& v)
  {
    cout << v.month () << '-' << v.day ();

    if (v.zone_present ())
      cout << (v.zone_hours () < 0 ? "" : "+") << v.zone_hours ()
           << ':' << v.zone_minutes ();

    cout << endl;
  }

  virtual void
  gyear_month (xml_schema::gyear_month const& v)
  {
    cout << v.year () << '-' << v.month ();

    if (v.zone_present ())
      cout << (v.zone_hours () < 0 ? "" : "+") << v.zone_hours ()
           << ':' << v.zone_minutes ();

    cout << endl;
  }

  virtual void
  date (xml_schema::date const& v)
  {
    cout << v.year () << '-' << v.month () << '-' << v.day ();

    if (v.zone_present ())
      cout << (v.zone_hours () < 0 ? "" : "+") << v.zone_hours ()
           << ':' << v.zone_minutes ();

    cout << endl;
  }

  virtual void
  time (xml_schema::time const& v)
  {
    cout << v.hours () << ':' << v.minutes () << ':' << v.seconds ();

    if (v.zone_present ())
      cout << (v.zone_hours () < 0 ? "" : "+") << v.zone_hours ()
           << ':' << v.zone_minutes ();

    cout << endl;
  }

  virtual void
  date_time (xml_schema::date_time const& v)
  {
    cout << v.year () << '-' << v.month () << '-' << v.day () << 'T'
         << v.hours () << ':' << v.minutes () << ':' << v.seconds ();

    if (v.zone_present ())
      cout << (v.zone_hours () < 0 ? "" : "+") << v.zone_hours ()
           << ':' << v.zone_minutes ();

    cout << endl;
  }

  virtual void
  duration (xml_schema::duration const& v)
  {
    cout << (v.negative () ? "-" : "") << 'P'
         << v.years () << 'Y'
         << v.months () << 'M'
         << v.days () << 'D'
         << 'T'
         << v.hours () << 'H'
         << v.minutes () << 'M'
         << v.seconds () << 'S'
         << endl;
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
    any_type_pimpl any_type_p;
    any_simple_type_pimpl any_simple_type_p;

    xml_schema::boolean_pimpl boolean_p;

    xml_schema::byte_pimpl byte_p;
    xml_schema::unsigned_byte_pimpl unsigned_byte_p;
    xml_schema::short_pimpl short_p;
    xml_schema::unsigned_short_pimpl unsigned_short_p;
    xml_schema::int_pimpl int_p;
    xml_schema::unsigned_int_pimpl unsigned_int_p;
    xml_schema::long_pimpl long_p;
    xml_schema::unsigned_long_pimpl unsigned_long_p;

    xml_schema::integer_pimpl integer_p;
    xml_schema::negative_integer_pimpl negative_integer_p;
    xml_schema::non_positive_integer_pimpl non_positive_integer_p;
    xml_schema::positive_integer_pimpl positive_integer_p;
    xml_schema::non_negative_integer_pimpl non_negative_integer_p;

    xml_schema::float_pimpl float_p;
    xml_schema::double_pimpl double_p;
    xml_schema::decimal_pimpl decimal_p;

    xml_schema::string_pimpl string_p;
    xml_schema::normalized_string_pimpl normalized_string_p;
    xml_schema::token_pimpl token_p;
    xml_schema::name_pimpl name_p;
    xml_schema::nmtoken_pimpl nmtoken_p;
    xml_schema::nmtokens_pimpl nmtokens_p;
    xml_schema::ncname_pimpl ncname_p;
    xml_schema::id_pimpl id_p;
    xml_schema::idref_pimpl idref_p;
    xml_schema::idrefs_pimpl idrefs_p;

    xml_schema::language_pimpl language_p;
    xml_schema::uri_pimpl uri_p;
    xml_schema::qname_pimpl qname_p;

    xml_schema::base64_binary_pimpl base64_binary_p;
    xml_schema::hex_binary_pimpl hex_binary_p;

    xml_schema::gday_pimpl gday_p;
    xml_schema::gmonth_pimpl gmonth_p;
    xml_schema::gyear_pimpl gyear_p;
    xml_schema::gmonth_day_pimpl gmonth_day_p;
    xml_schema::gyear_month_pimpl gyear_month_p;
    xml_schema::date_pimpl date_p;
    xml_schema::time_pimpl time_p;
    xml_schema::date_time_pimpl date_time_p;
    xml_schema::duration_pimpl duration_p;

    type_pimpl type_p;

    type_p.parsers (any_type_p,
                    any_simple_type_p,
                    boolean_p,
                    byte_p,
                    unsigned_byte_p,
                    short_p,
                    unsigned_short_p,
                    int_p,
                    unsigned_int_p,
                    long_p,
                    unsigned_long_p,
                    integer_p,
                    negative_integer_p,
                    non_positive_integer_p,
                    positive_integer_p,
                    non_negative_integer_p,
                    float_p,
                    double_p,
                    decimal_p,
                    string_p,
                    normalized_string_p,
                    token_p,
                    name_p,
                    nmtoken_p,
                    nmtokens_p,
                    ncname_p,
                    id_p,
                    idref_p,
                    idrefs_p,
                    language_p,
                    uri_p,
                    qname_p,
                    base64_binary_p,
                    hex_binary_p,
                    gday_p,
                    gmonth_p,
                    gyear_p,
                    gmonth_day_p,
                    gyear_month_p,
                    date_p,
                    time_p,
                    date_time_p,
                    duration_p);

    xml_schema::document doc_p (type_p, "test", "root");

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
