// file      : tests/cxx/tree/wildcard/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test wildcard (any & anyAttribute) mapping.
//

#include <memory> // std::auto_ptr
#include <sstream>
#include <iostream>

#include <xsd/cxx/xml/string.hxx>

#include "test.hxx"

using namespace std;
using namespace test;
using namespace xercesc;

namespace xml = xsd::cxx::xml;

void
print (type& t)
{
  if (t.att ())
    cout << *t.att () << endl;

  type::any_attribute_set& as (t.any_attribute ());

  for (type::any_attribute_iterator i (as.begin ()); i != as.end (); ++i)
  {
    cout << xml::transcode<char> (i->getTextContent ()) << endl;
  }

  cout << xml::transcode<char> (t.any ().getTextContent ()) << endl
       << t.foo () << endl;

  if (t.any1 ())
    cout << xml::transcode<char> (t.any1 ()->getTextContent ()) << endl;

  cout << t.bar () << endl;

  type::any2_sequence& es (t.any2 ());

  for (type::any2_iterator i (es.begin ()); i != es.end (); ++i)
  {
    cout << xml::transcode<char> (i->getTextContent ()) << endl;
  }

  cout << endl;
}

int
main (int argc, char* argv[])
{
  if (argc != 2)
  {
    cerr << "usage: " << argv[0] << " test.xml" << endl;
    return 1;
  }

  XMLPlatformUtils::Initialize ();

  try
  {
    // Test accessors/modifiers for various cardinalities.
    //
    type t;

    DOMDocument& doc (t.dom_document ());

    // one
    //
    {
      DOMElement* e (doc.createElement (xml::string ("a").c_str ()));
      t.any (*e);
      e->release ();
      assert (xml::transcode<char> (t.any ().getTagName ()) == "a");

      t.any (doc.createElement (xml::string ("b").c_str ()));
      assert (xml::transcode<char> (t.any ().getTagName ()) == "b");
    }

    // optional
    //
    {
      assert (!t.any1 ().present ());

      DOMElement* e (doc.createElement (xml::string ("a").c_str ()));
      t.any1 (*e);
      e->release ();
      assert (t.any1 ().present ());
      assert (xml::transcode<char> (t.any1 ().get ().getTagName ()) == "a");

      t.any1 (doc.createElement (xml::string ("b").c_str ()));
      assert (xml::transcode<char> (t.any1 ()->getTagName ()) == "b");

      type::any1_optional c (
        doc.createElement (xml::string ("c").c_str ()), doc);
      t.any1 (c);
      assert (xml::transcode<char> (t.any1 ()->getTagName ()) == "c");
    }


    // sequence
    //
    {
      type::any2_sequence& s (t.any2 ());

      DOMElement* e (doc.createElement (xml::string ("a").c_str ()));
      s.push_back (*e);
      e->release ();
      s.push_back (doc.createElement (xml::string ("b").c_str ()));
      assert (s.size () == 2);

      for (type::any2_iterator i (s.begin ()); i != s.end (); ++i)
      {
        if (i == s.begin ())
          assert (xml::transcode<char> (i->getTagName ()) == "a");
        else
          assert (xml::transcode<char> ((*i).getTagName ()) == "b");
      }

      // copy c-tor
      type::any2_sequence cs (s, doc);
      assert (cs.size () == 2);
      assert (xml::transcode<char> (cs[0].getTagName ()) == "a");
      assert (xml::transcode<char> (cs[1].getTagName ()) == "b");

      // assignment
      t.any2 (cs);
      assert (s.size () == 2);
      assert (xml::transcode<char> (s[0].getTagName ()) == "a");
      assert (xml::transcode<char> (s[1].getTagName ()) == "b");
    }

    // anyAttribute
    //
    {
      type::any_attribute_set& s (t.any_attribute ());

      DOMAttr* a (doc.createAttribute (xml::string ("a").c_str ()));
      s.insert (*a);
      a->release ();
      s.insert (doc.createAttribute (xml::string ("b").c_str ()));
      assert (s.size () == 2);

      assert (s.find ("a") != s.end ());
      assert (s.find ("b") != s.end ());

      for (type::any_attribute_iterator i (s.begin ()); i != s.end (); ++i)
      {
        assert (xml::transcode<char> (i->getName ()) == "a" ||
                xml::transcode<char> ((*i).getName ()) == "b");
      }

      // copy c-tor
      type::any_attribute_set cs (s, doc);
      assert (cs.size () == 2);
      assert (cs.count ("a"));
      assert (cs.count ("b"));

      // assignment
      t.any_attribute (cs);
      assert (s.size () == 2);
      assert (s.count ("a"));
      assert (s.count ("b"));
    }

    // Test parsing
    //
    auto_ptr<type> r (root (argv[1]));
    print (*r);

    // Test serialization.
    //
    xml_schema::namespace_infomap map;

    map["t"].name   = "test";
    map["t"].schema = "test.xsd";
    map["o"].name   = "other";

    stringstream iostr;
    root (iostr, *r, map);

    // cout << iostr.str () << endl
    //      << endl;

    auto_ptr<type> copy (root (iostr, argv[1]));
    assert (*copy == *r);

    print (*copy);
  }
  catch (xml_schema::exception const& e)
  {
    cerr << e << endl;
    return 1;
  }

  XMLPlatformUtils::Terminate ();
}
