// file      : tests/cxx/tree/containment/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test tree node containment.
//

#include <memory> // std::auto_ptr
#include <cassert>

#include "test.hxx"

using namespace std;
using namespace test;

int
main ()
{
  // Change of a container in a sub-tree without ID.
  //
  {
    auto_ptr<inner> i (new inner ());
    i->ref ("foo");

    outer o;
    o.i (i);
    o.ref ("foo");

    assert (o.i ()->ref ()->get () == 0);
    assert (o.ref ()->get () == 0);
  }

  // Change of container in a sub-tree with ID inside.
  //
  {
    auto_ptr<inner> i (new inner ());
    inner* p (i.get ());
    i->id ("foo");
    i->ref ("foo");
    assert (i->ref ()->get () == p);

    outer o;
    o.i (i);
    o.ref ("foo");

    assert (o.i ()->ref ()->get () == p);
    assert (o.ref ()->get () == p);
  }

  // Change of a container in ID.
  //
  {
    auto_ptr<xml_schema::id> id (new xml_schema::id ("foo"));

    inner i;
    i.id (id);
    i.ref ("foo");
    assert (i.ref ()->get () == &i);
  }

  // Change of a container in a type derived from ID with ID inside.
  //
  {
    auto_ptr<id_ex> id (new id_ex ("foo"));
    id_ex* p (id.get ());
    id->id ("bar");

    inner i;
    i.id_ex (id);

    i.ref ("foo");
    assert (i.ref ()->get () == &i);

    i.ref ("bar");
    assert (i.ref ()->get () == p);
  }

  // IDREF lists
  //
  {
    id i1 ("a"), i2 ("b");

    auto_ptr<ids> ic (new ids);
    ic->id ().push_back (i1);
    ic->id ().push_back (i2);

    auto_ptr<xml_schema::idrefs> r1 (new xml_schema::idrefs);
    r1->push_back (xml_schema::idref ("a"));
    r1->push_back (xml_schema::idref ("b"));

    auto_ptr<idref_list> r2 (new idref_list);
    r2->push_back (xml_schema::idref ("a"));
    r2->push_back (xml_schema::idref ("b"));

    auto_ptr<idrefs1> rc1 (new idrefs1);
    auto_ptr<idrefs2> rc2 (new idrefs2);

    rc1->idrefs (r1);
    rc2->idrefs (r2);

    model m;
    m.ids (ic);
    m.idrefs1 (rc1);
    m.idrefs2 (rc2);

    assert (m.idrefs1 ().idrefs ()[0].get () != 0);
    assert (m.idrefs1 ().idrefs ()[1].get () != 0);

    assert (m.idrefs2 ().idrefs ()[0].get () != 0);
    assert (m.idrefs2 ().idrefs ()[1].get () != 0);
  }
}
