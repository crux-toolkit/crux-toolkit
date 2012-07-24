// file      : tests/cxx/tree/detach/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test the detach functionality.
//

#include <memory> // std::auto_ptr
#include <cassert>

#include "test.hxx"

using namespace std;
using namespace test;

int
main ()
{
  // Construct the model.
  //
  object o1 ("o1");
  o1.data ().push_back ("1-1");
  o1.data ().push_back ("1-2");
  o1.data ().push_back ("1-3");

  object o2 ("o2");
  o1.data ().push_back ("2-1");
  o1.data ().push_back ("2-2");
  o1.data ().push_back ("2-3");

  object o3 ("o3");
  o1.data ().push_back ("3-1");
  o1.data ().push_back ("3-2");
  o1.data ().push_back ("3-3");

  object o4 ("o4");
  o1.data ().push_back ("4-1");
  o1.data ().push_back ("4-2");
  o1.data ().push_back ("4-3");

  subtree s1;
  s1.o ().push_back (o1);
  s1.o ().push_back (o2);
  s1.r ().push_back (ref ("o2"));
  s1.r ().push_back (ref ("o3"));

  subtree s2;
  s2.o ().push_back (o3);
  s2.o ().push_back (o4);
  s2.r ().push_back (ref ("o4"));
  s2.r ().push_back (ref ("o1"));

  model m;
  m.one (s1);
  m.opt (s2);

  // Detach one.
  //
  auto_ptr<subtree> p (m.detach_one ());
  assert (p->_container () == 0);
  assert (p->r ()[0].get () == &p->o ()[1]);
  assert (m.opt ()->r ()[1].get () == 0);

  m.one (p);
  assert (m.opt ()->r ()[1].get () == &m.one ().o ()[0]);
  p = m.detach_one ();

  model m1;
  m1.one (p);
  m1.opt (s2);
  assert (m1.opt ()->r ()[1].get () == &m1.one ().o ()[0]);

  p = m1.detach_one ();
  m.seq ().push_back (p);

  // Detach opt.
  //
  p = m.opt ().detach ();
  assert (!m.opt ());
  assert (p->_container () == 0);
  assert (p->r ()[0].get () == &p->o ()[1]);
  assert (m.seq ()[0].r ()[1].get () == 0);

  m.seq ().push_back (p);

  // Detach seq.
  //
  p = m.seq ().detach_back ();
  assert (p->_container () == 0);
  assert (p->r ()[0].get () == &p->o ()[1]);
  assert (m.seq ()[0].r ()[1].get () == 0);

  m.seq ().push_back (p);
  assert (m.seq ()[0].r ()[1].get () == &m.seq ()[1].o ()[0]);

  m.seq ().detach (m.seq ().begin (), p);
  assert (p->_container () == 0);
  assert (p->r ()[0].get () == &p->o ()[1]);
  assert (m.seq ()[0].r ()[1].get () == 0);
}
