// file      : xsd/xsd.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef XSD_HXX
#define XSD_HXX

#include <xsd-frontend/semantic-graph/elements.hxx> // Path

#include <cult/types.hxx>
#include <cult/containers/set.hxx>
#include <cult/containers/vector.hxx>


#include <cstdio> // std::remove

using namespace Cult::Types;

//
//
typedef Cult::Containers::Set<NarrowString> WarningSet;

//
//
typedef Cult::Containers::Vector<NarrowString> FileList;

//
//
struct AutoUnlink
{
  AutoUnlink (XSDFrontend::SemanticGraph::Path const& file)
      : file_ (file), canceled_ (false)
  {
  }

  ~AutoUnlink ()
  {
    if (!canceled_)
    {
#if !defined(BOOST_FILESYSTEM_VERSION) || BOOST_FILESYSTEM_VERSION == 2
      std::remove (file_.native_file_string ().c_str ());
#else
      std::remove (file_.string ().c_str ());
#endif
    }
  }

  void
  cancel ()
  {
    canceled_ = true;
  }

private:
  XSDFrontend::SemanticGraph::Path file_;
  Boolean canceled_;
};

//
//
struct AutoUnlinks
{
  Void
  add (XSDFrontend::SemanticGraph::Path const& file)
  {
    unlinks_.push_back (Evptr<AutoUnlink> (new AutoUnlink (file)));
  }

  Void
  cancel ()
  {
    for (Unlinks::Iterator i (unlinks_.begin ()); i != unlinks_.end (); ++i)
    {
      (*i)->cancel ();
    }
  }

private:
  typedef Cult::Containers::Vector<Evptr<AutoUnlink> > Unlinks;
  Unlinks unlinks_;
};

#endif // XSD_HXX
