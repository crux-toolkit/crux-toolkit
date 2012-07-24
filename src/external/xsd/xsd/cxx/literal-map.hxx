// file      : xsd/cxx/literal-map.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef CXX_LITERAL_MAP_HXX
#define CXX_LITERAL_MAP_HXX

#include <cult/types.hxx>
#include <cult/containers/map.hxx>

namespace CXX
{
  using namespace Cult;
  typedef WideString String;

  typedef Cult::Containers::Map<String, String> StringLiteralMap;

  bool
  read_literal_map (NarrowString const& file, StringLiteralMap& map);
}

#endif // CXX_LITERAL_MAP_HXX
