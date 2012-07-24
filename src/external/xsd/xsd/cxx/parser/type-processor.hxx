// file      : xsd/cxx/parser/type-processor.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef CXX_PARSER_TYPE_PROCESSOR_HXX
#define CXX_PARSER_TYPE_PROCESSOR_HXX

#include <cult/types.hxx>

#include <xsd-frontend/semantic-graph.hxx>

#include <type-map/type-map.hxx>

#include <cxx/parser/cli.hxx>

namespace CXX
{
  namespace Parser
  {
    using namespace Cult::Types;

    class TypeProcessor
    {
    public:
      TypeProcessor (); // Dummy ctor, helps with long symbols on HP-UX.

      Void
      process (CLI::Options const& options,
               XSDFrontend::SemanticGraph::Schema&,
               Boolean gen_driver,
               TypeMap::Namespaces&);
    };
  }
}

#endif // CXX_PARSER_TYPE_PROCESSOR_HXX
