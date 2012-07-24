// file      : xsd/cxx/parser/state-processor.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef CXX_PARSER_STATE_PROCESSOR_HXX
#define CXX_PARSER_STATE_PROCESSOR_HXX

#include <cult/types.hxx>
#include <xsd-frontend/semantic-graph.hxx>

namespace CXX
{
  namespace Parser
  {
    using namespace Cult::Types;

    class StateProcessor
    {
    public:
      Void
      process (XSDFrontend::SemanticGraph::Schema&,
               XSDFrontend::SemanticGraph::Path const& file);
    };
  }
}

#endif // CXX_PARSER_STATE_PROCESSOR_HXX
