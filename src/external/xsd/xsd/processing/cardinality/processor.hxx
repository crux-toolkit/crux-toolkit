// file      : processing/cardinality/processor.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef PROCESSING_CARDINALITY_PROCESSOR_HXX
#define PROCESSING_CARDINALITY_PROCESSOR_HXX

#include <cult/types.hxx>

#include <xsd-frontend/semantic-graph/elements.hxx> // Path
#include <xsd-frontend/semantic-graph/schema.hxx>

namespace Processing
{
  namespace Cardinality
  {
    using namespace Cult::Types;

    class Processor
    {
    public:
      struct Failed {};

      Void
      process (XSDFrontend::SemanticGraph::Schema&,
               XSDFrontend::SemanticGraph::Path const& file);
    };
  }
}

#endif // PROCESSING_CARDINALITY_PROCESSOR_HXX
