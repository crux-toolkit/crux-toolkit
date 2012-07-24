// file      : xsd/cxx/tree/default-value.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef CXX_TREE_DEFAULT_VALUE_HXX
#define CXX_TREE_DEFAULT_VALUE_HXX

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

#include <cxx/tree/elements.hxx>

namespace CXX
{
  namespace Tree
  {
    struct IsLiteralValue: IsFundamentalType, Traversal::Complex
    {
      IsLiteralValue (Boolean& r);

      virtual Void
      traverse (SemanticGraph::Complex&);

    private:
      Traversal::Inherits inherits_;
    };

    struct LiteralValue: Traversal::Fundamental::Byte,
                         Traversal::Fundamental::UnsignedByte,
                         Traversal::Fundamental::Short,
                         Traversal::Fundamental::UnsignedShort,
                         Traversal::Fundamental::Int,
                         Traversal::Fundamental::UnsignedInt,
                         Traversal::Fundamental::Long,
                         Traversal::Fundamental::UnsignedLong,
                         Traversal::Fundamental::Integer,
                         Traversal::Fundamental::NonPositiveInteger,
                         Traversal::Fundamental::NonNegativeInteger,
                         Traversal::Fundamental::PositiveInteger,
                         Traversal::Fundamental::NegativeInteger,

                         Traversal::Fundamental::Boolean,

                         Traversal::Fundamental::Float,
                         Traversal::Fundamental::Double,
                         Traversal::Fundamental::Decimal,

                         Traversal::Complex,

                         Context
    {
      LiteralValue (Context&);

      String
      dispatch (SemanticGraph::Node& type, String const& value);

      // Handle inheritance.
      //
      virtual Void
      traverse (SemanticGraph::Complex&);

      // Boolean.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Boolean&);

      // Integral types.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Byte&);

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedByte&);

      virtual Void
      traverse (SemanticGraph::Fundamental::Short&);

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedShort&);

      virtual Void
      traverse (SemanticGraph::Fundamental::Int&);

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedInt&);

      virtual Void
      traverse (SemanticGraph::Fundamental::Long&);

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedLong&);

      virtual Void
      traverse (SemanticGraph::Fundamental::Integer&);

      virtual Void
      traverse (SemanticGraph::Fundamental::NonPositiveInteger&);

      virtual Void
      traverse (SemanticGraph::Fundamental::NonNegativeInteger&);

      virtual Void
      traverse (SemanticGraph::Fundamental::PositiveInteger&);

      virtual Void
      traverse (SemanticGraph::Fundamental::NegativeInteger&);

      // Floats.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Float&);

      virtual Void
      traverse (SemanticGraph::Fundamental::Double&);

      virtual Void
      traverse (SemanticGraph::Fundamental::Decimal&);

    private:
      String value_;
      String literal_;

      Traversal::Inherits inherits_;
    };

    // Some initialization (e.g., list) need a function body while others
    // (e.g., *binary) require extra data.
    //
    struct InitKind: Traversal::List,
                     Traversal::Complex,

                     Traversal::Fundamental::Base64Binary,
                     Traversal::Fundamental::HexBinary,

                     Traversal::Fundamental::NameTokens,
                     Traversal::Fundamental::IdRefs,
                     Traversal::Fundamental::Entities
    {
      enum Kind
      {
        simple,
        data,
        function
      };

      // Should be simple initially.
      //
      InitKind (Kind& r);

      virtual Void
      traverse (SemanticGraph::List&);

      virtual Void
      traverse (SemanticGraph::Complex&);

      virtual Void
      traverse (SemanticGraph::Fundamental::Base64Binary&);

      virtual Void
      traverse (SemanticGraph::Fundamental::HexBinary&);

      virtual Void
      traverse (SemanticGraph::Fundamental::NameTokens&);

      virtual Void
      traverse (SemanticGraph::Fundamental::IdRefs&);

      virtual Void
      traverse (SemanticGraph::Fundamental::Entities&);

    private:
      Kind& r_;
      Traversal::Inherits inherits_;
    };

    struct InitValue: Traversal::List,
                      Traversal::Union,
                      Traversal::Complex,

                      Traversal::AnySimpleType,

                      Traversal::Fundamental::String,
                      Traversal::Fundamental::NormalizedString,
                      Traversal::Fundamental::Token,
                      Traversal::Fundamental::Name,
                      Traversal::Fundamental::NameToken,
                      Traversal::Fundamental::NameTokens,
                      Traversal::Fundamental::NCName,
                      Traversal::Fundamental::Language,

                      Traversal::Fundamental::QName,

                      Traversal::Fundamental::Id,
                      Traversal::Fundamental::IdRef,
                      Traversal::Fundamental::IdRefs,

                      Traversal::Fundamental::AnyURI,

                      Traversal::Fundamental::Base64Binary,
                      Traversal::Fundamental::HexBinary,

                      Traversal::Fundamental::Date,
                      Traversal::Fundamental::DateTime,
                      Traversal::Fundamental::Duration,
                      Traversal::Fundamental::Day,
                      Traversal::Fundamental::Month,
                      Traversal::Fundamental::MonthDay,
                      Traversal::Fundamental::Year,
                      Traversal::Fundamental::YearMonth,
                      Traversal::Fundamental::Time,

                      Traversal::Fundamental::Entity,
                      Traversal::Fundamental::Entities,

                      Context
    {
      InitValue (Context&);

      Void
      data (String const& data)
      {
        data_ = data;
        dispatch_count_ = 0;
      }

      Void
      dispatch (SemanticGraph::Node& type, String const& value);

      virtual Void
      traverse (SemanticGraph::List&);

      virtual Void
      traverse (SemanticGraph::Union&);

      virtual Void
      traverse (SemanticGraph::Complex&);

      // anySimpleType.
      //
      virtual Void
      traverse (SemanticGraph::AnySimpleType&);

      // Strings.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::String&);

      virtual Void
      traverse (SemanticGraph::Fundamental::NormalizedString&);

      virtual Void
      traverse (SemanticGraph::Fundamental::Token&);

      virtual Void
      traverse (SemanticGraph::Fundamental::NameToken&);

      virtual Void
      traverse (SemanticGraph::Fundamental::NameTokens&);

      virtual Void
      traverse (SemanticGraph::Fundamental::Name&);

      virtual Void
      traverse (SemanticGraph::Fundamental::NCName&);

      virtual Void
      traverse (SemanticGraph::Fundamental::Language&);

      // Qualified name.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::QName&);

      // ID/IDREF.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Id&);

      virtual Void
      traverse (SemanticGraph::Fundamental::IdRef&);

      virtual Void
      traverse (SemanticGraph::Fundamental::IdRefs&);

      // URI.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::AnyURI&);

      // Binary.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Base64Binary&);

      virtual Void
      traverse (SemanticGraph::Fundamental::HexBinary&);

      // Date/time.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Date&);

      virtual Void
      traverse (SemanticGraph::Fundamental::DateTime&);

      virtual Void
      traverse (SemanticGraph::Fundamental::Duration&);

      virtual Void
      traverse (SemanticGraph::Fundamental::Day&);

      virtual Void
      traverse (SemanticGraph::Fundamental::Month&);

      virtual Void
      traverse (SemanticGraph::Fundamental::MonthDay&);

      virtual Void
      traverse (SemanticGraph::Fundamental::Year&);

      virtual Void
      traverse (SemanticGraph::Fundamental::YearMonth&);

      virtual Void
      traverse (SemanticGraph::Fundamental::Time&);

      // Entity.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Entity&);

      virtual Void
      traverse (SemanticGraph::Fundamental::Entities&);

    private:
      Void
      string_sequence_type (SemanticGraph::Type& element_type);

      Void
      time_zone (Size pos);

    private:
      String value_;
      String data_;
      Size dispatch_count_;
      MemberTypeName type_name_;
      LiteralValue literal_value_;
    };
  }
}

#endif // CXX_TREE_DEFAULT_VALUE_HXX
