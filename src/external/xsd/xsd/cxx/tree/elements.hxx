// file      : xsd/cxx/tree/elements.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef CXX_TREE_ELEMENTS_HXX
#define CXX_TREE_ELEMENTS_HXX

#include <sstream>

#include <cult/containers/map.hxx>
#include <cult/containers/deque.hxx>
#include <cult/containers/vector.hxx>

#include <backend-elements/regex.hxx>

#include <cxx/elements.hxx>

#include <cxx/tree/cli.hxx>


namespace CXX
{
  namespace Tree
  {
    struct Counts
    {
      Counts ()
          : global_types (0),
            global_elements (0),
            generated_global_elements (0),
            complexity_total (0)
      {
      }

      UnsignedLong global_types;
      UnsignedLong global_elements;
      UnsignedLong generated_global_elements;

      // Complexity value for each global type and generated global
      // element, in order.
      //
      Cult::Containers::Vector<UnsignedLong> complexity;
      UnsignedLong complexity_total;
    };

    struct InvalidCustomTypeMapping
    {
      InvalidCustomTypeMapping (String const& mapping,
                                String const& reason)
          : mapping_ (mapping), reason_ (reason)
      {
      }

      String const&
      mapping () const
      {
        return mapping_;
      }

      String const&
      reason () const
      {
        return reason_;
      }

    private:
      String mapping_;
      String reason_;
    };

    //
    //
    class Context: public CXX::Context
    {
    public:
      typedef BackendElements::Regex::Expression<Char> Regex;
      typedef BackendElements::Regex::Pattern<WideChar> WideRegexPat;

      struct DirectCustomTypeMapInfo
      {
        DirectCustomTypeMapInfo (String const& t = L"",
                                 String const& b = L"")
            : type (t), base (b)
        {
        }

        String type;
        String base;
      };

      struct RegexCustomTypeMapInfo
      {
        RegexCustomTypeMapInfo (WideRegexPat const& p,
                                String const& t,
                                String const& b)
            : pat (p), type_sub (t), base_sub (b)
        {
        }

        WideRegexPat pat;
        String type_sub;
        String base_sub;
      };

      typedef
      Cult::Containers::Vector<RegexCustomTypeMapInfo>
      RegexCustomTypeMap;

      typedef
      Cult::Containers::Map<String, DirectCustomTypeMapInfo>
      DirectCustomTypeMap;

    public:
      Context (std::wostream& o,
               SemanticGraph::Schema& root,
               SemanticGraph::Path const& path,
               CLI::Options const& ops,
               Counts const& counts_,
               Boolean generate_xml_schema,
               StringLiteralMap const*,
               Regex const* fwd_expr,
               Regex const* hxx_expr,
               Regex const* ixx_expr);

    protected:
      Context (Context& c);
      Context (Context& c, std::wostream& o);

      // Custom type mapping.
      //
    public:
      // Returns empty string if custom mapping is not required.
      //
      String
      custom_type (SemanticGraph::Type const&) const;

      // Returns true if custom mapping is required. name is
      // populated with the custom type name or empty if the
      // original name should be used.
      //
      Boolean
      custom_type (SemanticGraph::Type const&, String& name) const;

      // Returns true if this type has been renamed as part of the
      // customization process. If the function returns true, the
      // name string is populated with the new name or empty if
      // the type should not be generated at all.
      //
      Boolean
      renamed_type (SemanticGraph::Type const&, String& name) const;

    public:
      // Performs a number of processing steps, including forcing a new
      // line after 80 characters as well as "commentizing" the text by
      // adding '* ' after each newline.
      //
      Void
      write_annotation (SemanticGraph::Annotation&);

      //
      //
    public:
      Boolean
      polymorphic_p (SemanticGraph::Type&);

      Boolean
      anonymous_p (SemanticGraph::Type const& t)
      {
        return t.context ().count ("anonymous");
      }

      // Escaped names.
      //
    public:
      // Accessor name.
      //
      static String const&
      eaname (SemanticGraph::Member const& m)
      {
        return m.context ().get<String> ("aname");
      }

      static String const&
      eaname (SemanticGraph::Any const& a)
      {
        return a.context ().get<String> ("aname");
      }

      static String const&
      eaname (SemanticGraph::AnyAttribute const& a)
      {
        return a.context ().get<String> ("aname");
      }

      // Modifier name.
      //
      static String const&
      emname (SemanticGraph::Member const& m)
      {
        return m.context ().get<String> ("mname");
      }

      static String const&
      emname (SemanticGraph::Any const& a)
      {
        return a.context ().get<String> ("mname");
      }

      static String const&
      emname (SemanticGraph::AnyAttribute const& a)
      {
        return a.context ().get<String> ("mname");
      }

      // Detach name.
      //
      static String const&
      edname (SemanticGraph::Member const& m)
      {
        return m.context ().get<String> ("dname");
      }

      //
      //
      static String const&
      etype (SemanticGraph::Member const& m)
      {
        return m.context ().get<String> ("type");
      }

      static String const&
      etraits (SemanticGraph::Member const& m)
      {
        return m.context ().get<String> ("traits");
      }

      static String const&
      econtainer (SemanticGraph::Member const& m)
      {
        return m.context ().get<String> ("container");
      }

      static String const&
      econtainer (SemanticGraph::Any const& a)
      {
        return a.context ().get<String> ("container");
      }

      static String const&
      econtainer (SemanticGraph::AnyAttribute const& a)
      {
        return a.context ().get<String> ("container");
      }

      static String const&
      eiterator (SemanticGraph::Member const& m)
      {
        return m.context ().get<String> ("iterator");
      }

      static String const&
      eiterator (SemanticGraph::Any const& a)
      {
        return a.context ().get<String> ("iterator");
      }

      static String const&
      eiterator (SemanticGraph::AnyAttribute const& a)
      {
        return a.context ().get<String> ("iterator");
      }

      static String const&
      econst_iterator (SemanticGraph::Member const& m)
      {
        return m.context ().get<String> ("const-iterator");
      }

      static String const&
      econst_iterator (SemanticGraph::Any const& a)
      {
        return a.context ().get<String> ("const-iterator");
      }

      static String const&
      econst_iterator (SemanticGraph::AnyAttribute const& a)
      {
        return a.context ().get<String> ("const-iterator");
      }

      static String const&
      emember (SemanticGraph::Member const& m)
      {
        return m.context ().get<String> ("member");
      }

      static String const&
      emember (SemanticGraph::Any const& a)
      {
        return a.context ().get<String> ("member");
      }

      static String const&
      emember (SemanticGraph::AnyAttribute const& a)
      {
        return a.context ().get<String> ("member");
      }

      static String const&
      edefault_value (SemanticGraph::Member const& m)
      {
        return m.context ().get<String> ("default-value");
      }

      static String const&
      edefault_value_member (SemanticGraph::Member const& m)
      {
        return m.context ().get<String> ("default-value-member");
      }

      // Underlying enum value type.
      //
      static String const&
      evalue (SemanticGraph::Enumeration const& e)
      {
        return e.context ().get<String> ("value");
      }

      // dom_document
      //
      static Boolean
      edom_document_p (SemanticGraph::Complex const& c)
      {
        return c.context ().count ("dom-document");
      }

      static String const&
      edom_document (SemanticGraph::Complex const& c)
      {
        return c.context ().get<String> ("dom-document");
      }

      static Boolean
      edom_document_member_p (SemanticGraph::Complex const& c)
      {
        return c.context ().count ("dom-document-member");
      }

      static String const&
      edom_document_member (SemanticGraph::Complex const& c)
      {
        return c.context ().get<String> ("dom-document-member");
      }

      // Parsing and serialization function names.
      //
      static String const&
      eparser (SemanticGraph::Element const& e)
      {
        return e.context ().get<String> ("parser");
      }

      static String const&
      eserializer (SemanticGraph::Element const& e)
      {
        return e.context ().get<String> ("serializer");
      }

    public:
      Void
      enter_ns_scope (String const& name)
      {
        ns_scope_stack.push_back (name);
        update_ns_scope ();
      }

      Void
      leave_ns_scope ()
      {
        ns_scope_stack.pop_back ();
        update_ns_scope ();
      }

    private:
      Void
      update_ns_scope ();

    private:
      // Write text that may contain characters that we will have
      // to escape (indicated by the rogue flag).
      //
      Void
      write_rogue_text (WideChar const* s, Size size, Boolean rogue);

    public:
      CLI::Options const& options;
      Counts const& counts;
      String& any_type;
      String& any_simple_type;
      String& element_type;
      String& container;
      String& flags_type;
      String& qname_type;
      String& xs_string_type;
      String& properties_type;
      String& error_handler_type;
      String& list_stream_type;
      String& namespace_infomap_type;
      String& parser_type;
      String& std_ostream_type;
      String& ostream_type;
      String& istream_type;
      String& xerces_ns;
      String& dom_auto_ptr;
      String& dom_node_key;
      String& as_double_type;
      String& as_decimal_type;

      Boolean& generate_xml_schema;
      Boolean& doxygen;
      Boolean polymorphic;
      Boolean polymorphic_all;
      Boolean detach;

      Regex const* fwd_expr;
      Regex const* hxx_expr;
      Regex const* ixx_expr;

      String& ns_scope;

      RegexCustomTypeMap& regex_custom_type_map;
      DirectCustomTypeMap& direct_custom_type_map;

    private:
      String any_type_;
      String any_simple_type_;
      String element_type_;
      String container_;
      String flags_type_;
      String qname_type_;
      String xs_string_type_;
      String properties_type_;
      String error_handler_type_;
      String list_stream_type_;
      String namespace_infomap_type_;
      String parser_type_;
      String std_ostream_type_;
      String ostream_type_;
      String istream_type_;
      String xerces_ns_;
      String dom_auto_ptr_;
      String dom_node_key_;
      String as_double_type_;
      String as_decimal_type_;

      Boolean generate_xml_schema_;
      Boolean doxygen_;

      typedef
      Cult::Containers::Deque<String>
      NamespaceStack;

      typedef
      Cult::Containers::Deque<String>
      ScopeStack;

      String ns_scope_;

      NamespaceStack& ns_scope_stack;
      NamespaceStack ns_scope_stack_;

      RegexCustomTypeMap regex_custom_type_map_;
      DirectCustomTypeMap direct_custom_type_map_;

    private:
      WideRegexPat const cxx_uq_id_expr_;
      WideRegexPat const& cxx_uq_id_expr;
    };

    // Check whether this Schema type maps to a fundamental C++ type.
    //
    struct IsFundamentalType: Traversal::Fundamental::Byte,
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
                              Traversal::Fundamental::Decimal

    {
      IsFundamentalType (Boolean& r)
          : r_ (r)
      {
      }

      // Integral types.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Byte&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedByte&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Short&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedShort&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Int&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedInt&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Long&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedLong&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Integer&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NonPositiveInteger&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NonNegativeInteger&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::PositiveInteger&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NegativeInteger&)
      {
        r_ = true;
      }

      // Boolean.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Boolean&)
      {
        r_ = true;
      }

      // Floats.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Float&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Double&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Decimal&)
      {
        r_ = true;
      }

    private:
      Boolean& r_;
    };

    // Check whether this is a string-based type.
    //
    struct IsStringBasedType: Traversal::Complex,
                              Traversal::Union,
                              Traversal::Fundamental::String,
                              Traversal::Fundamental::NormalizedString,
                              Traversal::Fundamental::Token,
                              Traversal::Fundamental::Name,
                              Traversal::Fundamental::NameToken,
                              Traversal::Fundamental::NCName,
                              Traversal::Fundamental::Language
    {
      IsStringBasedType (Boolean& r)
          : r_ (r)
      {
        *this >> inherits_ >> *this;
      }

      virtual Void
      traverse (SemanticGraph::Complex& c)
      {
        inherits (c, inherits_);
      }

      virtual Void
      traverse (SemanticGraph::Union&)
      {
        // Current mapping of union is string-based.
        //
        r_ = true;
      }

      // Strings.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::String&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NormalizedString&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Token&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NameToken&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Name&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NCName&)
      {
        r_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Language&)
      {
        r_ = true;
      }

    private:
      Boolean& r_;
      Traversal::Inherits inherits_;
    };


    // Check whether this is a enumeration-based type.
    //
    struct IsEnumBasedType: Traversal::Complex
    {
      IsEnumBasedType (SemanticGraph::Enumeration*& e)
          : enum_ (e)
      {
        *this >> inherits_;

        inherits_ >> *this;
        inherits_ >> enum_;
      }

      virtual Void
      traverse (SemanticGraph::Complex& c)
      {
        inherits (c, inherits_);
      }

    private:
      struct Enumeration: Traversal::Enumeration
      {
        Enumeration (SemanticGraph::Enumeration*& e)
            : e_ (e)
        {
        }

        virtual Void
        traverse (Type& e)
        {
          if (e_ == 0)
            e_ = &e;
        }

      private:
        SemanticGraph::Enumeration*& e_;
      };


    private:
      Enumeration enum_;
      Traversal::Inherits inherits_;
    };


    //
    //
    struct MemberTypeName : Context,
                            Traversal::Type,
                            Traversal::List,
                            Traversal::Union,
                            Traversal::Complex,

                            Traversal::AnyType,
                            Traversal::AnySimpleType,

                            Traversal::Fundamental::Byte,
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
                            Traversal::Fundamental::Entities

    {
      MemberTypeName (Context& c)
          : Context (c)
      {
      }

      MemberTypeName (Context& c, std::wostream& o)
          : Context (c, o)
      {
      }

      virtual Void
      traverse (SemanticGraph::Type&)
      {
        abort ();
      }

      virtual Void
      traverse (SemanticGraph::List& l)
      {
        os << fq_name (l);
      }

      virtual Void
      traverse (SemanticGraph::Union& u)
      {
        os << fq_name (u);
      }

      virtual Void
      traverse (SemanticGraph::Complex& c)
      {
        os << fq_name (c);
      }

      // anyType & anySimpleType.
      //
      virtual Void
      traverse (SemanticGraph::AnyType& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::AnySimpleType& t)
      {
        os << fq_name (t);
      }

      // Integral types.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Byte& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedByte& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Short& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedShort& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Int& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedInt& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Long& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedLong& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Integer& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NonPositiveInteger& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NonNegativeInteger& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::PositiveInteger& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NegativeInteger& t)
      {
        os << fq_name (t);
      }

      // Boolean.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Boolean& t)
      {
        os << fq_name (t);
      }

      // Floats.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Float& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Double& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Decimal& t)
      {
        os << fq_name (t);
      }

      // Strings.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::String& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NormalizedString& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Token& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NameToken& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NameTokens& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Name& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NCName& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Language& t)
      {
        os << fq_name (t);
      }


      // Qualified name.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::QName& t)
      {
        os << fq_name (t);
      }


      // ID/IDREF.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Id& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::IdRef& t)
      {
        if (t.named_p ())
        {
          // IDREF<anyType>
          //
          os << fq_name (t);
        }
        else
        {
          SemanticGraph::Nameable& ncname (
            xs_ns ().find ("NCName").first->named ());

          os << "::xsd::cxx::tree::idref< " << char_type << ", " <<
            fq_name (ncname) << ", "  <<
            type_name (t.argumented ().type ()) << " >";
        }
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::IdRefs& t)
      {
        if (t.named_p ())
        {
          // IDREFS<anyType>
          //
          os << fq_name (t);
        }
        else
        {
          SemanticGraph::Nameable& ncname (
            xs_ns ().find ("NCName").first->named ());

          os << "::xsd::cxx::tree::idrefs< " << char_type << ", " <<
            any_simple_type << ", ::xsd::cxx::tree::idref< " <<
            char_type << ", " << fq_name (ncname) << ", " <<
            type_name (t.argumented ().type ()) << " > >";
        }
      }

      // URI.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::AnyURI& t)
      {
        os << fq_name (t);
      }

      // Binary.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Base64Binary& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::HexBinary& t)
      {
        os << fq_name (t);
      }


      // Date/time.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Date& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::DateTime& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Duration& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Day& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Month& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::MonthDay& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Year& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::YearMonth& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Time& t)
      {
        os << fq_name (t);
      }

      // Entity.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Entity& t)
      {
        os << fq_name (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Entities& t)
      {
        os << fq_name (t);
      }

    private:
      // For idref/idrefs
      //
      String
      type_name (SemanticGraph::Type& t)
      {
        // This type is always named.
        //
        std::wostringstream o;

        MemberTypeName type (*this, o);
        type.dispatch (t);

        return o.str ();
      }
    };


    //
    //
    struct BaseTypeName : MemberTypeName
    {
      BaseTypeName (Context& c)
          : MemberTypeName (c)
      {
      }

      BaseTypeName (Context& c, std::wostream& o)
          : MemberTypeName (c, o)
      {
      }

      virtual Void
      fundamental_base (SemanticGraph::Type& t)
      {
        os << "::xsd::cxx::tree::fundamental_base< " <<
          fq_name (t) << ", " << char_type << ", " <<
          any_simple_type << " >";
      }

      // Integrals.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Byte& t)
      {
        fundamental_base (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedByte& t)
      {
        fundamental_base (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Short& t)
      {
        fundamental_base (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedShort& t)
      {
        fundamental_base (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Int& t)
      {
        fundamental_base (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedInt& t)
      {
        fundamental_base (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Long& t)
      {
        fundamental_base (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedLong& t)
      {
        fundamental_base (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Integer& t)
      {
        fundamental_base (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NonPositiveInteger& t)
      {
        fundamental_base (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NonNegativeInteger& t)
      {
        fundamental_base (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::PositiveInteger& t)
      {
        fundamental_base (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NegativeInteger& t)
      {
        fundamental_base (t);
      }

      // Boolean.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Boolean& t)
      {
        fundamental_base (t);
      }

      // Floats.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Float& t)
      {
        fundamental_base (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Double& t)
      {
        os << "::xsd::cxx::tree::fundamental_base< " <<
          fq_name (t) << ", " << char_type << ", " <<
          any_simple_type << ", " <<
          "::xsd::cxx::tree::schema_type::double_ >";
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Decimal& t)
      {
        os << "::xsd::cxx::tree::fundamental_base< " <<
          fq_name (t) << ", " << char_type << ", " <<
          any_simple_type << ", " <<
          "::xsd::cxx::tree::schema_type::decimal >";
      }
    };

    // Initial value should be true.
    //
    struct IsSimpleType: Traversal::Complex,
                         Traversal::Member,
                         Traversal::Any,
                         Traversal::AnyAttribute
    {
      IsSimpleType (Boolean& v)
          : v_ (v)
      {
        *this >> names_ >> *this;
        *this >> inherits_ >> *this;
      }

      virtual Void
      traverse (SemanticGraph::Complex& c)
      {
        names (c, names_);

        if (v_)
          inherits (c, inherits_);
      }

      virtual Void
      traverse (SemanticGraph::Member&)
      {
        v_ = false;
      }

      virtual Void
      traverse (SemanticGraph::Any&)
      {
        v_ = false;
      }

      virtual Void
      traverse (SemanticGraph::AnyAttribute&)
      {
        v_ = false;
      }

    private:
      Boolean& v_;
      Traversal::Names names_;
      Traversal::Inherits inherits_;
    };

    // Test whether we need to generate default c-tor. Note that we are not
    // interested in anyAttribute since it is always mapped to a sequence.
    //
    struct GenerateDefaultCtor: Traversal::Complex,
                                Traversal::Enumeration,
                                Traversal::Type,
                                Traversal::Element,
                                Traversal::Attribute,
                                Traversal::Any,
                                Context
    {
      // generate should initially be false.
      //
      GenerateDefaultCtor (Context&, Boolean& generate, Boolean no_base);

      virtual Void
      traverse (SemanticGraph::Complex&);

      virtual Void
      traverse (SemanticGraph::Type&);

      virtual Void
      traverse (SemanticGraph::Enumeration&);

      virtual Void
      traverse (SemanticGraph::Element&);

      virtual Void
      traverse (SemanticGraph::Attribute&);

      virtual Void
      traverse (SemanticGraph::Any&);

    private:
      Boolean& generate_;
      Boolean no_base_;

    private:
      Traversal::Inherits inherits_;
      Traversal::Names names_;
    };

    // Test whether we need to generate from-base c-tor.
    //
    struct GenerateFromBaseCtor: Traversal::Complex
    {
      // generate should initially be false.
      //
      GenerateFromBaseCtor (Context& c, Boolean& generate);

      virtual Void
      traverse (SemanticGraph::Complex& c);

    private:
      Boolean& generate_;
      Boolean custom_;

      // Note that we are not interested in anyAttribute since it is always
      // mapped to a sequence.
      //
      struct Traverser: Traversal::Type,
                        Traversal::Complex,
                        Traversal::Element,
                        Traversal::Attribute,
                        Traversal::Any,
                        Context
      {
        Traverser (Context& c, Boolean& generate, Boolean& custom);

        virtual Void
        traverse (SemanticGraph::Type&);

        virtual Void
        traverse (SemanticGraph::Complex&);

        virtual Void
        traverse (SemanticGraph::Attribute&);

        virtual Void
        traverse (SemanticGraph::Element&);

        virtual Void
        traverse (SemanticGraph::Any&);

      private:
        Boolean& generate_;
        Boolean& custom_;

      private:
        Traversal::Inherits inherits_;
        Traversal::Names names_;
      } traverser_;

      Traversal::Inherits inherits_;
    };

    // Test whether the type has any non-optional element of complex
    // (has attributes/elements) and polymorpjic types.
    //
    struct HasComplexPolyNonOptArgs: Traversal::Complex,
                                     Traversal::Element,
                                     Context
    {
      // complex and poly should initially be false. clash
      // should initially be true.
      //
      HasComplexPolyNonOptArgs (Context& c,
                                Boolean including_base,
                                Boolean& complex,
                                Boolean& poly,
                                Boolean& clash);

      virtual Void
      traverse (SemanticGraph::Complex&);

      virtual Void
      traverse (SemanticGraph::Element&);

    private:
      Boolean& complex_;
      Boolean& poly_;
      Boolean& clash_;

      Traversal::Inherits inherits_;
      Traversal::Names names_;
    };

    // Immediate non-optional member. Note that AnyAttribute is always
    // mapped to a sequence.
    //
    struct FromBaseCtorArg: Traversal::Any,
                            Traversal::Element,
                            Traversal::Attribute,
                            Context
    {
      enum ArgType
      {
        arg_type,
        arg_complex_auto_ptr,
        arg_poly_auto_ptr
      };

      FromBaseCtorArg (Context& c, ArgType, Boolean arg);

      virtual Void
      traverse (SemanticGraph::Any&);

      virtual Void
      traverse (SemanticGraph::Attribute&);

      virtual Void
      traverse (SemanticGraph::Element&);

    private:
      ArgType arg_type_;
      Boolean arg_;
    };

    // List of all non-optional members and a simple base. Note that
    // AnyAttribute is always mapped to a sequence.
    //
    struct CtorArgs: Traversal::Complex,
                     Traversal::Enumeration,
                     Traversal::Type,
                     Traversal::Any,
                     Traversal::Element,
                     Traversal::Attribute,
                     Context
    {
      enum ArgType
      {
        arg_type,
        arg_complex_auto_ptr,
        arg_poly_auto_ptr
      };

      // The second version outputs the argument name and stores
      // in in the base_arg string.
      //
      CtorArgs (Context&, ArgType);
      CtorArgs (Context&, ArgType, String& base_arg);

      virtual Void
      traverse (SemanticGraph::Type&);

      virtual Void
      traverse (SemanticGraph::Enumeration&);

      virtual Void
      traverse (SemanticGraph::Any&);

      virtual Void
      traverse (SemanticGraph::Attribute&);

      virtual Void
      traverse (SemanticGraph::Element&);

    private:
      String
      comma ();

    private:
      ArgType arg_type_;
      String base_;
      String* base_arg_;
      Boolean first_;

    private:
      Traversal::Inherits inherits_;
      Traversal::Names names_;

      MemberTypeName member_name_;
    };


    // Check whether we need to generate c-tor without the base argument.
    //
    struct GenerateWithoutBaseCtor: Traversal::List,
                                    Traversal::Union,
                                    Traversal::Complex,
                                    Traversal::Enumeration,

                                    Traversal::AnyType,
                                    Traversal::AnySimpleType,

                                    Traversal::Fundamental::String,
                                    Traversal::Fundamental::NormalizedString,
                                    Traversal::Fundamental::Token,

                                    Traversal::Fundamental::NameTokens,
                                    Traversal::Fundamental::IdRefs,

                                    Traversal::Fundamental::Base64Binary,
                                    Traversal::Fundamental::HexBinary
    {
      // generate should initially be false.
      //
      GenerateWithoutBaseCtor (Boolean& generate)
          : generate_ (generate)
      {
        *this >> inherits_ >> *this;
      }

      virtual Void
      traverse (SemanticGraph::List&)
      {
        generate_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Union&)
      {
        // No default initialization.
      }

      virtual Void
      traverse (SemanticGraph::Complex& c)
      {
        Complex::inherits (c);
      }

      virtual Void
      traverse (SemanticGraph::Enumeration&)
      {
        // No default initialization.
      }

      virtual Void
      traverse (SemanticGraph::AnyType&)
      {
        generate_ = true;
      }

      virtual Void
      traverse (SemanticGraph::AnySimpleType&)
      {
        generate_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::String&)
      {
        generate_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NormalizedString&)
      {
        generate_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Token&)
      {
        generate_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NameTokens&)
      {
        generate_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::IdRefs&)
      {
        generate_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Base64Binary&)
      {
        generate_ = true;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::HexBinary&)
      {
        generate_ = true;
      }

    private:
      Boolean& generate_;
      Traversal::Inherits inherits_;
    };


    // List of all non-optional members sans simple base. Note that
    // AnyAttribute is always mapped to a sequence.
    //
    struct CtorArgsWithoutBase: Traversal::Complex,
                                Traversal::Any,
                                Traversal::Element,
                                Traversal::Attribute,
                                Context
    {
      enum ArgType
      {
        arg_type,
        arg_complex_auto_ptr,
        arg_poly_auto_ptr
      };

      CtorArgsWithoutBase (Context& c, ArgType, Boolean arg, Boolean first);

      virtual Void
      traverse (SemanticGraph::Any&);

      virtual Void
      traverse (SemanticGraph::Element&);

      virtual Void
      traverse (SemanticGraph::Attribute&);

    private:
      String
      comma ();

    private:
      ArgType arg_type_;
      Boolean arg_;
      Boolean first_;

    private:
      Traversal::Inherits inherits_;
      Traversal::Names names_;
    };

    //
    //
    struct GlobalElementBase
    {
      GlobalElementBase (Context& c)
          : ctx_ (c)
      {
      }

      Boolean
      generate_p (SemanticGraph::Element&);

      Boolean
      doc_root_p (SemanticGraph::Element&);

    private:
      Context& ctx_;
    };


    //
    //
    struct Namespace: CXX::Namespace,
                      GlobalElementBase,
                      CXX::Namespace::ScopeTracker
    {
      Namespace (Context&,
                 UnsignedLong first = 1,
                 UnsignedLong last = 0);

      virtual Void
      traverse (Type&);

    protected:
      virtual Void
      enter (Type&, String const& name, Boolean last);

      virtual Void
      leave ();

    protected:
      Context& ctx_;

    private:
      UnsignedLong first_;
      UnsignedLong last_;
      UnsignedLong count_;
    };

    //
    //
    struct DocumentedNamespace: Namespace
    {
      DocumentedNamespace (Context& c)
          : Namespace (c)
      {
      }

      virtual Void
      enter (Type& ns, String const& name, Boolean last)
      {
        Namespace::enter (ns, name, last);

        // Only add documentation to the innermost namespace.
        //
        if (ctx_.doxygen && name && last)
        {
          ctx_.os << "/**" << endl
                  << " * @brief C++ namespace for the %" <<
		  ctx_.comment (ns.name ()) << endl
                  << " * schema namespace." << endl
                  << " */" << endl;
        }
      }
    };

    //
    //
    struct TypeForward: Traversal::Type, Context
    {
      TypeForward (Context& c)
          : Context (c)
      {
      }

      virtual Void
      traverse (SemanticGraph::Type& t);
    };

    struct Includes: Traversal::Imports,
                     Traversal::Includes
    {
      enum Type
      {
        forward,
        header,
        inline_,
        source
      };

      Includes (Context& c, Type type)
          : ctx_ (c),
            type_ (type),
            forward_ (c.options.value<CLI::generate_forward> ()),
            namespace_ (c),
            type_forward_ (c)
      {
        schema_ >> schema_names_ >> namespace_ >> names_ >> type_forward_;
      }

      virtual Void
      traverse (SemanticGraph::Imports& i)
      {
        traverse_ (i);
      }

      virtual Void
      traverse (SemanticGraph::Includes& i)
      {
        traverse_ (i);
      }

    private:
      Void
      traverse_ (SemanticGraph::Uses&);

    private:
      Context& ctx_;
      Type type_;
      Boolean forward_;

      Traversal::Schema schema_;
      Traversal::Names schema_names_;
      Namespace namespace_;
      Traversal::Names names_;
      TypeForward type_forward_;
    };

    //
    //
    struct FundIncludes: Traversal::Fundamental::Byte,
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
                         Context
    {
      FundIncludes (Context& c, String const& prefix)
          : Context (c), prefix_ (prefix),
            long_ (false), unsigned_long_ (false)
      {
      }

      // Integral types.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Byte& t)
      {
        gen_include (t, "byte.hxx");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedByte& t)
      {
        gen_include (t, "unsigned-byte.hxx");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Short& t)
      {
        gen_include (t, "short.hxx");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedShort& t)
      {
        gen_include (t, "unsigned-short.hxx");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Int& t)
      {
        gen_include (t, "int.hxx");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedInt& t)
      {
        gen_include (t, "unsigned-int.hxx");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Long& t)
      {
        if (!long_)
          long_ = gen_include (t, "long.hxx");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedLong& t)
      {
        if (!unsigned_long_)
          unsigned_long_ = gen_include (t, "unsigned-long.hxx");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Integer& t)
      {
        if (!long_)
          long_ = gen_include (t, "long.hxx");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NonPositiveInteger& t)
      {
        if (!long_)
          long_ = gen_include (t, "long.hxx");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NonNegativeInteger& t)
      {
        if (!unsigned_long_)
          unsigned_long_ = gen_include (t, "unsigned-long.hxx");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::PositiveInteger& t)
      {
        if (!unsigned_long_)
          unsigned_long_ = gen_include (t, "unsigned-long.hxx");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NegativeInteger& t)
      {
        if (!long_)
          long_ = gen_include (t, "long.hxx");
      }

      // Boolean.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Boolean& t)
      {
        gen_include (t, "boolean.hxx");
      }

      // Floats.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Float& t)
      {
        gen_include (t, "float.hxx");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Double& t)
      {
        gen_include (t, "double.hxx");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Decimal& t)
      {
        gen_include (t, "decimal.hxx");
      }

    private:
      Boolean
      gen_include (SemanticGraph::Type& t, String const& file)
      {
        String custom;

        // XML Schema built-in type customization is only possible when
        // we are generating separate header.
        //
        if (generate_xml_schema && custom_type (t, custom))
        {
          String new_name;
          renamed_type (t, new_name);

          if (!new_name)
            return false;
        }

        os << "#include <xsd/cxx/tree/" << prefix_ << "/" << file << ">"
           << endl;

        return true;
      }

    private:
      String prefix_;
      Boolean long_;
      Boolean unsigned_long_;
    };
  }
}

#endif  // CXX_TREE_ELEMENTS_HXX
