// file      : xsd/cxx/tree/fundamental-header.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef CXX_TREE_FUNDAMENTAL_HEADER_HXX
#define CXX_TREE_FUNDAMENTAL_HEADER_HXX

#include <cult/containers/set.hxx>
#include <cult/containers/vector.hxx>

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

#include <cxx/tree/elements.hxx>

namespace CXX
{
  namespace Tree
  {
    struct FundamentalNamespace : DocumentedNamespace,

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

                                  Traversal::Fundamental::Id,
                                  Traversal::Fundamental::IdRef,
                                  Traversal::Fundamental::IdRefs,

                                  Traversal::Fundamental::AnyURI,

                                  Traversal::Fundamental::QName,

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
      FundamentalNamespace (Context& c)
          : DocumentedNamespace (c),
            Context (c),
            export_ (c.options.value<CLI::export_xml_schema> () && type_exp)
      {
        *this >> names_ >> *this;

        if (export_)
          xs_ns_ = ns_name (xs_ns ());
      }

      Void
      gen_typedef (String const& name,
                   String const& type,
                   String const& arg1 = L"",
                   String const& arg2 = L"",
                   String const& arg3 = L"",
                   Boolean export_type = true)
      {
        os << "typedef " << type;

        // Use unqualified arguments since we are in the same
        // namespace.
        //
        if (arg1)
        {
          os << arg1;

          if (arg2)
          {
            os << ", " << arg2;

            if (arg3)
              os << ", " << arg3;
          }

          os << " >";
        }

        os << " " << name << ";";

        if (export_type && export_ && type.find (L'<') != String::npos)
        {
          String s (type);

          // Use qualified arguments.
          //
          if (arg1)
          {
            s += xs_ns_;
            s += L"::";
            s += arg1;

            if (arg2)
            {
              s += L", ";
              s += xs_ns_;
              s += L"::";
              s += arg2;

              if (arg3)
              {
                s += L", ";
                s += xs_ns_;
                s += L"::";
                s += arg3;
              }
            }

            s += " >";
          }

          if (exports_set_.count (s) == 0)
          {
            exports_.push_back (s);
            exports_set_.insert (s);
          }
        }
      }

      String
      built_in_type (SemanticGraph::Type& t,
                     String const& type,
                     String const& arg1 = L"",
                     String const& arg2 = L"",
                     String const& arg3 = L"")
      {
        String custom;

        String name (ename (t));

        // XML Schema built-in type customization is only possible when
        // we are generating separate header.
        //
        if (generate_xml_schema && custom_type (t, custom))
        {
          if (custom.empty ())
            custom = name;

          String new_name;
          renamed_type (t, new_name);

          if (new_name)
          {
            gen_typedef (new_name, type, arg1, arg2, arg3);

            if (doxygen)
              os << endl;
          }

          if (doxygen)
            os << "/**" << endl
               << " * @brief C++ type corresponding to the " <<
              comment (t.name ()) << " XML Schema" << endl
               << " * built-in type." << endl
               << " */" << endl;

          if (custom == name)
            os << "class " << name << ";";
          else
            os << "typedef " << custom << " " << name << ";";

          if (doxygen)
            os << endl;
        }
        else
        {
          // Otherwise generate simple typedef.
          //

          if (doxygen)
            os << "/**" << endl
               << " * @brief C++ type corresponding to the " <<
              comment (t.name ()) << " XML Schema" << endl
               << " * built-in type." << endl
               << " */" << endl;

          gen_typedef (name, type, arg1, arg2, arg3);

          if (doxygen)
            os << endl;
        }

        return name;
      }

      // anyType and anySimpleType
      //
      virtual Void
      traverse (SemanticGraph::AnyType& t)
      {
        os << "// anyType and anySimpleType." << endl
           << "//" << endl;

        if (doxygen)
          os << endl;

        type_ = built_in_type (t, "::xsd::cxx::tree::type");
      }

      virtual Void
      traverse (SemanticGraph::AnySimpleType& t)
      {
        simple_type_ = built_in_type (
          t, L"::xsd::cxx::tree::simple_type< ", type_);

        if (doxygen)
          os << "/**" << endl
             << " * @brief Alias for the anyType type." << endl
             << " */" << endl;

        gen_typedef (xs_ns ().context().get<String> ("container"),
                     "::xsd::cxx::tree::type");

        os << endl;

        if (doxygen)
          os << endl;
      }

      // Integrals.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Byte& t)
      {
        os << "// 8-bit" << endl
           << "//" << endl;

        if (doxygen)
          os << endl;

        built_in_type (t, "signed char");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedByte& t)
      {
        built_in_type (t, "unsigned char");
        os << endl;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Short& t)
      {
        os << "// 16-bit" << endl
           << "//" << endl;

        if (doxygen)
          os << endl;

        built_in_type (t, "short");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedShort& t)
      {
        built_in_type (t, "unsigned short");
        os << endl;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Int& t)
      {
        os << "// 32-bit" << endl
           << "//" << endl;

        if (doxygen)
          os << endl;

        built_in_type (t, "int");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedInt& t)
      {
        built_in_type (t, "unsigned int");
        os << endl;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Long& t)
      {
        os << "// 64-bit" << endl
           << "//" <<endl;

        if (doxygen)
          os << endl;

        built_in_type (t, "long long");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedLong& t)
      {
        built_in_type (t, "unsigned long long");
        os << endl;
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Integer& t)
      {
        os << "// Supposed to be arbitrary-length integral types." << endl
           << "//" << endl;

        if (doxygen)
          os << endl;

        built_in_type (t, "long long");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NonPositiveInteger& t)
      {
        built_in_type (t, "long long");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NonNegativeInteger& t)
      {
        built_in_type (t, "unsigned long long");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::PositiveInteger& t)
      {
        built_in_type (t, "unsigned long long");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NegativeInteger& t)
      {
        built_in_type (t, "long long");
        os << endl;
      }

      // Boolean.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Boolean& t)
      {
        os << "// Boolean." << endl
           << "//" << endl;

        if (doxygen)
          os << endl;

        built_in_type (t, "bool");

        os << endl;
      }

      // Floats.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Float& t)
      {
        os << "// Floating-point types." << endl
           << "//" << endl;

        if (doxygen)
          os << endl;

        built_in_type (t, "float");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Double& t)
      {
        double_ = built_in_type (t, "double");
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Decimal& t)
      {
        decimal_ = built_in_type (t, "double");
        os << endl;
      }


      // Strings.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::String& t)
      {
        os << "// String types." << endl
           << "//" << endl;

        if (doxygen)
          os << endl;

        string_ = built_in_type (
          t, L"::xsd::cxx::tree::string< " + char_type + L", ", simple_type_);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NormalizedString& t)
      {
        norm_string_ = built_in_type (
          t,
          L"::xsd::cxx::tree::normalized_string< " + char_type + L", ",
          string_);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Token& t)
      {
        token_ = built_in_type (
          t, L"::xsd::cxx::tree::token< " + char_type + L", ", norm_string_);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NameToken& t)
      {
        nmtoken_ = built_in_type (
          t, L"::xsd::cxx::tree::nmtoken< " + char_type + L", ", token_);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NameTokens& t)
      {
        built_in_type (
          t,
          L"::xsd::cxx::tree::nmtokens< " + char_type + L", ",
          simple_type_,
          nmtoken_);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Name& t)
      {
        name_ = built_in_type (
          t, L"::xsd::cxx::tree::name< " + char_type + L", ", token_);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NCName& t)
      {
        ncname_ = built_in_type (
          t, L"::xsd::cxx::tree::ncname< " + char_type + L", ", name_);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Language& t)
      {
        built_in_type (
          t, L"::xsd::cxx::tree::language< " + char_type + L", ", token_);

        os << endl;
      }

      // ID/IDREF.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Id& t)
      {
        os << "// ID/IDREF." << endl
           << "//" << endl;

        if (doxygen)
          os << endl;

        built_in_type (
          t, L"::xsd::cxx::tree::id< " + char_type + L", ", ncname_);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::IdRef& t)
      {
        idref_ = built_in_type (
          t, L"::xsd::cxx::tree::idref< " + char_type + L", ", ncname_, type_);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::IdRefs& t)
      {
        built_in_type (
          t,
          L"::xsd::cxx::tree::idrefs< " + char_type + L", ",
          simple_type_,
          idref_);

        os << endl;
      }


      // URI.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::AnyURI& t)
      {
        os << "// URI." << endl
           << "//" << endl;

        if (doxygen)
          os << endl;

        uri_ = built_in_type (
          t, L"::xsd::cxx::tree::uri< " + char_type + L", ", simple_type_);

        os << endl;
      }

      // Qualified name.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::QName& t)
      {
        os << "// Qualified name." << endl
           << "//" << endl;

        if (doxygen)
          os << endl;

        built_in_type (
          t,
          L"::xsd::cxx::tree::qname< " + char_type + L", ",
          simple_type_,
          uri_,
          ncname_);

        os << endl;
      }

      // Binary.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Base64Binary& t)
      {
        os << "// Binary." << endl
           << "//" << endl;

        if (doxygen)
          os << endl
             << "/**" << endl
             << " * @brief Binary buffer type." << endl
             << " */" << endl;

        gen_typedef (xs_ns ().context().get<String> ("buffer"),
                     L"::xsd::cxx::tree::buffer< " + char_type + L" >");

        if (doxygen)
          os << endl;

        built_in_type (
          t,
          L"::xsd::cxx::tree::base64_binary< " + char_type + L", ",
          simple_type_);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::HexBinary& t)
      {
        built_in_type (
          t,
          L"::xsd::cxx::tree::hex_binary< " + char_type + L", ",
          simple_type_);

        os << endl;
      }


      // Date/time.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Date& t)
      {
        os << "// Date/time." << endl
           << "//" << endl;

        if (doxygen)
          os << endl
             << "/**" << endl
             << " * @brief Time zone type." << endl
             << " */" << endl;

        gen_typedef (xs_ns ().context().get<String> ("time-zone"),
                     "::xsd::cxx::tree::time_zone");

        if (doxygen)
          os << endl;

        built_in_type (
          t, L"::xsd::cxx::tree::date< " + char_type + L", ", simple_type_);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::DateTime& t)
      {
        built_in_type (
          t,
          L"::xsd::cxx::tree::date_time< " + char_type + L", ",
          simple_type_);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Duration& t)
      {
        built_in_type (
          t,
          L"::xsd::cxx::tree::duration< " + char_type + L", ",
          simple_type_);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Day& t)
      {
        built_in_type (
          t, L"::xsd::cxx::tree::gday< " + char_type + L", ", simple_type_);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Month& t)
      {
        built_in_type (
          t, L"::xsd::cxx::tree::gmonth< " + char_type + L", ", simple_type_);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::MonthDay& t)
      {
        built_in_type (
          t,
          L"::xsd::cxx::tree::gmonth_day< " + char_type + L", ",
          simple_type_);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Year& t)
      {
        built_in_type (
          t, L"::xsd::cxx::tree::gyear< " + char_type + L", ", simple_type_);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::YearMonth& t)
      {
        built_in_type (
          t,
          L"::xsd::cxx::tree::gyear_month< " + char_type + L", ",
          simple_type_);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Time& t)
      {
        built_in_type (
          t, L"::xsd::cxx::tree::time< " + char_type + L", ", simple_type_);

        os << endl;
      }

      // Entity.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Entity& t)
      {
        os << "// Entity." << endl
           << "//" << endl;

        if (doxygen)
          os << endl;

        entity_ = built_in_type (
          t, L"::xsd::cxx::tree::entity< " + char_type + L", ", ncname_);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Entities& t)
      {
        built_in_type (
          t,
          L"::xsd::cxx::tree::entities< " + char_type + L", ",
          simple_type_,
          entity_);

        os << endl;
      }

      virtual Void
      post (SemanticGraph::Namespace& n)
      {
        SemanticGraph::Context& c (xs_ns ().context());

        Boolean parsing (!options.value<CLI::suppress_parsing> ());
        Boolean serialization (options.value<CLI::generate_serialization> ());
        Boolean element_map (options.value<CLI::generate_element_map> ());

        if (options.value<CLI::generate_element_type> ())
        {
          if (doxygen)
            os << "/**" << endl
               << " * @brief Base class for element types." << endl
               << " */" << endl;
          else
            os << "// Base class for element types." << endl
               << "//" << endl;

          gen_typedef (
            c.get<String> ("element-type"),
            L"::xsd::cxx::tree::element_type< " + char_type + L", ",
            type_);

          os << endl;
        }

        if (element_map)
        {
          if (doxygen)
            os << "/**" << endl
               << " * @brief Root element map." << endl
               << " */" << endl;
          else
            os << "// Root element map." << endl
               << "//" << endl;

          gen_typedef (
            c.get<String> ("element-map"),
            L"::xsd::cxx::tree::element_map< " + char_type + L", ",
            type_);

          os << endl;
        }

        if (serialization)
        {
          os << "// Namespace information and list stream. Used in" << endl
             << "// serialization functions." << endl
             << "//" << endl;

          if (doxygen)
            os << "/**" << endl
               << " * @brief Namespace serialization information." << endl
               << " */" << endl;

          gen_typedef (
            c.get<String> ("namespace-info"),
            L"::xsd::cxx::xml::dom::namespace_info< " + char_type + L" >");

            if (doxygen)
              os << endl
                 << "/**" << endl
                 << " * @brief Namespace serialization information map." << endl
                 << " */" << endl;

            gen_typedef (c.get<String> ("namespace-infomap"),
                         L"::xsd::cxx::xml::dom::namespace_infomap< " +
                         char_type + L" >");

            if (doxygen)
              os << endl
                 << "/**" << endl
                 << " * @brief List serialization stream." << endl
                 << " */" << endl;

            gen_typedef (
              c.get<String> ("list-stream"),
              L"::xsd::cxx::tree::list_stream< " + char_type + L" >");

          if (doxygen)
            os << endl
               << "/**" << endl
               << " * @brief Serialization wrapper for the %double type." << endl
               << " */" << endl;

          // Do not export as_double and as_decimal since they are already
          // instantiated.
          //
          gen_typedef (c.get<String> ("as-double"),
                       L"::xsd::cxx::tree::as_double< ",
                       double_,
                       "",
                       "",
                       false);

          if (doxygen)
            os << endl
               << "/**" << endl
               << " * @brief Serialization wrapper for the %decimal type." << endl
               << " */" << endl;

          gen_typedef (c.get<String> ("as-decimal"),
                       L"::xsd::cxx::tree::as_decimal< ",
                       decimal_,
                       "",
                       "",
                       false);

          if (doxygen)
            os << endl
               << "/**" << endl
               << " * @brief Simple type facet." << endl
               << " */" << endl;

          gen_typedef (c.get<String> ("facet"), "::xsd::cxx::tree::facet");

          os << endl;
        }

        //@@ Can't change names of ostream/istream since they are
        //   templates.
        //
        if (!options.value<CLI::generate_insertion> ().empty ())
        {
          if (doxygen)
            os << "/**" << endl
               << " * @brief Data representation output stream template." << endl
               << " */" << endl;
          else
            os << "// Data representation output stream template." << endl
               << "//" << endl;

          os << "using ::xsd::cxx::tree::ostream;"
             << endl;
        }

        if (!options.value<CLI::generate_extraction> ().empty ())
        {
          if (doxygen)
            os << "/**" << endl
               << " * @brief Data representation input stream template." << endl
               << " */" << endl;
          else
            os << "// Data representation input stream template." << endl
               << "//" << endl;

          os << "using ::xsd::cxx::tree::istream;"
             << endl;
        }

        os << "// Flags and properties." << endl
           << "//" << endl;

        if (doxygen)
          os << endl
             << "/**" << endl
             << " * @brief Parsing and serialization flags." << endl
             << " */" << endl;

        gen_typedef (c.get<String> ("flags"), "::xsd::cxx::tree::flags");

        if (doxygen)
          os << endl
             << "/**" << endl
             << " * @brief Parsing properties." << endl
             << " */" << endl;

        gen_typedef (c.get<String> ("properties"),
                     L"::xsd::cxx::tree::properties< " + char_type + L" >");
        os << endl;


        //
        //
        if (parsing || serialization)
        {
          os << "// Parsing/serialization diagnostics." << endl
             << "//" << endl;

          if (doxygen)
            os << endl
               << "/**" << endl
               << " * @brief Error severity." << endl
               << " */" << endl;

          gen_typedef (c.get<String> ("severity"),
                       "::xsd::cxx::tree::severity");

          if (doxygen)
            os << endl
               << "/**" << endl
               << " * @brief Error condition." << endl
               << " */" << endl;

          gen_typedef (c.get<String> ("error"),
                       L"::xsd::cxx::tree::error< " + char_type + L" >");

          if (doxygen)
            os << endl
               << "/**" << endl
               << " * @brief List of %error conditions." << endl
               << " */" << endl;

          gen_typedef (c.get<String> ("diagnostics"),
                       L"::xsd::cxx::tree::diagnostics< " + char_type + L" >");
          os << endl;
        }

        //
        //
        os << "// Exceptions." << endl
           << "//" << endl;

        if (doxygen)
          os << endl
             << "/**" << endl
             << " * @brief Root of the C++/Tree %exception hierarchy." << endl
             << " */" << endl;

        gen_typedef (c.get<String> ("exception"),
                     L"::xsd::cxx::tree::exception< " + char_type + L" >");

        if (doxygen)
          os << endl
             << "/**" << endl
             << " * @brief Exception indicating that the size argument exceeds" <<
            endl
             << " * the capacity argument." << endl
             << " */" << endl;

        gen_typedef (c.get<String> ("bounds"),
                     L"::xsd::cxx::tree::bounds< " + char_type + L" >");

        if (doxygen)
          os << endl
             << "/**" << endl
             << " * @brief Exception indicating that a duplicate ID value" <<
            endl
             << " * was encountered in the object model." << endl
             << " */" << endl;

        gen_typedef (c.get<String> ("duplicate-id"),
                     L"::xsd::cxx::tree::duplicate_id< " + char_type + L" >");

        if (parsing)
        {
          if (doxygen)
            os << endl
               << "/**" << endl
               << " * @brief Exception indicating a parsing failure." << endl
               << " */" << endl;

          gen_typedef (c.get<String> ("parsing"),
                       L"::xsd::cxx::tree::parsing< " + char_type + L" >");

          if (doxygen)
            os << endl
               << "/**" << endl
               << " * @brief Exception indicating that an expected element" <<
              endl
               << " * was not encountered." << endl
               << " */" << endl;

          gen_typedef (
            c.get<String> ("expected-element"),
            L"::xsd::cxx::tree::expected_element< " + char_type + L" >");
        }

        if (parsing || serialization)
        {
          if (doxygen)
            os << endl
               << "/**" << endl
               << " * @brief Exception indicating that an unexpected " <<
              "element" << endl
               << " * was encountered." << endl
               << " */" << endl;

          gen_typedef (
            c.get<String> ("unexpected-element"),
            L"::xsd::cxx::tree::unexpected_element< " + char_type + L" >");
        }

        if (parsing)
        {
          if (doxygen)
            os << endl
               << "/**" << endl
               << " * @brief Exception indicating that an expected " <<
              "attribute" << endl
               << " * was not encountered." << endl
               << " */" << endl;

          gen_typedef (
            c.get<String> ("expected-attribute"),
            L"::xsd::cxx::tree::expected_attribute< " + char_type + L" >");

          if (doxygen)
            os << endl
               << "/**" << endl
               << " * @brief Exception indicating that an unexpected " <<
              "enumerator" << endl
               << " * was encountered." << endl
               << " */" << endl;

          gen_typedef (
            c.get<String> ("unexpected-enumerator"),
            L"::xsd::cxx::tree::unexpected_enumerator< " + char_type + L" >");

          if (doxygen)
            os << endl
               << "/**" << endl
               << " * @brief Exception indicating that the text content " <<
              "was" << endl
               << " * expected for an element." << endl
               << " */" << endl;

          gen_typedef (
            c.get<String> ("expected-text-content"),
            L"::xsd::cxx::tree::expected_text_content< " + char_type + L" >");

          if (doxygen)
            os << endl
               << "/**" << endl
               << " * @brief Exception indicating that a prefix-namespace" <<
              endl
               << " * mapping was not provided." << endl
               << " */" << endl;

          gen_typedef (
            c.get<String> ("no-prefix-mapping"),
            L"::xsd::cxx::tree::no_prefix_mapping< " + char_type + L" >");
        }

        if (options.value<CLI::generate_polymorphic> ())
        {
          if (parsing || serialization)
          {
            if (doxygen)
              os << endl
                 << "/**" << endl
                 << " * @brief Exception indicating that the type " <<
                "information" << endl
                 << " * is not available for a type." << endl
                 << " */" << endl;

            gen_typedef (
              c.get<String> ("no-type-info"),
              L"::xsd::cxx::tree::no_type_info< " + char_type + L" >");
          }

          if (parsing)
          {
            if (doxygen)
              os << endl
                 << "/**" << endl
                 << " * @brief Exception indicating that the types are not" <<
                endl
                 << " * related by inheritance." << endl
                 << " */" << endl;

            gen_typedef (
              c.get<String> ("not-derived"),
              L"::xsd::cxx::tree::not_derived< " + char_type + L" >");
          }
        }

        if (element_map)
        {
          if (doxygen)
            os << endl
               << "/**" << endl
               << " * @brief Exception indicating that parsing or " <<
              "serialization" << endl
               << " * information is not available for an element." << endl
               << " */" << endl;

          gen_typedef (
            c.get<String> ("no-element-info"),
            L"::xsd::cxx::tree::no_element_info< " + char_type + L" >");
        }

        if (serialization)
        {
          if (doxygen)
            os << endl
               << "/**" << endl
               << " * @brief Exception indicating a serialization " <<
              "failure." << endl
               << " */" << endl;

          gen_typedef (
            c.get<String> ("serialization"),
            L"::xsd::cxx::tree::serialization< " + char_type + L" >");
        }

        os << endl;

        if (parsing || serialization)
        {
          if (doxygen)
            os << "/**" << endl
               << " * @brief Error handler callback interface." << endl
               << " */" << endl;
          else
            os << "// Error handler callback interface." << endl
               << "//" << endl;

          gen_typedef (
            c.get<String> ("error-handler"),
            L"::xsd::cxx::xml::error_handler< " + char_type + L" >");

          os << endl;
        }

        if (parsing || serialization)
        {
          if (doxygen)
            os << "/**" << endl
               << " * @brief DOM interaction." << endl
               << " */" << endl;
          else
            os << "// DOM interaction." << endl
               << "//" << endl;

          os << "namespace dom"
             << "{";

          // @@ Disregarding current naming convention by using the
          //    fixed name (no template typedef).
          //
          if (doxygen)
            os << "/**" << endl
               << " * @brief Automatic pointer for DOMDocument." << endl
               << " */" << endl;
          else
            os << "// Automatic pointer for DOMDocument." << endl
               << "//" << endl;

          os << "using ::xsd::cxx::xml::dom::auto_ptr;"
             << endl;

          if (parsing)
          {
            if (!generate_xml_schema)
            {
              String g (L"XSD_CXX_TREE_TREE_NODE_KEY" + ns_name (n));

              std::transform (g.begin (), g.end(), g.begin (), upcase);
              g = escape (g); // Make it a C++ id.

              os << "#ifndef " << g << endl
                 << "#define " << g << endl;
            }

            if (doxygen)
              os << "/**" << endl
                 << " * @brief DOM user data key for back pointers to " <<
                "tree nodes." << endl
                 << " */" << endl;
            else
              os << "// DOM user data key for back pointers to tree nodes." << endl
                 << "//" << endl;

            os << "const XMLCh* const " << c.get<String> ("tree-node-key") <<
              " = ::xsd::cxx::tree::user_data_keys::node;";

            if (!generate_xml_schema)
              os << "#endif" << endl;
          }

          os << "}"; // namespace dom
        }

        if (element_map)
        {
          if (doxygen)
            os << "//@cond" << endl
               << endl;

          if (!generate_xml_schema)
          {
            String g (L"XSD_CXX_TREE_ELEMENT_MAP_INIT" + ns_name (n));

            std::transform (g.begin (), g.end(), g.begin (), upcase);
            g = escape (g); // Make it a C++ id.

            os << "#ifndef " << g << endl
               << "#define " << g << endl;
          }

          os << "static" << endl
             << "const ::xsd::cxx::tree::element_map_init< " <<
            char_type << ", " << type_ << " >" << endl
             << "_xsd_element_map_init;";

          if (!generate_xml_schema)
            os << "#endif" << endl;

          if (doxygen)
            os << endl
               << "//@endcond" << endl;
        }

        Namespace::post (n);

        // Generate exports.
        //
        if (export_)
        {
          StringSet ns_set;

          for (StringList::ConstIterator i (exports_.begin ());
               i != exports_.end (); ++i)
          {
            String const& e (*i);

            // 12 is to skip ::xsd::cxx::
            //
            ns_set.insert (String (e, 12, e.rfind (':', e.find ('<')) - 13));
          }

          os << "#ifndef XSD_NO_EXPORT" << endl
             << endl
             << "namespace xsd"
             << "{"
             << "namespace cxx"
             << "{";

          for (StringSet::ConstIterator i (ns_set.begin ());
               i != ns_set.end (); ++i)
          {
            String const& ns (*i);
            String prefix (L"::xsd::cxx::" + ns);

            Size n (1);
            for (Size b (0), e (ns.find (':')); ; n++)
            {
              os << "namespace " << String (ns, b, e)
                 << "{";

              if (e == String::npos)
                break;

              b = e + 2;
              e = ns.find (':', b);
            }

            for (StringList::ConstIterator i (exports_.begin ());
                 i != exports_.end (); ++i)
            {
              String const& e (*i);
              String ens (e, 12, e.rfind (':', e.find ('<')) - 13);

              if (ns == ens)
              {
                String type (e, e.rfind (':', e.find ('<')) + 1);
                os << "template class " << type_exp << type << ";";
              }
            }

            while (n--)
              os << "}";
          }

          os << "}"  // cxx
             << "}"  // xsd
             << "#endif // XSD_NO_EXPORT" << endl
             << endl;
        }
      }

    private:
      typedef Cult::Containers::Set<String> StringSet;
      typedef Cult::Containers::Vector<String> StringList;

      Boolean export_;
      StringList exports_;
      StringSet exports_set_;
      String xs_ns_;

      Traversal::Names names_;

      String type_;
      String simple_type_;
      String string_;
      String norm_string_;
      String token_;
      String nmtoken_;
      String name_;
      String ncname_;
      String idref_;
      String uri_;
      String entity_;

      String double_;
      String decimal_;
    };
  }
}

#endif // CXX_TREE_FUNDAMENTAL_HEADER_HXX
