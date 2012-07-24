// file      : xsd/cxx/parser/parser-header.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/parser/parser-header.hxx>

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

namespace CXX
{
  namespace Parser
  {
    namespace
    {
      struct Enumeration: Traversal::Enumeration, Context
      {
        Enumeration (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& e)
        {
          String const& name (ename (e));
          SemanticGraph::Type& base (e.inherits ().base ());

          os << "class " << type_exp << name << ": public virtual " <<
            fq_name (base)
             << "{"
             << "public:" << endl
             << "// Parser callbacks. Override them in your " <<
            "implementation." << endl
             << "//" << endl;

          os << "// virtual void" << endl
             << "// pre ();" << endl
             << endl;

          String const& ret (ret_type (e));

          Boolean same (ret == ret_type (base));

          os << "virtual " << ret << endl
             << post_name (e) << " ()" <<
            (same || ret == L"void" ? ";" : " = 0;");

          if (polymorphic)
          {
            os << endl
               << "public:" << endl
               << "static const " << char_type << "*" << endl
               << "_static_type ();"
               << endl
               << "virtual const " << char_type << "*" << endl
               << "_dynamic_type () const;";
          }

          os << "};";
        }
      };

      //
      //
      struct List: Traversal::List, Context
      {
        List (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& l)
        {
          String const& name (ename (l));
          SemanticGraph::Type& t (l.argumented ().type ());

          String item (unclash (name, "item"));

          os << "class " << type_exp << name << ": public " << list_base
             << "{"
             << "public:" << endl
             << "// Parser callbacks. Override them in your " <<
            "implementation." << endl
             << "//" << endl;

          // pre
          //
          os << "// virtual void" << endl
             << "// pre ();" << endl
             << endl;

          // item
          //
          String const& arg (arg_type (t));

          os << "virtual void" << endl
             << item;

          if (arg == L"void")
            os << " ();";
          else
            os << " (" << arg << ");";

          os << endl;

          // post
          //
          String const& ret (ret_type (l));

          os << "virtual " << ret << endl
             << post_name (l) << " ()" << (ret == L"void" ? ";" : " = 0;")
             << endl;

          //
          //
          os << "// Parser construction API." << endl
             << "//" << endl;

          // item_parser
          //
          os << "void" << endl
             << unclash (name, "item_parser") << " (" << fq_name (t) << "&);"
             << endl;

          // parsers
          //
          os << "void" << endl
             << "parsers (" << fq_name (t) << "& /* item */);"
             << endl;

          // c-tor
          //
          os << "// Constructor." << endl
             << "//" << endl
             << name << " ();"
             << endl;


          if (polymorphic)
          {
            os << "public:" << endl
               << "static const " << char_type << "*" << endl
               << "_static_type ();"
               << endl
               << "virtual const " << char_type << "*" << endl
               << "_dynamic_type () const;"
               << endl;
          }

          //
          //
          os << "// Implementation." << endl
             << "//" << endl
             << "protected:" << endl;

          os << "virtual void" << endl
             << "_xsd_parse_item (const " << string_type << "&);"
             << endl;

          os << "protected:" << endl
             << fq_name (t) << "* _xsd_" << item << "_;"
             << "};";
        }
      };

      struct Union: Traversal::Union, Context
      {
        Union (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& u)
        {
          String const& name (ename (u));

          os << "class " << type_exp << name << ": public " << simple_base
             << "{"
             << "public:" << endl
             << "// Parser callbacks. Override them in your " <<
            "implementation." << endl
             << "//" << endl;

          os << "// virtual void" << endl
             << "// pre ();" << endl
             << "//" << endl
             << "// virtual void" << endl
             << "// _characters (const " << string_type << "&);" << endl
             << endl;

          String const& ret (ret_type (u));

          os << "virtual " << ret << endl
             << post_name (u) << " ()" << (ret == L"void" ? ";" : " = 0;");

          if (polymorphic)
          {
            os << endl
               << "public:" << endl
               << "static const " << char_type << "*" << endl
               << "_static_type ();"
               << endl
               << "virtual const " << char_type << "*" << endl
               << "_dynamic_type () const;";
          }

          os << "};";
        }
      };


      //
      //
      struct ParserCallback: Traversal::Member, Context
      {
        ParserCallback (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& m)
        {
	  if (skip (m)) return;

          String const& arg (arg_type (m.type ()));

          os << "virtual void" << endl
             << ename (m);

          if (arg == L"void")
            os << " ();";
          else
            os << " (" << arg << ");";

          os << endl;
        }
      };


      //
      //
      struct ParserModifier: Traversal::Member, Context
      {
        ParserModifier (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& m)
        {
	  if (skip (m))
            return;

          os << "void" << endl
             << eparser (m) << " (" << fq_name (m.type ()) << "&);"
             << endl;

          if (polymorphic &&
              m.is_a<SemanticGraph::Element> () &&
              !anonymous (m.type ()))
          {
            os << "void" << endl
               << eparser (m) << " (const " << parser_map << "&);"
               << endl;
          }
        }
      };

      //
      //
      struct ParserMember: Traversal::Member, Context
      {
        ParserMember (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& m)
        {
	  if (skip (m))
            return;

          String type (fq_name (m.type ()));

          os << type << "* " << emember (m) << ";";

          if (polymorphic &&
              m.is_a<SemanticGraph::Element> () &&
              !anonymous (m.type ()))
          {
            os << type << "* " << emember_cache (m) << ";"
               << "const " << parser_map << "* " << emember_map (m) << ";"
               << endl;
          }
        }
      };

      //
      //
      struct Particle: Traversal::All,
                       Traversal::Choice,
                       Traversal::Sequence,
                       Context
      {
        Particle (Context& c)
            : Context (c)
        {
          *this >> contains_particle_ >> *this;
        }


        virtual Void
        traverse (SemanticGraph::All& a)
        {
          if (!a.context().count ("comp-number"))
            return;

          UnsignedLong state_count (
            a.context().get<UnsignedLong> ("state-count"));

          os << "void" << endl
             << "all_0 (unsigned long& state," << endl
             << "unsigned char* count," << endl
             << "const " << string_type << "& ns," << endl
             << "const " << string_type << "& n," << endl
             << "const " << string_type << "* t," << endl
             << "bool start);"
             << endl
             << "unsigned char v_all_first_[" << state_count << "UL];"
             << "::xsd::cxx::parser::validating::all_stack v_all_count_;"
             << endl;
        }

        virtual Void
        traverse (SemanticGraph::Choice& c)
        {
          if (!c.context().count ("comp-number"))
            return;

          UnsignedLong n (c.context ().get<UnsignedLong> ("comp-number"));

          os << "void" << endl
             << "choice_" << n << " (unsigned long& state," << endl
             << "unsigned long& count," << endl
             << "const " << string_type << "& ns," << endl
             << "const " << string_type << "& n," << endl
             << "const " << string_type << "* t," << endl
             << "bool start);"
             << endl;

          Traversal::Choice::traverse (c);
        }

        virtual Void
        traverse (SemanticGraph::Sequence& s)
        {
          if (!s.context().count ("comp-number"))
            return;

          UnsignedLong n (s.context ().get<UnsignedLong> ("comp-number"));

          os << "void" << endl
             << "sequence_" << n << " (unsigned long& state," << endl
             << "unsigned long& count," << endl
             << "const " << string_type << "& ns," << endl
             << "const " << string_type << "& n," << endl
             << "const " << string_type << "* t," << endl
             << "bool start);"
             << endl;

          Traversal::Sequence::traverse (s);
        }

      private:
        Traversal::ContainsParticle contains_particle_;
      };


      //
      //
      struct AttributeValidationState: Traversal::Attribute, Context
      {
        AttributeValidationState (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& a)
        {
          if (!a.optional_p ())
          {
            os << "bool " << ename (a) << ";";
          }
        }
      };

      //
      //
      struct Complex: Traversal::Complex, Context
      {
        Complex (Context& c)
            : Context (c),
              parser_callback_ (c),
              parser_member_ (c),
              parser_modifier_ (c),
              attribute_validation_state_ (c)
        {
          names_parser_callback_ >> parser_callback_;
          names_parser_member_ >> parser_member_;
          names_parser_modifier_ >> parser_modifier_;
          names_attribute_validation_state_ >> attribute_validation_state_;
        }

        virtual Void
        traverse (Type& c)
        {
          String const& name (ename (c));

          // In case of an inheritance-by-restriction, we don't need to
          // generate parser callbacks, etc. since they are the same as in
          // the base. We only need the parsing/validation code.
          //
          Boolean restriction (restriction_p (c));

          Boolean he (has<Traversal::Element> (c));
          Boolean ha (has<Traversal::Attribute> (c));

          Boolean hae (has_particle<Traversal::Any> (c));
          Boolean haa (has<Traversal::AnyAttribute> (c));

          Boolean hra (false); // Has required attribute.
          if (ha)
          {
            RequiredAttributeTest test (hra);
            Traversal::Names names_test (test);
            names (c, names_test);
          }


          //
          //
          os << "class " << type_exp << name << ": public ";

          if (c.inherits_p ())
            os << "virtual " << fq_name (c.inherits ().base ());
          else
            os << complex_base;

          os << "{"
             << "public:" << endl
             << "// Parser callbacks. Override them in your " <<
            "implementation." << endl
             << "//" << endl;

          os << "// virtual void" << endl
             << "// pre ();" << endl
             << endl;


          if (!restriction && (ha || he))
          {
            names (c, names_parser_callback_);
          }

          String const& ret (ret_type (c));

          Boolean same (c.inherits_p () &&
                        ret == ret_type (c.inherits ().base ()));

          os << "virtual " << ret << endl
             << post_name (c) << " ()" <<
            (same || ret == L"void" ? ";" : " = 0;")
             << endl;

          //
          //
          if (!restriction && (he || ha))
          {
            os << "// Parser construction API." << endl
               << "//" << endl;

            names (c, names_parser_modifier_);

            os << "void" << endl
               << "parsers (";

            {
              ParserParamDecl decl (*this, false);
              decl.traverse (c);
            }

            os << ");"
               << endl;
          }

          // Default c-tor.
          //
          if ((!restriction && (he || ha)) ||
              (validation && (he || hae || hra)))
          {
            os << "// Constructor." << endl
               << "//" << endl
               << name << " ();"
               << endl;
          }

          if (polymorphic)
          {
            os << "public:" << endl
               << "static const " << char_type << "*" << endl
               << "_static_type ();"
               << endl
               << "virtual const " << char_type << "*" << endl
               << "_dynamic_type () const;"
               << endl;
          }

          // Implementation.
          //
          if (he || ha || (validation && (hae || haa)))
          {
            os << "// Implementation." << endl
               << "//" << endl
               << "protected:" << endl;
          }

          // element
          //
          if (he || (validation && hae))
          {
            // _start_element_impl
            //
            os << "virtual bool" << endl
               << "_start_element_impl (const " << string_type << "&," << endl
               << "const " << string_type << "&," << endl
               << "const " << string_type << "*);"
               << endl;

            // end_element
            //
            os << "virtual bool" << endl
               << "_end_element_impl (const " << string_type << "&," << endl
               << "const " << string_type << "&);"
               << endl;
          }

          // attribute
          //
          if (validation)
          {
            if (ha)
            {
              os << "virtual bool" << endl
                 << "_attribute_impl_phase_one (const " << string_type <<
                "&," << endl
                 << "const " << string_type << "&," << endl
                 << "const " << string_type << "&);" << endl
                 << endl;
            }

            if (haa)
            {
              os << "virtual bool" << endl
                 << "_attribute_impl_phase_two (const " << string_type <<
                "&," << endl
                 << "const " << string_type << "&," << endl
                 << "const " << string_type << "&);"
                 << endl;
            }
          }
          else
          {
            if (ha)
            {
              os << "virtual bool" << endl
                 << "_attribute_impl (const " << string_type << "&," << endl
                 << "const " << string_type << "&," << endl
                 << "const " << string_type << "&);"
                 << endl;
            }
          }

          // characters
          //
          if (validation && c.mixed_p ())
          {
            os << "virtual bool" << endl
               << "_characters_impl (const " << string_type << "&);"
               << endl;
          }

          if (!restriction && (he || ha))
          {
            os << "protected:" << endl;
            names (c, names_parser_member_);
            os << endl;
          }

          if (validation && (he || hae))
          {
            UnsignedLong depth (c.context ().get<UnsignedLong> ("depth"));

            os << "protected:" << endl;

            os << "struct v_state_descr_"
               << "{"
               << "void (" << fq_name (c) << "::*func) (" << endl
               << "unsigned long&," << endl
               << "unsigned long&," << endl
               << "const " << string_type << "&," << endl
               << "const " << string_type << "&," << endl
               << "const " << string_type << "*," << endl
               << "bool);"
               << "unsigned long state;"
               << "unsigned long count;"
               << "};";

            // Allocate one extra slot for the special state.
            //
            os << "struct v_state_"
               << "{"
               << "v_state_descr_ data[" << depth + 1 << "UL];"
               << "unsigned long size;"
               << "};";

            os << "v_state_ v_state_first_;"
               << "::xsd::cxx::parser::pod_stack v_state_stack_;"
               << endl;

            os << "virtual void" << endl
               << "_pre_e_validate ();"
               << endl;

            os << "virtual void" << endl
               << "_post_e_validate ();"
               << endl;

            Particle t (*this);
            t.dispatch (c.contains_compositor ().compositor ());
          }

          if (validation && hra)
          {
            os << "protected:" << endl;

            os << "struct v_state_attr_"
               << "{";

            names (c, names_attribute_validation_state_);

            os << "};";

            os << "v_state_attr_ v_state_attr_first_;"
               << "::xsd::cxx::parser::pod_stack v_state_attr_stack_;"
               << endl;

            os << "virtual void" << endl
               << "_pre_a_validate ();"
               << endl;

            os << "virtual void" << endl
               << "_post_a_validate ();"
               << endl;
          }

          os << "};";
        }

      private:
        //
        //
        ParserCallback parser_callback_;
        Traversal::Names names_parser_callback_;

        //
        //
        ParserMember parser_member_;
        Traversal::Names names_parser_member_;

        //
        //
        ParserModifier parser_modifier_;
        Traversal::Names names_parser_modifier_;

        //
        //
        AttributeValidationState attribute_validation_state_;
        Traversal::Names names_attribute_validation_state_;
      };

      struct FundType : Context,

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
        FundType (Context& c)
            : Context (c), xs_ns_ (xs_ns_name ())
        {
          impl_ns_ = "::xsd::cxx::parser::";
          impl_ns_ += (validation ? L"validating" : L"non_validating");

          if (char_type == L"char")
            string_type_ = L"::std::string";
          else if (char_type == L"wchar_t")
            string_type_ = L"::std::wstring";
          else
            string_type_ = L"::std::basic_string< " + char_type + L" >";
        }

        // anyType & anySimpleType.
        //
        virtual Void
        traverse (SemanticGraph::AnyType& t)
        {
          gen_typedef (t, "void", "any_type_pskel", "any_type_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::AnySimpleType& t)
        {
          gen_typedef (t, "void",
                       "any_simple_type_pskel", "any_simple_type_pimpl");
        }

        // Boolean.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Boolean& t)
        {
          gen_typedef (t, "bool", "boolean_pskel", "boolean_pimpl");
        }

        // Integral types.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Byte& t)
        {
          gen_typedef (t, "signed char", "byte_pskel", "byte_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::UnsignedByte& t)
        {
          gen_typedef (t, "unsigned char",
                       "unsigned_byte_pskel", "unsigned_byte_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Short& t)
        {
          gen_typedef (t, "short", "short_pskel", "short_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::UnsignedShort& t)
        {
          gen_typedef (t, "unsigned short",
                       "unsigned_short_pskel", "unsigned_short_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Int& t)
        {
          gen_typedef (t, "int", "int_pskel", "int_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::UnsignedInt& t)
        {
          gen_typedef (t, "unsigned int",
                       "unsigned_int_pskel", "unsigned_int_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Long& t)
        {
          gen_typedef (t, "long long", "long_pskel", "long_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::UnsignedLong& t)
        {
          gen_typedef (t, "unsigned long long",
                       "unsigned_long_pskel", "unsigned_long_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Integer& t)
        {
          gen_typedef (t, "long long", "integer_pskel", "integer_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NegativeInteger& t)
        {
          gen_typedef (t, "long long",
                       "negative_integer_pskel", "negative_integer_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NonPositiveInteger& t)
        {
          gen_typedef (t, "long long",
                       "non_positive_integer_pskel",
                       "non_positive_integer_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::PositiveInteger& t)
        {
          gen_typedef (t, "unsigned long long",
                       "positive_integer_pskel", "positive_integer_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NonNegativeInteger& t)
        {
          gen_typedef (t, "unsigned long long",
                       "non_negative_integer_pskel",
                       "non_negative_integer_pimpl");
        }

        // Floats.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Float& t)
        {
          gen_typedef (t, "float", "float_pskel", "float_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Double& t)
        {
          gen_typedef (t, "double", "double_pskel", "double_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Decimal& t)
        {
          gen_typedef (t, "double", "decimal_pskel", "decimal_pimpl");
        }

        // Strings.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::String& t)
        {
          gen_typedef (t, string_type_, "string_pskel", "string_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NormalizedString& t)
        {
          gen_typedef (t, string_type_,
                       "normalized_string_pskel", "normalized_string_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Token& t)
        {
          gen_typedef (t, string_type_, "token_pskel", "token_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NameToken& t)
        {
          nmtoken_ = gen_typedef (t, string_type_,
                                  "nmtoken_pskel", "nmtoken_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NameTokens& t)
        {
          // NMTOKENS uses NMTOKEN implementation to parse individual items.
          // As a result, we don't generate NMTOKENS if we didn't generate
          // NMTOKEN. Here we assume NMTOKEN is handled before NMTOKENS.
          //
          if(nmtoken_)
            gen_typedef (t, xs_ns_ + L"::string_sequence",
                         "nmtokens_pskel", "nmtokens_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Name& t)
        {
          gen_typedef (t, string_type_, "name_pskel", "name_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NCName& t)
        {
          gen_typedef (t, string_type_, "ncname_pskel", "ncname_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Language& t)
        {
          gen_typedef (t, string_type_, "language_pskel", "language_pimpl");
        }

        // Qualified name.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::QName& t)
        {
          gen_typedef (t, xs_ns_ + L"::qname", "qname_pskel", "qname_pimpl");
        }

        // ID/IDREF.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Id& t)
        {
          gen_typedef (t, string_type_, "id_pskel", "id_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::IdRef& t)
        {
          idref_ = gen_typedef (t, string_type_, "idref_pskel", "idref_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::IdRefs& t)
        {
          // IDREFS uses IDREF implementation to parse individual items.
          // As a result, we don't generate IDREFS if we didn't generate
          // IDREF. Here we assume IDREF is handled before IDREFS.
          //
          if (idref_)
            gen_typedef (t, xs_ns_ + L"::string_sequence",
                         "idrefs_pskel", "idrefs_pimpl");
        }

        // URI.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::AnyURI& t)
        {
          gen_typedef (t, string_type_, "uri_pskel", "uri_pimpl");
        }

        // Binary.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Base64Binary& t)
        {
          String buffer (L"::std::auto_ptr< " + xs_ns_ + L"::buffer >");
          gen_typedef (t, buffer,
                       "base64_binary_pskel", "base64_binary_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::HexBinary& t)
        {
          String buffer (L"::std::auto_ptr< " + xs_ns_ + L"::buffer >");
          gen_typedef (t, buffer, "hex_binary_pskel", "hex_binary_pimpl");
        }


        // Date/time.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Date& t)
        {
          gen_typedef (t, xs_ns_ + L"::date", "date_pskel", "date_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::DateTime& t)
        {
          gen_typedef (t, xs_ns_ + L"::date_time",
                       "date_time_pskel", "date_time_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Duration& t)
        {
          gen_typedef (t, xs_ns_ + L"::duration",
                       "duration_pskel", "duration_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Day& t)
        {
          gen_typedef (t, xs_ns_ + L"::gday", "gday_pskel", "gday_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Month& t)
        {
          gen_typedef (t, xs_ns_ + L"::gmonth",
                       "gmonth_pskel", "gmonth_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::MonthDay& t)
        {
          gen_typedef (t, xs_ns_ + L"::gmonth_day",
                       "gmonth_day_pskel", "gmonth_day_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Year& t)
        {
          gen_typedef (t, xs_ns_ + L"::gyear", "gyear_pskel", "gyear_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::YearMonth& t)
        {
          gen_typedef (t, xs_ns_ + L"::gyear_month",
                       "gyear_month_pskel", "gyear_month_pimpl");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Time& t)
        {
          gen_typedef (t, xs_ns_ + L"::time", "time_pskel", "time_pimpl");
        }

        // Entity.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Entity&)
        {
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Entities&)
        {
        }

      private:
        Boolean
        gen_typedef (SemanticGraph::Type& t,
                     String const& type,
                     String const& pskel,
                     String const& pimpl)
        {
          if (ret_type (t) == type)
          {
            os << "typedef " << impl_ns_ << "::" << pskel << "< " <<
              char_type << " > " << ename (t) << ";";

            String const& pimpl_name (t.context ().get<String> ("impl"));

            os << "typedef " << impl_ns_ << "::" << pimpl << "< " <<
              char_type << " > " << pimpl_name << ";"
               << endl;

            return true;
          }

          return false;
        }

        String xs_ns_;
        String impl_ns_;
        String string_type_;

        Boolean idref_;
        Boolean nmtoken_;
      };

      struct FundNamespace: Namespace, Context
      {
        FundNamespace (Context& c)
            : Namespace (c), Context (c)
        {
        }

        void
        traverse (Type& ns)
        {
          pre (ns);

          String impl ("::xsd::cxx::parser::");
          impl += (validation ? L"validating" : L"non_validating");

          String const c (char_type);

          os << "// Built-in XML Schema types mapping." << endl
             << "//" << endl
             << "typedef ::xsd::cxx::parser::string_sequence< " << c <<
            " > string_sequence;"
             << "typedef ::xsd::cxx::parser::qname< " << c << " > qname;"
             << "typedef ::xsd::cxx::parser::buffer buffer;"
             << "typedef ::xsd::cxx::parser::time_zone time_zone;"
             << "typedef ::xsd::cxx::parser::gday gday;"
             << "typedef ::xsd::cxx::parser::gmonth gmonth;"
             << "typedef ::xsd::cxx::parser::gyear gyear;"
             << "typedef ::xsd::cxx::parser::gmonth_day gmonth_day;"
             << "typedef ::xsd::cxx::parser::gyear_month gyear_month;"
             << "typedef ::xsd::cxx::parser::date date;"
             << "typedef ::xsd::cxx::parser::time time;"
             << "typedef ::xsd::cxx::parser::date_time date_time;"
             << "typedef ::xsd::cxx::parser::duration duration;"
             << endl;

          os << "// Base parser skeletons." << endl
             << "//" << endl
             << "typedef ::xsd::cxx::parser::parser_base< " << c <<
            " > parser_base;"
             << "typedef " << impl << "::empty_content< " << c <<
            " > empty_content;"
             << "typedef " << impl << "::simple_content< " << c <<
            " > simple_content;"
             << "typedef " << impl << "::complex_content< " << c <<
            " > complex_content;"
             << "typedef " << impl << "::list_base< " << c << " > list_base;"
             << endl;

          if (polymorphic)
          {
            os << "// Parser map interface and default implementation." << endl
               << "//" << endl
               << "typedef ::xsd::cxx::parser::parser_map< " << c <<
              " > parser_map;"
               << "typedef ::xsd::cxx::parser::parser_map_impl< " << c <<
              " > parser_map_impl;"
               << endl;
          }

          os << "// Parser skeletons and implementations for the XML Schema" << endl
             << "// built-in types." << endl
             << "//" << endl;

          names (ns);

          os << "// Exceptions. See libxsd/xsd/cxx/parser/exceptions.hxx " <<
            "for details." << endl
             << "//" << endl
             << "typedef ::xsd::cxx::parser::exception< " <<
            char_type << " > exception;"
             << endl
             << "// Parsing diagnostics." << endl
             << "//" << endl
             << "typedef ::xsd::cxx::parser::severity severity;"
             << "typedef ::xsd::cxx::parser::error< " << c << " > error;"
             << "typedef ::xsd::cxx::parser::diagnostics< " << c <<
            " > diagnostics;"
             << "typedef ::xsd::cxx::parser::parsing< " << c << " > parsing;"
             << endl;

          os << "// Error handler. See " <<
            "libxsd/xsd/cxx/xml/error-handler.hxx for details." << endl
             << "//" << endl
             << "typedef ::xsd::cxx::xml::error_handler< " << c <<
            " > error_handler;"
             << endl;

          os << "// Read-only string." << endl
             << "//" << endl
             << "typedef ::xsd::cxx::ro_string< " << c << " > ro_string;"
             << endl;

          if (xml_parser == L"xerces")
          {
            os << "// Parsing flags. See " <<
              "libxsd/xsd/cxx/parser/xerces/elements.hxx" << endl
               << "// for details." << endl
               << "//" << endl
               << "typedef ::xsd::cxx::parser::xerces::flags flags;"
               << endl;

            os << "// Parsing properties. See " <<
              "libxsd/xsd/cxx/parser/xerces/elements.hxx" << endl
               << "// for details." << endl
               << "//" << endl
               << "typedef ::xsd::cxx::parser::xerces::properties< " << c <<
              " > properties;"
               << endl;

            os << "// Document type. See " <<
              "libxsd/xsd/cxx/parser/xerces/elements.hxx" << endl
               << "// for details." << endl
               << "//" << endl
               << "typedef ::xsd::cxx::parser::xerces::document< " << c <<
              " > document;"
               << endl;

          }
          else if (xml_parser == L"expat")
          {
            os << "// Document type. See " <<
              "libxsd/xsd/cxx/parser/expat/elements.hxx" << endl
               << "// for details." << endl
               << "//" << endl
               << "typedef ::xsd::cxx::parser::expat::document< " << c <<
              " > document;"
               << endl;
          }

          post (ns);
        }
      };
    }

    Void
    generate_parser_header (Context& ctx, Boolean generate_xml_schema)
    {
      String c (ctx.char_type);

      //
      //
      if (c == L"char")
      {
        ctx.os << "#ifndef XSD_USE_CHAR" << endl
               << "#define XSD_USE_CHAR" << endl
               << "#endif" << endl
               << endl;

        ctx.os << "#ifndef XSD_CXX_PARSER_USE_CHAR" << endl
               << "#define XSD_CXX_PARSER_USE_CHAR" << endl
               << "#endif" << endl
               << endl;
      }
      else if (c == L"wchar_t")
      {
        ctx.os << "#ifndef XSD_USE_WCHAR" << endl
               << "#define XSD_USE_WCHAR" << endl
               << "#endif" << endl
               << endl;

        ctx.os << "#ifndef XSD_CXX_PARSER_USE_WCHAR" << endl
               << "#define XSD_CXX_PARSER_USE_WCHAR" << endl
               << "#endif" << endl
               << endl;
      }

      //
      //
      NarrowString extern_xml_schema;

      if (!generate_xml_schema)
        extern_xml_schema = ctx.options.value<CLI::extern_xml_schema> ();

      if (extern_xml_schema)
      {
        String name (ctx.hxx_expr->merge (extern_xml_schema));

        ctx.os << "#include " << ctx.process_include_path (name) << endl
               << endl;

        // Generate includes that came from the type map.
        //
        if (ctx.schema_root.context ().count ("includes"))
        {
          typedef Cult::Containers::Set<String> Includes;

          Includes const& is (
            ctx.schema_root.context ().get<Includes> ("includes"));

          for (Includes::ConstReverseIterator i (is.rbegin ());
               i != is.rend (); ++i)
          {
            ctx.os << "#include " << *i << endl;
          }

          ctx.os << endl;
        }
      }
      else
      {
        if (ctx.char_type == L"char" &&
            ctx.xml_parser == L"xerces" &&
            ctx.char_encoding != L"custom")
        {
          ctx.os << "#include <xsd/cxx/xml/char-" << ctx.char_encoding << ".hxx>" << endl;
        }

        ctx.os << "#include <xsd/cxx/xml/error-handler.hxx>" << endl
               << "#include <xsd/cxx/parser/exceptions.hxx>" << endl
               << "#include <xsd/cxx/parser/elements.hxx>" << endl
               << "#include <xsd/cxx/parser/xml-schema.hxx>" << endl;

        if (ctx.polymorphic)
          ctx.os << "#include <xsd/cxx/parser/map.hxx>" << endl;

        if (ctx.validation)
          ctx.os << "#include <xsd/cxx/parser/validating/parser.hxx>" << endl
                 << "#include <xsd/cxx/parser/validating/exceptions.hxx>" << endl
                 << "#include <xsd/cxx/parser/validating/xml-schema-pskel.hxx>" << endl
                 << "#include <xsd/cxx/parser/validating/xml-schema-pimpl.hxx>" << endl;
        else
          ctx.os << "#include <xsd/cxx/parser/non-validating/parser.hxx>" << endl
                 << "#include <xsd/cxx/parser/non-validating/xml-schema-pskel.hxx>" << endl
                 << "#include <xsd/cxx/parser/non-validating/xml-schema-pimpl.hxx>" << endl;

        ctx.os << "#include <xsd/cxx/parser/" << ctx.xml_parser <<
          "/elements.hxx>" << endl
               << endl;

        // Generate includes that came from the type map.
        //
        if (ctx.schema_root.context ().count ("includes"))
        {
          typedef Cult::Containers::Set<String> Includes;

          Includes const& is (
            ctx.schema_root.context ().get<Includes> ("includes"));

          for (Includes::ConstReverseIterator i (is.rbegin ());
               i != is.rend (); ++i)
          {
            ctx.os << "#include " << *i << endl;
          }

          ctx.os << endl;
        }

        // Generate fundamental types.
        //
        if (generate_xml_schema)
        {
          Traversal::Schema schema;
          Traversal::Names names;
          FundNamespace ns (ctx);

          schema >> names >> ns;

          Traversal::Names ns_names;
          FundType type (ctx);

          ns >> ns_names >> type;

          schema.dispatch (ctx.schema_root);
        }
        else
        {
          Traversal::Schema schema, xsd;
          Traversal::Implies implies;
          Traversal::Names names;
          FundNamespace ns (ctx);

          schema >> implies >> xsd >> names >> ns;

          Traversal::Names ns_names;
          FundType type (ctx);

          ns >> ns_names >> type;

          schema.dispatch (ctx.schema_root);
        }
      }

      // Generate user type mapping.
      //
      if (!generate_xml_schema)
      {
        Traversal::Schema schema;

        Traversal::Sources sources;
        Includes includes (ctx, Includes::header);
        Traversal::Names schema_names;

        Namespace ns (ctx);
        Traversal::Names names;

        schema >> includes;
        schema >> sources >> schema;
        schema >> schema_names >> ns >> names;

        List list (ctx);
        Union union_ (ctx);
        Complex complex (ctx);
        Enumeration enumeration (ctx);

        names >> list;
        names >> union_;
        names >> complex;
        names >> enumeration;

        schema.dispatch (ctx.schema_root);
      }
    }
  }
}
