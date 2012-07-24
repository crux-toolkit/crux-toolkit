// file      : xsd/cxx/parser/name-processor.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/parser/name-processor.hxx>

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

#include <cult/containers/set.hxx>

#include <sstream>
#include <iostream>

namespace CXX
{
  namespace Parser
  {
    NameProcessor::
    NameProcessor ()
    {
      // Dummy ctor, helps with long symbols on HP-UX.
    }

    namespace
    {
      //
      //
      typedef Cult::Containers::Set<String> NameSet;

      class Context: public CXX::Context
      {
      public:
        Context (CLI::Options const& ops,
                 SemanticGraph::Schema& root,
                 SemanticGraph::Path const& path,
                 StringLiteralMap const* map)
            : CXX::Context (std::wcerr,
                            root,
                            path,
                            map,
                            ops.value<CLI::char_type> (),
                            ops.value<CLI::char_encoding> (),
                            ops.value<CLI::include_with_brackets> (),
                            ops.value<CLI::include_prefix> (),
                            ops.value<CLI::export_symbol> (),
                            ops.value<CLI::namespace_map> (),
                            ops.value<CLI::namespace_regex> (),
                            ops.value<CLI::namespace_regex_trace> (),
                            ops.value<CLI::include_regex> (),
                            ops.value<CLI::include_regex_trace> (),
                            ops.value<CLI::generate_inline> (),
                            ops.value<CLI::reserved_name> ()),
              skel_suffix_ (ops.value<CLI::skel_type_suffix> ()),
              impl_suffix_ (ops.value<CLI::impl_type_suffix> ()),
              impl (ops.value<CLI::generate_noop_impl> () ||
                    ops.value<CLI::generate_print_impl> () ||
                    ops.value<CLI::generate_test_driver> ()),
              skel_suffix (skel_suffix_),
              impl_suffix (impl_suffix_),
              global_type_names (global_type_names_),
              polymorphic (ops.value<CLI::generate_polymorphic> ())
        {
        }

      protected:
        Context (Context& c)
            : CXX::Context (c),
              impl (c.impl),
              skel_suffix (c.skel_suffix),
              impl_suffix (c.impl_suffix),
              global_type_names (c.global_type_names),
              polymorphic (c.polymorphic)
        {
        }

      public:
        String
        find_name (String const& n, NameSet& set)
        {
          String base_name (escape (n));
          String name (base_name);

          for (UnsignedLong i (1); set.find (name) != set.end (); ++i)
          {
            std::wostringstream os;
            os << i;
            name = base_name + os.str ();
          }

          set.insert (name);
          return name;
        }

      private:
        String const skel_suffix_;
        String const impl_suffix_;

        Cult::Containers::Map<String, NameSet> global_type_names_;

      public:
        Boolean const impl;
        String const& skel_suffix;
        String const& impl_suffix;

        Cult::Containers::Map<String, NameSet>& global_type_names;

        Boolean polymorphic;
      };


      //
      //
      struct PrimaryMember: Traversal::Member, Context
      {
        PrimaryMember (Context& c, NameSet& set)
            : Context (c), set_ (set)
        {
        }

        virtual Void
        traverse (Type& m)
        {
          if (Parser::Context::skip (m))
            return;

          m.context ().set ("name", find_name (m.name (), set_));
        }

      private:
        NameSet& set_;
      };

      struct DerivedMember: Traversal::Member, Context
      {
        DerivedMember (Context& c, NameSet& set)
            : Context (c), set_ (set)
        {
        }

        virtual Void
        traverse (Type& m)
        {
          if (Parser::Context::skip (m))
            return;

          m.context ().set ("parser",
                            find_name (m.name () + L"_parser", set_));

          String const& base (m.context ().get<String> ("name"));
          m.context ().set ("member", find_name (base + L"_parser_", set_));

          if (polymorphic &&
              m.is_a<SemanticGraph::Element> () &&
              !m.type ().context ().count ("anonymous"))
          {
            m.context ().set (
              "member-cache", find_name (base + L"_parser_cache_", set_));

            m.context ().set (
              "member-map", find_name (base + L"_parser_map_", set_));

            m.context ().set (
              "member-map-impl",
              find_name (base + L"_parser_map_impl_", set_));
          }
        }

      private:
        NameSet& set_;
      };


      //
      //
      struct MemberInRestrictionBase: Traversal::Member
      {
      protected:
        MemberInRestrictionBase (NameSet& set, SemanticGraph::Complex& base)
            : set_ (set), base_ (base)
        {
        }

        struct NotFound {};

        Type&
        find_member (SemanticGraph::Complex& c, Type& m)
        {
          using SemanticGraph::Complex;

          Complex::NamesIteratorPair r (c.find (m.name ()));

          for (; r.first != r.second; ++r.first)
          {
            if (r.first->named ().is_a<Type> ())
            {
              Type& bm (dynamic_cast<Type&> (r.first->named ()));

              if (typeid (bm) != typeid (m))
                continue;

              if (m.qualified_p ())
              {
                if (bm.qualified_p () &&
                    m.name () == bm.name () &&
                    m.namespace_ ().name () == bm.namespace_ ().name ())
                  return bm;
              }
              else
              {
                if (!bm.qualified_p () && m.name () == bm.name ())
                  return bm;
              }
            }
          }

          // If we didn't find anything, try our base.
          //
          if (c.inherits_p ())
          {
            SemanticGraph::Type& base (c.inherits ().base ());

            if (base.is_a<Complex> ())
              return find_member (dynamic_cast<Complex&> (base), m);
          }

          //std::wcerr << "unable to find member " << m.name () << " in "
          //           << c.name () << std::endl;

          throw NotFound ();
        }

      protected:
        NameSet& set_;
        SemanticGraph::Complex& base_;
      };

      struct PrimaryMemberInRestriction: MemberInRestrictionBase, Context
      {
        PrimaryMemberInRestriction (Context& c,
                                    NameSet& set,
                                    SemanticGraph::Complex& base)
            : MemberInRestrictionBase (set, base), Context (c)
        {
        }

        virtual Void
        traverse (Type& m)
        {
          if (Parser::Context::skip (m))
            return;

          String name;

          try
          {
            // Try to find corresponding member in one of our bases.
            // This may fail if we use an element that substitutes
            // one in our base.
            //
            Type& bm (find_member (base_, m));
            name = bm.context ().get<String> ("name");
          }
          catch (NotFound const&)
          {
            // Fall back to the standard name assignment.
            //
            name = find_name (m.name (), set_);
          }

          m.context ().set ("name", name);
        }
      };

      struct DerivedMemberInRestriction: MemberInRestrictionBase, Context
      {
        DerivedMemberInRestriction (Context& c,
                                    NameSet& set,
                                    SemanticGraph::Complex& base)
            : MemberInRestrictionBase (set, base), Context (c)
        {
        }

        virtual Void
        traverse (Type& m)
        {
          if (Parser::Context::skip (m))
            return;

          Boolean poly (polymorphic &&
                        m.is_a<SemanticGraph::Element> () &&
                        !m.type ().context ().count ("anonymous"));

          String parser, member, member_cache, member_map, member_map_impl;

          try
          {
            // Try to find corresponding member in one of our bases.
            // This may fail if we use an element that substitutes
            // one in our base.
            //
            Type& bm (find_member (base_, m));
            parser = bm.context ().get<String> ("parser");
            member = bm.context ().get<String> ("member");

            if (poly)
            {
              member_cache = bm.context ().get<String> ("member-cache");
              member_map = bm.context ().get<String> ("member-map");
              member_map_impl = bm.context ().get<String> ("member-map-impl");
            }
          }
          catch (NotFound const&)
          {
            // Fall back to the standard name assignment.
            //
            String const& base (m.context ().get<String> ("name"));

            parser = find_name (m.name () + L"_parser", set_);
            member = find_name (base + L"_parser_", set_);

            if (poly)
            {
              member_cache = find_name (base + L"_parser_cache_", set_);
              member_map = find_name (base + L"_parser_map_", set_);
              member_map_impl = find_name (base + L"_parser_map_impl_", set_);
            }
          }

          m.context ().set ("parser", parser);
          m.context ().set ("member", member);

          if (poly)
          {
            m.context ().set ("member-cache", member_cache);
            m.context ().set ("member-map", member_map);
            m.context ().set ("member-map-impl", member_map_impl);
          }
        }
      };

      //
      //
      struct Complex: Traversal::Complex, Context
      {
        Complex (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& c)
        {
          SemanticGraph::Context& cc (c.context ());

          // Use processed name.
          //
          String const& name (cc.get<String> ("name"));

          // We leave this set around to allow other mappings to use
          // this information.
          //
          cc.set ("cxx-parser-name-processor-member-set", NameSet ());
          NameSet& member_set (
            cc.get<NameSet> ("cxx-parser-name-processor-member-set"));

          member_set.insert (name);

          // Add our base's members to the initial list unless we are
          // inheriting by restriction in which case we need to have
          // the same names as our base.
          //
          Boolean restriction (false);

          if (c.inherits_p ())
          {
            // @@ What if this types name is the same as one of base's
            //    members?
            //
            SemanticGraph::Type& base (c.inherits ().base ());

            if (base.is_a<SemanticGraph::Complex> ())
            {
              if (!base.context ().count (
                    "cxx-parser-name-processor-member-set"))
              {
                dispatch (base);
              }

              NameSet const& base_set (
                base.context ().get<NameSet> (
                  "cxx-parser-name-processor-member-set"));

              member_set.insert (base_set.begin (), base_set.end ());

              // Inheritance by restriction from anyType is a special case.
	      //
              restriction = c.inherits ().is_a<SemanticGraph::Restricts> () &&
	        !c.inherits ().base ().is_a<SemanticGraph::AnyType> ();
            }
          }

          if (restriction)
          {
            // First assign the "primary" names.
            //
            {
              PrimaryMemberInRestriction member (
                *this,
                member_set,
                dynamic_cast<SemanticGraph::Complex&> (
                  c.inherits ().base ()));

              Traversal::Names names (member);

              Complex::names (c, names);
            }

            // Assign "derived" names.
            //
            {
              DerivedMemberInRestriction member (
                *this,
                member_set,
                dynamic_cast<SemanticGraph::Complex&> (
                  c.inherits ().base ()));

              Traversal::Names names (member);

              Complex::names (c, names);
            }
          }
          else
          {
            // First assign the "primary" names.
            //
            {
              PrimaryMember member (*this, member_set);
              Traversal::Names names (member);

              Complex::names (c, names);
            }

            // Assign "derived" names.
            //
            {
              DerivedMember member (*this, member_set);
              Traversal::Names names (member);

              Complex::names (c, names);
            }
          }
        }
      };


      //
      //
      struct GlobalType: Traversal::Type, Context
      {
        GlobalType (Context& c, NameSet& set)
            : Context (c), set_ (set)
        {
        }

        virtual Void
        traverse (SemanticGraph::Type& t)
        {
          SemanticGraph::Context& c (t.context ());
          String const& n (t.name ());

          String name (find_name (n + skel_suffix, set_));
          c.set ("name", name);

          // Assign the post_* name.
          //
          c.set ("post", find_post_name (t));

          // Note that we do not add this name to the set so that it
          // does not influence other names.
          //
          if (impl)
            c.set ("impl", escape (n + impl_suffix));
        }

      private:
        String
        find_post_name (SemanticGraph::Type& t)
        {
          String const& n (t.name ());

          // It is possible that our base has the same type name (just
          // in a different namespaces). Avoid name clash in this case.
          //
          using SemanticGraph::Complex;

          Complex* c = dynamic_cast<Complex*> (&t);

          if (c == 0 || !c->inherits_p ())
          {
            return escape (L"post_" + n);
          }
          else
          {
            NameSet set;

            // Collect all base's post_*. In some mutual inclusion cases it
            // is possible that our base won't have the post name assigned
            // yet. In this situation will will have to figure it out
            // ourselves (we can do it since we use the "raw" type name).
            //
            SemanticGraph::Type* b (&c->inherits ().base ());

            while (true)
            {
              if (b->context ().count ("post"))
                set.insert (b->context ().get<String> ("post"));
              else
                set.insert (find_post_name (*b));

              Complex* cb (dynamic_cast<Complex*> (b));

              if (cb != 0 && cb->inherits_p ())
              {
                b = &cb->inherits ().base ();
                continue;
              }

              break;
            }

            String base_name (escape (L"post_" + n));
            String post (base_name);

            for (UnsignedLong i (1); set.find (post) != set.end (); ++i)
            {
              std::wostringstream os;
              os << i;
              post = base_name + os.str ();
            }

            return post;
          }
        }

      private:
        NameSet& set_;
      };


      struct Namespace: Traversal::Namespace, Context
      {
        Namespace (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& ns)
        {
          NameSet& type_set (global_type_names[ns.name ()]);

          GlobalType type (*this, type_set);
          Traversal::Names names (type);

          Traversal::Namespace::names (ns, names);
        }
      };


      struct FundType: Traversal::AnyType,
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
                       Traversal::Fundamental::Entities,

                       Context
      {
        FundType (Context& c)
            : Context (c)
        {
        }

        // anyType & anySimpleType.
        //
        virtual Void
        traverse (SemanticGraph::AnyType& t)
        {
          t.context ().set ("name", make_skel_name ("any_type"));
          t.context ().set ("impl", make_impl_name ("any_type"));
          t.context ().set ("post", String ("post_any_type"));
        }

        virtual Void
        traverse (SemanticGraph::AnySimpleType& t)
        {
          t.context ().set ("name", make_skel_name ("any_simple_type"));
          t.context ().set ("impl", make_impl_name ("any_simple_type"));
          t.context ().set ("post", String ("post_any_simple_type"));
        }

        // Boolean.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Boolean& t)
        {
          t.context ().set ("name", make_skel_name ("boolean"));
          t.context ().set ("impl", make_impl_name ("boolean"));
          t.context ().set ("post", String ("post_boolean"));
        }

        // Integral types.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Byte& t)
        {
          t.context ().set ("name", make_skel_name ("byte"));
          t.context ().set ("impl", make_impl_name ("byte"));
          t.context ().set ("post", String ("post_byte"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::UnsignedByte& t)
        {
          t.context ().set ("name", make_skel_name ("unsigned_byte"));
          t.context ().set ("impl", make_impl_name ("unsigned_byte"));
          t.context ().set ("post", String ("post_unsigned_byte"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Short& t)
        {
          t.context ().set ("name", make_skel_name ("short"));
          t.context ().set ("impl", make_impl_name ("short"));
          t.context ().set ("post", String ("post_short"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::UnsignedShort& t)
        {
          t.context ().set ("name", make_skel_name ("unsigned_short"));
          t.context ().set ("impl", make_impl_name ("unsigned_short"));
          t.context ().set ("post", String ("post_unsigned_short"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Int& t)
        {
          t.context ().set ("name", make_skel_name ("int"));
          t.context ().set ("impl", make_impl_name ("int"));
          t.context ().set ("post", String ("post_int"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::UnsignedInt& t)
        {
          t.context ().set ("name", make_skel_name ("unsigned_int"));
          t.context ().set ("impl", make_impl_name ("unsigned_int"));
          t.context ().set ("post", String ("post_unsigned_int"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Long& t)
        {
          t.context ().set ("name", make_skel_name ("long"));
          t.context ().set ("impl", make_impl_name ("long"));
          t.context ().set ("post", String ("post_long"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::UnsignedLong& t)
        {
          t.context ().set ("name", make_skel_name ("unsigned_long"));
          t.context ().set ("impl", make_impl_name ("unsigned_long"));
          t.context ().set ("post", String ("post_unsigned_long"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Integer& t)
        {
          t.context ().set ("name", make_skel_name ("integer"));
          t.context ().set ("impl", make_impl_name ("integer"));
          t.context ().set ("post", String ("post_integer"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NonPositiveInteger& t)
        {
          t.context ().set ("name", make_skel_name ("non_positive_integer"));
          t.context ().set ("impl", make_impl_name ("non_positive_integer"));
          t.context ().set ("post", String ("post_non_positive_integer"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NonNegativeInteger& t)
        {
          t.context ().set ("name", make_skel_name ("non_negative_integer"));
          t.context ().set ("impl", make_impl_name ("non_negative_integer"));
          t.context ().set ("post", String ("post_non_negative_integer"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::PositiveInteger& t)
        {
          t.context ().set ("name", make_skel_name ("positive_integer"));
          t.context ().set ("impl", make_impl_name ("positive_integer"));
          t.context ().set ("post", String ("post_positive_integer"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NegativeInteger& t)
        {
          t.context ().set ("name", make_skel_name ("negative_integer"));
          t.context ().set ("impl", make_impl_name ("negative_integer"));
          t.context ().set ("post", String ("post_negative_integer"));
        }

        // Floats.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Float& t)
        {
          t.context ().set ("name", make_skel_name ("float"));
          t.context ().set ("impl", make_impl_name ("float"));
          t.context ().set ("post", String ("post_float"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Double& t)
        {
          t.context ().set ("name", make_skel_name ("double"));
          t.context ().set ("impl", make_impl_name ("double"));
          t.context ().set ("post", String ("post_double"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Decimal& t)
        {
          t.context ().set ("name", make_skel_name ("decimal"));
          t.context ().set ("impl", make_impl_name ("decimal"));
          t.context ().set ("post", String ("post_decimal"));
        }

        // Strings.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::String& t)
        {
          t.context ().set ("name", make_skel_name ("string"));
          t.context ().set ("impl", make_impl_name ("string"));
          t.context ().set ("post", String ("post_string"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NormalizedString& t)
        {
          t.context ().set ("name", make_skel_name ("normalized_string"));
          t.context ().set ("impl", make_impl_name ("normalized_string"));
          t.context ().set ("post", String ("post_normalized_string"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Token& t)
        {
          t.context ().set ("name", make_skel_name ("token"));
          t.context ().set ("impl", make_impl_name ("token"));
          t.context ().set ("post", String ("post_token"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NameToken& t)
        {
          t.context ().set ("name", make_skel_name ("nmtoken"));
          t.context ().set ("impl", make_impl_name ("nmtoken"));
          t.context ().set ("post", String ("post_nmtoken"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NameTokens& t)
        {
          t.context ().set ("name", make_skel_name ("nmtokens"));
          t.context ().set ("impl", make_impl_name ("nmtokens"));
          t.context ().set ("post", String ("post_nmtokens"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Name& t)
        {
          t.context ().set ("name", make_skel_name ("name"));
          t.context ().set ("impl", make_impl_name ("name"));
          t.context ().set ("post", String ("post_name"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NCName& t)
        {
          t.context ().set ("name", make_skel_name ("ncname"));
          t.context ().set ("impl", make_impl_name ("ncname"));
          t.context ().set ("post", String ("post_ncname"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Language& t)
        {
          t.context ().set ("name", make_skel_name ("language"));
          t.context ().set ("impl", make_impl_name ("language"));
          t.context ().set ("post", String ("post_language"));
        }


        // Qualified name.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::QName& t)
        {
          t.context ().set ("name", make_skel_name ("qname"));
          t.context ().set ("impl", make_impl_name ("qname"));
          t.context ().set ("post", String ("post_qname"));
        }


        // ID/IDREF.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Id& t)
        {
          t.context ().set ("name", make_skel_name ("id"));
          t.context ().set ("impl", make_impl_name ("id"));
          t.context ().set ("post", String ("post_id"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::IdRef& t)
        {
          t.context ().set ("name", make_skel_name ("idref"));
          t.context ().set ("impl", make_impl_name ("idref"));
          t.context ().set ("post", String ("post_idref"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::IdRefs& t)
        {
          t.context ().set ("name", make_skel_name ("idrefs"));
          t.context ().set ("impl", make_impl_name ("idrefs"));
          t.context ().set ("post", String ("post_idrefs"));
        }

        // URI.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::AnyURI& t)
        {
          t.context ().set ("name", make_skel_name ("uri"));
          t.context ().set ("impl", make_impl_name ("uri"));
          t.context ().set ("post", String ("post_uri"));
        }

        // Binary.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Base64Binary& t)
        {
          t.context ().set ("name", make_skel_name ("base64_binary"));
          t.context ().set ("impl", make_impl_name ("base64_binary"));
          t.context ().set ("post", String ("post_base64_binary"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::HexBinary& t)
        {
          t.context ().set ("name", make_skel_name ("hex_binary"));
          t.context ().set ("impl", make_impl_name ("hex_binary"));
          t.context ().set ("post", String ("post_hex_binary"));
        }


        // Date/time.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Date& t)
        {
          t.context ().set ("name", make_skel_name ("date"));
          t.context ().set ("impl", make_impl_name ("date"));
          t.context ().set ("post", String ("post_date"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::DateTime& t)
        {
          t.context ().set ("name", make_skel_name ("date_time"));
          t.context ().set ("impl", make_impl_name ("date_time"));
          t.context ().set ("post", String ("post_date_time"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Duration& t)
        {
          t.context ().set ("name", make_skel_name ("duration"));
          t.context ().set ("impl", make_impl_name ("duration"));
          t.context ().set ("post", String ("post_duration"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Day& t)
        {
          t.context ().set ("name", make_skel_name ("gday"));
          t.context ().set ("impl", make_impl_name ("gday"));
          t.context ().set ("post", String ("post_gday"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Month& t)
        {
          t.context ().set ("name", make_skel_name ("gmonth"));
          t.context ().set ("impl", make_impl_name ("gmonth"));
          t.context ().set ("post", String ("post_gmonth"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::MonthDay& t)
        {
          t.context ().set ("name", make_skel_name ("gmonth_day"));
          t.context ().set ("impl", make_impl_name ("gmonth_day"));
          t.context ().set ("post", String ("post_gmonth_day"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Year& t)
        {
          t.context ().set ("name", make_skel_name ("gyear"));
          t.context ().set ("impl", make_impl_name ("gyear"));
          t.context ().set ("post", String ("post_gyear"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::YearMonth& t)
        {
          t.context ().set ("name", make_skel_name ("gyear_month"));
          t.context ().set ("impl", make_impl_name ("gyear_month"));
          t.context ().set ("post", String ("post_gyear_month"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Time& t)
        {
          t.context ().set ("name", make_skel_name ("time"));
          t.context ().set ("impl", make_impl_name ("time"));
          t.context ().set ("post", String ("post_time"));
        }

        // Entity.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Entity& t)
        {
          t.context ().set ("name", make_skel_name ("entity"));
          t.context ().set ("impl", make_impl_name ("entity"));
          t.context ().set ("post", String ("post_entity"));
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Entities& t)
        {
          t.context ().set ("name", make_skel_name ("entities"));
          t.context ().set ("impl", make_impl_name ("entities"));
          t.context ().set ("post", String ("post_entities"));
        }

      private:
        String
        make_skel_name (String const& base)
        {
          return escape (base + skel_suffix);
        }

        String
        make_impl_name (String const& base)
        {
          return escape (base + impl_suffix);
        }
      };

      // Go into sourced/included/imported schemas while making sure
      // we don't process the same stuff more than once.
      //
      struct Uses: Traversal::Sources,
                   Traversal::Includes,
                   Traversal::Imports
      {
        virtual Void
        traverse (SemanticGraph::Sources& sr)
        {
          SemanticGraph::Schema& s (sr.schema ());

          if (!s.context ().count ("cxx-parser-name-processor-seen"))
          {
            s.context ().set ("cxx-parser-name-processor-seen", true);
            Traversal::Sources::traverse (sr);
          }
        }

        virtual Void
        traverse (SemanticGraph::Includes& i)
        {
          SemanticGraph::Schema& s (i.schema ());

          if (!s.context ().count ("cxx-parser-name-processor-seen"))
          {
            s.context ().set ("cxx-parser-name-processor-seen", true);
            Traversal::Includes::traverse (i);
          }
        }

        virtual Void
        traverse (SemanticGraph::Imports& i)
        {
          SemanticGraph::Schema& s (i.schema ());

          if (!s.context ().count ("cxx-parser-name-processor-seen"))
          {
            s.context ().set ("cxx-parser-name-processor-seen", true);
            Traversal::Imports::traverse (i);
          }
        }
      };

      // Go into implied schemas while making sure we don't process
      // the same stuff more than once.
      //
      struct Implies: Traversal::Implies
      {
        virtual Void
        traverse (SemanticGraph::Implies& i)
        {
          SemanticGraph::Schema& s (i.schema ());

          if (!s.context ().count ("cxx-parser-name-processor-seen"))
          {
            s.context ().set ("cxx-parser-name-processor-seen", true);
            Traversal::Implies::traverse (i);
          }
        }
      };

      Void
      process_impl (CLI::Options const& ops,
                    SemanticGraph::Schema& tu,
                    SemanticGraph::Path const& file,
                    StringLiteralMap const& map)
      {
        Context ctx (ops, tu, file, &map);

        if (tu.names_begin ()->named ().name () ==
            L"http://www.w3.org/2001/XMLSchema")
        {
          // XML Schema namespace.
          //
          Traversal::Schema schema;

          Traversal::Names schema_names;
          Traversal::Namespace ns;
          Traversal::Names ns_names;
          FundType fund_type (ctx);

          schema >> schema_names >> ns >> ns_names >> fund_type;

          schema.dispatch (tu);
        }
        else
        {

          // Pass one - assign names to fundamental types.
          //
          {
            Traversal::Schema schema;
            Implies implies;
            Traversal::Schema xs_schema;

            schema >> implies >> xs_schema;

            Traversal::Names xs_schema_names;
            Traversal::Namespace ns;
            Traversal::Names ns_names;
            FundType fund_type (ctx);

            xs_schema >> xs_schema_names >> ns >> ns_names >> fund_type;

            schema.dispatch (tu);
          }

          // Pass two - assign names to global types. This pass cannot
          // be combined with pass three because of possible recursive
          // schema inclusions. Also note that we check first if this
          // schema has already been processed which may happen in the
          // file-per-type compilation mode.
          //
          if (!tu.context ().count ("cxx-parser-name-processor-seen"))
          {
            Traversal::Schema schema;
            Uses uses;

            schema >> uses >> schema;

            Traversal::Names schema_names;
            Namespace ns (ctx);

            schema >> schema_names >> ns;

            // Some twisted schemas do recusive self-inclusion.
            //
            tu.context ().set ("cxx-parser-name-processor-seen", true);

            schema.dispatch (tu);
          }

          // Pass three - assign names inside complex types. Here we don't
          // need to go into included/imported schemas.
          //
          {
            Traversal::Schema schema;
            Traversal::Sources sources;

            schema >> sources >> schema;

            Traversal::Names schema_names;
            Traversal::Namespace ns;
            Traversal::Names ns_names;

            schema >> schema_names >> ns >> ns_names;

            Complex complex (ctx);

            ns_names >> complex;

            schema.dispatch (tu);
          }
        }
      }
    }

    Void NameProcessor::
    process (CLI::Options const& ops,
             SemanticGraph::Schema& tu,
             SemanticGraph::Path const& file,
             StringLiteralMap const& map)
    {
      process_impl (ops, tu, file, map);
    }
  }
}
