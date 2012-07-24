// file      : xsd/cxx/parser/driver-source.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/parser/driver-source.hxx>
#include <cxx/parser/print-impl-common.hxx>

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

#include <cult/containers/map.hxx>
#include <cult/containers/set.hxx>

#include <sstream>

namespace CXX
{
  namespace Parser
  {
    namespace
    {
      typedef
      Cult::Containers::Map<SemanticGraph::Type*, String>
      TypeInstanceMap;

      typedef Cult::Containers::Set<String> InstanceSet;

      // For base types we only want member's types, but not the
      // base itself.
      //
      struct BaseType: Traversal::Complex, Context
      {
        BaseType (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (SemanticGraph::Complex& c)
        {
          inherits (c);

          if (!restriction_p (c))
            names (c);
        }
      };

      struct ParserDef: Traversal::Type,
                        Traversal::List,
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
                        Traversal::Fundamental::Entities,

                        Context
      {
        ParserDef (Context& c, TypeInstanceMap& map, InstanceSet& set)
            : Context (c), map_ (map), set_ (set), base_ (c)
        {
          *this >> inherits_ >> base_ >> inherits_;

          *this >> names_;
          base_ >> names_;

          names_ >> member_ >> belongs_ >> *this;
        }

        virtual Void
        traverse (SemanticGraph::Type& t)
        {
          if (map_.find (&t) == map_.end ())
          {
            String inst (find_instance_name (t));
            map_[&t] = inst;

            os << fq_name (t, "impl") << " " << inst << ";";
          }
        }

        virtual Void
        traverse (SemanticGraph::List& l)
        {
          if (map_.find (&l) == map_.end ())
          {
            String inst (find_instance_name (l));
            map_[&l] = inst;

            os << fq_name (l, "impl") << " " << inst << ";";

            dispatch (l.argumented ().type ());
          }
        }

        virtual Void
        traverse (SemanticGraph::Complex& c)
        {
          if (map_.find (&c) == map_.end ())
          {
            String inst (find_instance_name (c));
            map_[&c] = inst;

            os << fq_name (c, "impl") << " " << inst << ";";

            inherits (c);

            if (!restriction_p (c))
              names (c);
          }
        }

        // anyType & anySimpleType.
        //
        virtual Void
        traverse (SemanticGraph::AnyType& t)
        {
          fund_type (t, "any_type");
        }

        virtual Void
        traverse (SemanticGraph::AnySimpleType& t)
        {
          fund_type (t, "any_simple_type");
        }

        // Boolean.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Boolean& t)
        {
          fund_type (t, "boolean");
        }

        // Integral types.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Byte& t)
        {
          fund_type (t, "byte");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::UnsignedByte& t)
        {
          fund_type (t, "unsigned_byte");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Short& t)
        {
          fund_type (t, "short");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::UnsignedShort& t)
        {
          fund_type (t, "unsigned_short");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Int& t)
        {
          fund_type (t, "int");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::UnsignedInt& t)
        {
          fund_type (t, "unsigned_int");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Long& t)
        {
          fund_type (t, "long");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::UnsignedLong& t)
        {
          fund_type (t, "unsigned_long");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Integer& t)
        {
          fund_type (t, "integer");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NonPositiveInteger& t)
        {
          fund_type (t, "non_positive_integer");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NonNegativeInteger& t)
        {
          fund_type (t, "non_negative_integer");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::PositiveInteger& t)
        {
          fund_type (t, "positive_integer");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NegativeInteger& t)
        {
          fund_type (t, "negative_integer");
        }

        // Floats.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Float& t)
        {
          fund_type (t, "float");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Double& t)
        {
          fund_type (t, "double");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Decimal& t)
        {
          fund_type (t, "decimal");
        }

        // Strings.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::String& t)
        {
          fund_type (t, "string");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NormalizedString& t)
        {
          fund_type (t, "normalized_string");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Token& t)
        {
          fund_type (t, "token");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NameToken& t)
        {
          fund_type (t, "nmtoken");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NameTokens& t)
        {
          fund_type (t, "nmtokens");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Name& t)
        {
          fund_type (t, "name");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NCName& t)
        {
          fund_type (t, "ncname");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Language& t)
        {
          fund_type (t, "language");
        }


        // Qualified name.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::QName& t)
        {
          fund_type (t, "qname");
        }


        // ID/IDREF.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Id& t)
        {
          fund_type (t, "id");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::IdRef& t)
        {
          fund_type (t, "idref");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::IdRefs& t)
        {
          fund_type (t, "idrefs");
        }

        // URI.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::AnyURI& t)
        {
          fund_type (t, "uri");
        }

        // Binary.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Base64Binary& t)
        {
          fund_type (t, "base64_binary");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::HexBinary& t)
        {
          fund_type (t, "hex_binary");
        }


        // Date/time.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Date& t)
        {
          fund_type (t, "date");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::DateTime& t)
        {
          fund_type (t, "date_time");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Duration& t)
        {
          fund_type (t, "duration");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Day& t)
        {
          fund_type (t, "day");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Month& t)
        {
          fund_type (t, "month");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::MonthDay& t)
        {
          fund_type (t, "month_day");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Year& t)
        {
          fund_type (t, "year");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::YearMonth& t)
        {
          fund_type (t, "year_month");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Time& t)
        {
          fund_type (t, "time");
        }

        // Entity.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Entity& t)
        {
          fund_type (t, "entity");
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Entities& t)
        {
          fund_type (t, "entities");
        }

      private:
        virtual Void
        fund_type (SemanticGraph::Type& t, String const& name)
        {
          if (map_.find (&t) == map_.end ())
          {
            String inst (find_instance_name (name));
            map_[&t] = inst;

            os << fq_name (t, "impl") << " " << inst << ";";
          }
        }

        String
        find_instance_name (String const& raw_name)
        {
          String base_name (escape (raw_name + L"_p"));
          String name (base_name);

          for (UnsignedLong i (1); set_.find (name) != set_.end (); ++i)
          {
            std::wostringstream os;
            os << i;
            name = base_name + os.str ();
          }

          set_.insert (name);
          return name;
        }

        String
        find_instance_name (SemanticGraph::Type& t)
        {
          return find_instance_name (t.name ());
        }

        TypeInstanceMap& map_;
        InstanceSet& set_;

        BaseType base_;
        Traversal::Inherits inherits_;

        Traversal::Names names_;
        Traversal::Member member_;
        Traversal::Belongs belongs_;
      };

      struct ArgList: Traversal::Complex,
                      Traversal::List,
                      Traversal::Member,
                      Context
      {
        ArgList (Context& c, TypeInstanceMap& map)
            : Context (c), map_ (map), first_ (true)
        {
          inherits_ >> *this;
          names_ >> *this;
        }

        virtual Void
        traverse (SemanticGraph::Complex& c)
        {
          inherits (c, inherits_);

          if (!restriction_p (c))
            names (c, names_);
        }

        virtual Void
        traverse (SemanticGraph::List& l)
        {
          if (!first_)
            os << "," << endl;
          else
            first_ = false;

          os << map_[&l.argumented ().type ()];
        }

        virtual Void
        traverse (SemanticGraph::Member& m)
        {
          if (skip (m))
            return;

          if (!first_)
            os << "," << endl;
          else
            first_ = false;

          os << map_[&m.type ()];
        }

      private:
        TypeInstanceMap& map_;

        Traversal::Inherits inherits_;
        Traversal::Names names_;

        Boolean first_;
      };

      struct ParserConnect: Traversal::List,
                            Traversal::Complex,
                            Context
      {
        ParserConnect (Context& c, TypeInstanceMap& map)
            : Context (c), map_ (map), base_ (c)
        {
          *this >> inherits_ >> base_ >> inherits_;

          *this >> names_;
          base_ >> names_;

          names_ >> member_ >> belongs_ >> *this;
        }

        virtual Void
        traverse (SemanticGraph::List& l)
        {
          if (type_set_.find (&l) == type_set_.end ())
          {
            os << map_[&l] << ".parsers (" <<
              map_[&l.argumented ().type ()] << ");"
               << endl;

            type_set_.insert (&l);
          }
        }

        virtual Void
        traverse (SemanticGraph::Complex& c)
        {
          if (type_set_.find (&c) == type_set_.end ())
          {
            if (has_members (c))
            {
              os << map_[&c] << ".parsers (";

              ArgList args (*this, map_);
              args.dispatch (c);

              os << ");"
                 << endl;
            }

            type_set_.insert (&c);

            inherits (c);

            if (!restriction_p (c))
              names (c);
          }
        }

      private:
        Boolean
        has_members (SemanticGraph::Complex& c)
        {
          using SemanticGraph::Complex;

          if (has<Traversal::Member> (c))
            return true;

          if (c.inherits_p ())
          {
            SemanticGraph::Type& b (c.inherits ().base ());

            if (Complex* cb = dynamic_cast<Complex*> (&b))
              return has_members (*cb);

            return b.is_a<SemanticGraph::List> ();
          }

          return false;
        }

      private:
        TypeInstanceMap& map_;
        Cult::Containers::Set<SemanticGraph::Type*> type_set_;

        BaseType base_;
        Traversal::Inherits inherits_;

        Traversal::Names names_;
        Traversal::Member member_;
        Traversal::Belongs belongs_;
      };
    }

    Void
    generate_driver_source (Context& ctx)
    {
      // Figure out the root element. Validator should have made sure
      // it is unique.
      //
      SemanticGraph::Element* root (0);
      {
        Traversal::Schema schema;
        Traversal::Sources sources;

        schema >> sources >> schema;

        Traversal::Names schema_names;
        Traversal::Namespace ns;
        Traversal::Names ns_names;
        RootElement root_element (ctx.options, root);

        schema >> schema_names >> ns >> ns_names >> root_element;

        schema.dispatch (ctx.schema_root);
      }

      std::wostream& os (ctx.os);
      String const& L (ctx.L);
      String const& cerr (ctx.cerr_inst);

      InstanceSet set;
      TypeInstanceMap map;
      SemanticGraph::Type& root_type (root->type ());

      set.insert ("doc_p");

      os << "#include <iostream>" << endl
         << endl
         << "int" << endl
         << "main (int argc, char* argv[])"
         << "{"
         << "if (argc != 2)"
         << "{"
         << cerr << " << " << L << "\"usage: \" << argv[0] << " <<
        L << "\" file.xml\" << std::endl;"
         << "return 1;"
         << "}"
         << "try"
         << "{"
         << "// Instantiate individual parsers." << endl
         << "//" << endl;

      {
        ParserDef def (ctx, map, set);
        def.dispatch (root_type);
      }

      os << endl
         << "// Connect the parsers together." << endl
         << "//" << endl;

      {
        ParserConnect connect (ctx, map);
        connect.dispatch (root_type);
      }

      String const& root_p (map[&root_type]);

      os << "// Parse the XML document." << endl
         << "//" << endl;

      if (root->namespace_().name ())
        os << ctx.xs_ns_name () << "::document doc_p (" << endl
           << root_p << "," << endl
           << ctx.strlit (root->namespace_().name ()) << "," << endl
           << ctx.strlit (root->name ()) << ");"
           << endl;
      else
        os << ctx.xs_ns_name () << "::document doc_p (" << root_p << ", " <<
          ctx.strlit (root->name ()) << ");"
           << endl;

      os << root_p << ".pre ();"
         << "doc_p.parse (argv[1]);";

      String const& ret (Context::ret_type (root_type));
      String const& post (Context::post_name (root_type));

      if (ret == L"void")
        os << root_p << "." << post << " ();";
      else
      {
        os << Context::arg_type (root_type) << " v (" <<
          root_p << "." << post << " ());"
           << endl;

        if (ctx.options.value<CLI::generate_print_impl> ())
        {
          PrintCall t (ctx, root->name (), "v");
          t.dispatch (root_type);
        }
        else
          os << "// TODO" << endl
             << "//" << endl;
      }

      os << "}" // try
         << "catch (const " << ctx.xs_ns_name () << "::exception& e)"
         << "{"
         << cerr << " << e << std::endl;"
         << "return 1;"
         << "}"
         << "catch (const std::ios_base::failure&)"
         << "{"
         << cerr << " << argv[1] << " <<
        L << "\": error: io failure\" << std::endl;"
         << "return 1;"
         << "}"
         << "}";
    }
  }
}
