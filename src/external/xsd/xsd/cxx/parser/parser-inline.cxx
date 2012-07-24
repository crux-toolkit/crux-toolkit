// file      : xsd/cxx/parser/parser-inline.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/parser/parser-inline.hxx>

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

namespace CXX
{
  namespace Parser
  {
    namespace
    {
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

          os << "// " << name << endl
             << "//" << endl
             << endl;

          // item_parser
          //
          os << inl
             << "void " << name << "::" << endl
             << unclash (name, "item_parser") << " (" <<
            fq_name (t) << "& " << item << ")"
             << "{"
             << "this->_xsd_" << item << "_ = &" << item << ";"
             << "}";

          // parsers
          //
          os << inl
             << "void " << name << "::" << endl
             << "parsers (" << fq_name (t) << "& " << item << ")"
             << "{"
             << "this->_xsd_" << item << "_ = &" << item << ";"
             << "}";

          // c-tor
          //
          os << inl
             << name << "::" << endl
             << name << " ()" << endl
             << ": _xsd_" << item << "_ (0)"
             << "{"
             << "}";
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

          String const& scope (ename (m.scope ()));
          String const& parser (eparser (m));

          Boolean poly (polymorphic &&
                        m.is_a<SemanticGraph::Element> () &&
                        !anonymous (m.type ()));

          os << inl
             << "void " << scope << "::" << endl
             << parser << " (" << fq_name (m.type ()) << "& p)"
             << "{"
             << "this->" << emember (m) << " = &p;"
             << "}";

          if (poly)
          {
            os << inl
               << "void " << scope << "::" << endl
               << parser << " (const " << parser_map << "& m)"
               << "{"
               << "this->" << emember_map (m) << " = &m;"
               << "}";
          }
        }
      };


      //
      //
      struct ParserMemberSet: Traversal::Member, Context
      {
        ParserMemberSet (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& m)
        {
	  if (skip (m)) return;

          String const& name (ename (m));

          os << "this->" << emember (m) << " = &" << name << ";";
        }
      };


      //
      //
      struct ParserMemberInit: Traversal::Member, Context
      {
        ParserMemberInit (Context& c)
            : Context (c), first_ (true)
        {
        }

        virtual Void
        traverse (Type& m)
        {
	  if (skip (m)) return;

          if (first_)
            first_ = false;
          else
            os << "," << endl << "  ";

          os << emember (m) << " (0)";

          if (polymorphic &&
              m.is_a<SemanticGraph::Element> () &&
              !anonymous (m.type ()))
          {
            os << "," << endl
               << "  " << emember_map (m) << " (0)";
          }
        }

        Boolean
        comma () const
        {
          return !first_;
        }

      private:
        Boolean first_;
      };

      struct ParserBaseSet: Traversal::Complex,
                            Traversal::List,
                            Context
      {
        ParserBaseSet (Context& c)
            : Context (c), member_ (c)
        {
          inherits_ >> *this;
          names_ >> member_;
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
          String const& name (ename (l));
          String item (unclash (name, "item"));

          os << "this->_xsd_" << item << "_ = &" << name << "_item;";
        }

      private:
        Traversal::Inherits inherits_;

        ParserMemberSet member_;
        Traversal::Names names_;
      };

      struct Particle: Traversal::All, Context
      {
        Particle (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (SemanticGraph::All& a)
        {
          if (!a.context().count ("comp-number"))
            return;

          UnsignedLong state_count (
            a.context().get<UnsignedLong> ("state-count"));

          os << "," << endl
             << "  v_all_count_ (" << state_count << "UL, v_all_first_)";
        }
      };

      //
      //
      struct Complex: Traversal::Complex, Context
      {
        Complex (Context& c)
            : Context (c),
              parser_modifier_ (c),
              parser_base_set_ (c),
              parser_member_set_ (c),
              particle_ (c)
        {
          names_parser_modifier_ >> parser_modifier_;
          inherits_parser_base_set_ >> parser_base_set_;
          names_parser_member_set_ >> parser_member_set_;
        }

        virtual Void
        traverse (Type& c)
        {
          Boolean he (has<Traversal::Element> (c));
          Boolean ha (has<Traversal::Attribute> (c));

          Boolean hae (has_particle<Traversal::Any> (c));

          Boolean hra (false); // Has required attribute.
          if (ha)
          {
            RequiredAttributeTest test (hra);
            Traversal::Names names_test (test);
            names (c, names_test);
          }

          Boolean restriction (restriction_p (c));

          if (!((!restriction && (he || ha)) ||
                (validation && (he || hae || hra))))
            return;

          String const& name (ename (c));

          os << "// " << name << endl
             << "//" << endl
             << endl;

          if (!restriction && (he || ha))
          {
            // <name>_parser ()
            //
            names (c, names_parser_modifier_);


            // parsers ()
            //

            os << inl
               << "void " << name << "::" << endl
               << "parsers (";

            {
              ParserParamDecl decl (*this, true);
              decl.traverse (c);
            }

            os << ")"
               << "{";

            inherits (c, inherits_parser_base_set_);
            names (c, names_parser_member_set_);

            os << "}";
          }

          // Default c-tor.
          //
          os << inl
             << name << "::" << endl
             << name << " ()" << endl
             << ": ";

          Boolean comma (false);

          if (!restriction && (he || ha))
          {
            ParserMemberInit member_init (*this);
            Traversal::Names names_member_init (member_init);

            names (c, names_member_init);

            comma = member_init.comma ();
          }

          if (validation && (he || hae))
          {
            if (comma)
              os << "," << endl << "  ";

            os << "v_state_stack_ (sizeof (v_state_), &v_state_first_)";

            particle_.dispatch (c.contains_compositor ().compositor ());

            comma = true;
          }

          if (validation && (hra))
          {
            if (comma)
              os << "," << endl << "  ";

            os << "v_state_attr_stack_ (sizeof (v_state_attr_), " <<
              "&v_state_attr_first_)";
          }

          os << "{"
             << "}";
        }

      private:
        //
        //
        ParserModifier parser_modifier_;
        Traversal::Names names_parser_modifier_;

        //
        //
        ParserBaseSet parser_base_set_;
        Traversal::Inherits inherits_parser_base_set_;

        //
        //
        ParserMemberSet parser_member_set_;
        Traversal::Names names_parser_member_set_;

        //
        //
        Particle particle_;
      };
    }

    Void
    generate_parser_inline (Context& ctx)
    {
      // Emit "weak" header includes that are used in the file-per-type
      // compilation model.
      //
      if (!ctx.options.value<CLI::generate_inline> ())
      {
        Traversal::Schema schema;
        Includes includes (ctx, Includes::source);

        schema >> includes;
        schema.dispatch (ctx.schema_root);
      }

      Traversal::Schema schema;

      Traversal::Sources sources;
      Traversal::Names schema_names;

      Namespace ns (ctx);
      Traversal::Names names;

      schema >> sources >> schema;
      schema >> schema_names >> ns >> names;

      List list (ctx);
      Complex complex (ctx);

      names >> list;
      names >> complex;

      schema.dispatch (ctx.schema_root);
    }
  }
}
