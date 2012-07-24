// file      : xsd/cxx/parser/impl-header.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/parser/impl-header.hxx>

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
          String const& name (eimpl (e));
          String const& ret (ret_type (e));
          SemanticGraph::Type& base (e.inherits ().base ());

          os << "class " << type_exp << name << ": " <<
            "public virtual " << ename (e) << "," << endl
             << "  public " << fq_name (base, "impl")
             << "{"
             << "public:" << endl
             << "virtual void" << endl
             << "pre ();"
             << endl
             << "virtual " << ret << endl
             << post_name (e) << " ();"
             << "};";
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
          String const& name (eimpl (l));
          SemanticGraph::Type& t (l.argumented ().type ());

          String item (unclash (ename (l), "item"));

          os << "class " << type_exp << name << ": public virtual " <<
            ename (l)
             << "{"
             << "public:" << endl
             << "virtual void" << endl
             << "pre ();"
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
             << post_name (l) << " ();"
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
          String const& name (eimpl (u));
          String const& ret (ret_type (u));

          os << "class " << type_exp << name << ": public virtual " <<
            ename (u)
             << "{"
             << "public:" << endl
             << "virtual void" << endl
             << "pre ();"
             << endl
             << "virtual void" << endl
             << "_characters (const " << string_type << "&);"
             << endl
             << "virtual " << ret << endl
             << post_name (u) << " ();"
             << "};";
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
	  if (skip (m))
            return;

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
      struct Complex: Traversal::Complex, Context
      {
        Complex (Context& c)
            : Context (c),
              parser_callback_ (c)
        {
          names_parser_callback_ >> parser_callback_;
        }

        virtual Void
        traverse (Type& c)
        {
          String const& name (eimpl (c));
          String const& ret (ret_type (c));

          os << "class " << type_exp << name << ": public virtual " <<
            ename (c);

          if (c.inherits_p ())
            os << "," << endl
               << "  public " << fq_name (c.inherits ().base (), "impl");

          os << "{"
             << "public:" << endl
             << "virtual void" << endl
             << "pre ();"
             << endl;

          // In case of an inheritance-by-restriction, we don't need to
          // generate parser callbacks, etc. since they are the same as in
          // the base.
          //
          if (!restriction_p (c))
          {
            names (c, names_parser_callback_);
          }

          os << "virtual " << ret << endl
             << post_name (c) << " ();"
             << "};";
        }

      private:
        //
        //
        ParserCallback parser_callback_;
        Traversal::Names names_parser_callback_;
      };
    }

    Void
    generate_impl_header (Context& ctx)
    {
      Traversal::Schema schema;

      Traversal::Sources sources;
      Includes includes (ctx, Includes::impl_header);
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
