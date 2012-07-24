// file      : xsd/cxx/parser/impl-source.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/parser/impl-source.hxx>
#include <cxx/parser/print-impl-common.hxx>

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
          String const& base_ret (ret_type (base));

          os << "// " << name << endl
             << "//" << endl
             << endl;

          // pre
          //
          os << "void " << name << "::" << endl
             << "pre ()"
             << "{"
             << "}";

          // post
          //
          os << ret << " " << name << "::" << endl
             << post_name (e) << " ()"
             << "{";

          if (ret == base_ret)
          {
            os << (ret != L"void" ? "return " : "") <<
              post_name (base) << " ();";
          }
          else if (ret == L"void")
          {
            os << arg_type (base) << " v (" << post_name (base) << " ());"
               << endl;

            if (options.value<CLI::generate_print_impl> ())
            {
              PrintCall t (*this, e.name (), "v");
              t.dispatch (base);
            }
            else
              os << "// TODO" << endl
                 << "//" << endl;
          }
          else
          {
            if (base_ret == L"void")
              os << post_name (base) << " ();";
            else
              os << arg_type (base) << " v (" << post_name (base) << " ());"
                 << endl
                 << "// TODO" << endl
                 << "//" << endl
                 << "// return ... ;" << endl;
          }

          os << "}";
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
          SemanticGraph::Type& type (l.argumented ().type ());

          String item (unclash (ename (l), "item"));

          os << "// " << name << endl
             << "//" << endl
             << endl;

          // pre
          //
          os << "void " << name << "::" << endl
             << "pre ()"
             << "{"
             << "}";

          // item
          //
          String const& arg (arg_type (type));

          os << "void " << name << "::" << endl
             << item;

          if (arg == L"void")
            os << " ()";
          else
            os << " (" << arg << " " << item << ")";

          os << "{";

          if (arg != L"void")
          {
            if (options.value<CLI::generate_print_impl> ())
            {
              PrintCall t (*this, type.name (), item);
              t.dispatch (type);
            }
            else
              os << "// TODO" << endl
                 << "//" << endl;
          }

          os << "}";

          // post
          //
          String const& ret (ret_type (l));

          os << ret << " " << name << "::" << endl
             << post_name (l) << " ()"
             << "{";

          if (ret != L"void")
            os << "// TODO" << endl
               << "//" << endl
               << "// return ... ;" << endl;

          os << "}";
        }
      };

      //
      //
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

          os << "// " << name << endl
             << "//" << endl
             << endl;

          // pre
          //
          os << "void " << name << "::" << endl
             << "pre ()"
             << "{"
             << "}";

          // _characters
          //
          os << "void " << name << "::" << endl
             << "_characters (const " << string_type << "& s)"
             << "{";

          if (options.value<CLI::generate_print_impl> ())
            os << cout_inst << " << " << strlit (u.name () + L": ") <<
              " << s << std::endl;";
          else
            os << "// TODO" << endl
               << "//" << endl;

          os << "}";

          // post
          //
          String const& ret (ret_type (u));

          os << ret << " " << name << "::" << endl
             << post_name (u) << " ()"
             << "{";

          if (ret != L"void")
            os << "// TODO" << endl
               << "//" << endl
               << "// return ... ;" << endl;

          os << "}";
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

          String const& name (ename (m));
          String const& arg (arg_type (m.type ()));

          os << "void " <<
            eimpl (dynamic_cast<SemanticGraph::Complex&> (m.scope ())) <<
            "::" << endl
             << name;

          if (arg == L"void")
            os << " ()";
          else
            os << " (" << arg << " " << name << ")";

          os << "{";

          if (arg != L"void")
          {
            if (options.value<CLI::generate_print_impl> ())
            {
              PrintCall t (*this, m.name (), name);
              t.dispatch (m.type ());
            }
            else
              os << "// TODO" << endl
                 << "//" << endl;
          }

          os << "}";
        }
      };

      //
      //
      struct Complex: Traversal::Complex, Context
      {
        Complex (Context& c)
            : Context (c), parser_callback_ (c)
        {
          names_parser_callback_ >> parser_callback_;
        }

        virtual Void
        traverse (Type& c)
        {
          String const& name (eimpl (c));

          Boolean restriction (restriction_p (c));

          os << "// " << name << endl
             << "//" << endl
             << endl;

          // pre
          //
          os << "void " << name << "::" << endl
             << "pre ()"
             << "{"
             << "}";

          // Parser callbacks.
          //
          if (!restriction)
            names (c, names_parser_callback_);

          // post
          //
          String const& ret (ret_type (c));

          os << ret << " " << name << "::" << endl
             << post_name (c) << " ()"
             << "{";

          if (c.inherits_p ())
          {
            SemanticGraph::Type& base (c.inherits ().base ());
            String const& base_ret (ret_type (base));

            if (ret == base_ret)
            {
              os << (ret != L"void" ? "return " : "") <<
                post_name (base) << " ();";
            }
            else if (ret == L"void")
            {
              os << arg_type (base) << " v (" << post_name (base) << " ());"
                 << endl;

              if (options.value<CLI::generate_print_impl> ())
              {
                PrintCall t (*this, c.name (), "v");
                t.dispatch (base);
              }
              else
                os << "// TODO" << endl
                   << "//" << endl;
            }
            else
            {
              if (base_ret == L"void")
                os << post_name (base) << " ();";
              else
                os << arg_type (base) << " v (" << post_name (base) << " ());"
                   << endl
                   << "// TODO" << endl
                   << "//" << endl
                   << "// return ... ;" << endl;
            }
          }
          else
          {
            if (ret != L"void")
              os << "// TODO" << endl
                 << "//" << endl
                 << "// return ... ;" << endl;
          }

          os << "}";
        }

      private:
        //
        //
        ParserCallback parser_callback_;
        Traversal::Names names_parser_callback_;
      };
    }

    Void
    generate_impl_source (Context& ctx)
    {
      if (ctx.options.value<CLI::generate_print_impl> ())
        ctx.os << "#include <iostream>" << endl
               << endl;

      Traversal::Schema schema;
      Traversal::Sources sources;
      Traversal::Names schema_names;
      Namespace ns (ctx);
      Traversal::Names names;

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
