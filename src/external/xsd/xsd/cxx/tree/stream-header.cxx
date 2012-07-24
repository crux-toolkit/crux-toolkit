// file      : xsd/cxx/tree/stream-header.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/tree/stream-header.hxx>

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

namespace CXX
{
  namespace Tree
  {
    namespace
    {
      struct List: Traversal::List, Context
      {
        List (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& l)
        {
          String name (ename (l));

          // If renamed name is empty then we do not need to generate
          // anything for this type.
          //
          if (renamed_type (l, name) && !name)
            return;

          os << inst_exp
             << std_ostream_type << "&" << endl
             << "operator<< (" << std_ostream_type << "&, const " <<
            name << "&);"
             << endl;
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
          String name (ename (u));

          // If renamed name is empty then we do not need to generate
          // anything for this type.
          //
          if (renamed_type (u, name) && !name)
            return;

          os << inst_exp
             << std_ostream_type << "&" << endl
             << "operator<< (" << std_ostream_type << "&, const " <<
            name << "&);"
             << endl;
        }
      };


      struct Enumeration: Traversal::Enumeration, Context
      {
        Enumeration (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& e)
        {
          String name (ename (e));

          // If renamed name is empty then we do not need to generate
          // anything for this type.
          //
          if (renamed_type (e, name) && !name)
            return;

          Boolean string_based (false);
          {
            IsStringBasedType t (string_based);
            t.dispatch (e);
          }

          Boolean enum_based (false);
          if (string_based)
          {
            SemanticGraph::Enumeration* be (0);
            IsEnumBasedType t (be);
            t.dispatch (e);

            enum_based = (be != 0);
          }

          // If we are based on an enum then the value type is just an
          // alias and we don't need to generate this operator again.
          //
          if (string_based && !enum_based)
          {
            os << inst_exp
               << std_ostream_type << "&" << endl
               << "operator<< (" << std_ostream_type << "&, " <<
              name << "::" << evalue (e) << ");"
               << endl;
          }

          os << inst_exp
             << std_ostream_type << "&" << endl
             << "operator<< (" << std_ostream_type << "&, const " <<
            name << "&);"
             << endl;
        }
      };

      struct Complex: Traversal::Complex, Context
      {
        Complex (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& c)
        {
          String name (ename (c));

          // If renamed name is empty then we do not need to generate
          // anything for this type.
          //
          if (renamed_type (c, name) && !name)
            return;

          os << inst_exp
             << std_ostream_type << "&" << endl
             << "operator<< (" << std_ostream_type << "&, const " <<
            name << "&);"
             << endl;
        }
      };
    }

    Void
    generate_stream_header (Context& ctx)
    {
      String c (ctx.char_type);

      ctx.os << "#include <iosfwd>" << endl
             << endl;

      Traversal::Schema schema;

      Traversal::Sources sources;
      Traversal::Names names_ns, names;

      Namespace ns (ctx);

      List list (ctx);
      Union union_ (ctx);
      Complex complex (ctx);
      Enumeration enumeration (ctx);

      schema >> sources >> schema;
      schema >> names_ns >> ns >> names;

      names >> list;
      names >> union_;
      names >> complex;
      names >> enumeration;

      schema.dispatch (ctx.schema_root);
    }
  }
}
