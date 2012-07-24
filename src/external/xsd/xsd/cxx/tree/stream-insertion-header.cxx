// file      : xsd/cxx/tree/stream-insertion-header.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/tree/stream-insertion-header.hxx>

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

namespace CXX
{
  namespace Tree
  {
    namespace
    {
      typedef Containers::Vector<NarrowString> Streams;

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

          Streams const& st (options.value<CLI::generate_insertion> ());
          for (Streams::ConstIterator i (st.begin ()); i != st.end (); ++i)
          {
            String stream_type (ostream_type + L"< " + String (*i) + L" >");

            os << inst_exp
               << stream_type << "&" << endl
               << "operator<< (" << stream_type << "&," << endl
               << "const " << name << "&);"
               << endl;
          }
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

          Streams const& st (options.value<CLI::generate_insertion> ());
          for (Streams::ConstIterator i (st.begin ()); i != st.end (); ++i)
          {
            String stream_type (ostream_type + L"< " + String (*i) + L" >");

            os << inst_exp
               << stream_type << "&" << endl
               << "operator<< (" << stream_type << "&," << endl
               << "const " << name << "&);"
               << endl;
          }
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

          Streams const& st (options.value<CLI::generate_insertion> ());
          for (Streams::ConstIterator i (st.begin ()); i != st.end (); ++i)
          {
            String stream_type (ostream_type + L"< " + String (*i) + L" >");

            os << inst_exp
               << stream_type << "&" << endl
               << "operator<< (" << stream_type << "&," << endl
               << "const " << name << "&);"
               << endl;
          }
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

          Streams const& st (options.value<CLI::generate_insertion> ());
          for (Streams::ConstIterator i (st.begin ()); i != st.end (); ++i)
          {
            String stream_type (ostream_type + L"< " + String (*i) + L" >");

            os << inst_exp
               << stream_type << "&" << endl
               << "operator<< (" << stream_type << "&," << endl
               << "const " << name << "&);"
               << endl;
          }
        }
      };
    }

    Void
    generate_stream_insertion_header (Context& ctx)
    {
      String c (ctx.char_type);

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
