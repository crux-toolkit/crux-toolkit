// file      : xsd/cxx/parser/characters-validation-source.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/parser/characters-validation-source.hxx>

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
      struct Complex: Traversal::Complex, Context
      {
        Complex (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& c)
        {
          if (!c.mixed_p ())
            return;

          String const& name (ename (c));

          os <<"// Character validation functions for " << name << "." << endl
             <<"//" << endl;

            // _characters_impl
            //
          os << "bool " << name << "::" << endl
             << "_characters_impl (const " << string_type << "& s)"
             << "{"
             << "this->_any_characters (s);"
             << "return true;"
             << "}";
        }
      };
    }

    Void
    generate_characters_validation_source (Context& ctx)
    {
      //@@ Most of the time there is no mixed content type so
      //   we generate an empty namespace which looks ugly. Will
      //   need to implement smart namespace to handle this at
      //   some point.
      //
      Traversal::Schema schema;

      Traversal::Sources sources;
      Traversal::Names schema_names;

      Namespace ns (ctx);
      Traversal::Names names;

      schema >> sources >> schema;
      schema >> schema_names >> ns >> names;

      Complex complex (ctx);

      names >> complex;

      schema.dispatch (ctx.schema_root);
    }
  }
}
