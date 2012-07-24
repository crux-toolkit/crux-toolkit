// file      : xsd/cxx/parser/parser-source.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/parser/parser-source.hxx>

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
          String const& ret (ret_type (e));

          SemanticGraph::Type& base (e.inherits ().base ());

          Boolean same (ret == ret_type (base));

          if (same || ret == L"void" || polymorphic)
          {
            os << "// " << name << endl
               << "//" << endl
               << endl;
          }

          if (same || ret == L"void")
          {
            os << ret << " " << name << "::" << endl
               << post_name (e) << " ()"
               << "{";

            if (same)
            {
              if (ret == L"void")
                os << post_name (base) << " ();";
              else
                os << "return " << post_name (base) << " ();";
            }

            os << "}";
          }

          if (polymorphic)
          {
            String id (e.name ());

            if (String ns = xml_ns_name (e))
            {
              id += L' ';
              id += ns;
            }

            os << "const " << char_type << "* " << name << "::" << endl
               << "_static_type ()"
               << "{"
               << "return " << strlit (id) << ";"
               << "}";

            os << "const " << char_type << "* " << name << "::" << endl
               << "_dynamic_type () const"
               << "{"
               << "return _static_type ();"
               << "}";

            if (validation)
            {
              Boolean gen (!anonymous (e));

              // We normally don't need to enter anonymous types into
              // the inheritance map. The only exception is when an
              // anonymous types is defined inside an element that
              // is a member of a substitution group.
              //
              if (!gen)
              {
                // The first instance that this anonymous type classifies
                // is the prototype for others if any.
                //
                SemanticGraph::Instance& i (
                  e.classifies_begin ()->instance ());

                if (SemanticGraph::Element* e =
                    dynamic_cast<SemanticGraph::Element*> (&i))
                {
                  if (e->substitutes_p ())
                    gen = true;
                }
              }

              if (gen)
              {
                os << "static" << endl
                   << "const ::xsd::cxx::parser::validating::inheritance_map_entry< " <<
                  char_type << " >" << endl
                   << "_xsd_" << name << "_inheritance_map_entry_ (" << endl
                   << name << "::_static_type ()," << endl
                   << fq_name (base) << "::_static_type ());"
                   << endl;
              }
            }
          }
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

          os << "// " << name << endl
             << "//" << endl
             << endl;

          // item
          //
          String const& arg (arg_type (t));

          os << "void " << name << "::" << endl
             << item;

          if (arg == L"void")
            os << " ()";
          else
            os << " (" << arg << ")";

          os << "{"
             << "}";

          // post
          //
          if (ret_type (l) == L"void")
            os << "void " << name << "::" << endl
               << post_name (l) << " ()"
               << "{"
               << "}";

          // parse_item
          //
          String inst (L"_xsd_" + item + L"_");
          String const& post (post_name (t));

          os << "void " << name << "::" << endl
             << "_xsd_parse_item (const " << string_type << "& v)"
             << "{"
             << "if (this->" << inst << ")"
             << "{"
             << "this->" << inst << "->pre ();"
             << "this->" << inst << "->_pre_impl ();"
             << "this->" << inst << "->_characters (v);"
             << "this->" << inst << "->_post_impl ();";

          if (ret_type (t) == L"void")
            os << "this->" << inst << "->" << post << " ();"
               << "this->" << item << " ();";
          else
            os << "this->" << item << " (this->" << inst << "->" <<
              post << " ());";

          os << "}"
             << "}";

          //
          //
          if (polymorphic)
          {
            String id (l.name ());

            if (String ns = xml_ns_name (l))
            {
              id += L' ';
              id += ns;
            }

            os << "const " << char_type << "* " << name << "::" << endl
               << "_static_type ()"
               << "{"
               << "return " << strlit (id) << ";"
               << "}";

            os << "const " << char_type << "* " << name << "::" << endl
               << "_dynamic_type () const"
               << "{"
               << "return _static_type ();"
               << "}";
          }
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
          String const& name (ename (u));
          String const& ret (ret_type (u));

          if (ret == L"void" || polymorphic)
          {
            os << "// " << name << endl
               << "//" << endl
               << endl;
          }

          if (ret == L"void")
          {
            os << "void " << name << "::" << endl
               << post_name (u) << " ()"
               << "{"
               << "}";
          }

          if (polymorphic)
          {
            String id (u.name ());

            if (String ns = xml_ns_name (u))
            {
              id += L' ';
              id += ns;
            }

            os << "const " << char_type << "* " << name << "::" << endl
               << "_static_type ()"
               << "{"
               << "return " << strlit (id) << ";"
               << "}";

            os << "const " << char_type << "* " << name << "::" << endl
               << "_dynamic_type () const"
               << "{"
               << "return _static_type ();"
               << "}";
          }
        }
      };

      //
      //
      struct StartElement: Traversal::Element, Context
      {
        StartElement (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& e)
        {
	  if (skip (e))
            return;

          Boolean poly (polymorphic && !anonymous (e.type ()));

          String const& inst (poly ? emember_cache (e) : emember (e));

          os << "if (";

          if (poly && e.global_p ())
            os << "(";

          if (e.qualified_p () && e.namespace_ ().name ())
          {
            os << "n == " << strlit (e.name ()) << " && " <<
              "ns == " << strlit (e.namespace_ ().name ());
          }
          else
          {
            os << "n == " << strlit (e.name ()) << " && ns.empty ()";
          }

          // Only a globally-defined element can be a subst-group root.
          //
          if (poly && e.global_p ())
          {
            os << ") ||" << endl
               << "::xsd::cxx::parser::substitution_map_instance< " <<
              char_type << " > ().check (" << endl
               << "ns, n, " << strlit (e.namespace_ ().name ()) <<
              ", " << strlit (e.name ()) << ", t)";
          }

          os << ")"
             << "{";

          if (poly)
          {
            SemanticGraph::Type& t (e.type ());

            // For pre-computing length.
            //
            String type_id (t.name ());

            if (String type_ns = xml_ns_name (t))
            {
              type_id += L' ';
              type_id += type_ns;
            }

            String fq_type (fq_name (t));
            String const& member (emember (e));
            String const& member_map (emember_map (e));

            os << "if (t == 0 && this->" << member << " != 0)" << endl
               << "this->" << inst << " = this->" << member << ";"
               << "else"
               << "{"
               << string_type << " ts (" << fq_type <<
              "::_static_type (), " << type_id.size () << "UL);"
               << endl
               << "if (t == 0)" << endl
               << "t = &ts;"
               << endl
               << "if (this->" << member << " != 0 && *t == ts)" << endl
               << "this->" << inst << " = this->" << member << ";"
               << "else if (this->" << member_map << " != 0)" << endl
               << "this->" << inst << " = dynamic_cast< " << fq_type <<
              "* > (" << endl
               << "this->" << member_map << "->find (*t));"
               << "else" << endl
               << "this->" << inst << " = 0;"
               << "}";
          }

          os << "this->" << complex_base << "::context_.top ().parser_ = " <<
            "this->" << inst << ";"
             << endl
             << "if (this->" << inst << ")" << endl
             << "this->" << inst << "->pre ();" // _start_element calls _pre
             << endl
             << "return true;"
             << "}";
        }
      };


      //
      //
      struct EndElement: Traversal::Element, Context
      {
        EndElement (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& e)
        {
	  if (skip (e))
            return;

          Boolean poly (polymorphic && !anonymous (e.type ()));

          String const& name (ename (e));
          String const& inst (poly ? emember_cache (e) : emember (e));

          os << "if (";

          if (poly && e.global_p ())
            os << "(";

          if (e.qualified_p () && e.namespace_ ().name ())
          {
            os << "n == " << strlit (e.name ()) << " && " <<
              "ns == " << strlit (e.namespace_ ().name ());
          }
          else
          {
            os << "n == " << strlit (e.name ()) << " && ns.empty ()";
          }

          // Only a globally-defined element can be a subst-group root.
          //
          if (poly && e.global_p ())
          {
            os << ") ||" << endl
               << "::xsd::cxx::parser::substitution_map_instance< " <<
              char_type << " > ().check (" << endl
               << "ns, n, " << strlit (e.namespace_ ().name ()) <<
              ", " << strlit (e.name ()) << ")";
          }

          os << ")"
             << "{";

          // _end_element calls post
          //

          SemanticGraph::Type& type (e.type ());
          String const& post (post_name (type));

          os << "if (this->" << inst << ")";

          if (ret_type (type) == L"void")
            os << "{"
               << "this->" << inst << "->" << post << " ();"
               << "this->" << name << " ();"
               << "}";
          else
            os << endl
               << "this->" << name << " (this->" << inst << "->" <<
              post << " ());"
               << endl;

          os << "return true;"
             << "}";
        }
      };

      //
      //
      struct Attribute: Traversal::Attribute, Context
      {
        Attribute (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& a)
        {
          String const& name (ename (a));
          String const& inst (emember (a));

          if (a.qualified_p () && a.namespace_ ().name ())
          {
            os << "if (n == " << strlit (a.name ()) << " && " <<
              "ns == " << strlit (a.namespace_ ().name ()) << ")"
               << "{";
          }
          else
          {
            os << "if (n == " << strlit (a.name ()) << " && ns.empty ())"
               << "{";
          }

          SemanticGraph::Type& type (a.type ());
          String const& post (post_name (type));
          String const& ret (ret_type (type));

          os << "if (this->" << inst << ")"
             << "{"
             << "this->" << inst << "->pre ();"
             << "this->" << inst << "->_pre_impl ();"
             << "this->" << inst << "->_characters (v);"
             << "this->" << inst << "->_post_impl ();";

          if (ret == L"void")
            os << "this->" << inst << "->" << post << " ();"
               << "this->" << name << " ();";
          else
            os << "this->" << name << " (this->" << inst << "->" <<
              post << " ());";

          os << "}"
             << "return true;"
             << "}";
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

          os << "void " << ename (m.scope ()) << "::" << endl
             << ename (m);

          if (arg == L"void")
            os << " ()";
          else
            os << " (" << arg << ")";

          os << "{"
             << "}";
        }
      };

      //
      //
      struct Complex: Traversal::Complex, Context
      {
        Complex (Context& c)
            : Context (c),
              parser_callback_ (c),
              start_element_ (c),
              end_element_ (c),
              attribute_ (c)
        {
          names_parser_callback_ >> parser_callback_;
          names_start_element_ >> start_element_;
          names_end_element_ >> end_element_;
          names_attribute_ >> attribute_;
        }

        virtual Void
        traverse (Type& c)
        {
          Boolean he (has<Traversal::Element> (c));
          Boolean ha (has<Traversal::Attribute> (c));

          String const& ret (ret_type (c));
          Boolean same (c.inherits_p () &&
                        ret == ret_type (c.inherits ().base ()));

          String const& name (ename (c));

          if ((he || ha || same || ret == L"void") || polymorphic)
          {
            os << "// " << name << endl
               << "//" << endl
               << endl;
          }

          if (polymorphic)
          {
            String id (c.name ());

            if (String ns = xml_ns_name (c))
            {
              id += L' ';
              id += ns;
            }

            os << "const " << char_type << "* " << name << "::" << endl
               << "_static_type ()"
               << "{"
               << "return " << strlit (id) << ";"
               << "}";

            os << "const " << char_type << "* " << name << "::" << endl
               << "_dynamic_type () const"
               << "{"
               << "return _static_type ();"
               << "}";

            if (c.inherits_p () && validation)
            {
              Boolean gen (!anonymous (c));

              // We normally don't need to enter anonymous types into
              // the inheritance map. The only exception is when an
              // anonymous types is defined inside an element that
              // is a member of a substitution group.
              //
              if (!gen)
              {
                // The first instance that this anonymous type classifies
                // is the prototype for others if any.
                //
                SemanticGraph::Instance& i (
                  c.classifies_begin ()->instance ());

                if (SemanticGraph::Element* e =
                    dynamic_cast<SemanticGraph::Element*> (&i))
                {
                  if (e->substitutes_p ())
                    gen = true;
                }
              }

              if (gen)
              {
                SemanticGraph::Type& base (c.inherits ().base ());

                os << "static" << endl
                   << "const ::xsd::cxx::parser::validating::inheritance_map_entry< " <<
                  char_type << " >" << endl
                   << "_xsd_" << name << "_inheritance_map_entry_ (" << endl
                   << name << "::_static_type ()," << endl
                   << fq_name (base) << "::_static_type ());"
                   << endl;
              }
            }
          }

          if (!(he || ha || same || ret == L"void"))
            return;

          // Parser callbacks.
          //
          if (!restriction_p (c))
            names (c, names_parser_callback_);

          if (same || ret == L"void")
          {
            os << ret << " " << name << "::" << endl
               << post_name (c) << " ()"
               << "{";

            if (same)
            {
              SemanticGraph::Type& base (c.inherits ().base ());

              if (ret == L"void")
                os << post_name (base) << " ();";
              else
                os << "return " << post_name (base) << " ();";
            }

            os << "}";
          }

          // The rest is parsing/validation code which is generated in
          // *-validation-source.cxx.
          //
          if (validation)
            return;

          // Don't use restriction_p here since we don't want special
          // treatment of anyType.
          //
          Boolean restriction (
            c.inherits_p () &&
            c.inherits ().is_a<SemanticGraph::Restricts> ());

          // _start_element_impl & _end_element_impl
          //
          if (he)
          {
            os << "bool " << name << "::" << endl
               << "_start_element_impl (const " << string_type << "& ns," << endl
               << "const " << string_type << "& n," << endl
               << "const " << string_type << "* t)"
               << "{"
               << "XSD_UNUSED (t);"
               << endl;

            if (!restriction)
            {
              os << "if (this->";

              if (c.inherits_p ())
                os << fq_name (c.inherits ().base ());
              else
                os << complex_base;

              os << "::_start_element_impl (ns, n, t))" << endl
                 << "return true;"
                 << endl;
            }

            names (c, names_start_element_);

            os << "return false;"
               << "}";


            // _end_element_impl
            //
            os << "bool " << name << "::" << endl
               << "_end_element_impl (const " << string_type << "& ns," << endl
               << "const " << string_type << "& n)"
               << "{";

            if (!restriction)
            {
              os << "if (this->";

              if (c.inherits_p () && !restriction)
                os << fq_name (c.inherits ().base ());
              else
                os << complex_base;

              os << "::_end_element_impl (ns, n))" << endl
                 << "return true;"
                 << endl;
            }

            names (c, names_end_element_);

            os << "return false;"
               << "}";
          }


          if (ha)
          {
            // _attribute_impl
            //
            os << "bool " << name << "::" << endl
               << "_attribute_impl (const " << string_type << "& ns," << endl
               << "const " << string_type << "& n," << endl
               << "const " << string_type << "& v)"
               << "{";

            if (!restriction)
            {
              os << "if (this->";

              if (c.inherits_p ())
                os << fq_name (c.inherits ().base ());
              else
                os << complex_base;

              os << "::_attribute_impl (ns, n, v))" << endl
                 << "return true;"
                 << endl;
            }

            names (c, names_attribute_);

            os << "return false;"
               << "}";
          }
        }

      private:
        //
        //
        ParserCallback parser_callback_;
        Traversal::Names names_parser_callback_;

        //
        //
        StartElement start_element_;
        Traversal::Names names_start_element_;

        //
        //
        EndElement end_element_;
        Traversal::Names names_end_element_;

        //
        //
        Attribute attribute_;
        Traversal::Names names_attribute_;
      };


      // Generate substitution group map entries.
      //
      struct GlobalElement: Traversal::Element, Context
      {
        GlobalElement (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& e)
        {
          if (e.substitutes_p ())
          {
            String name (escape (e.name ()));
            Type& r (e.substitutes ().root ());

            SemanticGraph::Type& type (e.type ());

            os << "// Substitution map entry for " << comment (e.name ()) << "." << endl
               << "//" << endl
               << "static" << endl
               << "const ::xsd::cxx::parser::substitution_map_entry< " <<
              char_type << " >" << endl
               << "_xsd_" << name << "_substitution_map_entry_ (" << endl
               << strlit (e.namespace_ ().name ()) << "," << endl
               << strlit (e.name ()) << "," << endl
               << strlit (r.namespace_ ().name ()) << "," << endl
               << strlit (r.name ()) << "," << endl
               << fq_name (type) << "::_static_type ());"
               << endl;
          }
        }
      };
    }

    Void
    generate_parser_source (Context& ctx)
    {
      if (ctx.polymorphic)
      {
        ctx.os << "#include <xsd/cxx/parser/substitution-map.hxx>" << endl;

        if (ctx.validation)
          ctx.os << "#include <xsd/cxx/parser/validating/inheritance-map.hxx>" << endl
                 << endl;
        else
          ctx.os << endl;

        Boolean import_maps (ctx.options.value<CLI::import_maps> ());
        Boolean export_maps (ctx.options.value<CLI::export_maps> ());

        if (import_maps || export_maps)
        {
          ctx.os << "#ifndef XSD_NO_EXPORT" << endl
                 << endl
                 << "namespace xsd"
                 << "{"
                 << "namespace cxx"
                 << "{"
                 << "namespace parser"
                 << "{"
                 << "#ifdef _MSC_VER" << endl;

          if (export_maps)
            ctx.os << "template struct __declspec (dllexport) " <<
              "substitution_map_init< " << ctx.char_type << " >;";

          if (import_maps)
            ctx.os << "template struct __declspec (dllimport) " <<
              "substitution_map_init< " << ctx.char_type << " >;";

          if (ctx.validation && export_maps)
            ctx.os << "template struct __declspec (dllexport) " <<
              "inheritance_map_init< " << ctx.char_type << " >;";

          if (ctx.validation && import_maps)
            ctx.os << "template struct __declspec (dllimport) " <<
              "inheritance_map_init< " << ctx.char_type << " >;";

          ctx.os << "#elif defined(__GNUC__) && __GNUC__ >= 4" << endl
                 << "template struct __attribute__ ((visibility(\"default\"))) " <<
            "substitution_map_init< " << ctx.char_type << " >;";

          if (ctx.validation)
            ctx.os << "template struct __attribute__ ((visibility(\"default\"))) " <<
              "inheritance_map_init< " << ctx.char_type << " >;";

          ctx.os << "#elif defined(XSD_MAP_VISIBILITY)" << endl
                 << "template struct XSD_MAP_VISIBILITY " <<
            "substitution_map_init< " << ctx.char_type << " >;";

          if (ctx.validation)
            ctx.os << "template struct XSD_MAP_VISIBILITY " <<
              "inheritance_map_init< " << ctx.char_type << " >;";

          ctx.os << "#endif" << endl
                 << "}"  // parser
                 << "}"  // cxx
                 << "}"  // xsd
                 << "#endif // XSD_NO_EXPORT" << endl
                 << endl;
        }

        ctx.os << "static" << endl
               << "const ::xsd::cxx::parser::substitution_map_init< " <<
          ctx.char_type << " >" << endl
               << "_xsd_substitution_map_init_;"
               << endl;

        if (ctx.validation)
        {
          ctx.os << "static" << endl
                 << "const ::xsd::cxx::parser::validating::inheritance_map_init< " <<
            ctx.char_type << " >" << endl
                 << "_xsd_inheritance_map_init_;"
                 << endl;
        }
      }

      // Emit "weak" header includes that are used in the file-per-type
      // compilation model.
      //
      if (ctx.options.value<CLI::generate_inline> ())
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
      Union union_ (ctx);
      Complex complex (ctx);
      Enumeration enumeration (ctx);
      GlobalElement global_element (ctx);

      names >> list;
      names >> union_;
      names >> complex;
      names >> enumeration;

      if (ctx.polymorphic)
        names >> global_element;

      schema.dispatch (ctx.schema_root);
    }
  }
}
