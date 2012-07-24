// file      : xsd/cxx/parser/validator.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/parser/validator.hxx>

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

#include <cxx/parser/elements.hxx>

#include <iostream>

using std::wcerr;

namespace CXX
{
  namespace Parser
  {
    namespace
    {
      class ValidationContext: public Context
      {
      public:
        ValidationContext (SemanticGraph::Schema& root,
                           SemanticGraph::Path const& path,
                           CLI::Options const& options,
                           const WarningSet& disabled_warnings,
                           Boolean& valid_)
            : Context (std::wcerr, root, path, options, 0, 0, 0, 0),
              disabled_warnings_ (disabled_warnings),
              disabled_warnings_all_ (false),
              valid (valid_),
              subst_group_warning_issued (subst_group_warning_issued_),
              subst_group_warning_issued_ (false)
        {
        }

      public:
        Boolean
        is_disabled (Char const* w)
        {
          return disabled_warnings_all_ ||
            disabled_warnings_.find (w) != disabled_warnings_.end ();
        }

      public:
        String
        xpath (SemanticGraph::Nameable& n)
        {
          if (n.is_a<SemanticGraph::Namespace> ())
            return L"<namespace-level>"; // There is a bug if you see this.

          assert (n.named_p ());

          SemanticGraph::Scope& scope (n.scope ());

          if (scope.is_a<SemanticGraph::Namespace> ())
            return n.name ();

          return xpath (scope) + L"/" + n.name ();
        }

      protected:
        ValidationContext (ValidationContext& c)
            :  Context (c),
               disabled_warnings_ (c.disabled_warnings_),
               disabled_warnings_all_ (c.disabled_warnings_all_),
               valid (c.valid),
               subst_group_warning_issued (c.subst_group_warning_issued)
        {
        }

      protected:
        const WarningSet& disabled_warnings_;
        Boolean disabled_warnings_all_;
        Boolean& valid;
        Boolean& subst_group_warning_issued;
        Boolean subst_group_warning_issued_;
      };

      //
      //
      struct Any: Traversal::Any, ValidationContext
      {
        Any (ValidationContext& c)
            : ValidationContext (c)
        {
        }

        struct Element: Traversal::Element, ValidationContext
        {
          Element (ValidationContext& c, SemanticGraph::Any& any)
              : ValidationContext (c),
                any_ (any),
                ns_ (any.definition_namespace ().name ())
          {
          }

          virtual Void
          traverse (SemanticGraph::Element& e)
          {
	    if (skip (e)) return;

            using SemanticGraph::Any;

            Boolean q (e.qualified_p ());
            String ns (q ? e.namespace_ ().name () : "");

            for (Any::NamespaceIterator i (any_.namespace_begin ());
                 i != any_.namespace_end (); ++i)
            {
              Boolean failed (false);

              if (*i == L"##any")
              {
                failed = true;
              }
              else if (*i == L"##other")
              {
                if (ns_)
                {
                  // Note that here I assume that ##other does not
                  // include names without target namespace. This
                  // is not what the spec says but that seems to be
                  // the consensus.
                  //
                  failed = q && ns != ns_;
                }
                else
                {
                  // No target namespace.
                  //
                  failed = q && ns != L"";
                }
              }
              else if (*i == L"##local")
              {
                failed = !q || ns == L"";
              }
              else if (*i == L"##targetNamespace")
              {
                failed = (q && ns_ == ns) || (!q && ns_ == L"");
              }
              else
              {
                failed = q && *i == ns;
              }

              if (failed)
              {
                Any& a (any_);

                os << a.file () << ":" << a.line () << ":" << a.column ()
                   << ": warning P001: namespace '" << *i << "' allows for "
                   << "element '" << e.name () << "'" << endl;

                os << a.file () << ":" << a.line () << ":" << a.column ()
                   << ": warning P001: generated code may not associate element '"
                   << e.name () << "' correctly if it appears in place of "
                   << "this wildcard" << endl;

                os << e.file () << ":" << e.line () << ":" << e.column ()
                   << ": info: element '" << e.name () << "' is defined "
                   << "here" << endl;

                os << a.file () << ":" << a.line () << ":" << a.column ()
                   << ": info: turn on validation to ensure correct "
                   << "association" << endl;
              }
            }
          }

        private:
          SemanticGraph::Any& any_;
          String ns_;
        };

        struct Complex: Traversal::Complex
        {
          Complex ()
              : up_ (true), down_ (true)
          {
          }

          virtual Void
          post (Type& c)
          {
            // Go down the inheritance hierarchy.
            //
            if (down_)
            {
              Boolean up = up_;
              up_ = false;

              if (c.inherits_p ())
                dispatch (c.inherits ().base ());

              up_ = up;
            }

            // Go up the inheritance hierarchy.
            //
            if (up_)
            {
              Boolean down = down_;
              down_ = false;

              for (Type::BegetsIterator i (c.begets_begin ());
                   i != c.begets_end (); ++i)
              {
                dispatch (i->derived ());
              }

              down_ = down;
            }
          }

        private:
          Boolean up_, down_;
        };

        virtual Void
        traverse (SemanticGraph::Any& a)
        {
          using SemanticGraph::Compositor;

          // Find our complex type.
          //
          Compositor* c (&a.contained_particle ().compositor ());

          while(!c->contained_compositor_p ())
            c = &c->contained_particle ().compositor ();

          SemanticGraph::Complex& type (
            dynamic_cast<SemanticGraph::Complex&> (
              c->contained_compositor ().container ()));

          Complex complex;
          Traversal::Names names;
          Element element (*this, a);

          complex >> names >> element;

          complex.dispatch (type);
        }
      };


      //
      //
      struct Traverser: Traversal::Schema,
                        Traversal::Complex,
                        Traversal::Type,
                        Traversal::Element,
                        ValidationContext
      {
        Traverser (ValidationContext& c)
            : ValidationContext (c),
              any_ (c)
        {
          *this >> sources_ >> *this;
          *this >> schema_names_ >> ns_ >> names_ >> *this;

          // Any
          //
          if (!validation && !is_disabled ("P001"))
          {
            *this >> contains_compositor_ >> compositor_ >> contains_particle_;
            contains_particle_ >> compositor_;
            contains_particle_ >> any_;
          }
        }

        virtual Void
        traverse (SemanticGraph::Complex& c)
        {
          using SemanticGraph::Schema;

          traverse (static_cast<SemanticGraph::Type&> (c));

          if (c.inherits_p ())
          {
            SemanticGraph::Type& t (c.inherits ().base ());

            if (t.named_p () &&
                types_.find (
                  t.scope ().name () + L"#" + t.name ()) == types_.end ())
            {
              // Don't worry about types that are in included/imported
              // schemas.
              //
              Schema& s (dynamic_cast<Schema&> (t.scope ().scope ()));

              if (&s == &schema_root || sources_p (schema_root, s))
              {
                valid = false;

                wcerr << c.file () << ":" << c.line () << ":" << c.column ()
                      << ": error: type '" << xpath (c) << "' inherits from "
                      << "yet undefined type '" << xpath (t) << "'" << endl;

                wcerr << t.file () << ":" << t.line () << ":" << t.column ()
                      << ": info: '" << xpath (t) << "' is defined here"
                      << endl;

                wcerr << c.file () << ":" << c.line () << ":" << c.column ()
                      << ": info: inheritance from a yet-undefined type is "
                      << "not supported" << endl;

                wcerr << c.file () << ":" << c.line () << ":" << c.column ()
                      << ": info: re-arrange your schema and try again"
                      << endl;
              }
            }
          }

          Complex::traverse (c);
        }

        virtual Void
        traverse (SemanticGraph::Type& t)
        {
          if (t.named_p ())
          {
            types_.insert (t.scope ().name () + L"#" + t.name ());
          }
        }

        virtual Void
        traverse (SemanticGraph::Element& e)
        {
          if (is_disabled ("P002"))
            return;

          if (e.substitutes_p () &&
              !options.value<CLI::generate_polymorphic> () &&
              !subst_group_warning_issued)
          {
            subst_group_warning_issued = true;

            os << e.file () << ":" << e.line () << ":" << e.column ()
               << ": warning P002: substitution groups are used but "
               << "--generate-polymorphic was not specified" << endl;

            os << e.file () << ":" << e.line () << ":" << e.column ()
               << ": info: generated code may not be able to parse "
               << "some conforming instances" << endl;
          }
        }

        // Return true if root sources s.
        //
        Boolean
        sources_p (SemanticGraph::Schema& root, SemanticGraph::Schema& s)
        {
          using SemanticGraph::Schema;
          using SemanticGraph::Sources;

          for (Schema::UsesIterator i (root.uses_begin ());
               i != root.uses_end (); ++i)
          {
            if (i->is_a<Sources> ())
            {
              if (&i->schema () == &s || sources_p (i->schema (), s))
                return true;
            }
          }

          return false;
        }

      private:
        Containers::Set<String> types_;

        Traversal::Sources sources_;

        Traversal::Names schema_names_;
        Traversal::Namespace ns_;

        Traversal::Names names_;

        // Any.
        //
        Any any_;
        Traversal::Compositor compositor_;
        Traversal::ContainsParticle contains_particle_;
        Traversal::ContainsCompositor contains_compositor_;
      };


      struct AnonymousType: Traversal::Schema,
                            Traversal::Complex,
                            Traversal::Element,
                            Traversal::Attribute,
                            ValidationContext
      {
        AnonymousType (ValidationContext& c)
            : ValidationContext (c),
              anonymous_error_issued_ (false)
        {
          *this >> sources_ >> *this;
          *this >> schema_names_ >> ns_ >> names_ >> *this;
          *this >> names_;
        }

        Boolean
        traverse_common (SemanticGraph::Member& m)
        {
          SemanticGraph::Type& t (m.type ());

          if (!t.named_p ()
              && !t.is_a<SemanticGraph::Fundamental::IdRef> ()
              && !t.is_a<SemanticGraph::Fundamental::IdRefs> ())
          {
            if (!anonymous_error_issued_)
            {
              valid = false;
              anonymous_error_issued_ = true;

              wcerr << t.file ()
                    << ": error: anonymous types detected"
                    << endl;

              wcerr << t.file ()
                    << ": info: "
                    << "anonymous types are not supported in this mapping"
                    << endl;

              wcerr << t.file ()
                    << ": info: consider explicitly naming these types or "
                    << "remove the --preserve-anonymous option to "
                    << "automatically name them"
                    << endl;

              if (!options.value<CLI::show_anonymous> ())
                wcerr << t.file ()
                      << ": info: use --show-anonymous option to see these "
                      << "types" << endl;
            }

            return true;
          }

          return false;
        }

        virtual Void
        traverse (SemanticGraph::Element& e)
        {
	  if (skip (e)) return;

          if (traverse_common (e))
          {
            if (options.value<CLI::show_anonymous> ())
            {
              wcerr << e.file () << ":" << e.line () << ":" << e.column ()
                    << ": error: element '" << xpath (e) << "' "
                    << "is of anonymous type" << endl;
            }
          }
          else
            Traversal::Element::traverse (e);
        }

        virtual Void
        traverse (SemanticGraph::Attribute& a)
        {
          if (traverse_common (a))
          {
            if (options.value<CLI::show_anonymous> ())
            {
              wcerr << a.file () << ":" << a.line () << ":" << a.column ()
                    << ": error: attribute '" << xpath (a) << "' "
                    << "is of anonymous type" << endl;
            }
          }
          else
            Traversal::Attribute::traverse (a);
        }

      private:
        Boolean anonymous_error_issued_;

        Containers::Set<String> types_;

        Traversal::Sources sources_;

        Traversal::Names schema_names_;
        Traversal::Namespace ns_;

        Traversal::Names names_;
      };

      struct GlobalElement: Traversal::Element, ValidationContext
      {
        GlobalElement (ValidationContext& c, SemanticGraph::Element*& element)
            : ValidationContext (c), element_ (element)
        {
        }

        virtual Void
        traverse (Type& e)
        {
          if (!valid)
            return;

          if (options.value<CLI::root_element_first> ())
          {
            if (element_ == 0)
              element_ = &e;
          }
          else if (options.value<CLI::root_element_last> ())
          {
            element_ = &e;
          }
          else if (String name = options.value<CLI::root_element> ())
          {
            if (e.name () == name)
              element_ = &e;
          }
          else
          {
            if (element_ == 0)
              element_ = &e;
            else
            {
              wcerr << schema_root.file () << ": error: unable to generate "
                    << "the test driver without a unique document root"
                    << endl;

              wcerr << schema_root.file () << ": info: use --root-element-* "
                    << "options to specify the document root" << endl;

              valid = false;
            }
          }
        }

      private:
        SemanticGraph::Element*& element_;
      };

    }

    Validator::
    Validator ()
    {
      // Dummy ctor, helps with long symbols on HP-UX.
    }

    Boolean Validator::
    validate (CLI::Options const& options,
              SemanticGraph::Schema& root,
              SemanticGraph::Path const& path,
              Boolean gen_driver,
              const WarningSet& disabled_warnings)
    {
      Boolean valid (true);
      ValidationContext ctx (root, path, options, disabled_warnings, valid);

      //
      //
      if (options.value<CLI::char_type> () != "char" &&
          options.value<CLI::char_type> () != "wchar_t" &&
          !ctx.is_disabled ("P003"))
      {
        wcerr << "warning P003: unknown base character type '" <<
          options.value<CLI::char_type> ().c_str () << "'" << endl;
      }

      //
      //
      if (options.value<CLI::xml_parser> () != "xerces" &&
          options.value<CLI::xml_parser> () != "expat" &&
          !ctx.is_disabled ("P004"))
      {
        wcerr << "warning P004: unknown underlying XML parser '" <<
          options.value<CLI::xml_parser> ().c_str () << "'" << endl;
      }

      //
      //
      if (options.value<CLI::xml_parser> () == "expat" &&
          options.value<CLI::char_type> () == "wchar_t")
      {
        wcerr << "error: using expat with wchar_t is not supported"
              << endl;

        return false;
      }

      //
      //
      if (options.value<CLI::xml_parser> () == "expat" &&
          !options.value<CLI::char_encoding> ().empty () &&
          options.value<CLI::char_encoding> () != "utf8")
      {
        wcerr << "error: using expat with character encoding other than "
              << "utf8 is not supported"
              << endl;

        return false;
      }

      //
      //
      if (options.value<CLI::generate_validation> () &&
          options.value<CLI::suppress_validation> ())
      {
        wcerr << "error: mutually exclusive options specified: "
              << "--generate-validation and --suppress-validation"
              << endl;

        return false;
      }

      //
      //
      if (options.value<CLI::generate_noop_impl> () &&
          options.value<CLI::generate_print_impl> ())
      {
        wcerr << "error: mutually exclusive options specified: "
              << "--generate-noop-impl and --generate-print-impl"
              << endl;

        return false;
      }

      //
      //
      {
        Boolean ref (options.value<CLI::root_element_first> ());
        Boolean rel (options.value<CLI::root_element_last> ());
        Boolean re (options.value<CLI::root_element> ());

        if ((ref && rel) || (ref && re) || (rel && re))
        {
          wcerr << "error: mutually exclusive options specified: "
                << "--root-element-last, --root-element-first, and "
                << "--root-element"
                << endl;

          return false;
        }
      }

      //
      //
      Boolean import_maps (options.value<CLI::import_maps> ());
      Boolean export_maps (options.value<CLI::export_maps> ());

      if (import_maps && export_maps)
      {
        wcerr << "error: --import-maps and --export-maps are "
              << "mutually exclusive" << endl;

        return false;
      }

      if (import_maps && !ctx.polymorphic)
      {
        wcerr << "error: --import-maps can only be specified together with "
              << "--generate-polymorphic" << endl;

        return false;
      }

      if (export_maps && !ctx.polymorphic)
      {
        wcerr << "error: --export-maps can only be specified together with "
              << "--generate-polymorphic" << endl;

        return false;
      }

      // Test for anonymout types.
      //
      {
        AnonymousType traverser (ctx);
        traverser.dispatch (root);
      }


      // Test the rest.
      //
      if (valid)
      {
        Traverser traverser (ctx);
        traverser.dispatch (root);
      }

      // Test that the document root is unique.
      //
      if (valid && gen_driver)
      {
        SemanticGraph::Element* element (0);

        Traversal::Schema schema;
        Traversal::Sources sources;

        schema >> sources >> schema;

        Traversal::Names schema_names;
        Traversal::Namespace ns;
        Traversal::Names ns_names;
        GlobalElement global_element (ctx, element);

        schema >> schema_names >> ns >> ns_names >> global_element;

        schema.dispatch (root);

        if (valid && element == 0)
        {
          wcerr << root.file () << ": error: unable to generate the "
                << "test driver without a global element (document root)"
                << endl;

          valid = false;
        }
      }

      return valid;
    }
  }
}
