// file      : xsde/cxx/tree/polymorphism-processor.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/tree/elements.hxx>
#include <cxx/tree/polymorphism-processor.hxx>

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

#include <cult/containers/set.hxx>

#include <iostream>

using std::wcerr;

namespace CXX
{
  namespace Tree
  {
    namespace
    {
      struct TypeSet
      {
        template <typename I>
        TypeSet (I begin, I end)
        {
          for (; begin != end; ++begin)
            insert (*begin);
        }

        Void
        insert (String const& name)
        {
          Size p (name.rfind ('#'));

          if (p == String::npos)
            unames_.insert (name);
          else
            qnames_.insert (name);
        }

        Boolean
        find (SemanticGraph::Type& t)
        {
          if (!unames_.empty ())
          {
            if (unames_.find (t.name ()) != unames_.end ())
              return true;
          }

          if (!qnames_.empty ())
          {
            if (qnames_.find (t.scope ().name () + L"#" + t.name ()) !=
                qnames_.end ())
              return true;
          }

          return false;
        }

      private:
        typedef Cult::Containers::Set<String> StringSet;

        StringSet unames_;
        StringSet qnames_;
      };


      //
      //
      struct Type: Traversal::Type,
                   Traversal::Complex
      {
        Type (TypeSet& poly_types)
            : poly_types_ (poly_types)
        {
        }

        virtual Void
        traverse (SemanticGraph::Type& t)
        {
          SemanticGraph::Context& ctx (t.context ());

          if (!ctx.count ("polymorphic"))
            ctx.set ("polymorphic", poly_types_.find (t));
        }

        virtual Void
        traverse (SemanticGraph::Complex& c)
        {
          SemanticGraph::Context& ctx (c.context ());

          if (!ctx.count ("polymorphic"))
          {
            // First check our base.
            //
            Boolean pb (false);
            if (c.inherits_p ())
            {
              SemanticGraph::Type& b (c.inherits ().base ());

              if (!b.context ().count ("polymorphic"))
                dispatch (b);

              pb = b.context ().get<Boolean> ("polymorphic");
            }

            ctx.set ("polymorphic", pb || poly_types_.find (c));
          }
        }

      private:
        TypeSet& poly_types_;
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
                       Traversal::Fundamental::Entities
      {
        FundType (TypeSet& poly_types, Boolean& valid)
            : poly_types_ (poly_types), valid_ (valid)
        {
        }

        Void
        check (SemanticGraph::Type& t, Boolean fund)
        {
          SemanticGraph::Context& ctx (t.context ());

          if (poly_types_.find (t))
          {
            if (!fund)
              ctx.set ("polymorphic", true);
            else
            {
              wcerr << "error: built-in type '" << t.name () << "' cannot "
                    << "be polymorphic because it is mapped to a fundamental "
                    << "C++ type"
                    << endl;

              valid_ = false;
            }
          }
          else
            ctx.set ("polymorphic", false);
        }

        // anyType & anySimpleType.
        //
        virtual Void
        traverse (SemanticGraph::AnyType& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::AnySimpleType& t)
        {
          check (t, false);
        }

        // Boolean.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Boolean& t)
        {
          check (t, true);
        }

        // Integral types.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Byte& t)
        {
          check (t, true);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::UnsignedByte& t)
        {
          check (t, true);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Short& t)
        {
          check (t, true);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::UnsignedShort& t)
        {
          check (t, true);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Int& t)
        {
          check (t, true);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::UnsignedInt& t)
        {
          check (t, true);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Long& t)
        {
          check (t, true);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::UnsignedLong& t)
        {
          check (t, true);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Integer& t)
        {
          check (t, true);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NonPositiveInteger& t)
        {
          check (t, true);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NonNegativeInteger& t)
        {
          check (t, true);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::PositiveInteger& t)
        {
          check (t, true);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NegativeInteger& t)
        {
          check (t, true);
        }

        // Floats.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Float& t)
        {
          check (t, true);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Double& t)
        {
          check (t, true);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Decimal& t)
        {
          check (t, true);
        }

        // Strings.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::String& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NormalizedString& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Token& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NameToken& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NameTokens& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Name& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::NCName& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Language& t)
        {
          check (t, false);
        }


        // Qualified name.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::QName& t)
        {
          check (t, false);
        }


        // ID/IDREF.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Id& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::IdRef& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::IdRefs& t)
        {
          check (t, false);
        }

        // URI.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::AnyURI& t)
        {
          check (t, false);
        }

        // Binary.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Base64Binary& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::HexBinary& t)
        {
          check (t, false);
        }


        // Date/time.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Date& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::DateTime& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Duration& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Day& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Month& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::MonthDay& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Year& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::YearMonth& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Time& t)
        {
          check (t, false);
        }

        // Entity.
        //
        virtual Void
        traverse (SemanticGraph::Fundamental::Entity& t)
        {
          check (t, false);
        }

        virtual Void
        traverse (SemanticGraph::Fundamental::Entities& t)
        {
          check (t, false);
        }

      private:
        TypeSet& poly_types_;
        Boolean& valid_;
      };

      struct GlobalElement: Traversal::Element
      {
        GlobalElement (TypeSet& poly_types,
                       Boolean& valid,
                       const WarningSet& disabled_warnings)
            : poly_types_ (poly_types), valid_ (valid), warning_ (true)
        {
          if (disabled_warnings.find ("all") != disabled_warnings.end () ||
              disabled_warnings.find ("T005") != disabled_warnings.end ())
            warning_ = false;
        }

        virtual Void
        traverse (Type& e)
        {
          using SemanticGraph::Schema;

          if (!e.substitutes_p ())
            return;

          // If we are a substitution for some element, then mark
          // that element's type as polymorphic.
          //
          Type& r (e.substitutes ().root ());
          SemanticGraph::Type& rt (r.type ());
          SemanticGraph::Context& ctx (rt.context ());

          // We may need to override the previous value.
          //
          if (ctx.count ("polymorphic") && ctx.get<Boolean> ("polymorphic"))
            return;

          // Built-in types that are mapped to fundamental types cannot
          // be declared polymorphic.
          //
          Boolean fund (false);
          IsFundamentalType test (fund);
          test.dispatch (rt);

          if (fund)
          {
            wcerr << r.file () << ":" << r.line () << ":" << r.column ()
                  << ": error: built-in type '" << rt.name () << "' "
                  << "is mapped to a fundamental C++ type and is expected "
                  << "to be polymorphic" << endl;

            wcerr << e.file () << ":" << e.line () << ":" << e.column ()
                  << ": info: because type '" << rt.name () << "' is "
                  << "used in a substitution group declared here" << endl;

            valid_ = false;
            return;
          }

          ctx.set ("polymorphic", true);

          if (!warning_)
            return;

          Schema& es (dynamic_cast<Schema&> (e.scope ().scope ()));
          Schema& rts (dynamic_cast<Schema&> (rt.scope ().scope ()));

          // If the root type and this element are in different schemas
          // and the root type is not explicitly marked as polymorphic,
          // then issue a warning.
          //
          if (&es != &rts && !sources_p (es, rts) && !poly_types_.find (rt))
          {
            wcerr << rt.file () << ":" << rt.line () << ":" << rt.column ()
                  << ": warning T005: assuming type '" << rt.name () << "' "
                  << "is polymorphic" << endl;

            wcerr << e.file () << ":" << e.line () << ":" << e.column ()
                  << ": info: because type '" << rt.name () << "' is "
                  << "used in a substitution group declared here" << endl;

            wcerr << rt.file () << ":" << rt.line () << ":" << rt.column ()
                  << ": info: use --polymorphic-type to indicate this type "
                  << "is polymorphic when compiling schemas that "
                  << "reference it" << endl;
          }
        }

      private:
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
        TypeSet& poly_types_;
        Boolean& valid_;
        Boolean warning_;
      };

      // Go into sourced/included/imported schemas while making sure
      // we don't process the same stuff more than once.
      //
      struct Uses: Traversal::Sources,
                   Traversal::Includes,
                   Traversal::Imports
      {
        Uses (Char const* seen_key)
            : seen_key_ (seen_key)
        {
        }

        virtual Void
        traverse (SemanticGraph::Sources& sr)
        {
          SemanticGraph::Schema& s (sr.schema ());

          if (!s.context ().count (seen_key_))
          {
            s.context ().set (seen_key_, true);
            Traversal::Sources::traverse (sr);
          }
        }

        virtual Void
        traverse (SemanticGraph::Includes& i)
        {
          SemanticGraph::Schema& s (i.schema ());

          if (!s.context ().count (seen_key_))
          {
            s.context ().set (seen_key_, true);
            Traversal::Includes::traverse (i);
          }
        }

        virtual Void
        traverse (SemanticGraph::Imports& i)
        {
          SemanticGraph::Schema& s (i.schema ());

          if (!s.context ().count (seen_key_))
          {
            s.context ().set (seen_key_, true);
            Traversal::Imports::traverse (i);
          }
        }

      private:
        Char const* seen_key_;
      };

      Char const* pass_one_key = "cxx-tree-polymorphism-processor-seen-one";
      Char const* pass_two_key = "cxx-tree-polymorphism-processor-seen-two";

      Boolean
      process_impl (CLI::Options const& ops,
                    SemanticGraph::Schema& tu,
                    SemanticGraph::Path const&,
                    const WarningSet& disabled_warnings)
      {
        Boolean valid (true);

        // Prepare a set of polymorphic types.
        //

        TypeSet poly_types (ops.value<CLI::polymorphic_type> ().begin (),
                            ops.value<CLI::polymorphic_type> ().end ());

        // Root schema in the file-per-type mode is just a bunch
        // of includes without a namespace.
        //
        SemanticGraph::Schema::NamesIterator i (tu.names_begin ());

        if (i != tu.names_end () &&
            i->named ().name () == L"http://www.w3.org/2001/XMLSchema")
        {
          // XML Schema namespace.
          //
          Traversal::Schema schema;

          Traversal::Names schema_names;
          Traversal::Namespace ns;
          Traversal::Names ns_names;
          FundType fund_type (poly_types, valid);

          schema >> schema_names >> ns >> ns_names >> fund_type;

          schema.dispatch (tu);
        }
        else
        {
          // First handle fundamental types.
          //
          {
            Traversal::Schema schema;
            Traversal::Implies implies;
            Traversal::Schema xs_schema;

            schema >> implies >> xs_schema;

            Traversal::Names xs_schema_names;
            Traversal::Namespace ns;
            Traversal::Names ns_names;
            FundType fund_type (poly_types, valid);

            xs_schema >> xs_schema_names >> ns >> ns_names >> fund_type;

            schema.dispatch (tu);
          }

          // Note that we check first if this schema has already been
          // processed which may happen in the file-per-type compilation
          // mode.
          //
          if (!tu.context ().count (pass_two_key))
          {
            // Pass one - check substitution groups.
            //
            if (valid)
            {
              Traversal::Schema schema;
              Uses uses (pass_one_key);

              schema >> uses >> schema;

              Traversal::Names schema_names;
              Traversal::Namespace ns;
              Traversal::Names ns_names;
              GlobalElement element (poly_types, valid, disabled_warnings);

              schema >> schema_names >> ns >> ns_names >> element;

              // Some twisted schemas do recusive self-inclusion.
              //
              tu.context ().set (pass_one_key, true);

              schema.dispatch (tu);
            }

            // Pass two - process types.
            //
            if (valid)
            {
              Traversal::Schema schema;
              Uses uses (pass_two_key);

              schema >> uses >> schema;

              Traversal::Names schema_names;
              Traversal::Namespace ns;
              Traversal::Names ns_names;
              Type type (poly_types);

              schema >> schema_names >> ns >> ns_names >> type;

              // Some twisted schemas do recusive self-inclusion.
              //
              tu.context ().set (pass_two_key, true);

              schema.dispatch (tu);
            }
          }
        }

        return valid;
      }
    }

    Boolean PolymorphismProcessor::
    process (CLI::Options const& ops,
             SemanticGraph::Schema& tu,
             SemanticGraph::Path const& file,
             const WarningSet& disabled_warnings)
    {
      return process_impl (ops, tu, file, disabled_warnings);
    }
  }
}
