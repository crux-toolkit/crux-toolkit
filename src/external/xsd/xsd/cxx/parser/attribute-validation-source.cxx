// file      : xsd/cxx/parser/attribute-validation-source.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/parser/attribute-validation-source.hxx>

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

namespace CXX
{
  namespace Parser
  {
    namespace
    {
      struct Test: Traversal::Attribute,
                   Traversal::AnyAttribute,
                   Context
      {
        Test (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (SemanticGraph::Attribute& a)
        {
          String const& name (a.name ());

          if (a.qualified_p () && a.namespace_ ().name ())
          {
            String const& ns (a.namespace_ ().name ());

            os << "n == " << strlit (name) << " &&" << endl
               << "ns == " << strlit (ns);
          }
          else
            os << "n == " << strlit (name) << " && ns.empty ()";
        }

        virtual Void
        traverse (SemanticGraph::AnyAttribute& a)
        {
          String const& ns (a.definition_namespace ().name ());

          for (SemanticGraph::AnyAttribute::NamespaceIterator
                 i (a.namespace_begin ()), e (a.namespace_end ()); i != e;)
          {
            if (*i == L"##any")
            {
              os << "!n.empty ()";
            }
            else if (*i == L"##other")
            {
              if (ns)
              {
                // Note that here I assume that ##other does not include
                // unqualified names in a schema with target namespace.
                // This is not what the spec says but that seems to be
                // the consensus.
                //
                os << "(!ns.empty () && ns != " << strlit (ns) << ")";
              }
              else
                os << "!ns.empty ()";
            }
            else if (*i == L"##local")
            {
              os << "(ns.empty () && !n.empty ())";
            }
            else if (*i == L"##targetNamespace")
            {
              os << "ns == " << strlit (ns);
            }
            else
            {
              os << "ns == " << strlit (*i);
            }

            if (++i != e)
              os << " ||" << endl;
          }
        }
      };

      //
      //
      struct PhaseOne: Traversal::Attribute, Context
      {
        PhaseOne (Context& c)
            : Context (c), test_ (c)
        {
        }

        virtual Void
        traverse (Type& a)
        {
          String const& name (ename (a));
          String const& inst (emember (a));

          SemanticGraph::Type& type (a.type ());
          String const& post (post_name (type));
          String const& ret (ret_type (type));

          os << "if (";

          test_.traverse (a);

          os << ")"
             << "{"
             << "if (this->" << inst << ")"
             << "{"
             << "this->" << inst << "->pre ();"
             << "this->" << inst << "->_pre_impl ();"
             << "this->" << inst << "->_characters (s);"
             << "this->" << inst << "->_post_impl ();";

          if (ret == L"void")
            os << "this->" << inst << "->" << post << " ();"
               << "this->" << name << " ();";
          else
            os << arg_type (type) << " tmp (this->" << inst << "->" <<
              post << " ());"
               << "this->" << name << " (tmp);"
               << endl;

          os << "}";

          if (!a.optional_p ())
            os << "static_cast< v_state_attr_* > (" <<
              "this->v_state_attr_stack_.top ())->" << name << " = true;";

          os << "return true;"
             << "}";
        }

      private:
        Test test_;
      };


      //
      //
      struct PhaseTwo: Traversal::AnyAttribute, Context
      {
        PhaseTwo (Context& c)
            : Context (c), test_ (c)
        {
        }

        virtual Void
        traverse (Type& a)
        {
          os << "if (";

          test_.traverse (a);

          os << ")" << endl
             << "{"
             << "this->_any_attribute (ns, n, s);"
             << "return true;"
             << "}";
        }

      private:
        Test test_;
      };


      //
      //
      struct AttributeStateInit: Traversal::Attribute, Context
      {
        AttributeStateInit (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& a)
        {
          if (!a.optional_p ())
            os << "as." << ename (a) << " = false;";
        }
      };


      //
      //
      struct AttributeStateCheck: Traversal::Attribute, Context
      {
        AttributeStateCheck (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& a)
        {
          if (!a.optional_p ())
          {
            String ns (a.qualified_p () ? a.namespace_ ().name () : String ());

            os << "if (!as." << ename (a) << ")" << endl
               << "this->_expected_attribute (" << endl
               << strlit (ns) << ", " << strlit (a.name ()) << ");";
          }
        }
      };

      //
      //
      struct Complex: Traversal::Complex, Context
      {
        Complex (Context& c)
            : Context (c),
              phase_one_ (c),
              phase_two_ (c),
              attribute_state_init_ (c),
              attribute_state_check_ (c)
        {
          names_phase_one_ >> phase_one_;
          names_phase_two_ >> phase_two_;

          names_attribute_state_init_ >> attribute_state_init_;
          names_attribute_state_check_ >> attribute_state_check_;
        }

        virtual Void
        traverse (Type& c)
        {
          Boolean has_att (has<Traversal::Attribute> (c));
          Boolean has_any (has<Traversal::AnyAttribute> (c));

          if (!has_att && !has_any)
            return;

          Boolean has_req_att (false);
          if (has_att)
          {
            RequiredAttributeTest test (has_req_att);
            Traversal::Names names_test (test);
            names (c, names_test);
          }

          String const& name (ename (c));

          os <<"// Attribute validation and dispatch functions for " <<
            name << "." << endl
             <<"//" << endl;

          if (has_att)
          {
            // _attribute_impl_phase_one
            //
            os << "bool " << name << "::" << endl
               << "_attribute_impl_phase_one (const " << string_type <<
              "& ns," << endl
               << "const " << string_type << "& n," << endl
               << "const " << string_type << "& s)" << endl
               << "{";

            names (c, names_phase_one_);

            // Nothing matched - call our base (extension) or return false
            // if there is no base (or restriction (even from anyType)).
            //
            os << "return ";

            if (c.inherits_p () &&
                !c.inherits ().is_a<SemanticGraph::Restricts> ())
            {
              os << "this->" << fq_name (c.inherits ().base ()) <<
                "::_attribute_impl_phase_one (ns, n, s);";
            }
            else
              os << "false;";

            os << "}";
          }


          if (has_any)
          {
            // _attribute_impl_phase_two
            //
            os << "bool " << name << "::" << endl
               << "_attribute_impl_phase_two (const " << string_type <<
              "& ns," << endl
               << "const " << string_type << "& n," << endl
               << "const " << string_type << "& s)"
               << "{";

            names (c, names_phase_two_);

            // Nothing matched - call our base (extension) or return false
            // if there is no base (or restriction (even from anyType)).
            //
            os << "return ";

            if (c.inherits_p () &&
                !c.inherits ().is_a<SemanticGraph::Restricts> ())
            {
              os << "this->" << fq_name (c.inherits ().base ()) <<
                "::_attribute_impl_phase_two (ns, n, s);";
            }
            else
              os << "false;";

            os << "}";
          }

          if (has_req_att)
          {
            // _pre_a_validate
            //
            os << "void " << name << "::" << endl
               << "_pre_a_validate ()"
               << "{"
               << "this->v_state_attr_stack_.push ();"
               << "v_state_attr_& as = *static_cast< v_state_attr_* > (" <<
              "this->v_state_attr_stack_.top ());"
               << endl;

            names (c, names_attribute_state_init_);

            // Call our base (extension) last.
            //
            if (c.inherits_p () &&
                !c.inherits ().is_a<SemanticGraph::Restricts> ())
            {
              os << "this->" << fq_name (c.inherits ().base ()) <<
                "::_pre_a_validate ();";
            }

            os << "}";


            // _post_a_validate
            //
            os << "void " << name << "::" << endl
               << "_post_a_validate ()"
               << "{";

            // Call our base (extension) first.
            //
            if (c.inherits_p () &&
                !c.inherits ().is_a<SemanticGraph::Restricts> ())
            {
              os << "this->" << fq_name (c.inherits ().base ()) <<
                "::_post_a_validate ();"
                 << endl;
            }

            os << "v_state_attr_& as = *static_cast< v_state_attr_* > (" <<
              "this->v_state_attr_stack_.top ());"
               << endl;

            names (c, names_attribute_state_check_);

            os << endl
               << "this->v_state_attr_stack_.pop ();"
               << "}";
          }
        }

      private:
        PhaseOne phase_one_;
        Traversal::Names names_phase_one_;

        PhaseTwo phase_two_;
        Traversal::Names names_phase_two_;

        AttributeStateInit attribute_state_init_;
        Traversal::Names names_attribute_state_init_;

        AttributeStateCheck attribute_state_check_;
        Traversal::Names names_attribute_state_check_;
      };
    }

    Void
    generate_attribute_validation_source (Context& ctx)
    {
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
