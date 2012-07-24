// file      : xsd/cxx/parser/element-validation-source.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <cxx/parser/element-validation-source.hxx>

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

#include <cult/containers/vector.hxx>

namespace CXX
{
  namespace Parser
  {
    namespace
    {
      typedef Cult::Containers::Vector<SemanticGraph::Particle*> Particles;


      //
      //
      struct ParticleTest: Traversal::Compositor,
                           Traversal::Element,
                           Traversal::Any,
                           Context
      {
        ParticleTest (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (SemanticGraph::Element& e)
        {
          String const& name (e.name ());

          if (polymorphic && e.global_p ())
            os << "(";

          if (e.qualified_p () && e.namespace_ ().name ())
          {
            String const& ns (e.namespace_ ().name ());

            os << "n == " << strlit (name) << " &&" << endl
               << "ns == " << strlit (ns);
          }
          else
            os << "n == " << strlit (name) << " && ns.empty ()";


          // Only a globally-defined element can be a subst-group root.
          //
          if (polymorphic && e.global_p ())
          {
            os << ") ||" << endl
               << "::xsd::cxx::parser::substitution_map_instance< " <<
              char_type << " > ().check (" << endl
               << "ns, n, " << strlit (e.namespace_ ().name ()) <<
              ", " << strlit (name) << ", t)";
          }
        }

        virtual Void
        traverse (SemanticGraph::Any& a)
        {
          String const& ns (a.definition_namespace ().name ());

          // Note that we need to make sure the "flush" element (both name
          // and namespace are empty) does not match any compositor.
          //
          for (SemanticGraph::Any::NamespaceIterator i (a.namespace_begin ()),
                 e (a.namespace_end ()); i != e;)
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

        virtual Void
        traverse (SemanticGraph::Compositor& c)
        {
          // This compositor should already have been tested for
          // triviality (empty).
          //
          Particles const& p (c.context ().get<Particles> ("prefixes"));

          Boolean paren (p.size () != 1);

          for (Particles::ConstIterator i (p.begin ()), e (p.end ());
               i != e;)
          {
            if (paren)
              os << "(";

            dispatch (**i);

            if (paren)
              os << ")";

            if (++i != e)
              os << " ||" << endl;
          }
        }
      };


      // Generates particle namespace-name pair. Used to generate
      // the _expected_element call.
      //
      struct ParticleName: Traversal::Compositor,
                           Traversal::Element,
                           Traversal::Any,
                           Context
      {
        ParticleName (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (SemanticGraph::Element& e)
        {
          String ns (e.qualified_p () ? e.namespace_ ().name () : String ());
          os << strlit (ns) << ", " << strlit (e.name ());
        }

        virtual Void
        traverse (SemanticGraph::Any& a)
        {
          String const& ns (*a.namespace_begin ());
          os << strlit (ns) << ", " << L << "\"*\"";
        }

        virtual Void
        traverse (SemanticGraph::Compositor& c)
        {
          Particles const& p (c.context ().get<Particles> ("prefixes"));

          dispatch (**p.begin ());
        }
      };


      // Common base for the ParticleIn{All, Choice, Sequence} treversers.
      //
      struct ParticleInCompositor: Context
      {
      protected:
        ParticleInCompositor (Context& c, SemanticGraph::Complex& type)
            : Context (c), type_ (type), particle_test_ (c)
        {
        }


        // Generate sub-parser setup code as well as the pre/post calls.
        //
        Void
        pre_post_calls (SemanticGraph::Particle& p)
        {
          using SemanticGraph::Element;
          using SemanticGraph::Complex;

          if (Element* e = dynamic_cast<Element*> (&p))
          {
            SemanticGraph::Type& type (e->type ());
            Boolean poly (polymorphic && !anonymous (type));

            String name, inst, def_parser, map;

            if (e->context ().count("name"))
            {
              name = ename (*e);
              inst = poly ? emember_cache (*e) : emember (*e);

              if (poly)
              {
                def_parser = emember (*e);
                map = emember_map (*e);
              }
            }
            else
            {
              // This is the subsequent mentioning of this element in the
              // content. We need to find the first one in order to get
              // to the escaped names.
              //
              Complex::NamesIteratorPair ip (type_.find (e->name ()));
              assert (ip.first != ip.second);
              Element& fe (dynamic_cast<Element&>(ip.first->named ()));

              name = ename (fe);
              inst = poly ? emember_cache (fe) : emember (fe);

              if (poly)
              {
                def_parser = emember (fe);
                map = emember_map (fe);
              }
            }

            if (poly)
            {
              // For pre-computing length.
              //
              String type_id (type.name ());

              if (String type_ns = xml_ns_name (type))
              {
                type_id += L' ';
                type_id += type_ns;
              }

              os << "if (t == 0 && this->" << def_parser << " != 0)" << endl
                 << "this->" << inst << " = this->" << def_parser << ";"
                 << "else"
                 << "{"
                 << string_type << " ts (" << fq_name (type) <<
                "::_static_type (), " << type_id.size () << "UL);"
                 << endl
                 << "if (t == 0)" << endl
                 << "t = &ts;"
                 << endl
                 << "if (this->" << def_parser << " != 0 && *t == ts)" << endl
                 << "this->" << inst << " = this->" << def_parser << ";"
                 << "else"
                 << "{";

              // Check that the types are related by inheritance.
              //
              os << "if (t != &ts &&" << endl
                 << "!::xsd::cxx::parser::validating::" <<
                "inheritance_map_instance< " << char_type <<
                " > ().check (" << endl
                 << "t->data (), ts))" << endl
                 << "throw ::xsd::cxx::parser::dynamic_type< " << char_type <<
                " > (*t);"
                 << endl
                 << "if (this->" << map << " != 0)" << endl
                 << "this->" << inst << " = dynamic_cast< " <<
                fq_name (type) << "* > (" << endl
                 << "this->" << map << "->find (*t));"
                 << "else" << endl
                 << "this->" << inst << " = 0;"
                 << "}"
                 << "}";
            }

            os << "this->" << complex_base << "::context_.top ()." <<
              "parser_ = this->" << inst << ";"
               << endl
               << "if (this->" << inst << ")" << endl
               << "this->" << inst << "->pre ();"
               << "}"
               << "else" // start
               << "{"
               << "if (this->" << inst << ")"
               << "{";

            String const& ret (ret_type (type));
            String const& post (post_name (type));

            if (ret == L"void")
              os << "this->" << inst << "->" << post << " ();"
                 << "this->" << name << " ();";
            else
              os << arg_type (type) << " tmp (this->" << inst << "->" <<
                post << " ());"
                 << "this->" << name << " (tmp);";

            os << "}";
          }
          else
          {
            os << "this->_start_any_element (ns, n, t);"
               << "this->" << complex_base << "::context_.top ().any_ = true;"
               << "}"
               << "else" // start
               << "{"
               << "this->" << complex_base << "::context_.top ().any_ = false;"
               << "this->_end_any_element (ns, n);";
          }
        }

      protected:
        SemanticGraph::Complex& type_;
        ParticleTest particle_test_;
      };



      // The 'all' compositor can only contain elements with min={0,1}, max=1.
      //
      struct ParticleInAll: Traversal::Element,
                            ParticleInCompositor
      {
        ParticleInAll (Context& c, SemanticGraph::Complex& type)
            : ParticleInCompositor (c, type)
        {
        }

        virtual Void
        traverse (SemanticGraph::Element& e)
        {
          UnsignedLong state (e.context ().get<UnsignedLong> ("state"));

          if (state != 0)
            os << "else ";

          os << "if (";

          particle_test_.traverse (e);

          os << ")"
             << "{"
             << "if (count[" << state << "UL] == 0)"
             << "{"
             << "if (start)"
             << "{";

          pre_post_calls (e);

          os << "count[" << state << "UL] = 1;"
             << "}"
             << "}"
             << "else" // count != 0
             << "{"
             << "assert (start);" // Assuming well-formed XML.

            // Since there is never more content after 'all', we could have
            // as well thrown here. But instead we will let the code in
            // start_element handle this along with other unexpected
            // elements.
            //
             << "state = ~0UL;"
             << "}"
             << "}";
        }
      };


      //
      //
      struct ParticleInChoice: Traversal::Particle,
                               Traversal::Compositor,
                               ParticleInCompositor
      {
        ParticleInChoice (Context& c, SemanticGraph::Complex& type)
            : ParticleInCompositor (c, type), particle_name_ (c)
        {
        }

        virtual Void
        traverse (SemanticGraph::Particle& p)
        {
          using SemanticGraph::Element;

          UnsignedLong state (p.context ().get<UnsignedLong> ("state"));

          UnsignedLong min (p.min ()), max (p.max ());

          os << "case " << state << "UL:" << endl
             << "{";

          if (max != 1) // We don't need the test if max == 1.
          {
            os << "if (";

            particle_test_.dispatch (p);

            os << ")"
               << "{";
          }

          os << "if (start)"
             << "{";

          pre_post_calls (p);

          switch (max)
          {
          case 0:
            {
              os << "count++;";
              break;
            }
          case 1:
            {
              // We do not need to increment count because min <= max and
              // we do not generate min check for min <= 1 (see below).
              //
              os << "state = ~0UL;";
              break;
            }
          default:
            {
              os << "if (++count == " << max << "UL)" << endl
                 << "state = ~0UL;";
            }
          };

          os << "}"; // start

          // We've already moved to the final state if max == 1.
          //
          if (max != 1)
          {
            os << "}"
               << "else"
               << "{"
               << "assert (start);"; // Assuming well-formed XML

            // Check if min cardinality requirements have been met. Since
            // count is always >= 1, don't generate dead code if min <= 1.
            //
            if (min > 1)
            {
              os << "if (count < " << min << "UL)" << endl
                 << "this->_expected_element (" << endl;

              particle_name_.dispatch (p);

              os << "," << endl
                 << "ns, n);";
            }


            os << "state = ~0UL;"
               << "}";
          }

          os << "break;"
             << "}"; // case
        }

        virtual Void
        traverse (SemanticGraph::Compositor& c)
        {
          using SemanticGraph::Compositor;

          UnsignedLong max (c.max ());
          UnsignedLong min (c.context ().get<UnsignedLong> ("effective-min"));
          UnsignedLong n (c.context ().get<UnsignedLong> ("comp-number"));
          UnsignedLong state (c.context ().get<UnsignedLong> ("state"));

          String func (c.is_a<SemanticGraph::Choice> () ?
                       "choice_" : "sequence_");

          os << "case " << state << "UL:" << endl
             << "{"
             << "unsigned long s (~0UL);"
             << endl;

          Boolean first (true);

          for (Compositor::ContainsIterator ci (c.contains_begin ());
               ci != c.contains_end (); ++ci)
          {
            SemanticGraph::Particle& p (ci->particle ());

            if (p.is_a<Compositor> () && !c.context ().count ("comp-number"))
              continue; // Empty compositor.

            if (!p.context ().count ("prefix"))
              break;

            UnsignedLong state (p.context ().get<UnsignedLong> ("state"));

            if (first)
              first = false;
            else
              os << "else ";

            os << "if (";

            particle_test_.dispatch (p);

            os << ")" << endl
               << "s = " << state << "UL;";
          }

          // This compositor.
          //
          os << endl
             << "if (s != ~0UL)"
             << "{"
             << "assert (start);"; // End is handled by the sub-machine.

          switch (max)
          {
          case 0:
            {
              os << "count++;";
              break;
            }
          case 1:
            {
              // We do not need to increment count because min <= max and
              // we do not generate min check for min <= 1 (see below).
              //
              os << "state = ~0UL;";
              break;
            }
          default:
            {
              os << "if (++count == " << max << "UL)" << endl
                 << "state = ~0UL;";
            }
          };

          // Delegate to the sub-machine.
          //

          os << endl
             << "v_state_& vs = *static_cast< v_state_* > (" <<
            "this->v_state_stack_.top ());"
             << "v_state_descr_& vd = vs.data[vs.size++];" // push
             << endl
             << "vd.func = &" << ename (type_) << "::" << func << n << ";"
             << "vd.state = s;"
             << "vd.count = 0;"
             << endl
             << "this->" << func << n << " (vd.state, vd.count, ns, n, t, true);"
             << "}";


          // Not this compositor. We've elready moved to the final state
          // if max == 1.
          //
          if (max != 1)
          {
            os << "else"
               << "{"
               << "assert (start);"; // Assuming well-formed XML

            // Check if min cardinality requirements have been met. Since
            // count is always >= 1, don't generate dead code if min <= 1.
            //
            if (min > 1)
            {
              os << "if (count < " << min << "UL)" << endl
                 << "this->_expected_element (" << endl;

              particle_name_.dispatch (c);

              os << "," << endl
                 << "ns, n);";
            }

            os << "state = ~0UL;"
               << "}";
          }

          os << "break;"
             << "}"; // case
        }

      private:
        ParticleName particle_name_;
      };


      //
      //
      struct ParticleInSequence: Traversal::Particle,
                                 Traversal::Compositor,
                                 ParticleInCompositor
      {
        ParticleInSequence (Context& c,
                            UnsignedLong state,
                            UnsignedLong next_state,
                            SemanticGraph::Complex& type)
            : ParticleInCompositor (c, type),
              state_ (state), particle_name_ (c)
        {
          // next_state == 0 indicates the terminal state (~0UL).
          //
          if (next_state != 0)
          {
            std::wostringstream ostr;
            ostr << next_state;
            next_state_ = ostr.str ();
          }
          else
            next_state_ = L"~0";
        }

        virtual Void
        traverse (SemanticGraph::Particle& p)
        {
          UnsignedLong min (p.min ()), max (p.max ());

          os << "case " << state_ << "UL:" << endl
             << "{"
             << "if (";

          particle_test_.dispatch (p);

          os << ")"
             << "{";

          // This element.
          //

          os << "if (start)"
             << "{";

          pre_post_calls (p);

          switch (max)
          {
          case 0:
            {
              os << "count++;";
              break;
            }
          case 1:
            {
              os << "count = 0;"
                 << "state = " << next_state_ << "UL;";
              break;
            }
          default:
            {
              os << "if (++count == " << max << "UL)"
                 << "{"
                 << "count = 0;"
                 << "state = " << next_state_ << "UL;"
                 << "}";
            }
          };

          os << "}" // start
             << "break;"
             << "}";

          // Not this element.
          //

          os << "else"
             << "{"
             << "assert (start);"; // Assuming well-formed XML.

          // Check if min cardinality requirements have been met. Since
          // count is always >= 0, don't generate dead code if min == 0.
          //
          if (min != 0)
          {
            os << "if (count < " << min << "UL)" << endl
               << "this->_expected_element (" << endl;

            particle_name_.dispatch (p);

            os << "," << endl
               << "ns, n);";
          }

          os << "count = 0;"
             << "state = " << next_state_ << "UL;"
             << "// Fall through." << endl
             << "}"  // else
             << "}"; // case
        }

        virtual Void
        traverse (SemanticGraph::Compositor& c)
        {
          using SemanticGraph::Compositor;

          UnsignedLong max (c.max ());
          UnsignedLong min (c.context ().get<UnsignedLong> ("effective-min"));
          UnsignedLong n (c.context ().get<UnsignedLong> ("comp-number"));

          String func (c.is_a<SemanticGraph::Choice> () ?
                       "choice_" : "sequence_");

          os << "case " << state_ << "UL:" << endl
             << "{"
             << "unsigned long s (~0UL);"
             << endl;

          Boolean first (true);

          for (Compositor::ContainsIterator ci (c.contains_begin ());
               ci != c.contains_end (); ++ci)
          {
            SemanticGraph::Particle& p (ci->particle ());

            if (p.is_a<Compositor> () && !c.context ().count ("comp-number"))
              continue; // Empty compositor.

            if (!p.context ().count ("prefix"))
              break;

            UnsignedLong state (p.context ().get<UnsignedLong> ("state"));

            if (first)
              first = false;
            else
              os << "else ";

            os << "if (";

            particle_test_.dispatch (p);

            os << ")" << endl
               << "s = " << state << "UL;";
          }

          // This element.
          //

          os << endl
             << "if (s != ~0UL)"
             << "{"
             << "assert (start);"; // End is handled by the sub-machine.

          switch (max)
          {
          case 0:
            {
              os << "count++;"
                 << endl;
              break;
            }
          case 1:
            {
              os << "count = 0;"
                 << "state = " << next_state_ << "UL;"
                 << endl;
              break;
            }
          default:
            {
              os << "if (++count == " << max << "UL)"
                 << "{"
                 << "count = 0;"
                 << "state = " << next_state_ << "UL;"
                 << "}";
            }
          };

          // Delegate to the sub-machine.
          //

          os << "v_state_& vs = *static_cast< v_state_* > (" <<
            "this->v_state_stack_.top ());"
             << "v_state_descr_& vd = vs.data[vs.size++];" // push
             << endl
             << "vd.func = &" << ename (type_) << "::" << func << n << ";"
             << "vd.state = s;"
             << "vd.count = 0;"
             << endl
             << "this->" << func << n << " (vd.state, vd.count, ns, n, t, true);"
             << "break;"
             << "}";

          // Not this compositor.
          //

          os << "else"
             << "{"
             << "assert (start);"; // Assuming well-formed XML

          // Check if min cardinality requirements have been met. Since
          // count is always >= 0, don't generate dead code if min == 0.
          //
          if (min != 0)
          {
            os << "if (count < " << min << "UL)" << endl
               << "this->_expected_element (" << endl;

            particle_name_.dispatch (c);

            os << "," << endl
               << "ns, n);";
          }

          os << "count = 0;"
             << "state = " << next_state_ << "UL;"
             << "// Fall through." << endl
             << "}"  // else
             << "}"; // case
        }

      private:
        UnsignedLong state_;
        String next_state_;

        ParticleName particle_name_;
      };


      //
      //
      struct ParticleFunction: Traversal::All,
                               Traversal::Choice,
                               Traversal::Sequence,
                               Context
      {
        ParticleFunction (Context& c, SemanticGraph::Complex& type)
            : Context (c), type_ (type)
        {
          *this >> contains_particle_ >> *this;
        }


        virtual Void
        traverse (SemanticGraph::All& a)
        {
          if (!a.context().count ("comp-number")) // Empty compositor.
            return;

          using SemanticGraph::Element;
          using SemanticGraph::Compositor;


          os << "void " << ename (type_) << "::" << endl
             << "all_0 (unsigned long& state," << endl
             << "unsigned char* count," << endl
             << "const " << string_type << "& ns," << endl
             << "const " << string_type << "& n," << endl
             << "const " << string_type << "* t," << endl
             << "bool start)"
             << "{"
             << "XSD_UNUSED (t);"
             << endl;

          for (Compositor::ContainsIterator ci (a.contains_begin ()),
                 ce (a.contains_end ()); ci != ce; ++ci)
          {
            ParticleInAll t (*this, type_);
            t.dispatch (ci->particle ());
          }

          // Handle the flush.
          //
          os << "else if (n.empty () && ns.empty ())"
             << "{";

          for (Compositor::ContainsIterator ci (a.contains_begin ()),
                 ce (a.contains_end ()); ci != ce; ++ci)
          {
            if (ci->min () == 0)
              continue;

            Element& e (dynamic_cast<Element&> (ci->particle ()));
            String ns (e.qualified_p () ? e.namespace_ ().name () : String ());
            UnsignedLong state (e.context ().get<UnsignedLong> ("state"));

            os << "if (count[" << state << "UL] == 0)" << endl
               << "this->_expected_element (" << endl
               << strlit (ns) << ", " <<
              strlit (e.name ()) << ");"
               << endl;
          }

          os << "state = ~0UL;"
             << "}"
             << "else" << endl
             << "state = ~0UL;"
             << "}";
        }

        virtual Void
        traverse (SemanticGraph::Choice& c)
        {
          if (!c.context().count ("comp-number")) // Empty compositor.
            return;

          using SemanticGraph::Compositor;

          UnsignedLong n (c.context ().get<UnsignedLong> ("comp-number"));

          os << "void " << ename (type_) << "::" << endl
             << "choice_" << n << " (unsigned long& state," << endl
             << "unsigned long& count," << endl
             << "const " << string_type << "& ns," << endl
             << "const " << string_type << "& n," << endl
             << "const " << string_type << "* t," << endl
             << "bool start)"
             << "{"
             << "XSD_UNUSED (count);"
             << "XSD_UNUSED (ns);"
             << "XSD_UNUSED (n);"
             << "XSD_UNUSED (t);"
             << endl
             << "switch (state)"
             << "{";

          for (Compositor::ContainsIterator ci (c.contains_begin ()),
                 ce (c.contains_end ()); ci != ce; ++ci)
          {
            SemanticGraph::Particle& p (ci->particle ());

            if (p.is_a<Compositor> () && !p.context().count ("comp-number"))
              continue; // Empty compositor.

            ParticleInChoice t (*this, type_);
            t.dispatch (p);
          }

          os << "}" // switch
             << "}";

          // Generate nested compositor functions.
          //
          Traversal::Choice::traverse (c);
        }

        virtual Void
        traverse (SemanticGraph::Sequence& s)
        {
          if (!s.context().count ("comp-number")) // Empty compositor.
            return;

          using SemanticGraph::Compositor;

          UnsignedLong n (s.context ().get<UnsignedLong> ("comp-number"));

          os << "void " << ename (type_) << "::" << endl
             << "sequence_" << n << " (unsigned long& state," << endl
             << "unsigned long& count," << endl
             << "const " << string_type << "& ns," << endl
             << "const " << string_type << "& n," << endl
             << "const " << string_type << "* t," << endl
             << "bool start)"
             << "{"
             << "XSD_UNUSED (t);"
             << endl
             << "switch (state)"
             << "{";

          UnsignedLong state (0);

          for (Compositor::ContainsIterator ci (s.contains_begin ()),
                 ce (s.contains_end ()); ci != ce;)
          {
            SemanticGraph::Particle& p (ci->particle ());

            if (p.is_a<Compositor> () && !p.context().count ("comp-number"))
            {
              // Empty compositor.
              //
              ++ci;
              continue;
            }

            // Find the next state.
            //
            do
              ++ci;
            while (ci != ce &&
                   ci->particle ().is_a<Compositor> () &&
                   !ci->particle ().context().count ("comp-number"));

            UnsignedLong next (ci == ce ? 0 : state + 1);

            ParticleInSequence t (*this, state++, next, type_);
            t.dispatch (p);
          }

          os << "case ~0UL:" << endl
             << "break;"
             << "}" // switch
             << "}";

          // Generate nested compositor functions.
          //
          Traversal::Sequence::traverse (s);
        }

      private:
        SemanticGraph::Complex& type_;
        Traversal::ContainsParticle contains_particle_;
      };


      //
      //
      struct CompositorPre: Traversal::All,
                            Traversal::Compositor,
                            Context
      {
        CompositorPre (Context& c, SemanticGraph::Complex& type)
            : Context (c), type_ (type)
        {
        }

        virtual Void
        traverse (SemanticGraph::All& a)
        {
          // Clear the counts and push the initial state.
          //
          os << "v_all_count_.push ();"
             << endl;

          SemanticGraph::Compositor& c (a);
          traverse (c);
        }

        virtual Void
        traverse (SemanticGraph::Compositor&) // Choice and sequence.
        {
          os << "v_state_& vs = *static_cast< v_state_* > (" <<
            "this->v_state_stack_.top ());"
             << "v_state_descr_& vd = vs.data[vs.size++];" // push
             << endl
             << "vd.func = 0;"
             << "vd.state = 0;"
             << "vd.count = 0;";
        }

      private:
        SemanticGraph::Complex& type_;
      };


      //
      //
      struct CompositorStartElement: Traversal::All,
                                     Traversal::Compositor,
                                     Context
      {
        CompositorStartElement (Context& c, SemanticGraph::Complex& type)
            : Context (c), type_ (type),
              particle_test_ (c), particle_name_ (c)
        {
        }

        virtual Void
        traverse (SemanticGraph::All&)
        {
          // The 'all' state machine reaches the final state only
          // on an unknown element, in which case we won't get here
          // again (it would be a validation error). Note that 'all'
          // compositor cannot contain nested compositors so we don't
          // need to re-set vd.
          //
          os << "all_0 (vd->state, v_all_count_.top (), ns, n, t, true);"
             << endl
             << "if (vd->state != ~0UL)" << endl
             << "vd->count++;"
             << "else" << endl
             << "return false;" // Let our parent handle this.
             << endl;
        }

        virtual Void
        traverse (SemanticGraph::Compositor& c) // Choice and sequence.
        {
          using SemanticGraph::Compositor;

          UnsignedLong max (c.max ());
          UnsignedLong min (c.context ().get<UnsignedLong> ("effective-min"));
          UnsignedLong n (c.context ().get<UnsignedLong> ("comp-number"));

          String func (c.is_a<SemanticGraph::Choice> () ?
                       "choice_" : "sequence_");

          // Invoke the current state machine. If it reaches its
          // terminal state, pop it and invoke the next one until
          // we reach the top, which requires special handling.
          //
          os << "while (vd->func != 0)"
             << "{"
             << "(this->*vd->func) (vd->state, vd->count, ns, n, t, true);"
             << endl
             << "vd = vs.data + (vs.size - 1);" // re-acquire
             << endl
             << "if (vd->state == ~0UL)" << endl
             << "vd = vs.data + (--vs.size - 1);" // pop
             << "else" << endl
             << "break;"
             << "}";


          // Check if we got to the top. This code is pretty much the
          // same as the one found in ParticleInSequence.
          //
          os << "if (vd->func == 0)"
             << "{"
             << "if (vd->state != ~0UL)"
             << "{"
             << "unsigned long s = ~0UL;"
             << endl;

          Boolean first (true);

          // Note that we don't need to worry about the compositor
          // being empty - this case is handled by our caller.
          //
          for (Compositor::ContainsIterator ci (c.contains_begin ());
               ci != c.contains_end (); ++ci)
          {
            SemanticGraph::Particle& p (ci->particle ());

            if (p.is_a<Compositor> () && !c.context ().count ("comp-number"))
              continue; // Empty compositor.

            if (!p.context ().count ("prefix"))
              break;

            UnsignedLong state (p.context ().get<UnsignedLong> ("state"));

            if (first)
              first = false;
            else
              os << "else ";

            os << "if (";

            particle_test_.dispatch (p);

            os << ")" << endl
               << "s = " << state << "UL;";
          }

          os << endl
             << "if (s != ~0UL)"
             << "{";

          // This element is a prefix of the root compositor.
          //

          switch (max)
          {
          case 0:
            {
              os << "vd->count++;";
              break;
            }
          case 1:
            {
              os << "vd->count++;"
                 << "vd->state = ~0UL;";
              break;
            }
          default:
            {
              os << "if (++vd->count == " << max << "UL)" << endl
                 << "vd->state = ~0UL;";
            }
          };

          // Delegate to the sub-machine.
          //

          os << endl
             << "vd = vs.data + vs.size++;" // push
             << "vd->func = &" << ename (type_) << "::" << func << n << ";"
             << "vd->state = s;"
             << "vd->count = 0;"
             << endl
             << "this->" << func << n << " (vd->state, vd->count, ns, n, t, true);"
             << "}";

          // This element is not our prefix.
          //

          os << "else"
             << "{";

          // Check if min cardinality requirements have been met. Since
          // count is always >= 0, don't generate dead code if min == 0.
          //
          if (min != 0)
          {
            os << "if (vd->count < " << min << "UL)" << endl
               << "this->_expected_element (" << endl;

            particle_name_.dispatch (c);

            os << "," << endl
               << "ns, n);";
          }

          // Return false to indicate that we are not handling this element.
          //
          os << "return false;"
             << "}"
             << "}" // if (state != ~0)
             << "else" << endl
             << "return false;"
             << "}"; // if (function == 0)
        }

      private:
        SemanticGraph::Complex& type_;
        ParticleTest particle_test_;
        ParticleName particle_name_;
      };


      //
      //
      struct CompositorEndElement: Traversal::All,
                                   Traversal::Compositor,
                                   Context
      {
        CompositorEndElement (Context& c, SemanticGraph::Complex& type)
            : Context (c), type_ (type)
        {
        }

        virtual Void
        traverse (SemanticGraph::All&)
        {
          os << "all_0 (vd.state, v_all_count_.top (), " <<
            "ns, n, 0, false);"
             << endl;
        }

        virtual Void
        traverse (SemanticGraph::Compositor&) // Choice and sequence.
        {
          os << "assert (vd.func != 0);"
             << "(this->*vd.func) (vd.state, vd.count, ns, n, 0, false);"
             << endl
             << "if (vd.state == ~0UL)" << endl
             << "vs.size--;" // pop
             << endl;
        }

      private:
        SemanticGraph::Complex& type_;
      };


      //
      //
      struct CompositorPost: Traversal::All,
                             Traversal::Compositor,
                             Context
      {
        CompositorPost (Context& c, SemanticGraph::Complex& type)
            : Context (c), type_ (type), particle_name_ (c)
        {
        }

        virtual Void
        traverse (SemanticGraph::All& a)
        {
          using SemanticGraph::Element;

          os << "v_state_& vs = *static_cast< v_state_* > (" <<
            "this->v_state_stack_.top ());"
             << "v_state_descr_& vd = vs.data[vs.size - 1];"
             << endl;

          // Flush the state machine with the empty element name. This
          // allows us to detect missing content.
          //
          os << "if (vd.count != 0)"
             << "{"
             << string_type << " empty;"
             << "all_0 (vd.state, v_all_count_.top (), empty, empty, 0, true);"
             << "}";

          if (a.context ().get<UnsignedLong> ("effective-min") != 0)
          {
            os << "else" << endl
               << "this->_expected_element (" << endl;

            particle_name_.dispatch (a);

            os << ");";
          }

          os << endl
             << "vs.size--;" // pop
             << "v_all_count_.pop ();";
        }

        virtual Void
        traverse (SemanticGraph::Compositor& c) // Choice and sequence.
        {
          UnsignedLong min (c.context ().get<UnsignedLong> ("effective-min"));

          os << "v_state_& vs = *static_cast< v_state_* > (" <<
            "this->v_state_stack_.top ());"
             << "v_state_descr_* vd = vs.data + (vs.size - 1);"
             << endl;


          // Flush unfinished state machines with the empty element name.
          // This allows us to detect missing content. Note that I am
          // not re-setting vd since no new compositors are pushed on
          // flush.
          //
          os << string_type << " empty;"
             << "while (vd->func != 0)"
             << "{"
             << "(this->*vd->func) (vd->state, vd->count, empty, empty, 0, true);"
             << "assert (vd->state == ~0UL);"
             << "vd = vs.data + (--vs.size - 1);" // pop
             << "}";

          // Check if min cardinality requirements have been met. Since
          // count is always >= 0, don't generate dead code if min == 0.
          //
          if (min != 0)
          {
            os << "if (vd->count < " << min << "UL)" << endl
               << "this->_expected_element (" << endl;

            particle_name_.dispatch (c);

            os << ");";
          }
        }

      private:
        SemanticGraph::Complex& type_;
        ParticleName particle_name_;
      };


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
          // Nothing to generate if we don't have any elements and wildcards.
          //
          if (!has<Traversal::Element> (c) &&
              !has_particle<Traversal::Any> (c))
            return;

          using SemanticGraph::Compositor;

          String const& name (ename (c));
          Compositor& comp (c.contains_compositor ().compositor ());

          // Don't use restriction_p here since we don't want special
          // treatment of anyType.
          //
          Boolean restriction (
            c.inherits_p () &&
            c.inherits ().is_a<SemanticGraph::Restricts> ());

          os <<"// Element validation and dispatch functions for " <<
            name << "." << endl
             <<"//" << endl;

          // _start_element_impl
          //

          os << "bool " << name << "::" << endl
             << "_start_element_impl (const " << string_type << "& ns," << endl
             << "const " << string_type << "& n," << endl
             << "const " << string_type << "* t)"
             << "{"
             << "XSD_UNUSED (t);"
             << endl;

          os << "v_state_& vs = *static_cast< v_state_* > (" <<
            "this->v_state_stack_.top ());"
             << "v_state_descr_* vd = vs.data + (vs.size - 1);"
             << endl;

          //@@ OPT: I don't really need to call parser_base since it always
          // returns false.
          //
          // In case of an inheritance-by-extension, call our base first.
          // We don't need to generate this code for the 'all' compositor
          // because it can only inherit from the empty content model.
          // States of the root machine for sequence and choice:
          //
          //  0 - calling base
          //  1 - base returned false
          // ~0 - terminal state
          //
          if (!restriction && !comp.is_a<SemanticGraph::All> ())
          {
            os << "if (vd->func == 0 && vd->state == 0)"
               << "{"
               << "if (this->";

            if (c.inherits_p ())
              os << fq_name (c.inherits ().base ());
            else
              os << complex_base;


            os << "::_start_element_impl (ns, n, t))" << endl
               << "return true;"
               << "else" << endl
               << "vd->state = 1;"
               << "}";
          }

          {
            CompositorStartElement t (*this, c);
            t.dispatch (comp);
          }

          os << "return true;"
             << "}";


          // _end_element_impl
          //

          os << "bool " << name << "::" << endl
             << "_end_element_impl (const " << string_type << "& ns," << endl
             << "const " << string_type << "& n)"
             << "{";

          os << "v_state_& vs = *static_cast< v_state_* > (" <<
            "this->v_state_stack_.top ());"
             << "v_state_descr_& vd = vs.data[vs.size - 1];"
             << endl;

          //@@ OPT: I don't really need to call parser_base since it always
          // returns false.
          //
          // In case of an inheritance-by-extension, call our base first.
          // We don't need to generate this code for the 'all' compositor
          // because it can only inherit from the empty content model.
          //
          if (!restriction && !comp.is_a<SemanticGraph::All> ())
          {
            os << "if (vd.func == 0 && vd.state == 0)"
               << "{"
               << "if (!";

            if (c.inherits_p ())
              os << fq_name (c.inherits ().base ());
            else
              os << complex_base;

            os << "::_end_element_impl (ns, n))" << endl
               << "assert (false);" // Start and end should match.
               << "return true;"
               << "}";
          }

          {
            CompositorEndElement t (*this, c);
            t.dispatch (comp);
          }

          os << "return true;"
             << "}";


          // _pre_e_validate
          //
          os << "void " << name << "::" << endl
             << "_pre_e_validate ()"
             << "{"
             << "this->v_state_stack_.push ();"
             << "static_cast< v_state_* > (this->v_state_stack_.top ())->" <<
            "size = 0;"
             << endl;

          {
            CompositorPre t (*this, c);
            t.dispatch (comp);
          }

          // In case of an inheritance-by-extension, call our base
	  // _pre_e_validate. We don't need to generate this code for the
	  // 'all' compositor because it can only inherit from the empty
	  // content model.
          //
          if (!restriction && !comp.is_a<SemanticGraph::All> ())
          {
            // We don't need to call parser_base's implementation
            // since it does nothing.
            //
            if (c.inherits_p ())
            {
              os << endl
                 << fq_name (c.inherits ().base ()) << "::_pre_e_validate ();";
            }
          }

          os << "}";


          // _post_e_validate
          //
          os << "void " << name << "::" << endl
             << "_post_e_validate ()"
             << "{";

          // In case of an inheritance-by-extension, call our base
	  // _post_e_validate. We don't need to generate this code for
	  // the 'all' compositor because it can only inherit from
	  // the empty content model.
          //
          if (!restriction && !comp.is_a<SemanticGraph::All> ())
          {
            // We don't need to call parser_base's implementation
            // since it does nothing.
            //
            if (c.inherits_p ())
            {
              os << fq_name (c.inherits ().base ()) << "::_post_e_validate ();"
                 << endl;
            }
          }

          {
            CompositorPost t (*this, c);
            t.dispatch (c.contains_compositor ().compositor ());
          }

          os << endl
             << "this->v_state_stack_.pop ();"
             << "}";

          //
          //
          ParticleFunction t (*this, c);
          t.dispatch (c.contains_compositor ().compositor ());
        }
      };
    }

    Void
    generate_element_validation_source (Context& ctx)
    {
      ctx.os << "#include <cassert>" << endl
             << endl;

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
