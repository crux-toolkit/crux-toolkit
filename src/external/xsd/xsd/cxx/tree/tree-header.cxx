// file      : xsd/cxx/tree/tree-header.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

#include <cxx/tree/tree-header.hxx>
#include <cxx/tree/default-value.hxx>
#include <cxx/tree/fundamental-header.hxx>

namespace CXX
{
  namespace Tree
  {
    namespace
    {
      typedef Containers::Vector<NarrowString> Streams;

      // List mapping.
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
          String name (ename (l));

          // If renamed name is empty then we do not need to generate
          // anything for this type.
          //
          if (renamed_type (l, name) && !name)
            return;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief List class corresponding to the %" <<
              comment (l.name ()) << endl
               << " * schema type." << endl
               << " *" << endl
               << " * This class has an interface of a standard C++ " <<
              "sequence (e.g.," << endl
               << " * std::vector)." << endl;

            if (l.annotated_p ())
            {
              os << " *" << endl;
              write_annotation (l.annotation ());
            }

            os << " */" << endl;
          }

          SemanticGraph::Type& item_type (l.argumented ().type ());
          String item_name (item_type_name (item_type));
          String base_type (L"::xsd::cxx::tree::list< " + item_name + L", " +
                            char_type);

          if (item_type.is_a<SemanticGraph::Fundamental::Double> ())
            base_type += L", ::xsd::cxx::tree::schema_type::double_";
          else if (item_type.is_a<SemanticGraph::Fundamental::Decimal> ())
            base_type += L", ::xsd::cxx::tree::schema_type::decimal";

          base_type += L" >";

          os << "class " << type_exp << name <<
            ": public " << any_simple_type << "," << endl
             << "  public " << base_type
             << "{"
             << "public:" << endl;

          // c-tor ()
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Default constructor." << endl
               << " *" << endl
               << " * Creates an empty list." << endl
               << " */" << endl;
          }
          os << name << " ();"
             << endl;

          // c-tor (size_type, const X& x)
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Create a list with copies of the specified " <<
               "element." << endl
               << " *" << endl
               << " * @param n A number of elements to copy." << endl
               << " * @param x An element to copy." << endl
               << " *" << endl
               << " * This constructor creates a list with @a n copies " <<
              "of @a x." << endl
               << " */" << endl;
          }

          String size_type (name != L"size_type"
                            ? String (L"size_type")
                            : base_type + L"::size_type");

          os << name << " (" << size_type << " n, const " << item_name <<
            "& x);"
             << endl;

          // c-tor (const I& begin, const I& end)
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Create a list from an iterator range." << endl
               << " *" << endl
               << " * @param begin An iterator pointing to the first " <<
              "element." << endl
               << " * @param end An iterator pointing to the one past " <<
              "the last element." << endl
               << " *" << endl
               << " * This constructor creates a list consisting of " <<
              "copies of the" << endl
               << " * elements in the range [begin,end)." << endl
               << " */" << endl;
          }

          String iter_type (unclash (name, "I"));

          os << "template < typename " << iter_type << " >" << endl
             << name << " (const " << iter_type << "& begin, const " <<
            iter_type << "& end)" << endl
             << ": " << base_type << " (begin, end, this)"
             << "{"
             << "}";

          // c-tor (istream&)
          //
          Streams const& st (options.value<CLI::generate_extraction> ());
          for (Streams::ConstIterator i (st.begin ()); i != st.end (); ++i)
          {
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a data " <<
                "representation" << endl
                 << " * stream." << endl
                 << " *" << endl
                 << " * @param s A stream to extract the data from." << endl
                 << " * @param f Flags to create the new instance with." << endl
                 << " * @param c A pointer to the object that will " <<
                "contain the new" << endl
                 << " * instance." << endl
                 << " */" << endl;
            }

            os << name << " (" << istream_type << "< " << i->c_str () <<
              " >& s," << endl
               << flags_type << " f = 0," << endl
               << container << "* c = 0);"
               << endl;
          }

          if (!options.value<CLI::suppress_parsing> ())
          {
            // c-tor (xercesc::DOMElement)
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a DOM element." << endl
                 << " *" << endl
                 << " * @param e A DOM element to extract the data from." << endl
                 << " * @param f Flags to create the new instance with." << endl
                 << " * @param c A pointer to the object that will " <<
                "contain the new" << endl
                 << " * instance." << endl
                 << " */" << endl;
            }

            os << name << " (const " << xerces_ns << "::DOMElement& e," << endl
               << flags_type << " f = 0," << endl
               << container << "* c = 0);"
               << endl;

            // c-tor (xercesc::DOMAttr)
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a DOM attribute." << endl
                 << " *" << endl
                 << " * @param a A DOM attribute to extract the data from." << endl
                 << " * @param f Flags to create the new instance with." << endl
                 << " * @param c A pointer to the object that will " <<
                "contain the new" << endl
                 << " * instance." << endl
                 << " */" << endl;
            }

            os << name << " (const " << xerces_ns << "::DOMAttr& a," << endl
               << flags_type << " f = 0," << endl
               << container << "* c = 0);"
               << endl;

            // c-tor (std::basic_string const&, xercesc::DOMElement)
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a string fragment." << endl
                 << " *" << endl
                 << " * @param s A string fragment to extract the data from." << endl
                 << " * @param e A pointer to DOM element containing the " <<
                "string fragment." << endl
                 << " * @param f Flags to create the new instance with." << endl
                 << " * @param c A pointer to the object that will " <<
                "contain the new" << endl
                 << " * instance." << endl
                 << " */" << endl;
            }

            os << name << " (const " << string_type << "& s," << endl
               << "const " << xerces_ns << "::DOMElement* e," << endl
               << flags_type << " f = 0," << endl
               << container << "* c = 0);"
               << endl;
          }

          // copy c-tor ()
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Copy constructor." << endl
               << " *" << endl
               << " * @param x An instance to make a copy of." << endl
               << " * @param f Flags to create the copy with." << endl
               << " * @param c A pointer to the object that will contain " <<
              "the copy." << endl
               << " *" << endl
               << " * For polymorphic object models use the @c _clone " <<
              "function instead." << endl
               << " */" << endl;
          }

          os << name << " (const " << name << "& x," << endl
             << flags_type << " f = 0," << endl
             << container << "* c = 0);"
             << endl;

          // clone
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Copy the instance polymorphically." << endl
               << " *" << endl
               << " * @param f Flags to create the copy with." << endl
               << " * @param c A pointer to the object that will contain " <<
              "the copy." << endl
               << " * @return A pointer to the dynamically allocated copy." << endl
               << " *" << endl
               << " * This function ensures that the dynamic type of the " <<
              "instance is" << endl
               << " * used for copying and should be used for polymorphic " <<
              "object" << endl
               << " * models instead of the copy constructor." << endl
               << " */" << endl;
          }

          os << "virtual " << name << "*" << endl
             << "_clone (" << flags_type << " f = 0," << endl
             << container << "* c = 0) const;"
             << endl;

          // d-tor
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Destructor." << endl
               << " */" << endl;
          }

          os << "virtual " << endl
             << "~" << name << " ();";

          os << "};";
        }

      private:
        String
        item_type_name (SemanticGraph::Type& t)
        {
          std::wostringstream o;

          MemberTypeName type (*this, o);
          type.dispatch (t);

          return o.str ();
        }
      };


      // Union mapping.
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
          String name (ename (u));

          // If renamed name is empty then we do not need to generate
          // anything for this type.
          //
          if (renamed_type (u, name) && !name)
            return;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Union class corresponding to the %" <<
              comment (u.name ()) << endl
               << " * schema type." << endl
               << " *" << endl
               << " * The mapping represents unions as strings." << endl;

            if (u.annotated_p ())
            {
              os << " *" << endl;
              write_annotation (u.annotation ());
            }

            os << " */" << endl;
          }

          os << "class " << type_exp << name <<
	    ": public " <<  xs_string_type
             << "{"
             << "public:" << endl
             << endl;

          if (options.value<CLI::generate_default_ctor> ())
          {
            // c-tor ()
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Default constructor." << endl
                 << " *" << endl
                 << " * Note that this constructor may leave the " <<
                "instance in an" << endl
                 << " * invalid state." << endl
                 << " */" << endl;
            }

            os << name << " ();"
               << endl;
          }

          // c-tor (const char*)
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Create an instance from a C string." << endl
               << " *" << endl
               << " * @param v A string value." << endl
               << " */" << endl;
          }
          os << name << " (const " << char_type << "* v);"
             << endl;

          // c-tor (string const&)
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Create an instance from a string." << endl
               << " *" << endl
               << " * @param v A string value." << endl
               << " */" << endl;
          }
          os << name << " (const " << string_type << "& v);"
             << endl;

          // c-tor (istream&)
          //
          Streams const& st (options.value<CLI::generate_extraction> ());
          for (Streams::ConstIterator i (st.begin ()); i != st.end (); ++i)
          {
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a data " <<
                "representation" << endl
                 << " * stream." << endl
                 << " *" << endl
                 << " * @param s A stream to extract the data from." << endl
                 << " * @param f Flags to create the new instance with." << endl
                 << " * @param c A pointer to the object that will " <<
                "contain the new" << endl
                 << " * instance." << endl
                 << " */" << endl;
            }

            os << name << " (" << istream_type << "< " << i->c_str () <<
              " >& s," << endl
               << flags_type << " f = 0," << endl
               << container << "* c = 0);"
               << endl;
          }

          if (!options.value<CLI::suppress_parsing> ())
          {
            // c-tor (xercesc::DOMElement)
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a DOM element." << endl
                 << " *" << endl
                 << " * @param e A DOM element to extract the data from." << endl
                 << " * @param f Flags to create the new instance with." << endl
                 << " * @param c A pointer to the object that will " <<
                "contain the new" << endl
                 << " * instance." << endl
                 << " */" << endl;
            }

            os << name << " (const " << xerces_ns << "::DOMElement& e," << endl
               << flags_type << " f = 0," << endl
               << container << "* c = 0);"
               << endl;

            // c-tor (xercesc::DOMAttr)
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a DOM attribute." << endl
                 << " *" << endl
                 << " * @param a A DOM attribute to extract the data from." << endl
                 << " * @param f Flags to create the new instance with." << endl
                 << " * @param c A pointer to the object that will " <<
                "contain the new" << endl
                 << " * instance." << endl
                 << " */" << endl;
            }

            os << name << " (const " << xerces_ns << "::DOMAttr& a," << endl
               << flags_type << " f = 0," << endl
               << container << "* c = 0);"
               << endl;

            // c-tor (std::basic_string const&, xercesc::DOMElement)
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a string fragment." << endl
                 << " *" << endl
                 << " * @param s A string fragment to extract the data from." << endl
                 << " * @param e A pointer to DOM element containing the " <<
                "string fragment." << endl
                 << " * @param f Flags to create the new instance with." << endl
                 << " * @param c A pointer to the object that will " <<
                "contain the new" << endl
                 << " * instance." << endl
                 << " */" << endl;
            }

            os << name << " (const " << string_type << "& s," << endl
               << "const " << xerces_ns << "::DOMElement* e," << endl
               << flags_type << " f = 0," << endl
               << container << "* c = 0);"
               << endl;
          }

          // copy c-tor ()
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Copy constructor." << endl
               << " *" << endl
               << " * @param x An instance to make a copy of." << endl
               << " * @param f Flags to create the copy with." << endl
               << " * @param c A pointer to the object that will contain " <<
              "the copy." << endl
               << " *" << endl
               << " * For polymorphic object models use the @c _clone " <<
              "function instead." << endl
               << " */" << endl;
          }

          os << name << " (const " << name << "& x," << endl
             << flags_type << " f = 0," << endl
             << container << "* c = 0);"
             << endl;

          // clone
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Copy the instance polymorphically." << endl
               << " *" << endl
               << " * @param f Flags to create the copy with." << endl
               << " * @param c A pointer to the object that will contain " <<
              "the copy." << endl
               << " * @return A pointer to the dynamically allocated copy." << endl
               << " *" << endl
               << " * This function ensures that the dynamic type of the " <<
              "instance is" << endl
               << " * used for copying and should be used for polymorphic " <<
              "object" << endl
               << " * models instead of the copy constructor." << endl
               << " */" << endl;
          }

          os << "virtual " << name << "*" << endl
             << "_clone (" << flags_type << " f = 0," << endl
             << container << "* c = 0) const;"
             << endl;

          os << "};";
        }
      };

      // Enum mapping.
      //
      struct Enumerator: Traversal::Enumerator, Context
      {
        Enumerator (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& e)
        {
          if (doxygen && e.annotated_p ())
          {
            os << "/**" << endl;
            write_annotation (e.annotation ());
            os << " */" << endl;
          }

          os << ename (e);
        }
      };

      struct Enumeration: Traversal::Enumeration, Context
      {
        Enumeration (Context& c)
            : Context (c),
              base_ (c),
              member_ (c),
              enumerator_ (c)
        {
          inherits_base_ >> base_;
          inherits_member_ >> member_;

          names_ >> enumerator_;
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
          SemanticGraph::Enumeration* base_enum (0);

          if (string_based)
          {
            IsEnumBasedType t (base_enum);
            t.dispatch (e);

            if (base_enum != 0)
              enum_based = true;
          }

          String value;
          if (string_based)
            value = evalue (e);

          // Get to the ultimate base and see if is a fundamental type.
          //
          Boolean fund_based (false);
          SemanticGraph::Type& ult_base (ultimate_base (e));
          {
            IsFundamentalType t (fund_based);
            t.dispatch (ult_base);
          }

          // Count enumerators.
          //
          UnsignedLong enum_count (0);

          for (Type::NamesIterator i (e.names_begin ()), end (e.names_end ());
               i != end; ++i)
            ++enum_count;

          //
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Enumeration class corresponding to the %" <<
              comment (e.name ()) << endl
               << " * schema type." << endl;

            if (e.annotated_p ())
            {
              os << " *" << endl;
              write_annotation (e.annotation ());
            }

            os << " */" << endl;
          }

          os << "class " << type_exp << name << ": public ";

          // Enumeration always has a base.
          //
          inherits (e, inherits_base_);

          os << "{"
             << "public:" << endl;

          if (string_based)
          {
            if (doxygen)
            {
              os << endl
                 << "/**" << endl
                 << " * @brief Underlying enum type." << endl
                 << " */" << endl;
            }

            if (enum_based)
            {
              os << "typedef ";

              inherits (e, inherits_base_);

              os << "::" << evalue (*base_enum) << " " << value << ";"
                 << endl;
            }
            else
            {
              os << "enum " << value
                 << "{";

              names<Enumeration> (e, names_, 0, 0, 0, &Enumeration::comma);

              os << "};";
            }
          }

          // default c-tor
          //
          if (options.value<CLI::generate_default_ctor> ())
          {
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Default constructor." << endl
                 << " *" << endl
                 << " * Note that this constructor may leave the " <<
                "instance in an" << endl
                 << " * invalid state." << endl
                 << " */" << endl;
            }

            os << name << " ();"
               << endl;
          }

          // c-tor (value)
          //
          if (string_based)
          {
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from the " <<
                "underlying enum value." << endl
                 << " *" << endl
                 << " * @param v A enum value." << endl
                 << " */" << endl;
            }

            os << name << " (" << value << " v);"
               << endl;
          }

          // c-tor (const char*)
          //
          if (string_based)
          {
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a C string." << endl
                 << " *" << endl
                 << " * @param v A string value." << endl
                 << " */" << endl;
            }

            os << name << " (const " << char_type << "* v);"
               << endl;
          }

          // c-tor (const std::string&)
          //
          if (string_based)
          {
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a string." << endl
                 << " *" << endl
                 << " * @param v A string value." << endl
                 << " */" << endl;
            }

            os << name << " (const " << string_type << "& v);"
               << endl;
          }

          // c-tor (fundamental)
          //
          if (fund_based)
          {
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a fundamental " <<
                "type value." << endl
                 << " *" << endl
                 << " * @param v A fundamental type value." << endl
                 << " */" << endl;
            }

            os << name << " (";

            member_.dispatch (ult_base);

            os << " v);"
               << endl;
          }

          // c-tor (base)
          //
          // If the ultimate is also our immediate base and it is a
          // fundamental type then this c-tor clashes with c-tor
          // (fundamental) above.
          //
          if (!fund_based || &ult_base != &e.inherits ().base ())
          {
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from the " <<
                "base value." << endl
                 << " *" << endl
                 << " * @param v A base value." << endl
                 << " */" << endl;
            }

            os << name << " (const ";

            inherits (e, inherits_member_);

            os << "& v);"
               << endl;
          }


          // c-tor (istream&)
          //
          Streams const& st (options.value<CLI::generate_extraction> ());
          for (Streams::ConstIterator i (st.begin ()); i != st.end (); ++i)
          {
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a data " <<
                "representation" << endl
                 << " * stream." << endl
                 << " *" << endl
                 << " * @param s A stream to extract the data from." << endl
                 << " * @param f Flags to create the new instance with." << endl
                 << " * @param c A pointer to the object that will " <<
                "contain the new" << endl
                 << " * instance." << endl
                 << " */" << endl;
            }

            os << name << " (" << istream_type << "< " << i->c_str () <<
              " >& s," << endl
               << flags_type << " f = 0," << endl
               << container << "* c = 0);"
               << endl;
          }

          if (!options.value<CLI::suppress_parsing> ())
          {
            // c-tor (xercesc::DOMElement)
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a DOM element." << endl
                 << " *" << endl
                 << " * @param e A DOM element to extract the data from." << endl
                 << " * @param f Flags to create the new instance with." << endl
                 << " * @param c A pointer to the object that will " <<
                "contain the new" << endl
                 << " * instance." << endl
                 << " */" << endl;
            }

            os << name << " (const " << xerces_ns << "::DOMElement& e," << endl
               << flags_type << " f = 0," << endl
               << container << "* c = 0);"
               << endl;

            // c-tor (xercesc::DOMAttr)
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a DOM attribute." << endl
                 << " *" << endl
                 << " * @param a A DOM attribute to extract the data from." << endl
                 << " * @param f Flags to create the new instance with." << endl
                 << " * @param c A pointer to the object that will " <<
                "contain the new" << endl
                 << " * instance." << endl
                 << " */" << endl;
            }

            os << name << " (const " << xerces_ns << "::DOMAttr& a," << endl
               << flags_type << " f = 0," << endl
               << container << "* c = 0);"
               << endl;

            // c-tor (std::basic_string const&, xercesc::DOMElement)
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a string fragment." << endl
                 << " *" << endl
                 << " * @param s A string fragment to extract the data from." << endl
                 << " * @param e A pointer to DOM element containing the " <<
                "string fragment." << endl
                 << " * @param f Flags to create the new instance with." << endl
                 << " * @param c A pointer to the object that will " <<
                "contain the new" << endl
                 << " * instance." << endl
                 << " */" << endl;
            }

            os << name << " (const " << string_type << "& s," << endl
               << "const " << xerces_ns << "::DOMElement* e," << endl
               << flags_type << " f = 0," << endl
               << container << "* c = 0);"
               << endl;
          }

          // copy c-tor
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Copy constructor." << endl
               << " *" << endl
               << " * @param x An instance to make a copy of." << endl
               << " * @param f Flags to create the copy with." << endl
               << " * @param c A pointer to the object that will contain " <<
              "the copy." << endl
               << " *" << endl
               << " * For polymorphic object models use the @c _clone " <<
              "function instead." << endl
               << " */" << endl;
          }

          os << name << " (const " << name << "& x," << endl
             << flags_type << " f = 0," << endl
             << container << "* c = 0);"
             << endl;

          // clone
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Copy the instance polymorphically." << endl
               << " *" << endl
               << " * @param f Flags to create the copy with." << endl
               << " * @param c A pointer to the object that will contain " <<
              "the copy." << endl
               << " * @return A pointer to the dynamically allocated copy." << endl
               << " *" << endl
               << " * This function ensures that the dynamic type of the " <<
              "instance is" << endl
               << " * used for copying and should be used for polymorphic " <<
              "object" << endl
               << " * models instead of the copy constructor." << endl
               << " */" << endl;
          }

          os << "virtual " << name << "*" << endl
             << "_clone (" << flags_type << " f = 0," << endl
             << container << "* c = 0) const;"
             << endl;

          // operator= (value)
          //
          if (string_based)
          {
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Assign the underlying enum value." << endl
                 << " *" << endl
                 << " * @param v A enum value." << endl
                 << " * @return A refernce to the instance." << endl
                 << " */" << endl;
            }

            os << name << "&" << endl
               << "operator= (" << value << " v);"
               << endl;
          }

          // operator value ()
          //
          // Name lookup differences in various compilers make generation
          // of this operator outside of the class a really hard task. So
          // we are going to make it "always inline".
          //
          if (string_based)
          {
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Implicit conversion operator to the " <<
                "underlying" << endl
                 << " * enum value." << endl
                 << " *" << endl
                 << " * @return A enum value." << endl
                 << " */" << endl;
            }

            os << "virtual" << endl
               << "operator " << value << " () const"
               << "{"
               << "return _xsd_" << name << "_convert ();"
               << "}";
          }

          //
          //
          if (string_based)
          {
            if (doxygen)
              os << "//@cond" << endl
                 << endl;

            os << "protected:" << endl
               << value << endl
               << "_xsd_" << name << "_convert () const;"
               << endl;

            os << "public:" << endl;

            if (enum_based)
            {
              // We are going to reuse our base's literals.
              //
              os << "static const " << char_type << "* const* " <<
                "_xsd_" << name << "_literals_;";
            }
            else
            {
              os << "static const " << char_type << "* const " <<
                "_xsd_" << name << "_literals_[" << enum_count << "];";
            }

            os << "static const " << value <<
              " _xsd_" << name << "_indexes_[" << enum_count << "];";

            if (doxygen)
              os << endl
                 << "//@endcond" << endl
                 << endl;
          }

          os << "};";
        }

        virtual Void
        comma (Type&)
        {
          os << "," << endl;
        }

      private:
        Traversal::Inherits inherits_base_;
        BaseTypeName base_;

        Traversal::Inherits inherits_member_;
        MemberTypeName member_;

        Traversal::Names names_;
        Enumerator enumerator_;
      };


      //
      //
      struct MemberFunction: Traversal::Member, Context
      {
        MemberFunction (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& m)
        {
          if (skip (m))
            return;

          String const& aname (eaname (m));
          String const& mname (emname (m));
          String kind (m.is_a<SemanticGraph::Element> ()
                       ? "element" : "attribute");

          Boolean fund (false);
          {
            IsFundamentalType t (fund);
            t.dispatch (m.type ());
          }

          Boolean def_attr (m.default_p () &&
                            m.is_a<SemanticGraph::Attribute> ());

          if (max (m) != 1)
          {
            // sequence
            //
            String container (econtainer (m));

            // container const&
            // name () const;
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Return a read-only (constant) reference " <<
                "to the element" << endl
                 << " * sequence." << endl
                 << " *" << endl
                 << " * @return A constant reference to the sequence " <<
                "container." << endl
                 << " */" << endl;
            }

            os << "const " << container << "&" << endl
               << aname << " () const;"
               << endl;

            // container&
            // name ();
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Return a read-write reference to the " <<
                "element sequence." << endl
                 << " *" << endl
                 << " * @return A reference to the sequence container." << endl
                 << " */" << endl;
            }

            os << container << "&" << endl
               << aname << " ();"
               << endl;

            // void
            // name (container const&);
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Copy elements from a given sequence." << endl
                 << " *" << endl
                 << " * @param s A sequence to copy elements from." << endl
                 << " *" << endl
                 << " * For each element in @a s this function " <<
                "makes a copy and adds it " << endl
                 << " * to the sequence. Note that this operation " <<
                "completely changes the " << endl
                 << " * sequence and all old elements will be lost." << endl
                 << " */" << endl;
            }

            os << "void" << endl
               << mname << " (const " << container << "& s);"
               << endl;
          }
          else if (min (m) == 0 && !def_attr)
          {
            // optional
            //
            String const& type (etype (m));
            String container (econtainer (m));

            // container const&
            // name () const;
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Return a read-only (constant) reference " <<
                "to the " << kind << endl
                 << " * container." << endl
                 << " *" << endl
                 << " * @return A constant reference to the optional " <<
                "container." << endl
                 << " */" << endl;
            }

            os << "const " << container << "&" << endl
               << aname << " () const;"
               << endl;

            // container&
            // name ();
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Return a read-write reference to the " <<
                kind << " container." << endl
                 << " *" << endl
                 << " * @return A reference to the optional container." << endl
                 << " */" << endl;
            }

            os << container << "&" << endl
               << aname << " ();"
               << endl;

            // void
            // name (type const&);
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Set the " << kind << " value." << endl
                 << " *" << endl
                 << " * @param x A new value to set." << endl
                 << " *" << endl
                 << " * This function makes a copy of its argument " <<
                "and sets it as" << endl
                 << " * the new value of the " << kind << "." << endl
                 << " */" << endl;
            }

            os << "void" << endl
               << mname << " (const " << type << "& x);"
               << endl;

            // void
            // name (container const&);
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Set the " << kind << " value." << endl
                 << " *" << endl
                 << " * @param x An optional container with the new value " <<
                "to set." << endl
                 << " *" << endl
                 << " * If the value is present in @a x then this function " <<
                "makes a copy " << endl
                 << " * of this value and sets it as the new value of the " <<
                kind << "." << endl
                 << " * Otherwise the " << kind << " container is set " <<
                "the 'not present' state." << endl
                 << " */" << endl;
            }

            os << "void" << endl
               << mname << " (const " << container << "& x);"
               << endl;

            // void
            // name (auto_ptr<type>);
            //
            if (!fund)
            {
              if (doxygen)
              {
                os << "/**" << endl
                   << " * @brief Set the " << kind << " value without " <<
                  "copying." << endl
                   << " *" << endl
                   << " * @param p A new value to use." << endl
                   << " *" << endl
                   << " * This function will try to use the passed value " <<
                  "directly instead" << endl
                   << " * of making a copy." << endl
                   << " */" << endl;
              }

              os << "void" << endl
                 << mname << " (::std::auto_ptr< " << type << " > p);"
                 << endl;
            }
          }
          else
          {
            // one
            //
            String const& type (etype (m));

            // type const&
            // name () const;
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Return a read-only (constant) reference " <<
                "to the " << kind << "." << endl
                 << " *" << endl
                 << " * @return A constant reference to the " << kind <<
                "." << endl
                 << " */" << endl;
            }

            os << "const " << type << "&" << endl
               << aname << " () const;"
               << endl;

            // Do not generate modifiers for fixed attributes.
            //
            if (!(def_attr && m.fixed_p ()))
            {
              // type&
              // name ();
              //
              if (doxygen)
              {
                os << "/**" << endl
                   << " * @brief Return a read-write reference to the " <<
                  kind << "." << endl
                   << " *" << endl
                   << " * @return A reference to the " << kind << "." << endl
                   << " */" << endl;
              }

              os << type << "&" << endl
                 << aname << " ();"
                 << endl;

              // void
              // name (type const&);
              //
              if (doxygen)
              {
                os << "/**" << endl
                   << " * @brief Set the " << kind << " value." << endl
                   << " *" << endl
                   << " * @param x A new value to set." << endl
                   << " *" << endl
                   << " * This function makes a copy of its argument " <<
                  "and sets it as" << endl
                   << " * the new value of the " << kind << "." << endl
                   << " */" << endl;
              }

              os << "void" << endl
                 << mname << " (const " << type << "& x);"
                 << endl;

              // void
              // name (auto_ptr<type>);
              //
              if (!fund)
              {
                if (doxygen)
                {
                  os << "/**" << endl
                     << " * @brief Set the " << kind << " value without " <<
                    "copying." << endl
                     << " *" << endl
                     << " * @param p A new value to use." << endl
                     << " *" << endl
                     << " * This function will try to use the passed value " <<
                    "directly" << endl
                     << " * instead of making a copy." << endl
                     << " */" << endl;
                }

                os << "void" << endl
                   << mname << " (::std::auto_ptr< " << type << " > p);"
                   << endl;

              }

              // auto_ptr<type>
              // detach_name ();
              //
              if (detach && !fund)
              {
                if (doxygen)
                {
                  os << "/**" << endl
                     << " * @brief Detach the " << kind << " value from " <<
                    "the object model." << endl
                     << " *" << endl
                     << " * @return A pointer to the " << kind << " value." << endl
                     << " *" << endl
                     << " * Note that this function leaves the required " <<
                    kind << " in " << endl
                     << " * the original object model uninitialized." << endl
                     << " */" << endl;
                }

                os << "::std::auto_ptr< " << type << " >" << endl
                   << edname (m) << " ();"
                   << endl;
              }
            }
          }

          // default_value
          //
          if (m.default_p ())
          {
            Boolean simple (true);

            if (m.is_a<SemanticGraph::Element> ())
            {
              IsSimpleType test (simple);
              test.dispatch (m.type ());
            }

            if (simple)
            {
              Boolean lit (false);
              {
                IsLiteralValue test (lit);
                test.dispatch (m.type ());
              }

              if (doxygen)
              {
                os << "/**" << endl
                   << " * @brief Return the default value for the " <<
                  kind << "." << endl
                   << " *" << endl;

                if (lit)
                  os << " * @return The " << kind << "'s default value." << endl;
                else
                  os << " * @return A read-only (constant) reference to the "
                     << kind << "'s" << endl
                     << " * default value." << endl;

                os << " */" << endl;
              }

              if (lit)
                os << "static " << etype (m) << endl;
              else
                os << "static const " << etype (m) << "&" << endl;

              os << edefault_value (m) << " ();"
                 << endl;
            }
          }
        }
      };

      struct AnyFunction: Traversal::Any, Context
      {
        AnyFunction (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (SemanticGraph::Any& a)
        {
          String const& aname (eaname (a));
          String const& mname (emname (a));

          SemanticGraph::Complex& c (
            dynamic_cast<SemanticGraph::Complex&> (a.scope ()));

          if (max (a) != 1)
          {
            // sequence
            //
            String container (econtainer (a));

            // container const&
            // name () const;
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Return a read-only (constant) reference " <<
                "to the wildcard" << endl
                 << " * element sequence." << endl
                 << " *" << endl
                 << " * @return A constant reference to the sequence " <<
                "container." << endl
                 << " */" << endl;
            }

            os << "const " << container << "&" << endl
               << aname << " () const;"
               << endl;

            // container&
            // name ();
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Return a read-write reference to the " <<
                "wildcard element" << endl
                 << " * sequence." << endl
                 << " *" << endl
                 << " * @return A reference to the sequence container." << endl
                 << " */" << endl;
            }

            os << container << "&" << endl
               << aname << " ();"
               << endl;

            // void
            // name (container const&);
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Copy elements from a given sequence." << endl
                 << " *" << endl
                 << " * @param s A sequence to copy elements from." << endl
                 << " *" << endl
                 << " * For each element in @a s this function " <<
                "makes a copy and adds" << endl
                 << " * it to the wildcard element sequence. Note that " <<
                "this operation" << endl
                 << " * completely changes the sequence and all old " <<
                "elements will be" << endl
                 << " * lost." << endl
                 << " */" << endl;
            }

            os << "void" << endl
               << mname << " (const " << container << "& s);"
               << endl;
          }
          else if (min (a) == 0)
          {
            // optional
            //
            String container (econtainer (a));

            // container const&
            // name () const;
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Return a read-only (constant) reference " <<
                "to the wildcard" << endl
                 << " * element container." << endl
                 << " *" << endl
                 << " * @return A constant reference to the optional " <<
                "container." << endl
                 << " */" << endl;
            }

            os << "const " << container << "&" << endl
               << aname << " () const;"
               << endl;

            // container&
            // name ();
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Return a read-write reference to the " <<
                "wildcard element" << endl
                 << " * container." << endl
                 << " *" << endl
                 << " * @return A reference to the optional container." << endl
                 << " */" << endl;
            }

            os << container << "&" << endl
               << aname << " ();"
               << endl;

            // void
            // name (type const&);
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Set the wildcard content." << endl
                 << " *" << endl
                 << " * @param e A new element to set." << endl
                 << " *" << endl
                 << " * This function makes a copy of its argument " <<
                "and sets it as" << endl
                 << " * the new wildcard content." << endl
                 << " */" << endl;
            }

            os << "void" << endl
               << mname << " (const " << xerces_ns << "::DOMElement& e);"
               << endl;

            // void
            // name (type*);
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Set the wildcard content without copying." << endl
                 << " *" << endl
                 << " * @param p A new element to use." << endl
                 << " *" << endl
                 << " * This function will use the passed element " <<
                "directly instead" << endl
                 << " * of making a copy. For this to work the element " <<
                "should belong" << endl
                 << " * to the DOM document associated with this instance." << endl
                 << " *" << endl
                 << " * @see " << edom_document (c) << endl
                 << " */" << endl;
            }

            os << "void" << endl
               << mname << " (" << xerces_ns << "::DOMElement* p);"
               << endl;

            // void
            // name (container const&);
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Set the wildcard content." << endl
                 << " *" << endl
                 << " * @param x An optional container with the new " <<
                "element to set." << endl
                 << " *" << endl
                 << " * If the element is present in @a x then this function " <<
                "makes a " << endl
                 << " * copy of this element and sets it as the new wildcard " <<
                "content." << endl
                 << " * Otherwise the element container is set the 'not " <<
                "present' state." << endl
                 << " */" << endl;
            }

            os << "void" << endl
               << mname << " (const " << container << "& x);"
               << endl;
          }
          else
          {
            // one
            //

            // type const&
            // name () const;
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Return a read-only (constant) reference " <<
                "to the wildcard" << endl
                 << " * element." << endl
                 << " *" << endl
                 << " * @return A constant reference to the DOM element." << endl
                 << " */" << endl;
            }

            os << "const " << xerces_ns << "::DOMElement&" << endl
               << aname << " () const;"
               << endl;

            // type&
            // name ();
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Return a read-write reference to the " <<
                "wildcard element." << endl
                 << " *" << endl
                 << " * @return A reference to the DOM element." << endl
                 << " */" << endl;
            }

            os << xerces_ns << "::DOMElement&" << endl
               << aname << " ();"
               << endl;

            // void
            // name (type const&);
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Set the wildcard content." << endl
                 << " *" << endl
                 << " * @param e A new element to set." << endl
                 << " *" << endl
                 << " * This function makes a copy of its argument " <<
                "and sets it as" << endl
                 << " * the new wildcard content." << endl
                 << " */" << endl;
            }

            os << "void" << endl
               << mname << " (const " << xerces_ns << "::DOMElement& e);"
               << endl;

            // void
            // name (const*);
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Set the wildcard content without copying." << endl
                 << " *" << endl
                 << " * @param p A new element to use." << endl
                 << " *" << endl
                 << " * This function will use the passed element " <<
                "directly instead" << endl
                 << " * of making a copy. For this to work the element " <<
                "should belong" << endl
                 << " * to the DOM document associated with this instance." << endl
                 << " *" << endl
                 << " * @see " << edom_document (c) << endl
                 << " */" << endl;
            }

            os << "void" << endl
               << mname << " (" << xerces_ns << "::DOMElement* p);"
               << endl;
          }
        }

        virtual Void
        traverse (SemanticGraph::AnyAttribute& a)
        {
          String const& aname (eaname (a));
          String const& mname (emname (a));

          String container (econtainer (a));

          // container const&
          // name () const;
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Return a read-only (constant) reference " <<
              "to the" << endl
               << " * attribute set." << endl
               << " *" << endl
               << " * @return A constant reference to the set " <<
              "container." << endl
               << " */" << endl;
          }

          os << "const " << container << "&" << endl
             << aname << " () const;"
             << endl;

          // container&
          // name ();
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Return a read-write reference to the " <<
              "attribute set." << endl
               << " *" << endl
               << " * @return A reference to the set container." << endl
               << " */" << endl;
          }

          os << container << "&" << endl
             << aname << " ();"
             << endl;

          // void
          // name (container const&);
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Copy attributes from a given set." << endl
               << " *" << endl
               << " * @param s A set to copy elements from." << endl
               << " *" << endl
               << " * For each attribute in @a s this function " <<
              "makes a copy and adds" << endl
               << " * it to the set. Note that this operation " <<
              "completely changes the " << endl
               << " * set and all old attributes will be lost." << endl
               << " */" << endl;
          }

          os << "void" << endl
             << mname << " (const " << container << "& s);"
             << endl;
        }
      };

      //
      //
      struct Member: Traversal::Member, Context
      {
        Member (Context& c)
            : Context (c),
              type_name_ (c),
              member_function_ (c)
        {
          belongs_ >> type_name_;
        }

        virtual Void
        traverse (Type& m)
        {
          if (skip (m))
            return;

          String const& type (etype (m));
          Boolean el (m.is_a<SemanticGraph::Element> ());

          Boolean def_attr (m.default_p () && !el);

          if (doxygen)
          {
            os << "/**" << endl
               << " * @name " << comment (m.name ()) << endl
               << " *" << endl
               << " * @brief Accessor and modifier functions for the %" <<
              comment (m.name ()) << endl
               << " * ";

            if (max (m) != 1)
            {
              os << "sequence element." << endl;
            }
            else if (min (m) == 0)
            {
              if (def_attr)
                os << "optional attribute with a default value." << endl;
              else
                os << "optional " << (el ? "element." : "attribute.") << endl;
            }
            else
            {
              os << "required " << (el ? "element." : "attribute.") << endl;
            }

            if (m.annotated_p ())
            {
              os << " *" << endl;
              write_annotation (m.annotation ());
            }

            os << " */" << endl
               << "//@{" << endl;
          }
          else
          {
            os << "// " << comment (m.name ()) << endl
               << "// " << endl;
          }

          // Typedefs.
          //
          if (doxygen)
          {
            os << endl
               << "/**" << endl
               << " * @brief " << (el ? "Element" : "Attribute") <<
              " type." << endl
               << " */" << endl;
          }

          os << "typedef ";

          belongs (m, belongs_);

          os << " " << type << ";";

          if (max (m) != 1)
          {
            String const& container (econtainer (m));
            Boolean isense (options.value<CLI::generate_intellisense> ());

            // sequence
            //
            if (doxygen)
            {
              os << endl
                 << "/**" << endl
                 << " * @brief Element sequence container type." << endl
                 << " */" << endl;
            }

            os << "typedef ::xsd::cxx::tree::sequence< " << type << " > " <<
              container << ";";

            if (doxygen)
            {
              os << endl
                 << "/**" << endl
                 << " * @brief Element iterator type." << endl
                 << " */" << endl;
            }

            // IntelliSense does not not like aliases and fully-qualified
            // names here.
            //
            if (!isense)
              os << "typedef " << container << "::iterator " <<
                eiterator (m) << ";";
            else
              os << "typedef xsd::cxx::tree::sequence< " << type <<
                " >::iterator " << eiterator (m) << ";";

            if (doxygen)
            {
              os << endl
                 << "/**" << endl
                 << " * @brief Element constant iterator type." << endl
                 << " */" << endl;
            }

            if (!isense)
              os << "typedef " << container << "::const_iterator " <<
                econst_iterator (m) << ";";
            else
              os << "typedef xsd::cxx::tree::sequence< " << type <<
                " >::const_iterator " << econst_iterator (m) << ";";

          }
          else if (min (m) == 0 && !def_attr)
          {
            // optional
            //
            if (doxygen)
            {
              os << endl
                 << "/**" << endl
                 << " * @brief " << (el ? "Element" : "Attribute") <<
                " optional container type." << endl
                 << " */" << endl;
            }

            os << "typedef ::xsd::cxx::tree::optional< " << type << " > " <<
              econtainer (m) << ";";
          }
          else
          {
            // one
            //
          }

          if (doxygen)
          {
            os << endl
               << "/**" << endl
               << " * @brief " << (el ? "Element" : "Attribute") <<
              " traits type." << endl
               << " */" << endl;
          }
          os << "typedef ::xsd::cxx::tree::traits< " << type << ", " <<
            char_type;

          SemanticGraph::Type& t (m.type ());

          if (t.is_a<SemanticGraph::Fundamental::Double> ())
            os << ", ::xsd::cxx::tree::schema_type::double_";
          else if (t.is_a<SemanticGraph::Fundamental::Decimal> ())
            os << ", ::xsd::cxx::tree::schema_type::decimal";

          os << " > " << etraits (m) << ";"
             << endl;

          member_function_.traverse (m);

          if (doxygen)
          {
            os << "//@}" << endl
               << endl;
          }
        }

      private:
        MemberTypeName type_name_;
        Traversal::Belongs belongs_;

        MemberFunction member_function_;
      };


      struct Any: Traversal::Any,
                  Traversal::AnyAttribute,
                  Context
      {
        Any (Context& c)
            : Context (c), any_function_ (c)
        {
        }

        virtual Void
        traverse (SemanticGraph::Any& a)
        {
          if (doxygen)
          {
            os << "/**" << endl
               << " * @name " << ename (a) << endl
               << " *" << endl
               << " * @brief Accessor and modifier functions for the " <<
              "any wildcard." << endl;

            if (a.annotated_p ())
            {
              os << " *" << endl;
              write_annotation (a.annotation ());
            }

            os << " */" << endl
               << "//@{" << endl;
          }
          else
          {
            os << "// " << ename (a) << endl
               << "// " << endl;
          }

          // Typedefs.
          //
          if (max (a) != 1)
          {
            String const& container (econtainer (a));

            // sequence
            //
            if (doxygen)
            {
              os << endl
                 << "/**" << endl
                 << " * @brief DOM element sequence container type." << endl
                 << " */" << endl;
            }

            os << "typedef ::xsd::cxx::tree::element_sequence " <<
              container << ";";

            if (doxygen)
            {
              os << endl
                 << "/**" << endl
                 << " * @brief DOM element iterator type." << endl
                 << " */" << endl;
            }

            os << "typedef " << container << "::iterator " <<
              eiterator (a) << ";";

            if (doxygen)
            {
              os << endl
                 << "/**" << endl
                 << " * @brief DOM element constant iterator type." << endl
                 << " */" << endl;
            }

            os << "typedef " << container << "::const_iterator " <<
              econst_iterator (a) << ";"
               << endl;

          }
          else if (min (a) == 0)
          {
            // optional
            //
            if (doxygen)
            {
              os << endl
                 << "/**" << endl
                 << " * @brief DOM element optional container type." << endl
                 << " */" << endl;
            }

            os << "typedef ::xsd::cxx::tree::element_optional " <<
              econtainer (a) << ";"
               << endl;
          }
          else
          {
            // one
            //
            if (doxygen)
              os << endl;
          }

          any_function_.traverse (a);

          if (doxygen)
          {
            os << "//@}" << endl
               << endl;
          }
        }

        virtual Void
        traverse (SemanticGraph::AnyAttribute& a)
        {
          String const& container (econtainer (a));

          if (doxygen)
          {
            os << "/**" << endl
               << " * @name " << ename (a) << endl
               << " *" << endl
               << " * @brief Accessor and modifier functions for the " <<
              "anyAttribute" << endl
               << " * wildcard." << endl;

            if (a.annotated_p ())
            {
              os << " *" << endl;
              write_annotation (a.annotation ());
            }

            os << " */" << endl
               << "//@{" << endl;
          }
          else
          {
            os << "// " << ename (a) << endl
               << "// " << endl;
          }

          if (doxygen)
          {
            os << endl
               << "/**" << endl
               << " * @brief DOM attribute set container type." << endl
               << " */" << endl;
          }

          os << "typedef ::xsd::cxx::tree::attribute_set< " << char_type <<
            " > " << container << ";";

          if (doxygen)
          {
            os << endl
               << "/**" << endl
               << " * @brief DOM attribute iterator type." << endl
               << " */" << endl;
          }

          os << "typedef " << container << "::iterator " <<
            eiterator (a) << ";";

          if (doxygen)
          {
            os << endl
               << "/**" << endl
               << " * @brief DOM attribute constant iterator type." << endl
               << " */" << endl;
          }

          os << "typedef " << container << "::const_iterator " <<
            econst_iterator (a) << ";"
             << endl;

          any_function_.traverse (a);

          if (doxygen)
          {
            os << "//@}" << endl
               << endl;
          }
        }

      private:
        AnyFunction any_function_;
      };

      struct DataMember: Traversal::Member, Context
      {
        DataMember (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (Type& m)
        {
          if (skip (m)) return;

          String const& member (emember (m));

          Boolean def_attr (m.default_p () &&
                            m.is_a<SemanticGraph::Attribute> ());

          if (max (m) != 1)
          {
            // sequence
            //
            os << econtainer (m) << " " << member << ";";
          }
          else if (min (m) == 0 && !def_attr)
          {
            // optional
            //
            os << econtainer (m) << " " << member << ";";
          }
          else
          {
            // one
            //
            os << "::xsd::cxx::tree::one< " << etype (m) << " > " <<
              member << ";";
          }

          // default_value
          //
          if (m.default_p ())
          {
            Boolean simple (true);

            if (m.is_a<SemanticGraph::Element> ())
            {
              IsSimpleType test (simple);
              test.dispatch (m.type ());
            }

            if (simple)
            {
              Boolean lit (false);
              {
                IsLiteralValue test (lit);
                test.dispatch (m.type ());
              }

              if (!lit)
              {
                os << "static const " << etype (m) << " " <<
                  edefault_value_member (m) << ";";
              }
            }
          }
        }
      };

      struct DataAny: Traversal::Any,
                      Traversal::AnyAttribute,
                      Context
      {
        DataAny (Context& c)
            : Context (c)
        {
        }

        virtual Void
        traverse (SemanticGraph::Any& a)
        {
          String const& member (emember (a));

          if (max (a) != 1)
          {
            // sequence
            //
            os << econtainer (a) << " " << member << ";";
          }
          else if (min (a) == 0)
          {
            // optional
            //
            os << econtainer (a) << " " << member << ";";
          }
          else
          {
            // one
            //
            os << "::xsd::cxx::tree::element_one " << member << ";";
          }
        }

        virtual Void
        traverse (SemanticGraph::AnyAttribute& a)
        {
          os << econtainer (a) << " " << emember (a) << ";";
        }
      };


      struct Complex: Traversal::Complex, Context
      {
        Complex (Context& c)
            : Context (c),
              base_name_ (c),
              member_name_ (c),
              any_ (c),
              member_ (c),
              data_any_ (c),
              data_member_ (c)
        {
          inherits_base_ >> base_name_;
          inherits_member_ >> member_name_;

          names_ >> member_;
          if (options.value<CLI::generate_wildcard> ())
            names_ >> any_;

          names_data_ >> data_member_;
          if (options.value<CLI::generate_wildcard> ())
            names_data_ >> data_any_;
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

          Boolean has_members (has<Traversal::Member> (c));

          Boolean hae (has<Traversal::Any> (c));
          Boolean haa (has<Traversal::AnyAttribute> (c));

          Boolean gen_wildcard (options.value<CLI::generate_wildcard> ());

          Boolean simple (true);
          {
            IsSimpleType t (simple);
            t.dispatch (c);
          }

          Boolean string_based (false);
          {
            IsStringBasedType t (string_based);
            t.dispatch (c);
          }

          SemanticGraph::Enumeration* enum_base (0);
          {
            IsEnumBasedType t (enum_base);
            t.dispatch (c);
          }

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Class corresponding to the %" <<
	       comment (c.name ()) << " schema type." << endl;

            if (c.annotated_p ())
            {
              os << " *" << endl;
              write_annotation (c.annotation ());
            }

            os << " *" << endl
               << " * @nosubgrouping" << endl
               << " */" << endl;
          }

          os << "class " << type_exp << name << ": public ";

          if (c.inherits_p ())
            inherits (c, inherits_base_);
          else
            os << any_type;

          os << "{"
             << "public:" << endl;

          // Members.
          //
          names (c, names_);

          // dom_document accessors.
          //
          if (edom_document_member_p (c))
          {

            if (!doxygen)
            {
              os << "// DOMDocument for wildcard content." << endl
                 << "//" << endl;
            }

            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Return a read-only (constant) reference " <<
                "to the DOM" << endl
                 << " * document associated with this instance." << endl
                 << " *" << endl
                 << " * @return A constant reference to the DOM document." << endl
                 << " *" << endl
                 << " * The DOM document returned by this function is " <<
                "used to store" << endl
                 << " * the raw XML content corresponding to wildcards." << endl
                 << " */" << endl;
            }

            os << "const " << xerces_ns << "::DOMDocument&" << endl
               << edom_document (c) << " () const;"
               << endl;

            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Return a read-write reference to the DOM " <<
                "document" << endl
                 << " * associated with this instance." << endl
                 << " *" << endl
                 << " * @return A reference to the DOM document." << endl
                 << " *" << endl
                 << " * The DOM document returned by this function is " <<
                "used to store" << endl
                 << " * the raw XML content corresponding to wildcards." << endl
                 << " */" << endl;
            }
            os << xerces_ns << "::DOMDocument&" << endl
               << edom_document (c) << " ();"
               << endl;
          }

          if (doxygen)
          {
            os << "/**" << endl
               << " * @name Constructors" << endl
               << " */" << endl
               << "//@{" << endl
               << endl;
          }
          else
          {
            os << "// Constructors." << endl
               << "//" << endl;
          }

          Boolean generate_no_base_ctor (false);
          {
            GenerateWithoutBaseCtor t (generate_no_base_ctor);
            t.traverse (c);
          }

          Boolean has_complex_non_op_args (false);
          Boolean has_poly_non_op_args (false);
          Boolean complex_poly_args_clash (true);
          {
            HasComplexPolyNonOptArgs t (*this, true,
                                        has_complex_non_op_args,
                                        has_poly_non_op_args,
                                        complex_poly_args_clash);
            t.traverse (c);
          }

          // default c-tor
          //
          if (options.value<CLI::generate_default_ctor> ())
          {
            // c-tor (ultimate-base, all-non-optional-members) will become
            // default c-tor if our inheritance hierarchy has no required
            // members and no simple base. We can also collide with
            // c-tor (all-non-optional-members) if we have no required
            // members.
            //
            Boolean generate (false);
            {
              GenerateDefaultCtor t (*this, generate, generate_no_base_ctor);
              t.traverse (c);
            }

            if (generate)
            {
              if (doxygen)
              {
                os << "/**" << endl
                   << " * @brief Default constructor." << endl
                   << " *" << endl
                   << " * Note that this constructor leaves required " <<
                  "elements and" << endl
                   << " * attributes uninitialized." << endl
                   << " */" << endl;
              }

              os << name << " ();"
                 << endl;
            }
          }

          // c-tor (base, all-non-optional-members)
          //
          if (options.value<CLI::generate_from_base_ctor> ())
          {
            // c-tor (base, all-non-optional-members) will be equivalent to
            // c-tor (ultimate-base, all-non-optional-members) unless our
            // immediate base's hierarchy has some non-optional members.
            // We also need to generate this c-tor when one of the types
            // in our inheritance hierarchy was customized since the
            // customized version may not necessarily be convertible to
            // the base without loss of information.
            //
            Boolean generate (false);
            {
              GenerateFromBaseCtor t (*this, generate);
              t.traverse (c);
            }

            if (generate)
            {
              Boolean has_complex_non_op_args (false);
              Boolean has_poly_non_op_args (false);
              Boolean complex_poly_args_clash (true);
              {
                HasComplexPolyNonOptArgs t (*this, false,
                                            has_complex_non_op_args,
                                            has_poly_non_op_args,
                                            complex_poly_args_clash);
                t.traverse (c);
              }

              //
              //
              if (doxygen)
              {
                os << "/**" << endl
                   << " * @brief Create an instance from the immediate "
                  "base and" << endl
                   << " * initializers for required elements and "
                   << "attributes." << endl
                   << " */" << endl;
              }

              os << name << " (const ";
              inherits (c, inherits_member_);
              os << "&";
              {
                FromBaseCtorArg args (*this, FromBaseCtorArg::arg_type, false);
                Traversal::Names args_names (args);
                names (c, args_names);
              }
              os << ");"
                 << endl;

              // If we have any complex arguments in the previous c-tor
              // then also generate the auto_ptr version.
              //
              if (has_complex_non_op_args)
              {
                if (doxygen)
                {
                  os << "/**" << endl
                     << " * @brief Create an instance from the immediate "
                    "base and" << endl
                     << " * initializers for required elements and "
                     << "attributes" << endl
                     << " * (auto_ptr version)." << endl
                     << " *" << endl
                     << " * This constructor will try to use the passed " <<
                    "values directly" << endl
                     << " * instead of making copies." << endl
                     << " */" << endl;
                }

                os << name << " (const ";
                inherits (c, inherits_member_);
                os << "&";
                {
                  FromBaseCtorArg args (
                    *this, FromBaseCtorArg::arg_complex_auto_ptr, false);
                  Traversal::Names args_names (args);
                  names (c, args_names);
                }
                os << ");"
                   << endl;
              }

              // If we are generating polymorphic code then we also need to
              // provide auto_ptr version for every polymorphic type.
              //
              if (polymorphic &&
                  has_poly_non_op_args && !complex_poly_args_clash)
              {
                if (doxygen)
                {
                  os << "/**" << endl
                     << " * @brief Create an instance from the immediate "
                    "base and" << endl
                     << " * initializers for required elements and "
                     << "attributes" << endl
                     << " * (auto_ptr version)." << endl
                     << " *" << endl
                     << " * This constructor will try to use the passed " <<
                    "values directly" << endl
                     << " * instead of making copies." << endl
                     << " */" << endl;
                }

                os << name << " (const ";
                inherits (c, inherits_member_);
                os << "&";
                {
                  FromBaseCtorArg args (
                    *this, FromBaseCtorArg::arg_poly_auto_ptr, false);
                  Traversal::Names args_names (args);
                  names (c, args_names);
                }
                os << ");"
                   << endl;
              }
            }
          }

          // c-tor (all-non-optional-members)
          //
          if (generate_no_base_ctor)
          {
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from initializers " <<
                "for required " << endl
                 << " * elements and attributes." << endl
                 << " */" << endl;
            }

            os << name << " (";
            {
              CtorArgsWithoutBase ctor_args (
                *this, CtorArgsWithoutBase::arg_type, false, true);
              ctor_args.dispatch (c);
            }
            os << ");"
               << endl;


            // If we have any complex arguments in the previous c-tor
            // then also generate the auto_ptr version. One case where
            // this c-tor will be generated is restriction of anyType.
            //
            if (has_complex_non_op_args)
            {
              if (doxygen)
              {
                os << "/**" << endl
                   << " * @brief Create an instance from initializers " <<
                  "for required " << endl
                   << " * elements and attributes (auto_ptr version)." << endl
                   << " *" << endl
                   << " * This constructor will try to use the passed " <<
                  "values directly" << endl
                   << " * instead of making copies." << endl
                   << " */" << endl;
              }

              os << name << " (";
              {
                CtorArgsWithoutBase ctor_args (
                  *this, CtorArgsWithoutBase::arg_complex_auto_ptr, false, true);
                ctor_args.dispatch (c);
              }
              os << ");"
                 << endl;
            }

            // If we are generating polymorphic code then we also need to
            // provide auto_ptr version for every polymorphic type.
            //
            if (polymorphic &&
                has_poly_non_op_args && !complex_poly_args_clash)
            {
              if (doxygen)
              {
                os << "/**" << endl
                   << " * @brief Create an instance from initializers " <<
                  "for required " << endl
                   << " * elements and attributes (auto_ptr version)." << endl
                   << " *" << endl
                   << " * This constructor will try to use the passed " <<
                  "values directly" << endl
                   << " * instead of making copies." << endl
                   << " */" << endl;
              }

              os << name << " (";
              {
                CtorArgsWithoutBase ctor_args (
                  *this, CtorArgsWithoutBase::arg_poly_auto_ptr, false, true);
                ctor_args.dispatch (c);
              }
              os << ");"
                 << endl;
            }
          }

          if (string_based)
          {
            // We might not have the value type if this enum is customized.
            //
            if (enum_base != 0 && enum_base->context ().count ("value"))
            {
              // c-tor (enum-value, all-non-optional-members)
              //
              if (doxygen)
              {
                os << "/**" << endl
                   << " * @brief Create an instance from the " <<
                  "underlying enum value" << endl
                   << " * and initializers for required elements and " <<
                  "attributes." << endl
                   << " */" << endl;
              }

              os << name << " (" << fq_name (*enum_base) << "::" <<
                evalue (*enum_base);

              {
                CtorArgsWithoutBase ctor_args (
                  *this, CtorArgsWithoutBase::arg_type, false, false);
                ctor_args.dispatch (c);
              }

              os << ");"
                 << endl;
            }

            // c-tor (const char*, all-non-optional-members)
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a C string and " <<
                "initializers" << endl
                 << " * for required elements and attributes." << endl
                 << " */" << endl;
            }

            os << name << " (const " << char_type << "*";

            {
              CtorArgsWithoutBase ctor_args (
                *this, CtorArgsWithoutBase::arg_type, false, false);
              ctor_args.dispatch (c);
            }

            os << ");"
               << endl;

            // c-tor (const std::string&, all-non-optional-members)
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a string and" <<
                "initializers" << endl
                 << " * for required elements and attributes." << endl
                 << " */" << endl;
            }

            os << name << " (const " << string_type << "&";

            {
              CtorArgsWithoutBase ctor_args (
                *this, CtorArgsWithoutBase::arg_type, false, false);
              ctor_args.dispatch (c);
            }

            os << ");"
               << endl;
          }

          // c-tor (ultimate-base, all-non-optional-members)
          //

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Create an instance from the ultimate "
              "base and" << endl
               << " * initializers for required elements and " <<
              "attributes." << endl
               << " */" << endl;
          }

          os << name << " (";

          {
            CtorArgs ctor_args (*this, CtorArgs::arg_type);
            ctor_args.dispatch (c);
          }

          os << ");"
             << endl;

          // If we have any complex arguments in the previous c-tor
          // then also generate the auto_ptr version.
          //
          if (has_complex_non_op_args)
          {
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from the ultimate "
                "base and" << endl
                 << " * initializers for required elements and " <<
                "attributes" << endl
                 << " * (auto_ptr version)." << endl
                 << " *" << endl
                 << " * This constructor will try to use the passed " <<
                "values directly" << endl
                 << " * instead of making copies." << endl
                 << " */" << endl;
            }

            os << name << " (";

            {
              CtorArgs ctor_args (*this, CtorArgs::arg_complex_auto_ptr);
              ctor_args.dispatch (c);
            }

            os << ");"
               << endl;
          }

          // If we are generating polymorphic code then we also need to
          // provide auto_ptr version for every polymorphic type.
          //
          if (polymorphic && has_poly_non_op_args && !complex_poly_args_clash)
          {
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from the ultimate "
                "base and" << endl
                 << " * initializers for required elements and " <<
                "attributes" << endl
                 << " * (auto_ptr version)." << endl
                 << " *" << endl
                 << " * This constructor will try to use the passed " <<
                "values directly" << endl
                 << " * instead of making copies." << endl
                 << " */" << endl;
            }

            os << name << " (";

            {
              CtorArgs ctor_args (*this, CtorArgs::arg_poly_auto_ptr);
              ctor_args.dispatch (c);
            }

            os << ");"
               << endl;
          }

          // c-tor (istream&)
          //
          Streams const& st (options.value<CLI::generate_extraction> ());
          for (Streams::ConstIterator i (st.begin ()); i != st.end (); ++i)
          {
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a data " <<
                "representation" << endl
                 << " * stream." << endl
                 << " *" << endl
                 << " * @param s A stream to extract the data from." << endl
                 << " * @param f Flags to create the new instance with." << endl
                 << " * @param c A pointer to the object that will " <<
                "contain the new" << endl
                 << " * instance." << endl
                 << " */" << endl;
            }

            os << name << " (" << istream_type << "< " << i->c_str () <<
              " >& s," << endl
               << flags_type << " f = 0," << endl
               << container << "* c = 0);"
               << endl;
          }


          if (!options.value<CLI::suppress_parsing> ())
          {
            // c-tor (xercesc::DOMElement)
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a DOM element." << endl
                 << " *" << endl
                 << " * @param e A DOM element to extract the data from." << endl
                 << " * @param f Flags to create the new instance with." << endl
                 << " * @param c A pointer to the object that will " <<
                "contain the new" << endl
                 << " * instance." << endl
                 << " */" << endl;
            }

            os << name << " (const " << xerces_ns << "::DOMElement& e," << endl
               << flags_type << " f = 0," << endl
               << container << "* c = 0);"
               << endl;


            if (simple)
            {
              // c-tor (xercesc::DOMAttr)
              //
              if (doxygen)
              {
                os << "/**" << endl
                   << " * @brief Create an instance from a DOM attribute." << endl
                   << " *" << endl
                   << " * @param a A DOM attribute to extract the data from." << endl
                   << " * @param f Flags to create the new instance with." << endl
                   << " * @param c A pointer to the object that will " <<
                  "contain the new" << endl
                   << " * instance." << endl
                   << " */" << endl;
              }

              os << name << " (const " << xerces_ns << "::DOMAttr& a," << endl
                 << flags_type << " f = 0," << endl
                 << container << "* c = 0);"
                 << endl;

              // c-tor (std::basic_string const&, xercesc::DOMElement)
              //
              if (doxygen)
              {
                os << "/**" << endl
                   << " * @brief Create an instance from a string fragment." << endl
                   << " *" << endl
                   << " * @param s A string fragment to extract the data from." << endl
                   << " * @param e A pointer to DOM element containing the " <<
                  "string fragment." << endl
                   << " * @param f Flags to create the new instance with." << endl
                   << " * @param c A pointer to the object that will " <<
                  "contain the new" << endl
                   << " * instance." << endl
                   << " */" << endl;
              }

              os << name << " (const " << string_type << "& s," << endl
                 << "const " << xerces_ns << "::DOMElement* e," << endl
                 << flags_type << " f = 0," << endl
                 << container << "* c = 0);"
                 << endl;
            }
          }

          // copy c-tor
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Copy constructor." << endl
               << " *" << endl
               << " * @param x An instance to make a copy of." << endl
               << " * @param f Flags to create the copy with." << endl
               << " * @param c A pointer to the object that will contain " <<
              "the copy." << endl
               << " *" << endl
               << " * For polymorphic object models use the @c _clone " <<
              "function instead." << endl
               << " */" << endl;
          }

          os << name << " (const " << name << "& x," << endl
             << flags_type << " f = 0," << endl
             << container << "* c = 0);"
             << endl;

          // clone
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Copy the instance polymorphically." << endl
               << " *" << endl
               << " * @param f Flags to create the copy with." << endl
               << " * @param c A pointer to the object that will contain " <<
              "the copy." << endl
               << " * @return A pointer to the dynamically allocated copy." << endl
               << " *" << endl
               << " * This function ensures that the dynamic type of the " <<
              "instance is" << endl
               << " * used for copying and should be used for polymorphic " <<
              "object" << endl
               << " * models instead of the copy constructor." << endl
               << " */" << endl;
          }

          os << "virtual " << name << "*" << endl
             << "_clone (" << flags_type << " f = 0," << endl
             << container << "* c = 0) const;"
             << endl;

          if (doxygen)
          {
            os << "//@}" << endl
               << endl;
          }

          // d-tor
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Destructor." << endl
               << " */" << endl;
          }

          os << "virtual " << endl
             << "~" << name << " ();"
             << endl;

          // Data members and implementation functions.
          //
          if (has_members || hae || (haa && gen_wildcard))
          {
            os << "// Implementation." << endl
               << "//" << endl;

            if (doxygen)
              os << endl
                 << "//@cond" << endl
                 << endl;

            if (!options.value<CLI::suppress_parsing> ())
            {
              // parse (xercesc::DOMElement)
              //
              os << "protected:" << endl
                 << "void" << endl
                 << unclash (name, "parse") << " (" <<
                parser_type << "&," << endl
                 << flags_type << ");"
                 << endl;
            }

            os << "protected:"
               << endl;

            // parse (istream)
            //
            if (has_members)
            {
              for (Streams::ConstIterator i (st.begin ()); i != st.end (); ++i)
              {
                os << "void" << endl
                   << unclash (name, "parse") << " (" <<
                  istream_type << "< " << i->c_str () << " >&," << endl
                   << flags_type << ");"
                   << endl;
              }
            }

            //
            //
            if (edom_document_member_p (c))
            {
              os << dom_auto_ptr << "< " << xerces_ns <<
                "::DOMDocument > " << edom_document_member (c) << ";"
                 << endl;
            }

            //
            //
            names (c, names_data_);

            if (doxygen)
              os << endl
                 << "//@endcond" << endl;
          }

          os << "};";

          // Comparison operators.
          //
          if (options.value<CLI::generate_comparison> () &&
              (has_members || !c.inherits_p () ||
               ((hae || haa) && gen_wildcard)))
          {
            os << inst_exp
               << "bool" << endl
               << "operator== (const " << name << "&, const " << name << "&);"
               << endl;

            os << inst_exp
               << "bool" << endl
               << "operator!= (const " << name << "&, const " << name << "&);"
               << endl
               << endl;
          }
        }

      private:
        Traversal::Inherits inherits_base_;
        BaseTypeName base_name_;

        Traversal::Inherits inherits_member_;
        MemberTypeName member_name_;

        Traversal::Names names_;
        Any any_;
        Member member_;

        Traversal::Names names_data_;
        DataAny data_any_;
        DataMember data_member_;
      };


      struct GlobalElement: Traversal::Element,
                            GlobalElementBase,
                            Context
      {
        GlobalElement (Context& c)
            : GlobalElementBase (c), Context (c), type_name_ (c)
        {
          belongs_ >> type_name_;
        }

        virtual Void
        traverse (Type& e)
        {
          if (!doc_root_p (e))
            return;

          SemanticGraph::Type& t (e.type ());

          Boolean fund (false);
          {
            IsFundamentalType test (fund);
            test.dispatch (t);
          }

          Boolean simple (true);
          if (!fund)
          {
            IsSimpleType test (simple);
            test.dispatch (t);
          }

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Class corresponding to the %" <<
              comment (e.name ()) << " root element." << endl;

            if (e.annotated_p ())
            {
              os << " *" << endl;
              write_annotation (e.annotation ());
            }

            os << " *" << endl
               << " * @nosubgrouping" << endl
               << " */" << endl;
          }

          String const& name (ename (e));

          os << "class " << type_exp << name << ": public " << element_type
             << "{"
             << "public:" << endl
             << endl;

          String const& type (etype (e));

          if (doxygen)
          {
            os << "/**" << endl
               << " * @name Element value" << endl
               << " *" << endl
               << " * @brief Accessor and modifier functions for the " <<
              "element value." << endl
               << " */" << endl
               << "//@{" << endl
               << endl;
          }
          else
          {
            os << "// Element value." << endl
               << "//" << endl;
          }

          // Typedefs.
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Element value type." << endl
               << " */" << endl;
          }

          os << "typedef ";

          belongs (e, belongs_);

          os << " " << type << ";";

          if (doxygen)
          {
            os << endl
               << "/**" << endl
               << " * @brief Element value traits type." << endl
               << " */" << endl;
          }

          os << "typedef ::xsd::cxx::tree::traits< " << type << ", " <<
            char_type;

          if (t.is_a<SemanticGraph::Fundamental::Double> ())
            os << ", ::xsd::cxx::tree::schema_type::double_";
          else if (t.is_a<SemanticGraph::Fundamental::Decimal> ())
            os << ", ::xsd::cxx::tree::schema_type::decimal";

          os << " > " << etraits (e) << ";"
             << endl;

          // Accessors/modifiers.
          //
          String const& aname (eaname (e));
          String const& mname (emname (e));

          // type const&
          // name () const;
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Return a read-only (constant) reference " <<
              "to the element" << endl
               << " * value." << endl
               << " *" << endl
               << " * @return A constant reference to the element value." <<
              endl
               << " */" << endl;
          }

          os << "const " << type << "&" << endl
             << aname << " () const;"
             << endl;

          // type&
          // name ();
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Return a read-write reference to the " <<
              "element value." << endl
               << " *" << endl
               << " * @return A reference to the element value." << endl
               << " */" << endl;
          }

          os << type << "&" << endl
             << aname << " ();"
             << endl;

          // void
          // name (type const&);
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Set the element value." << endl
               << " *" << endl
               << " * @param x A new value to set." << endl
               << " *" << endl
               << " * This function makes a copy of its argument " <<
              "and sets it as" << endl
               << " * the new value of the element." << endl
               << " */" << endl;
          }

          os << "void" << endl
             << mname << " (const " << type << "& x);"
             << endl;

          // void
          // name (auto_ptr<type>);
          //
          if (!fund)
          {
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Set the element value without " <<
                "copying." << endl
                 << " *" << endl
                 << " * @param p A new value to use." << endl
                 << " *" << endl
                 << " * This function will try to use the passed value " <<
                "directly" << endl
                 << " * instead of making a copy." << endl
                 << " */" << endl;
            }

            os << "void" << endl
               << mname << " (::std::auto_ptr< " << type << " > p);"
               << endl;
          }

          // auto_ptr<type>
          // detach_name ();
          //
          if (detach && !fund)
          {
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Detach the element value from " <<
                "the object." << endl
                 << " *" << endl
                 << " * @return A pointer to the element value." << endl
                 << " *" << endl
                 << " * Note that this function leaves the element " <<
                "object uninitialized." << endl
                 << " */" << endl;
            }

            os << "::std::auto_ptr< " << type << " >" << endl
               << edname (e) << " ();"
               << endl;
          }

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Return a read-only (constant) pointer " <<
              "to the element" << endl
               << " * value." << endl
               << " *" << endl
               << " * @return A constant pointer to the element value " <<
              "or 0 if this" << endl
               << " * element is of a fundamental type." << endl
               << " */" << endl;
          }

          os << "virtual const " << any_type << "*" << endl
             << "_value () const;"
             << endl;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Return a pointer to the element value." << endl
               << " *" << endl
               << " * @return A pointer to the element value or 0 if this " <<
              "element is" << endl
               << " * of a fundamental type." << endl
               << " */" << endl;
          }
          os << "virtual " << any_type << "*" << endl
             << "_value ();"
             << endl;

          if (doxygen)
          {
            os << "//@}" << endl
               << endl;
          }

          // Constructor.
          //

          if (doxygen)
          {
            os << "/**" << endl
               << " * @name Constructors" << endl
               << " */" << endl
               << "//@{" << endl
               << endl;
          }
          else
          {
            os << "// Constructors." << endl
               << "//" << endl;
          }

          // default c-tor
          //
          if (options.value<CLI::generate_default_ctor> ())
          {
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Default constructor." << endl
                 << " *" << endl
                 << " * Note that this constructor leaves the element " <<
                "value" << endl
                 << " * uninitialized." << endl
                 << " */" << endl;
            }

            os << name << " ();"
               << endl;
          }

          // c-tor (value)
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Create an instance from an initializer " <<
              "for the element" << endl
               << " * value." << endl
               << " *" << endl
               << " * @param x Element value." << endl
               << " */" << endl;
          }

          os << name << " (const " << type << "& x);"
             << endl;


          // If the element value is a complex type (has elements,
          // attributes, or wildcards) then also generate the auto_ptr
          // version. If we are generating polymorphic code then we
          // also need to provide auto_ptr version for simple types.
          //
          if (!simple || (polymorphic && polymorphic_p (t)))
          {
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from an initializer " <<
                "for" << endl
                 << " * the element value (auto_ptr version)." << endl
                 << " *" << endl
                 << " * @param p Element value to use." << endl
                 << " *" << endl
                 << " * This constructor will try to use the passed " <<
                "value directly" << endl
                 << " * instead of making a copy." << endl
                 << " */" << endl;
            }

            os << name << " (::std::auto_ptr< " << type << " > p);"
               << endl;
          }

          if (!options.value<CLI::suppress_parsing> ())
          {
            // c-tor (xercesc::DOMElement)
            //
            if (doxygen)
            {
              os << "/**" << endl
                 << " * @brief Create an instance from a DOM element." << endl
                 << " *" << endl
                 << " * @param e A DOM element to extract the data from." << endl
                 << " * @param f Flags to create the new instance with." << endl
                 << " */" << endl;
            }

            os << name << " (const " << xerces_ns << "::DOMElement& e, " <<
              flags_type << " f = 0);"
               << endl;
          }

          // copy c-tor
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Copy constructor." << endl
               << " *" << endl
               << " * @param x An instance to make a copy of." << endl
               << " * @param f Flags to create the copy with." << endl
               << " *" << endl
               << " * For polymorphic object models use the @c _clone " <<
              "function instead." << endl
               << " */" << endl;
          }

          os << name << " (const " << name << "& x, " <<
            flags_type << " f = 0);"
             << endl;

          // _clone
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Copy the instance polymorphically." << endl
               << " *" << endl
               << " * @param f Flags to create the copy with." << endl
               << " * @return A pointer to the dynamically allocated copy." << endl
               << " *" << endl
               << " * This function ensures that the dynamic type of the " <<
              "instance is" << endl
               << " * used for copying and should be used for polymorphic " <<
              "object" << endl
               << " * models instead of the copy constructor." << endl
               << " */" << endl;
          }

          os << "virtual " << name << "*" << endl
             << "_clone (" << flags_type << " f = 0) const;"
             << endl;

          if (doxygen)
          {
            os << "//@}" << endl
               << endl;
          }

          // Element name and namespace accessors.
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @name Element name and namespace" << endl
               << " *" << endl
               << " * @brief Accessor functions for the element name " <<
              "and namespace." << endl
               << " */" << endl
               << "//@{" << endl
               << endl;
          }
          else
          {
            os << "// Element name and namespace." << endl
               << "//" << endl;
          }

          SemanticGraph::Context& ec (e.context ());

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Return the element name (static function)." << endl
               << " *" << endl
               << " * @return A read-only string reference containing " <<
              "the element" << endl
               << " * name." << endl
               << " */" << endl;
          }
          os << "static const " << string_type << "&" << endl
             << ec.get<String> ("element-name") << " ();"
             << endl;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Return the element namespace (static " <<
              "function)." << endl
               << " *" << endl
               << " * @return A read-only string reference containing " <<
              "the element" << endl
               << " * namespace." << endl
               << " */" << endl;
          }
          os << "static const " << string_type << "&" << endl
             << ec.get<String> ("element-ns") << " ();"
             << endl;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Return the element name." << endl
               << " *" << endl
               << " * @return A read-only string reference containing " <<
              "the element" << endl
               << " * name." << endl
               << " */" << endl;
          }
          os << "virtual const " << string_type << "&" << endl
             << "_name () const;"
             << endl;

          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Return the element namespace." << endl
               << " *" << endl
               << " * @return A read-only string reference containing " <<
              "the element" << endl
               << " * namespace." << endl
               << " */" << endl;
          }
          os << "virtual const " << string_type << "&" << endl
             << "_namespace () const;"
             << endl;

          if (doxygen)
          {
            os << "//@}" << endl
               << endl;
          }

          // d-tor
          //
          if (doxygen)
          {
            os << "/**" << endl
               << " * @brief Destructor." << endl
               << " */" << endl;
          }

          os << "virtual " << endl
             << "~" << name << " ();"
             << endl;

          // Data member.
          //

          if (doxygen)
            os << "//@cond" << endl
               << endl;

          os << "protected:" << endl
             << "::xsd::cxx::tree::one< " << type << " > " <<
            emember (e) << ";"
             << "static const " << string_type << " " <<
            ec.get<String> ("element-name-member") << ";"
             << "static const " << string_type << " " <<
            ec.get<String> ("element-ns-member") << ";";

          if (doxygen)
            os << endl
               << "//@endcond" << endl;

          os << "};";
        }

      private:
        Traversal::Belongs belongs_;
        MemberTypeName type_name_;
      };
    }

    Void
    generate_tree_header (Context& ctx)
    {
      if (ctx.generate_xml_schema)
      {
        if (ctx.char_type == L"char" && ctx.char_encoding != L"custom")
        {
          ctx.os << "#include <xsd/cxx/xml/char-" << ctx.char_encoding <<
            ".hxx>" << endl
                 << endl;
        }

        ctx.os << "#include <xsd/cxx/tree/exceptions.hxx>" << endl
               << "#include <xsd/cxx/tree/elements.hxx>" << endl
               << "#include <xsd/cxx/tree/types.hxx>" << endl
               << endl;

        if (!ctx.options.value<CLI::suppress_parsing> () ||
            ctx.options.value<CLI::generate_serialization> ())
        {
          ctx.os << "#include <xsd/cxx/xml/error-handler.hxx>" << endl
                 << endl;
        }

        if (!ctx.options.value<CLI::suppress_parsing> () ||
            ctx.options.value<CLI::generate_serialization> ())
        {
          ctx.os << "#include <xsd/cxx/xml/dom/auto-ptr.hxx>" << endl
                 << endl;
        }

        Boolean element_map (ctx.options.value<CLI::generate_element_map> ());

        if (element_map)
          ctx.os << "#include <xsd/cxx/tree/element-map.hxx>" << endl
                 << endl;

        // I need to include all the "optional" headers here (instead of
        // later in the individual generators for each feature because
        // those headers provide implementation for the fundamental types.
        //
        if (!ctx.options.value<CLI::suppress_parsing> ())
        {
          ctx.os << "#include <xsd/cxx/tree/parsing.hxx>" << endl;

          Traversal::Schema schema;
          Traversal::Names names;
          Traversal::Namespace ns;
          Traversal::Names ns_names;
          FundIncludes type (ctx, "parsing");

          schema >> names >> ns >> ns_names >> type;

          schema.dispatch (ctx.schema_root);

          if (element_map)
            ctx.os << "#include <xsd/cxx/tree/parsing/element-map.txx>" <<
              endl;

          ctx.os << endl;
        }

        if (ctx.options.value<CLI::generate_serialization> ())
        {
          ctx.os << "#include <xsd/cxx/xml/dom/serialization-header.hxx>" << endl
                 << "#include <xsd/cxx/tree/serialization.hxx>" << endl;

          Traversal::Schema schema;
          Traversal::Names names;
          Traversal::Namespace ns;
          Traversal::Names ns_names;
          FundIncludes type (ctx, "serialization");

          schema >> names >> ns >> ns_names >> type;

          schema.dispatch (ctx.schema_root);

          if (element_map)
            ctx.os << "#include <xsd/cxx/tree/serialization/element-map.txx>" <<
              endl;

          ctx.os << endl;
        }

        if (ctx.options.value<CLI::generate_ostream> ())
        {
          ctx.os << "#include <xsd/cxx/tree/std-ostream-operators.hxx>" << endl
                 << endl;
        }

        Streams const& ist (ctx.options.value<CLI::generate_insertion> ());
        if (!ist.empty ())
        {
          for (Streams::ConstIterator i (ist.begin ()); i != ist.end (); ++i)
          {
            if (*i == "ACE_OutputCDR")
              ctx.os << "#include <xsd/cxx/tree/ace-cdr-stream-insertion.hxx>"
                     << endl;
            else if (*i == "XDR")
              ctx.os << "#include <xsd/cxx/tree/xdr-stream-insertion.hxx>"
                     << endl;
          }

          ctx.os << "#include <xsd/cxx/tree/stream-insertion.hxx>" << endl
                 << endl;
        }

        Streams const& est (ctx.options.value<CLI::generate_extraction> ());
        if (!est.empty ())
        {
          for (Streams::ConstIterator i (est.begin ()); i != est.end (); ++i)
          {
            if (*i == "ACE_InputCDR")
              ctx.os << "#include <xsd/cxx/tree/ace-cdr-stream-extraction.hxx>"
                     << endl;
            else if (*i == "XDR")
              ctx.os << "#include <xsd/cxx/tree/xdr-stream-extraction.hxx>"
                     << endl;
          }

          ctx.os << "#include <xsd/cxx/tree/stream-extraction.hxx>" << endl
                 << endl;
        }

        // Emit fundamental types.
        //
        {
          Traversal::Schema schema;
          Traversal::Names names;
          FundamentalNamespace ns (ctx);

          schema >> names >> ns;

          schema.dispatch (ctx.schema_root);
        }
      }
      else
      {
        Boolean inline_ (ctx.options.value<CLI::generate_inline> ());

        ctx.os << "#include <memory>    // std::auto_ptr" << endl
               << "#include <limits>    // std::numeric_limits" << endl
               << "#include <algorithm> // std::binary_search" << endl
               << endl;

        if (ctx.char_type == L"char" && ctx.char_encoding != L"custom")
        {
          ctx.os << "#include <xsd/cxx/xml/char-" << ctx.char_encoding <<
            ".hxx>" << endl
                 << endl;
        }

        ctx.os << "#include <xsd/cxx/tree/exceptions.hxx>" << endl
               << "#include <xsd/cxx/tree/elements.hxx>" << endl
               << "#include <xsd/cxx/tree/containers.hxx>" << endl
               << "#include <xsd/cxx/tree/list.hxx>" << endl
               << endl;

        if (!ctx.options.value<CLI::suppress_parsing> ())
        {
          ctx.os << "#include <xsd/cxx/xml/dom/parsing-header.hxx>" << endl
                 << endl;
        }

        if (ctx.options.value<CLI::generate_wildcard> ())
        {
          if (ctx.options.value<CLI::suppress_parsing> () ||
              !ctx.options.value<CLI::generate_serialization> ())
            ctx.os << "#include <xsd/cxx/xml/dom/auto-ptr.hxx>" << endl;

          ctx.os << "#include <xsd/cxx/tree/containers-wildcard.hxx>" << endl
                 << endl;
        }

        if (!ctx.options.value<CLI::generate_extraction> ().empty ())
          ctx.os << "#include <xsd/cxx/tree/istream-fwd.hxx>" << endl
                 << endl;

        // Emit header includes.
        //
        {
          if (inline_)
          {
            ctx.os << "#ifndef XSD_DONT_INCLUDE_INLINE" << endl
                   << "#define XSD_DONT_INCLUDE_INLINE" << endl
                   << endl;
          }

          Traversal::Schema schema;
          Includes includes (ctx, Includes::header);

          schema >> includes;

          schema.dispatch (ctx.schema_root);

          if (inline_)
          {
            ctx.os << "#undef XSD_DONT_INCLUDE_INLINE" << endl
                   << "#else" << endl
                   << endl;

            schema.dispatch (ctx.schema_root);

            ctx.os << "#endif // XSD_DONT_INCLUDE_INLINE" << endl
                   << endl;
          }
        }


        {
          Traversal::Schema schema;

          Traversal::Sources sources;
          Traversal::Names names_ns, names;

          DocumentedNamespace ns (ctx);

          List list (ctx);
          Union union_ (ctx);
          Complex complex (ctx);
          Enumeration enumeration (ctx);
          GlobalElement element (ctx);

          schema >> sources >> schema;
          schema >> names_ns >> ns >> names;

          names >> list;
          names >> union_;
          names >> complex;
          names >> enumeration;

          if (ctx.options.value<CLI::generate_element_type> ())
            names >> element;

          schema.dispatch (ctx.schema_root);
        }

        // Emit inline includes.
        //
        if (inline_)
        {
          ctx.os << "#ifndef XSD_DONT_INCLUDE_INLINE" << endl
                 << endl;

          Traversal::Schema schema;
          Includes ixx_includes (ctx, Includes::inline_);
          schema >> ixx_includes;

          schema.dispatch (ctx.schema_root);

          ctx.os << "#endif // XSD_DONT_INCLUDE_INLINE" << endl
                 << endl;
        }
      }
    }
  }
}
