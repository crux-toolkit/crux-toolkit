// file      : xsd/cxx/parser/print-impl-common.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2005-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef CXX_PARSER_PRINT_IMPL_COMMON_HXX
#define CXX_PARSER_PRINT_IMPL_COMMON_HXX

#include <xsd-frontend/semantic-graph.hxx>
#include <xsd-frontend/traversal.hxx>

#include <cxx/parser/elements.hxx>

namespace CXX
{
  namespace Parser
  {
    struct PrintCall: Traversal::Type,

                      Traversal::Fundamental::Boolean,

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

                      Context
    {
      PrintCall (Context& c, String const& tag, String const& arg)
          : Context (c), tag_ (tag), arg_ (arg)
      {
      }

      virtual Void
      traverse (SemanticGraph::Type&)
      {
        gen_user_type ();
      }

      // Boolean.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Boolean& t)
      {
        if (default_type (t, "bool"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << " << " <<
            arg_ << " << std::endl;";
        }
        else
          gen_user_type ();
      }

      // Integral types.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Byte& t)
      {
        if (default_type (t, "signed char"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") <<
            " << static_cast<short> (" << arg_ << ") << std::endl;";
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedByte& t)
      {
        if (default_type (t, "unsigned char"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") <<
            " << static_cast<unsigned short> (" << arg_ << ") << std::endl;";
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Short& t)
      {
        if (default_type (t, "short"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << " << " <<
            arg_ << " << std::endl;";
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedShort& t)
      {
        if (default_type (t, "unsigned short"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << " << " <<
            arg_ << " << std::endl;";
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Int& t)
      {
        if (default_type (t, "int"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << " << " <<
            arg_ << " << std::endl;";
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedInt& t)
      {
        if (default_type (t, "unsigned int"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << " << " <<
            arg_ << " << std::endl;";
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Long& t)
      {
        if (default_type (t, "long long"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << " << " <<
            arg_ << " << std::endl;";
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::UnsignedLong& t)
      {
        if (default_type (t, "unsigned long long"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << " << " <<
            arg_ << " << std::endl;";
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Integer& t)
      {
        if (default_type (t, "long long"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << " << " <<
            arg_ << " << std::endl;";
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NegativeInteger& t)
      {
        if (default_type (t, "long long"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << " << " <<
            arg_ << " << std::endl;";
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NonPositiveInteger& t)
      {
        if (default_type (t, "long long"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << " << " <<
            arg_ << " << std::endl;";
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::PositiveInteger& t)
      {
        if (default_type (t, "unsigned long long"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << " << " <<
            arg_ << " << std::endl;";
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NonNegativeInteger& t)
      {
        if (default_type (t, "unsigned long long"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << " << " <<
            arg_ << " << std::endl;";
        }
        else
          gen_user_type ();
      }

      // Floats.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Float& t)
      {
        if (default_type (t, "float"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << " << " <<
            arg_ << " << std::endl;";
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Double& t)
      {
        if (default_type (t, "double"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << " << " <<
            arg_ << " << std::endl;";
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Decimal& t)
      {
        if (default_type (t, "double"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << " << " <<
            arg_ << " << std::endl;";
        }
        else
          gen_user_type ();
      }

      // Strings.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::String& t)
      {
        gen_string (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NormalizedString& t)
      {
        gen_string (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Token& t)
      {
        gen_string (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NameToken& t)
      {
        gen_string (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Name& t)
      {
        gen_string (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::NCName& t)
      {
        gen_string (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Language& t)
      {
        gen_string (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Id& t)
      {
        gen_string (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::IdRef& t)
      {
        gen_string (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::AnyURI& t)
      {
        gen_string (t);
      }

      // String sequences.
      //

      virtual Void
      traverse (SemanticGraph::Fundamental::NameTokens& t)
      {
        gen_sequence (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::IdRefs& t)
      {
        gen_sequence (t);
      }

      // QName
      //

      virtual Void
      traverse (SemanticGraph::Fundamental::QName& t)
      {
        if (default_type (t, xs_ns_name () + L"::qname"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << ";"
             << endl
             << "if (" << arg_ << ".prefix ().empty ())" << endl
             << cout_inst << " << " << arg_ << ".name ();"
             << "else" << endl
             << cout_inst << " << " << arg_ << ".prefix () << " << L <<
            "':' << " << arg_ << ".name ();"
             << endl
             << cout_inst << " << std::endl;";
        }
        else
          gen_user_type ();
      }

      // Binary.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Base64Binary& t)
      {
        gen_buffer (t);
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::HexBinary& t)
      {
        gen_buffer (t);
      }

      // Date/time.
      //
      virtual Void
      traverse (SemanticGraph::Fundamental::Date& t)
      {
        if (default_type (t, xs_ns_name () + L"::date"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << endl
             << " << " << arg_ << ".year () << '-'" << endl
             << " << " << arg_ << ".month () << '-'" << endl
             << " << " << arg_ << ".day ();";

          gen_time_zone ();
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::DateTime& t)
      {
        if (default_type (t, xs_ns_name () + L"::date_time"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << endl
             << " << " << arg_ << ".year () << '-'" << endl
             << " << " << arg_ << ".month () << '-'" << endl
             << " << " << arg_ << ".day () << 'T'" << endl
             << " << " << arg_ << ".hours () << ':'" << endl
             << " << " << arg_ << ".minutes () << ':'" << endl
             << " << " << arg_ << ".seconds ();";

          gen_time_zone ();
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Duration& t)
      {
        if (default_type (t, xs_ns_name () + L"::duration"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << ";"
             << endl
             << "if (" << arg_ << ".negative ())" << endl
             << cout_inst << " << '-';"
             << endl
             << cout_inst << " << 'P'" << endl
             << " << " << arg_ << ".years () << 'Y'" << endl
             << " << " << arg_ << ".months () << 'M'" << endl
             << " << " << arg_ << ".days () << " << L << "\"DT\"" << endl
             << " << " << arg_ << ".hours () << 'H'" << endl
             << " << " << arg_ << ".minutes () << 'M'" << endl
             << " << " << arg_ << ".seconds () << 'S'"
             << " << std::endl;";
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Day& t)
      {
        if (default_type (t, xs_ns_name () + L"::gday"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ---") <<
            " << " << arg_ << ".day ();";

          gen_time_zone ();
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Month& t)
      {
        if (default_type (t, xs_ns_name () + L"::gmonth"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": --") <<
            " << " << arg_ << ".month ();";

          gen_time_zone ();
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::MonthDay& t)
      {
        if (default_type (t, xs_ns_name () + L"::gmonth_day"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": --") << endl
             << " << " << arg_ << ".month () << '-'" << endl
             << " << " << arg_ << ".day ();";

          gen_time_zone ();
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Year& t)
      {
        if (default_type (t, xs_ns_name () + L"::gyear"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << " << " <<
            arg_ << ".year ();";

          gen_time_zone ();
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::YearMonth& t)
      {
        if (default_type (t, xs_ns_name () + L"::gyear_month"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << endl
             << " << " << arg_ << ".year () << '-'" << endl
             << " << " << arg_ << ".month ();";

          gen_time_zone ();
        }
        else
          gen_user_type ();
      }

      virtual Void
      traverse (SemanticGraph::Fundamental::Time& t)
      {
        if (default_type (t, xs_ns_name () + L"::time"))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << endl
             << " << " << arg_ << ".hours () << ':'" << endl
             << " << " << arg_ << ".minutes () << ':'" << endl
             << " << " << arg_ << ".seconds ();";

          gen_time_zone ();
        }
        else
          gen_user_type ();
      }

    private:
      bool
      default_type (SemanticGraph::Type& t, String const& def_type)
      {
        return ret_type (t) == def_type;
      }

      void
      gen_user_type ()
      {
        os << "// TODO" << endl
           << "//" << endl;
      }

      void
      gen_string (SemanticGraph::Type& t)
      {
        if ((char_type == L"char" && default_type (t, "::std::string")) ||
            (char_type == L"wchar_t" && default_type (t, "::std::wstring")))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << " << " <<
            arg_ << " << std::endl;";
        }
        else
          gen_user_type ();
      }

      void
      gen_sequence (SemanticGraph::Type& t)
      {
        String type (xs_ns_name () + L"::string_sequence");

        if (default_type (t, type))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << ";"
             << endl;

          os << "for (" << type << "::const_iterator i (" << arg_ <<
            ".begin ()), e (" << arg_ << ".end ());" << endl
             << "i != e;)"
             << "{"
             << cout_inst << " << *i++;"
             << "if (i != e)" << endl
             << cout_inst << " << ' ';"
             << "}"
             << cout_inst << " << std::endl;";
        }
        else
          gen_user_type ();
      }

      void
      gen_buffer (SemanticGraph::Type& t)
      {
        String type (L"::std::auto_ptr< " + xs_ns_name () + L"::buffer >");

        if (default_type (t, type))
        {
          os << cout_inst << " << " << strlit (tag_ + L": ") << " << "
             << arg_ << "->size () << " << L << "\" bytes\" << std::endl;";
        }
        else
          gen_user_type ();
      }

      void
      gen_time_zone ()
      {
        os << endl
           << "if (" << arg_ << ".zone_present ())"
           << "{"
           << "if (" << arg_ << ".zone_hours () < 0)" << endl
           << cout_inst << " << " << arg_ << ".zone_hours () << ':' << -" <<
          arg_ << ".zone_minutes ();"
           << "else" << endl
           << cout_inst << " << '+' << " << arg_ << ".zone_hours () << " <<
          "':' << " << arg_ << ".zone_minutes ();";

        os << "}"
           << cout_inst << " << std::endl;";
      }

    private:
      String tag_;
      String arg_;
    };
  }
}

#endif // CXX_PARSER_PRINT_IMPL_COMMON_HXX
