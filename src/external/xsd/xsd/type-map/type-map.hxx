// file      : xsd/type-map/type-map.hxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2007-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

#ifndef XSD_TYPE_MAP_TYPE_MAP_HXX
#define XSD_TYPE_MAP_TYPE_MAP_HXX

#include <cult/types.hxx>
#include <cult/containers/vector.hxx>

#include <backend-elements/regex.hxx>

namespace TypeMap
{
  using namespace Cult::Types;
  typedef WideString String;
  typedef BackendElements::Regex::Pattern<WideChar> Pattern;

  class Type
  {
  public:
    Type (Pattern const& xsd_name,
          String const& cxx_ret_name,
          String const& cxx_arg_name)
        : xsd_name_ (xsd_name),
          cxx_ret_name_ (cxx_ret_name),
          cxx_arg_name_ (cxx_arg_name)
    {
    }

    Pattern const&
    xsd_name () const
    {
      return xsd_name_;
    }

    String const&
    cxx_ret_name () const
    {
      return cxx_ret_name_;
    }

    String const&
    cxx_arg_name () const
    {
      return cxx_arg_name_;
    }

  private:
    Pattern xsd_name_;
    String cxx_ret_name_;
    String cxx_arg_name_;
  };

  class Namespace
  {
  public:
    Namespace (Pattern const& xsd_name)
        : xsd_name_ (xsd_name), has_cxx_name_ (false)
    {
    }

    Namespace (Pattern const& xsd_name, String const& cxx_name)
        : xsd_name_ (xsd_name), has_cxx_name_ (true), cxx_name_ (cxx_name)
    {
    }

    Namespace (Pattern const& xsd_name,
               Boolean has_cxx_name,
               String const& cxx_name)
        : xsd_name_ (xsd_name),
          has_cxx_name_ (has_cxx_name),
          cxx_name_ (cxx_name)
    {
    }

    //
    //
    typedef Cult::Containers::Vector<String> Includes;
    typedef Includes::ConstIterator IncludesIterator;

    IncludesIterator
    includes_begin () const
    {
      return includes_.begin ();
    }

    IncludesIterator
    includes_end () const
    {
      return includes_.end ();
    }

    Void
    includes_push_back (String const& i)
    {
      includes_.push_back (i);
    }

    //
    //
    typedef Cult::Containers::Vector<Type> Types;
    typedef Types::ConstIterator TypesIterator;

    TypesIterator
    types_begin () const
    {
      return types_.begin ();
    }

    TypesIterator
    types_end () const
    {
      return types_.end ();
    }

    Void
    types_push_back (Pattern const& xsd_type,
                     String const& cxx_ret_type,
                     String const& cxx_arg_type = L"")
    {
      types_.push_back (Type (xsd_type, cxx_ret_type, cxx_arg_type));
    }

    //
    //
    Pattern const&
    xsd_name () const
    {
      return xsd_name_;
    }

    //
    //
    Boolean
    has_cxx_name () const
    {
      return has_cxx_name_;
    }

    String const&
    cxx_name () const
    {
      return cxx_name_;
    }

  private:
    Includes includes_;
    Types types_;
    Pattern xsd_name_;
    Boolean has_cxx_name_;
    String cxx_name_;
  };

  typedef Cult::Containers::Vector<Namespace> Namespaces;
}

#endif // XSD_TYPE_MAP_TYPE_MAP_HXX

