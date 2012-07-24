// file      : tests/cxx/tree/compilation/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Make sure the runtime library compiles by explicitly instantiating
// all the types.
//

#include <memory> // std::auto_ptr
#include <iostream>

#include "test.hxx"

using namespace std;
using namespace test;

template class xsd::cxx::tree::simple_type<xml_schema::type>;

// String types.
//
template class xsd::cxx::tree::string< char, xml_schema::simple_type >;
template class xsd::cxx::tree::normalized_string< char, xml_schema::string >;
template class xsd::cxx::tree::token< char, xml_schema::normalized_string >;
template class xsd::cxx::tree::name< char, xml_schema::token >;
template class xsd::cxx::tree::nmtoken< char, xml_schema::token >;
template class xsd::cxx::tree::nmtokens< char, xml_schema::simple_type, xml_schema::nmtoken >;
template class xsd::cxx::tree::ncname< char, xml_schema::name >;
template class xsd::cxx::tree::language< char, xml_schema::token >;

// ID/IDREF.
//
template class xsd::cxx::tree::id< char, xml_schema::ncname >;
template class xsd::cxx::tree::idref< char, xml_schema::ncname, xml_schema::type >;
template class xsd::cxx::tree::idrefs< char, xml_schema::simple_type, xml_schema::idref >;

// URI.
//
template class xsd::cxx::tree::uri< char, xml_schema::simple_type >;

// Qualified name.
//
template class xsd::cxx::tree::qname< char, xml_schema::simple_type, xml_schema::uri, xml_schema::ncname >;

// Binary.
//
template class xsd::cxx::tree::buffer< char >;
template class xsd::cxx::tree::base64_binary< char, xml_schema::simple_type >;
template class xsd::cxx::tree::hex_binary< char, xml_schema::simple_type >;

// Date/time.
//
template class xsd::cxx::tree::date< char, xml_schema::simple_type >;
template class xsd::cxx::tree::date_time< char, xml_schema::simple_type >;
template class xsd::cxx::tree::duration< char, xml_schema::simple_type >;
template class xsd::cxx::tree::gday< char, xml_schema::simple_type >;
template class xsd::cxx::tree::gmonth< char, xml_schema::simple_type >;
template class xsd::cxx::tree::gmonth_day< char, xml_schema::simple_type >;
template class xsd::cxx::tree::gyear< char, xml_schema::simple_type >;
template class xsd::cxx::tree::gyear_month< char, xml_schema::simple_type >;
template class xsd::cxx::tree::time< char, xml_schema::simple_type >;

// Entity.
//
template class xsd::cxx::tree::entity< char, xml_schema::ncname >;
template class xsd::cxx::tree::entities< char, xml_schema::simple_type, xml_schema::entity >;

// Namespace information and list stream. Used in
// serialization functions.
//
template class xsd::cxx::xml::dom::namespace_info < char >;
template class xsd::cxx::xml::dom::namespace_infomap < char >;
template class xsd::cxx::tree::list_stream < char >;

// Flags and properties.
//
template class xsd::cxx::tree::properties< char >;

// Exceptions.
//
template class xsd::cxx::tree::exception< char >;
template class xsd::cxx::tree::parsing< char >;
template class xsd::cxx::tree::expected_element< char >;
template class xsd::cxx::tree::unexpected_element< char >;
template class xsd::cxx::tree::expected_attribute< char >;
template class xsd::cxx::tree::unexpected_enumerator< char >;
template class xsd::cxx::tree::expected_text_content< char >;
template class xsd::cxx::tree::no_type_info< char >;
template class xsd::cxx::tree::not_derived< char >;
template class xsd::cxx::tree::duplicate_id< char >;
template class xsd::cxx::tree::serialization< char >;
template class xsd::cxx::tree::no_prefix_mapping< char >;
template class xsd::cxx::tree::bounds< char >;

// Parsing/serialization diagnostics.
//
template class xsd::cxx::tree::error< char >;
template class xsd::cxx::tree::diagnostics< char >;

// Error handler interface.
//
template class xsd::cxx::xml::error_handler< char >;


//
//
template class xsd::cxx::tree::fundamental_base<int, char, xml_schema::type>;
template class xsd::cxx::tree::one<int>;
template class xsd::cxx::tree::one<xml_schema::string>;
template class xsd::cxx::tree::optional<int>;
template class xsd::cxx::tree::optional<xml_schema::string>;
template class xsd::cxx::tree::sequence<int>;
template class xsd::cxx::tree::sequence<xml_schema::string>;


int
main ()
{
}
