// file      : tests/cxx/tree/binary/xdr/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test non-polymorphic binary serialization to XDR.
//

#include <memory>  // std::auto_ptr
#include <cstring> // std::memcpy
#include <cassert>
#include <iostream>

#include "test.hxx"

using namespace std;
using namespace test;

extern "C" int
overflow (char* p, char* buf, int n)
{
  xml_schema::buffer* dst (reinterpret_cast<xml_schema::buffer*> (p));

  std::size_t size (dst->size ());
  dst->size (size + n);
  memcpy (dst->data () + size, buf, n);

  return n;
}

struct underflow_info
{
  xml_schema::buffer* buf;
  std::size_t pos;
};

extern "C" int
underflow (char* p, char* buf, int n)
{
  underflow_info* ui (reinterpret_cast<underflow_info*> (p));

  std::size_t size (ui->buf->size () - ui->pos);
  n = size > n ? n : size;

  memcpy (buf, ui->buf->data () + ui->pos, n);
  ui->pos += n;

  return n;
}

int
main (int argc, char* argv[])
{
  if (argc != 2)
  {
    cerr << "usage: " << argv[0] << " test.xml" << endl;
    return 1;
  }

  try
  {
    auto_ptr<type> r (root (argv[1]));

    // Save to an XDR stream.
    //
    XDR xdr;
    xml_schema::buffer buf;
    xdrrec_create (&xdr, 0, 0, reinterpret_cast<char*> (&buf), 0, &overflow);
    xdr.x_op = XDR_ENCODE;
    xsd::cxx::tree::ostream<XDR> oxdr (xdr);
    oxdr << *r;
    xdrrec_endofrecord (&xdr, true); // flush the data
    xdr_destroy (&xdr);

    // Load from an XDR stream.
    //
    underflow_info ui;
    ui.buf = &buf;
    ui.pos = 0;
    xdrrec_create (&xdr, 0, 0, reinterpret_cast<char*> (&ui), &underflow, 0);
    xdr.x_op = XDR_DECODE;
    xdrrec_skiprecord (&xdr);
    xsd::cxx::tree::istream<XDR> ixdr (xdr);
    auto_ptr<type> c (new type (ixdr));
    xdr_destroy (&xdr);

    // Compare the two.
    //
    assert (r->list () == c->list ());
    assert (r->union_ () == c->union_ ());
    assert (r->enumeration () == c->enumeration ());

    type::complex_sequence rs (r->complex ()), cs (c->complex ());

    for (type::complex_iterator ri (rs.begin ()), ci (cs.begin ());
         ri != rs.end () && ci != rs.end (); ++ri, ++ci)
    {
      assert (ri->a () == ci->a ());
      if (ri->b ())
        assert (ri->b () == ci->b ());
      assert (ri->c () == ci->c ());

      assert (ri->x () == ci->x ());
      if (ri->y ())
        assert (ri->y () == ci->y ());
    }

    // integers
    //
    assert (r->byte () == c->byte ());
    assert (r->unsigned_byte () == c->unsigned_byte ());
    assert (r->short_ () == c->short_ ());
    assert (r->unsigned_short () == c->unsigned_short ());
    assert (r->int_ () == c->int_ ());
    assert (r->unsigned_int () == c->unsigned_int ());
    assert (r->long_ () == c->long_ ());
    assert (r->unsigned_long () == c->unsigned_long ());
    assert (r->integer () == c->integer ());
    assert (r->non_positive_integer () == c->non_positive_integer ());
    assert (r->non_negative_integer () == c->non_negative_integer ());
    assert (r->positive_integer () == c->positive_integer ());
    assert (r->negative_integer () == c->negative_integer ());

    // boolean
    //
    assert (r->boolean () == c->boolean ());

    // floats
    //
    assert (r->float_ () == c->float_ ());
    assert (r->double_ () == c->double_ ());
    assert (r->decimal () == c->decimal ());

    // strings
    //
    assert (r->string () == c->string ());
    assert (r->normalized_string () == c->normalized_string ());
    assert (r->token () == c->token ());
    assert (r->name () == c->name ());
    assert (r->name_token () == c->name_token ());
    assert (r->name_tokens () == c->name_tokens ());
    assert (r->ncname () == c->ncname ());
    assert (r->language () == c->language ());

    // qualified name
    //
    assert (r->qname () == c->qname ());

    // ID/IDREF
    //
    assert (r->id () == c->id ());
    assert (r->id_ref () == c->id_ref ());
    assert (r->id_refs () == c->id_refs ());

    // URI
    //
    assert (r->any_uri () == c->any_uri ());

    // binary
    //
    assert (r->base64_binary () == c->base64_binary ());
    assert (r->hex_binary () == c->hex_binary ());

    // date/time
    //
    assert (r->date () == c->date ());
    assert (r->date_time () == c->date_time ());
    assert (r->duration () == c->duration ());
    assert (r->day () == c->day ());
    assert (r->month () == c->month ());
    assert (r->month_day () == c->month_day ());
    assert (r->year () == c->year ());
    assert (r->year_month () == c->year_month ());
    assert (r->time () == c->time ());
  }
  catch (xml_schema::exception const& e)
  {
    cerr << e << endl;
    return 1;
  }
}
