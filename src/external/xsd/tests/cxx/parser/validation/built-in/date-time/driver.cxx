// file      : tests/cxx/parser/validation/built-in/date-time/driver.cxx
// author    : Boris Kolpackov <boris@codesynthesis.com>
// copyright : Copyright (c) 2006-2010 Code Synthesis Tools CC
// license   : GNU GPL v2 + exceptions; see accompanying LICENSE file

// Test the built-in date and time types validation.
//
#include <cassert>

#include <xsd/cxx/parser/validating/exceptions.hxx>
#include <xsd/cxx/parser/validating/xml-schema-pimpl.hxx>

using namespace xsd::cxx::parser::validating;

template <typename T>
bool
test_post_fail (T& p)
{
  try
  {
    p._post ();
  }
  catch (invalid_value<char> const&)
  {
    return true;
  }

  return false;
}

int
main ()
{
  typedef xsd::cxx::parser::gday gday;
  typedef xsd::cxx::parser::gmonth gmonth;
  typedef xsd::cxx::parser::gyear gyear;
  typedef xsd::cxx::parser::gmonth_day gmonth_day;
  typedef xsd::cxx::parser::gyear_month gyear_month;
  typedef xsd::cxx::parser::date date;
  typedef xsd::cxx::parser::time time;
  typedef xsd::cxx::parser::date_time date_time;
  typedef xsd::cxx::parser::duration duration;

  // Good.
  //

  // gday & time zone parsing
  //
  {
    gday_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \t\n ");
    p._characters ("---1");
    p._characters ("2+12:00");
    p._post ();
    assert (p.post_gday () == gday (12, 12, 00));
  }

  {
    gday_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("---01");
    p._post ();
    assert (p.post_gday () == gday (1));
  }

  {
    gday_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("---31");
    p._post ();
    assert (p.post_gday () == gday (31));
  }

  {
    gday_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("---15Z");
    p._post ();
    assert (p.post_gday () == gday (15, 0, 0));
  }

  {
    gday_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("---15-14:00");
    p._post ();
    assert (p.post_gday () == gday (15, -14, 0));
  }

  // gmonth
  //
  {
    gmonth_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \t\n ");
    p._characters ("--1");
    p._characters ("0+12:00");
    p._post ();
    assert (p.post_gmonth () == gmonth (10, 12, 0));
  }

  {
    gmonth_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("--01");
    p._post ();
    assert (p.post_gmonth () == gmonth (1));
  }

  {
    gmonth_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("--12Z");
    p._post ();
    assert (p.post_gmonth () == gmonth (12, 0, 0));
  }

  // gyear
  //
  {
    gyear_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \t\n ");
    p._characters ("20");
    p._characters ("07+12:00");
    p._post ();
    assert (p.post_gyear () == gyear (2007, 12, 00));
  }

  {
    gyear_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("0001");
    p._post ();
    assert (p.post_gyear () == gyear (1));
  }

  {
    gyear_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-20000Z");
    p._post ();
    assert (p.post_gyear () == gyear (-20000, 0, 0));
  }

  // gmonth_day
  //
  {
    gmonth_day_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \t\n ");
    p._characters ("--1");
    p._characters ("0-28+12:00  ");
    p._post ();
    assert (p.post_gmonth_day () == gmonth_day (10, 28, 12, 00));
  }

  {
    gmonth_day_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("--12-31");
    p._post ();
    assert (p.post_gmonth_day () == gmonth_day (12, 31));
  }

  {
    gmonth_day_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("--01-01Z");
    p._post ();
    assert (p.post_gmonth_day () == gmonth_day (1, 1, 0, 0));
  }

  // gyear_month
  //
  {
    gyear_month_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \t\n ");
    p._characters ("200");
    p._characters ("7-12+12:00  ");
    p._post ();
    assert (p.post_gyear_month () == gyear_month (2007, 12, 12, 00));
  }

  {
    gyear_month_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-2007-10");
    p._post ();
    assert (p.post_gyear_month () == gyear_month (-2007, 10));
  }

  {
    gyear_month_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("20007-10Z");
    p._post ();
    assert (p.post_gyear_month () == gyear_month (20007, 10, 0, 0));
  }

  {
    gyear_month_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-20007-01");
    p._post ();
    assert (p.post_gyear_month () == gyear_month (-20007, 1));
  }

  // date
  //
  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \t\n ");
    p._characters ("200");
    p._characters ("7-12-26+12:00  ");
    p._post ();
    assert (p.post_date () == date (2007, 12, 26, 12, 0));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-2007-10-15");
    p._post ();
    assert (p.post_date () == date (-2007, 10, 15));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("20007-12-31Z");
    p._post ();
    assert (p.post_date () == date (20007, 12, 31, 0, 0));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-20007-01-01");
    p._post ();
    assert (p.post_date () == date (-20007, 1, 1));
  }

  // time
  //
  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \t\n ");
    p._characters ("12:");
    p._characters ("46:23.456+12:00  ");
    p._post ();
    assert (p.post_time () == time (12, 46, 23.456, 12, 0));
  }

  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("12:13:14");
    p._post ();
    assert (p.post_time () == time (12, 13, 14.0));
  }

  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("12:13:14Z");
    p._post ();
    assert (p.post_time () == time (12, 13, 14.0, 0, 0));
  }

  // date_time
  //
  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \t\n ");
    p._characters ("200");
    p._characters ("7-12-26T12:13:14.123+12:00  ");
    p._post ();
    assert (p.post_date_time () ==
            date_time (2007, 12, 26, 12, 13, 14.123, 12, 0));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-2007-10-15T12:13:14");
    p._post ();
    assert (p.post_date_time () == date_time (-2007, 10, 15, 12, 13, 14.0));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("20007-12-31T12:13:14Z");
    p._post ();
    assert (p.post_date_time () ==
            date_time (20007, 12, 31, 12, 13, 14.0, 0, 0));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-20007-01-01T12:13:14");
    p._post ();
    assert (p.post_date_time () == date_time (-20007, 1, 1, 12, 13, 14.0));
  }

  // duration
  //
  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters (" \t\n ");
    p._characters ("-P200");
    p._characters ("7Y13M32DT25H61M61.123S  ");
    p._post ();
    assert (p.post_duration () ==
            duration (true, 2007, 13, 32, 25, 61, 61.123));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("P1Y");
    p._post ();
    assert (p.post_duration () == duration (false, 1, 0, 0, 0, 0, 0.0));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("P1M");
    p._post ();
    assert (p.post_duration () == duration (false, 0, 1, 0, 0, 0, 0.0));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("P1D");
    p._post ();
    assert (p.post_duration () == duration (false, 0, 0, 1, 0, 0, 0.0));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("PT1H");
    p._post ();
    assert (p.post_duration () == duration (false, 0, 0, 0, 1, 0, 0.0));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("PT1M");
    p._post ();
    assert (p.post_duration () == duration (false, 0, 0, 0, 0, 1, 0.0));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("PT1.1S");
    p._post ();
    assert (p.post_duration () == duration (false, 0, 0, 0, 0, 0, 1.1));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("P1YT1S");
    p._post ();
    assert (p.post_duration () == duration (false, 1, 0, 0, 0, 0, 1.0));
  }

  // Bad
  //

  // gday & time zone parsing
  //
  {
    gday_pimpl<char> p;
    p.pre ();
    p._pre ();
    // p._characters ("");
    assert (test_post_fail (p));
  }

  {
    gday_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("--12");
    assert (test_post_fail (p));
  }

  {
    gday_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("---1");
    assert (test_post_fail (p));
  }

  {
    gday_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("---00");
    assert (test_post_fail (p));
  }

  {
    gday_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("---32");
    assert (test_post_fail (p));
  }

  {
    gday_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("---2X");
    assert (test_post_fail (p));
  }

  {
    gday_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("---12asd");
    assert (test_post_fail (p));
  }

  {
    gday_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("---12X");
    assert (test_post_fail (p));
  }

  {
    gday_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("---1212:00");
    assert (test_post_fail (p));
  }

  {
    gday_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("---12+2:00");
    assert (test_post_fail (p));
  }

  {
    gday_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("---12+1200");
    assert (test_post_fail (p));
  }

  {
    gday_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("---12+15:00");
    assert (test_post_fail (p));
  }

  {
    gday_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("---12+12:60");
    assert (test_post_fail (p));
  }

  {
    gday_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("---12+14:01");
    assert (test_post_fail (p));
  }

  // gmonth
  //
  {
    gmonth_pimpl<char> p;
    p.pre ();
    p._pre ();
    // p._characters ("");
    assert (test_post_fail (p));
  }

  {
    gmonth_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-12");
    assert (test_post_fail (p));
  }

  {
    gmonth_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("--00");
    assert (test_post_fail (p));
  }

  {
    gmonth_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("--13");
    assert (test_post_fail (p));
  }

  {
    gmonth_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("--1X");
    assert (test_post_fail (p));
  }

  {
    gmonth_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("--11+12:3o");
    assert (test_post_fail (p));
  }

  // gyear
  //
  {
    gyear_pimpl<char> p;
    p.pre ();
    p._pre ();
    // p._characters ("");
    assert (test_post_fail (p));
  }

  {
    gyear_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("207");
    assert (test_post_fail (p));
  }

  {
    gyear_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-207");
    assert (test_post_fail (p));
  }

  {
    gyear_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-0000");
    assert (test_post_fail (p));
  }

  {
    gyear_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("20X7");
    assert (test_post_fail (p));
  }

  {
    gyear_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007+12:3o");
    assert (test_post_fail (p));
  }

  // gmonth_day
  //
  {
    gmonth_day_pimpl<char> p;
    p.pre ();
    p._pre ();
    // p._characters ("");
    assert (test_post_fail (p));
  }

  {
    gmonth_day_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-12-12");
    assert (test_post_fail (p));
  }

  {
    gmonth_day_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("--1212");
    assert (test_post_fail (p));
  }

  {
    gmonth_day_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("--12?12");
    assert (test_post_fail (p));
  }

  {
    gmonth_day_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("--00-12");
    assert (test_post_fail (p));
  }

  {
    gmonth_day_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("--12-00");
    assert (test_post_fail (p));
  }

  {
    gmonth_day_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("--13-23");
    assert (test_post_fail (p));
  }

  {
    gmonth_day_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("--12-32");
    assert (test_post_fail (p));
  }

  {
    gmonth_day_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("--1X-12");
    assert (test_post_fail (p));
  }

  {
    gmonth_day_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("--12-2X");
    assert (test_post_fail (p));
  }

  {
    gmonth_day_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("--11-11+12:3o");
    assert (test_post_fail (p));
  }

  // gyear_month
  //
  {
    gyear_month_pimpl<char> p;
    p.pre ();
    p._pre ();
    // p._characters ("");
    assert (test_post_fail (p));
  }

  {
    gyear_month_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("207-01");
    assert (test_post_fail (p));
  }

  {
    gyear_month_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-207-01");
    assert (test_post_fail (p));
  }

  {
    gyear_month_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("0000-01");
    assert (test_post_fail (p));
  }

  {
    gyear_month_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("20X7-01");
    assert (test_post_fail (p));
  }

  {
    gyear_month_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007");
    assert (test_post_fail (p));
  }

  {
    gyear_month_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007?12");
    assert (test_post_fail (p));
  }

  {
    gyear_month_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-0");
    assert (test_post_fail (p));
  }

  {
    gyear_month_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-00");
    assert (test_post_fail (p));
  }

  {
    gyear_month_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-13");
    assert (test_post_fail (p));
  }

  {
    gyear_month_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-1X");
    assert (test_post_fail (p));
  }

  {
    gyear_month_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-01+12:3o");
    assert (test_post_fail (p));
  }

  // date
  //
  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    // p._characters ("");
    assert (test_post_fail (p));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("207-01-01");
    assert (test_post_fail (p));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-207-01-01");
    assert (test_post_fail (p));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("0000-01-01");
    assert (test_post_fail (p));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("20X7-01-01");
    assert (test_post_fail (p));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007");
    assert (test_post_fail (p));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007?01-01");
    assert (test_post_fail (p));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-0-01");
    assert (test_post_fail (p));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-00-01");
    assert (test_post_fail (p));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-13-01");
    assert (test_post_fail (p));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-1X-01");
    assert (test_post_fail (p));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-10");
    assert (test_post_fail (p));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-10?12");
    assert (test_post_fail (p));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-10-");
    assert (test_post_fail (p));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-10-0");
    assert (test_post_fail (p));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-10-00");
    assert (test_post_fail (p));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-10-32");
    assert (test_post_fail (p));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-10-2X");
    assert (test_post_fail (p));
  }

  {
    date_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-01-01+12:3o");
    assert (test_post_fail (p));
  }

  // time
  //
  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    // p._characters ("");
    assert (test_post_fail (p));
  }

  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("1:01:01");
    assert (test_post_fail (p));
  }

  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2X:01:01");
    assert (test_post_fail (p));
  }

  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("23");
    assert (test_post_fail (p));
  }

  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("23?01:01");
    assert (test_post_fail (p));
  }

  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("23:0:01");
    assert (test_post_fail (p));
  }

  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("23:60:01");
    assert (test_post_fail (p));
  }

  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("23:4X:01");
    assert (test_post_fail (p));
  }

  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("23:10");
    assert (test_post_fail (p));
  }

  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("23:10?12");
    assert (test_post_fail (p));
  }

  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("23:10:");
    assert (test_post_fail (p));
  }

  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("23:10:0");
    assert (test_post_fail (p));
  }

  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("23:10:01.");
    assert (test_post_fail (p));
  }

  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("23:10:60");
    assert (test_post_fail (p));
  }

  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("23:10:2X");
    assert (test_post_fail (p));
  }

  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("24:01:00");
    assert (test_post_fail (p));
  }

  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("24:00:01");
    assert (test_post_fail (p));
  }

  {
    time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("23:01:01+12:3o");
    assert (test_post_fail (p));
  }

  // date_time
  //
  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    // p._characters ("");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("207-01-01T12:13:14");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-207-01-01T12:13:14");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("0000-01-01T12:13:14");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("20X7-01-01T12:13:14");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007?01-01T12:13:14");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-0-01T12:13:14");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-00-01T12:13:14");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-13-01T12:13:14");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-1X-01T12:13:14");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-10");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-10?12T12:13:14");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-10-T12:13:14");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-10-0T12:13:14");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-10-00T12:13:14");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-10-32T12:13:14");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-10-2XT12:13:14");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-01-01T1:01:01");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-01-01T2X:01:01");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-01-01T23");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-01-01T23?01:01");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-01-01T23:0:01");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-01-01T23:60:01");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-01-01T23:4X:01");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-01-01T23:10");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-01-01T23:10?12");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-01-01T23:10:");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-01-01T23:10:0");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-01-01T23:10:01.");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-01-01T23:10:60");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-01-01T23:10:2X");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-01-01T24:01:00");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-01-01T24:00:01");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("23:01:01+12:3o");
    assert (test_post_fail (p));
  }

  {
    date_time_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007-01-01T12:13:14+12:3o");
    assert (test_post_fail (p));
  }

  // duration
  //
  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    // p._characters ("");
    assert (test_post_fail (p));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("2007Y");
    assert (test_post_fail (p));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("-2007Y");
    assert (test_post_fail (p));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("P-2007Y");
    assert (test_post_fail (p));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("P-1M");
    assert (test_post_fail (p));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("P-1D");
    assert (test_post_fail (p));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("PT-1H");
    assert (test_post_fail (p));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("PT-1M");
    assert (test_post_fail (p));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("PT-1.1S");
    assert (test_post_fail (p));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("P1H1M1S");
    assert (test_post_fail (p));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("P1M1Y");
    assert (test_post_fail (p));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("PT1S1H");
    assert (test_post_fail (p));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("PT1H1Y");
    assert (test_post_fail (p));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("P1Ygarbage");
    assert (test_post_fail (p));
  }

  {
    duration_pimpl<char> p;
    p.pre ();
    p._pre ();
    p._characters ("P1YT");
    assert (test_post_fail (p));
  }
}
