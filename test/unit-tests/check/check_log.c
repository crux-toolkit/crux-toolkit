/*
 * Check: a unit test framework for C
 * Copyright (C) 2001, 2002 Arien Malec
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <check.h>

#include "check_error.h"
#include "check_list.h"
#include "check_impl.h"
#include "check_log.h"
#include "check_print.h"


static void srunner_send_evt (SRunner *sr, void *obj, enum cl_event evt);

void srunner_set_log (SRunner *sr, const char *fname)
{
  if (sr->log_fname)
    return;
  sr->log_fname = fname;
}

int srunner_has_log (SRunner *sr)
{
  return sr->log_fname != NULL;
}

const char *srunner_log_fname (SRunner *sr)
{
  return sr->log_fname;
}


void srunner_set_xml (SRunner *sr, const char *fname)
{
  if (sr->xml_fname)
    return;
  sr->xml_fname = fname;
}

int srunner_has_xml (SRunner *sr)
{
  return sr->xml_fname != NULL;
}

const char *srunner_xml_fname (SRunner *sr)
{
  return sr->xml_fname;
}

void srunner_register_lfun (SRunner *sr, FILE *lfile, int close,
			    LFun lfun, enum print_output printmode)
{
  Log *l = emalloc (sizeof(Log));

  if (printmode == CK_ENV) {
    printmode = get_env_printmode();
  }

  l->lfile = lfile;
  l->lfun = lfun;
  l->close = close;
  l->mode = printmode;
  list_add_end (sr->loglst, l);
  return;
}

void log_srunner_start (SRunner *sr)
{
  srunner_send_evt (sr, NULL, CLSTART_SR);
}

void log_srunner_end (SRunner *sr)
{
  srunner_send_evt (sr, NULL, CLEND_SR);
}

void log_suite_start (SRunner *sr, Suite *s)
{
  srunner_send_evt (sr, s, CLSTART_S);
}

void log_suite_end (SRunner *sr, Suite *s)
{
  srunner_send_evt (sr, s, CLEND_S);
}

void log_test_end (SRunner *sr, TestResult *tr)
{
  srunner_send_evt (sr, tr, CLEND_T);
}

static void srunner_send_evt (SRunner *sr, void *obj, enum cl_event evt)
{
  List *l;
  Log *lg;
  l = sr->loglst;
  for (list_front(l); !list_at_end(l); list_advance(l)) {
    lg = list_val(l);
    fflush(lg->lfile);
    lg->lfun (sr, lg->lfile, lg->mode, obj, evt);
    fflush(lg->lfile);
  }
}

void stdout_lfun (SRunner *sr, FILE *file, enum print_output printmode,
		  void *obj, enum cl_event evt)
{
  TestResult *tr;
  Suite *s;
  
  if (printmode == CK_ENV) {
    printmode = get_env_printmode();
  }

  switch (evt) {
  case CLINITLOG_SR:
    break;
  case CLENDLOG_SR:
    break;
  case CLSTART_SR:
    if (printmode > CK_SILENT) {
      fprintf(file, "Running suite(s):");
    }
    break;
  case CLSTART_S:
    s = obj;
    if (printmode > CK_SILENT) {
      fprintf(file, " %s\n", s->name);
    }
    break;
  case CLEND_SR:
    if (printmode > CK_SILENT) {
      /* we don't want a newline before printing here, newlines should
	 come after printing a string, not before.  it's better to add
	 the newline above in CLSTART_S.
      */
      srunner_fprint (file, sr, printmode);
    }
    break;
  case CLEND_S:
    s = obj;
    break;
  case CLEND_T:
    tr = obj;
    break;
  default:
    eprintf("Bad event type received in stdout_lfun", __FILE__, __LINE__);
  }

  
}

void lfile_lfun (SRunner *sr, FILE *file, enum print_output printmode,
		 void *obj, enum cl_event evt)
{
  TestResult *tr;
  Suite *s;
  
  switch (evt) {
  case CLINITLOG_SR:
    break;
  case CLENDLOG_SR:
    break;
  case CLSTART_SR:
    break;
  case CLSTART_S:
    s = obj;
    fprintf(file, "Running suite %s\n", s->name);
    break;
  case CLEND_SR:
    fprintf (file, "Results for all suites run:\n");
    srunner_fprint (file, sr, CK_MINIMAL);
    break;
  case CLEND_S:
    s = obj;
    break;
  case CLEND_T:
    tr = obj;
    tr_fprint(file, tr, CK_VERBOSE);
    break;
  default:
    eprintf("Bad event type received in stdout_lfun", __FILE__, __LINE__);
  }

  
}

void xml_lfun (SRunner *sr, FILE *file, enum print_output printmode,
		  void *obj, enum cl_event evt)
{
  TestResult *tr;
  Suite *s;
  static struct timeval inittv, endtv;
  static char t[sizeof "yyyy-mm-dd hh:mm:ss"] = {0};

  if (t[0] == 0)
  {
    struct tm now;
    gettimeofday(&inittv, NULL);
    localtime_r(&(inittv.tv_sec), &now);
    strftime(t, sizeof("yyyy-mm-dd hh:mm:ss"), "%Y-%m-%d %H:%M:%S", &now);
  }

  switch (evt) {
  case CLINITLOG_SR:
    fprintf(file, "<?xml version=\"1.0\"?>\n");
    fprintf(file, "<testsuites xmlns=\"http://check.sourceforge.net/ns\">\n");
    fprintf(file, "  <datetime>%s</datetime>\n", t);
    break;
  case CLENDLOG_SR:
    gettimeofday(&endtv, NULL);
    fprintf(file, "  <duration>%f</duration>\n",
        (endtv.tv_sec + (float)(endtv.tv_usec)/1000000) - \
        (inittv.tv_sec + (float)(inittv.tv_usec/1000000)));
    fprintf(file, "</testsuites>\n");
    break;
  case CLSTART_SR:
    break;
  case CLSTART_S:
    s = obj;
    fprintf(file, "  <suite>\n");
    fprintf(file, "    <title>%s</title>\n", s->name);
    break;
  case CLEND_SR:
    break;
  case CLEND_S:
    fprintf(file, "  </suite>\n");
    s = obj;
    break;
  case CLEND_T:
    tr = obj;
    tr_xmlprint(file, tr, CK_VERBOSE);
    break;
  default:
    eprintf("Bad event type received in xml_lfun", __FILE__, __LINE__);
  }

}


FILE *srunner_open_lfile (SRunner *sr)
{
  FILE *f = NULL;
  if (srunner_has_log (sr)) {
    f = fopen(sr->log_fname, "w");
    if (f == NULL)
      eprintf ("Could not open log file %s:", __FILE__, __LINE__,
	       sr->log_fname);
  }
  return f;
}

FILE *srunner_open_xmlfile (SRunner *sr)
{
  FILE *f = NULL;
  if (srunner_has_xml (sr)) {
    f = fopen(sr->xml_fname, "w");
    if (f == NULL)
      eprintf ("Could not open xml file %s:", __FILE__, __LINE__,
	       sr->xml_fname);
  }
  return f;
}

void srunner_init_logging (SRunner *sr, enum print_output print_mode)
{
  FILE *f;
  sr->loglst = check_list_create();
  srunner_register_lfun (sr, stdout, 0, stdout_lfun, print_mode);
  f = srunner_open_lfile (sr);
  if (f) {
    srunner_register_lfun (sr, f, 1, lfile_lfun, print_mode);
  }
  f = srunner_open_xmlfile (sr);
  if (f) {
    srunner_register_lfun (sr, f, 2, xml_lfun, print_mode);
  }
  srunner_send_evt (sr, NULL, CLINITLOG_SR);
}

void srunner_end_logging (SRunner *sr)
{
  List *l;
  int rval;

  srunner_send_evt (sr, NULL, CLENDLOG_SR);

  l = sr->loglst;
  for (list_front(l); !list_at_end(l); list_advance(l)) {
    Log *lg = list_val(l);
    if (lg->close) {
      rval = fclose (lg->lfile);
      if (rval != 0)
	eprintf ("Error closing log file:", __FILE__, __LINE__);
    }
    free (lg);
  }
  list_free(l);
  sr->loglst = NULL;
}
