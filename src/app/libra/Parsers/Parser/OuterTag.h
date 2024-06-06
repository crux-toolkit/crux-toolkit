#ifndef OUTER_TAG_H
#define OUTER_TAG_H
/*
Program       : OuterTag
Author        : ?
Date          : ?
SVN info      : $Id: OuterTag.h 7714 2017-12-20 00:16:03Z real_procopio $


Copyright (C) 2003 Andrew Keller

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Andrew Keller
Institute for Systems Biology
401 Terry Avenue North
Seattle, WA  98109  USA
akeller@systemsbiology.org
*/

#include "Tag.h"
#include <assert.h>

class OuterTag : public Tag {

 public:

  OuterTag(const Tag* tag);
  ~OuterTag();
  void enter(const char* namespace_prefix, const char* namesp, const char* ref_schema);
  void enterRefs(const OuterTag* outer);
  int getNumRefs() const;
  int getNumNameSpaces() const;
  const char* getNamespace_prefix(int k) const;
  const char* getNameSpace(int k) const;
  const char* getRefSchema(int k) const;
  const char* getRefNameSpace(int k) const;


 protected:

  StringArray namespace_prefixes_;
  StringArray namespaces_;
  StringArray referenced_schemas_;
  StringArray referenced_namespaces_;

};
#endif
