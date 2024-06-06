/*
    Program: seeker
    Description: simple utility to display file at specific byte offsets
    Date: June 9, 2008

    Copyright (C) 2008  Natalie Tasman (original author), ISB Seattle

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
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

    Natalie Tasman
    Institute for Systems Biology
    401 Terry Avenue North
    Seattle, WA  98109  USA
    
    email (remove underscores):
    n_tasman at systems_biology dot org
*/

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
  FILE* fp;
  char buf[1024];
  int done = 0;
  int offset;
  int seekok = -1;
  size_t readok = -1;

  if (argc != 2) {
    printf("usage: seeker [filename]\n\ndumps text of file at given byte offset\n");
    exit(0);
  }
  printf("opening file %s\n", argv[1]);
  fflush(0);
  fp=fopen(argv[1], "r");
  if (fp == NULL) {
    printf("unable to open file %s\n", argv[1]);
    exit(1);
  }
  while (!done) {
    printf("enter offset to seek to or 'q' to quit: ");
    fflush(0);
    char* errCheckMe = NULL;
    errCheckMe = fgets(buf,1024,fp);
    if (buf[0] == 'q') {
      done = 1;
      continue;
    }
    offset = atoi(buf);
    printf("\nseeking to %d", offset);
    fflush(0);

    seekok = fseek(fp, offset, SEEK_SET);
    if (seekok == -1) {
      printf("(unable to seek to %d\n", offset);
      continue;
    }

    readok = fread(buf, sizeof(char), 1000, fp);
    if (readok < 0) { 
      readok = 0;
    }
    if (readok < 1000) {
      buf[readok] = 0;
    }
    printf("read:\n***from here:***");
    fflush(0);
    printf("%s\n", buf);
    fflush(0);
  }
  

  return 0;
}
