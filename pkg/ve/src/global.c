/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */

#define DEBUG
#include "global.h"

/* This is here to be a breakpoint in the debugger, and
for no other reason */

void myexit(rc)
     int rc;
{
  exit(rc);
}


