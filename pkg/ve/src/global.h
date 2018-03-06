/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */
#ifndef GLOBAL_H
#define GLOBAL_H 1
#include <stdio.h>
/* #include <malloc.h> */

typedef unsigned int bool;

#define true 1
#define false 0

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)>(b)?(b):(a))

#ifdef __STDC__
#include <stdlib.h>
#define PT(ARGS) ARGS
#define FuncPt &
typedef void * pointer;
#else
#define PT(ARGS) ()
#define FuncPt
typedef char * pointer;
#define inline
#define const
#endif



#ifndef DEBUG
#define myexit(rc) exit(rc)
#else
extern void myexit();
#endif

#endif
