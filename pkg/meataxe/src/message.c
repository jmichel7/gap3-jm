/* ========================== C MeatAxe =============================
   message.c - Messages.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: message.c,v 1.2 1997/09/11 15:43:00 gap Exp $
 *
 * $Log: message.c,v $
 * Revision 1.2  1997/09/11 15:43:00  gap
 * New version 2.2.3. AH
 *
 * Revision 2.5  1994/12/17  20:33:08  mringe
 * Neu: Message() und ErrorExit().
 *
 * Revision 2.4  1994/07/05  15:28:05  mringe
 * Neu: ERR_UNKNOWN
 *
 * Revision 2.3  1994/05/19  11:36:27  mringe
 * Neue Werte fuer mtxerraction.
 *
 * Revision 2.2  1994/05/18  10:08:01  mringe
 * return statement eingebaut, um gcc-Warnung zu vermeiden.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.14  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.13  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.12  1993/10/05  19:05:02  mringe
 * Neue Fehlermeldungen. yerrno eliminiert.
 *
 * Revision 1.11  1993/10/02  15:26:39  mringe
 * Fehlermeldungen fuer Schrieb-/Lesefehler.
 *
 * Revision 1.10  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.9  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 */

#include <stdarg.h>
#include <stdlib.h>
#include "meataxe.h"


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

int msg_level = 0;
int mtxerrno = 0;
int mtxerraction = 1;

static char *mtxerrfile;
static int mtxerrline;

static struct msg_struct { int errno; char *smsg, *lmsg; } msgs[] =
  {
    { ERR_NOMEM, "Not enough memory", NULL },
    { ERR_GAMEOVER, "Time limit exceeded", NULL },
    { ERR_DIV0, "Division by zero", NULL },

    { ERR_FILEOPEN, "Could not open file", NULL },
    { ERR_FILEREAD, "File read error", NULL },
    { ERR_FILEWRITE, "File write error", NULL },
    { ERR_FILEFMT, "Bad file format", NULL },
    { ERR_NFILES, "Too many files", NULL },

    { ERR_BADARG, "Bad argument", NULL },
    { ERR_BADTYPE, "Bad type", NULL },
    { ERR_RANGE, "Argument out of range", NULL },
    { ERR_BADARG, "Bad argument", NULL },
    { ERR_NOTECH, "Matrix not in echelon form", NULL },
    { ERR_NOTSQUARE, "Matrix not square", NULL },
    { ERR_INCOMPAT, "Incompatible objects", NULL },

    { ERR_BADUSAGE, "Bad syntax, try `-help'", NULL},
    { ERR_OPTION, "Bad usage of option, try `-help'", NULL},
    { ERR_NARGS, "Bad number of arguments, try `-help'", NULL},

    { ERR_NOTMATRIX, "Not a matrix", NULL},
    { ERR_NOTPERM, "Not a permutation", NULL},
    { ERR_NOTPOLY, "Not a polynomial", NULL},
    
    { ERR_UNKNOWN, "Unknown symbol", NULL},
    { -1, NULL, NULL }
  };


/* ------------------------------------------------------------------
   errmsg()
   ------------------------------------------------------------------ */

static char *errmsg(en)
int en;
{
    static char buf[30];
    struct msg_struct *x;

    for (x = msgs; x->errno >= 0; ++x)
        if (x->errno == en)
	    return x->smsg;

    sprintf(buf,"Unknown error %d",en);
    return buf;
}



/* ------------------------------------------------------------------
   errhandler() - MeatAxe error handler
   ------------------------------------------------------------------ */

int errhandler(errno,errfile,errline)
int errno;
char *errfile;
int errline;

{
    if (mtxerraction == 0) return 0;
    mtxerrno = errno;
    mtxerrfile = errfile;
    mtxerrline = errline;
    if (mtxerraction <= 1) return 0;
    fprintf(stderr,"%s(%d): %s.\n",mtxerrfile,mtxerrline,
	errmsg(errno));
    if (mtxerraction <= 2) return 1;
    fprintf(stderr,"Program halted\n");
    exit(1);
    return 0;
}


/* ------------------------------------------------------------------
   mtxerror() - Print error message with prefix.
   ------------------------------------------------------------------ */

void mtxerror(text)
char *text;

{
    fprintf(stderr,"%s: %s.\n",text,errmsg(mtxerrno));
}


/* ------------------------------------------------------------------
   errexit() - Terminate program with error message
   ------------------------------------------------------------------ */

void errexit(code,text)
int code;
char *text;

{
    if (code >= 0)
	mtxerrno = code;
    mtxerror(text);
    exit(EXIT_ERR);
}


/* ------------------------------------------------------------------
   printmessage() - Print a message (internal function)
   ------------------------------------------------------------------ */

static void printmessage(FILE *f,const char *fmt, va_list vl)

{
    long l;
    int i;
    char *c;

    while (*fmt != 0)
    {
	if (*fmt == '%')
	{
	    switch (*++fmt)
	    {
		case 'd':
		    l = va_arg(vl,long);
		    fprintf(f,"%ld",l);
		    break;
		case 's':
		    c = va_arg(vl,char *);
		    fprintf(f,"%s",c);
		    break;
		case 'E':
		    i = va_arg(vl,int);
		    fprintf(f,"%s",errmsg(i));
		    break;
	    }
	    ++fmt;
	}
	else
	    putc(*fmt++,f);
    }
}


/* ------------------------------------------------------------------
   Message() - Print a message
   ------------------------------------------------------------------ */

int Message(FILE *f, const char *fmt, ...)

{
    va_list vl;

    va_start(vl,fmt);
    printmessage(f,fmt,vl);
    va_end(vl);
    return 0;
}



/* ------------------------------------------------------------------
   ErrorExit() - Print an error message and exit
   ------------------------------------------------------------------ */

int ErrorExit(const char *fmt, ...)

{
    va_list vl;

    va_start(vl,fmt);
    printmessage(stderr,fmt,vl);
    putc('\n',stderr);
    va_end(vl);
    exit(EXIT_ERR);
}


