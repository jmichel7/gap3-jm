/* ========================== C MeatAxe =============================
   intstr.c - The integer and string data type.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: intstr.c,v 1.2 1997/09/11 15:42:55 gap Exp $
 *
 * $Log: intstr.c,v $
 * Revision 1.2  1997/09/11 15:42:55  gap
 * New version 2.2.3. AH
 *
 * Revision 1.3  1994/07/10  15:33:29  mringe
 * Compiler warnings.
 *
 * Revision 1.3  1994/07/10  15:33:29  mringe
 * Compiler warnings.
 *
 * Revision 1.2  1994/06/29  07:52:12  mringe
 * Neu: dup und print.
 *
 * Revision 1.1  1994/06/25  13:34:43  mringe
 * Initial revision
 *
 */

#include <string.h>
#include <stdlib.h>
#include "meataxe.h"



/* ------------------------------------------------------------------
   Integers
   ------------------------------------------------------------------ */

char *intdesc(s)
integer_t *s;

{
    static char buf[50];
    sprintf(buf,"integer (%ld)",s->l);
    return buf;
}


void intprint(name,s)
char *name;
integer_t *s;

{
    if (name != NULL) printf("%s=",name);
    printf("%ld",s->l);
    if (name != NULL) printf("\n");
    fflush(stdout);
}


integer_t *intalloc(l)
long l;

{
    integer_t *n;

    n = ALLOC(integer_t);
    if (n != NULL)
    {
        n->id = T_INTEGER;
	n->l = l;
	return n;
    }
    MTXFAIL(ERR_NOMEM,(integer_t *)NULL);
}

integer_t *intdup(s)
integer_t *s;

{
    return intalloc(s->l);
}


void intfree(s)
integer_t *s;

{
    free(s);
}


integer_t *intread(f)
FILE *f;

{
    long l;

    if (zreadlong(f,&l,1) != 1)
	MTXFAIL(ERR_FILEREAD,NULL);
    return intalloc(l);
}


int intwrite(f, s)
FILE *f;
integer_t *s;

{
    long l = s->l;

    if (zwritelong(f,&l,1) != 1)
	return -1;
    return 0;
}


integer_t *intadd(d,s)
integer_t *d, *s;
{
    d->l += s->l;
    return d;
}


integer_t *intsub(d,s)
integer_t *d, *s;
{
    d->l -= s->l;
    return d;
}


integer_t *intmul(d,s)
integer_t *d, *s;
{
    d->l *= s->l;
    return d;
}


integer_t *intdiv(d,s)
integer_t *d, *s;
{
    d->l /= s->l;
    return d;
}


integer_t *intneg(x)
integer_t *x;
{
    x->l = - x->l;
    return x;
}


/* ------------------------------------------------------------------
   Strings
   ------------------------------------------------------------------ */

char *stringdesc(s)
string_t *s;

{
    static char buf[50];
    sprintf(buf,"string (%s)",s->s);
    return buf;
}


void stringprint(name,s)
char *name;
string_t *s;

{
    if (name != NULL) printf("%s=",name);
    printf("%s",s->s);
    if (name != NULL) printf("\n");
    fflush(stdout);
}


string_t *stringalloc(c)
char *c;

{
    string_t *n;

    n = ALLOC(string_t);
    if (n != NULL)
    {
        n->id = T_STRING;
	if (c == NULL) return n;
	n->s = Smalloc(strlen(c)+1);
	if (n->s != NULL)
	{
	    strcpy(n->s,c);
	    return n;
	}
	free(n);
    }
    MTXFAIL(ERR_NOMEM,(string_t *)NULL);
}

string_t *stringdup(s)
string_t *s;

{
    return stringalloc(s->s);
}


void stringfree(s)
string_t *s;

{
    free(s->s);
    free(s);
}


string_t *stringread(f)
FILE *f;

{
    long l;
    char *c;
    string_t *n;

    if (zreadlong(f,&l,1) != 1)
	MTXFAIL(ERR_FILEREAD,NULL);
    if ((c = Smalloc((size_t)l)) == NULL)
	MTXFAIL(ERR_NOMEM,NULL);
    if (fread(c,1,(size_t)l,f) != (size_t)l)
	MTXFAIL(ERR_FILEREAD,NULL);
    if ((n = stringalloc(NULL)) == NULL)
	return NULL;
    n->s = c;
    return n;
}


int stringwrite(f, s)
FILE *f;
string_t *s;

{
    long l = strlen(s->s) + 1;

    if (zwritelong(f,&l,1) != 1) return -1;
    if (fwrite(s->s,1,(size_t)l,f) != (size_t) l)
	MTXFAIL(ERR_FILEWRITE,-1);
    return 0;
}


